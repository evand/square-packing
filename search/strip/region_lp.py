#!/usr/bin/env python3
"""Weighted point-cover LP for the strip-leaf region lemma (float measurement, 10-11).

Region: R = [x0, T] x [0, T] minus the open notch int N, N = [x0, xn] x [yn, T].
Poses: closed unit squares Q inside R (closed containment) with Q ∩ int N = ∅.
Cover:  min sum w  s.t.  w(Q) >= 1 for every pose Q (closed: a point on dQ counts).
A cover of weight < 8 proves "no 8 pairwise closed-disjoint unit squares in R".

Finite LP over a candidate point set P and a pose sample S, grown by
  * pose separation (rows): dense lattice sweep (centres on a grid over the free centre region,
    with the box boundary and the notch-obstacle boundary sampled explicitly, angles on a grid),
    plus random local refinement around the worst poses;
  * point pricing (columns): dual coverage sum_{Q ∋ p} mu_Q on a fine candidate grid + special lines.
Directions: restricting points raises the LP value; restricting poses lowers it.
Independent of both: the dual measure mu, normalised by its max closed coverage M over the whole
region (grid with dilated squares, so M >= true max), gives L = sum(mu)/M <= nu_f <= any cover.

Usage: python3 search/strip/region_lp.py TAG x0 xn yn [--nonotch] [--open] [--T 4] [--iters N]
Outputs: runs/strip_region/TAG.{log,json}
"""
import argparse, ctypes, json, math, os, subprocess, sys, time
import numpy as np
from scipy.optimize import linprog
import scipy.sparse as sp

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "runs", "strip_region")
CSRC = r"""
#include <math.h>
#include <omp.h>
/* out[k] = sum_j w[j] * [ point j in square k ], square k: centre (cx,cy), angle th, half side hs */
void cov_poses(int n, const double*px, const double*py, const double*w,
               int m, const double*cx, const double*cy, const double*th, double hs, double*out, int nt){
  #pragma omp parallel for num_threads(nt) schedule(dynamic,256)
  for(int k=0;k<m;k++){ double c=cos(th[k]), s=sin(th[k]), acc=0;
    for(int j=0;j<n;j++){ double dx=px[j]-cx[k], dy=py[j]-cy[k];
      double u=dx*c+dy*s, v=-dx*s+dy*c;
      if(fabs(u)<=hs && fabs(v)<=hs) acc+=w[j]; }
    out[k]=acc; } }
/* out[j] = sum_k mu[k] * [ point j in square k ] */
void cov_points(int n, const double*px, const double*py, double*out,
                int m, const double*cx, const double*cy, const double*th, const double*mu, double hs, int nt){
  #pragma omp parallel for num_threads(nt) schedule(dynamic,256)
  for(int j=0;j<n;j++){ double acc=0;
    for(int k=0;k<m;k++){ double dx=px[j]-cx[k], dy=py[j]-cy[k];
      if(fabs(dx)>0.75||fabs(dy)>0.75) continue;
      double c=cos(th[k]), s=sin(th[k]);
      double u=dx*c+dy*s, v=-dx*s+dy*c;
      if(fabs(u)<=hs && fabs(v)<=hs) acc+=mu[k]; }
    out[j]=acc; } }
"""
NT = 4


def kernel():
    so = os.path.join(OUT, "_region_kern.so")
    src = os.path.join(OUT, "_region_kern.c")
    if not os.path.exists(so):
        os.makedirs(OUT, exist_ok=True)
        open(src, "w").write(CSRC)
        subprocess.check_call(["gcc", "-O3", "-march=native", "-fopenmp", "-shared", "-fPIC", src, "-o", so, "-lm"])
    lib = ctypes.CDLL(so)
    D = ctypes.POINTER(ctypes.c_double)
    lib.cov_poses.argtypes = [ctypes.c_int, D, D, D, ctypes.c_int, D, D, D, ctypes.c_double, D, ctypes.c_int]
    lib.cov_points.argtypes = [ctypes.c_int, D, D, D, ctypes.c_int, D, D, D, D, ctypes.c_double, ctypes.c_int]
    return lib


LIB = None


def _p(a):
    return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def cov_poses(P, w, S, hs):
    P = np.ascontiguousarray(P, float); w = np.ascontiguousarray(w, float); S = np.ascontiguousarray(S, float)
    cx, cy, th = [np.ascontiguousarray(S[:, i]) for i in range(3)]
    px, py = np.ascontiguousarray(P[:, 0]), np.ascontiguousarray(P[:, 1])
    out = np.zeros(len(S))
    LIB.cov_poses(len(P), _p(px), _p(py), _p(w), len(S), _p(cx), _p(cy), _p(th), hs, _p(out), NT)
    return out


def cov_points(Q, S, mu, hs):
    S = np.ascontiguousarray(S, float); mu = np.ascontiguousarray(mu, float)
    cx, cy, th = [np.ascontiguousarray(S[:, i]) for i in range(3)]
    px, py = np.ascontiguousarray(Q[:, 0]), np.ascontiguousarray(Q[:, 1])
    out = np.zeros(len(Q))
    LIB.cov_points(len(Q), _p(px), _p(py), _p(out), len(S), _p(cx), _p(cy), _p(th), _p(mu), hs, NT)
    return out


class Region:
    def __init__(s, x0, xn, yn, T=4.0, notch=True):
        s.x0, s.xn, s.yn, s.T, s.notch = x0, xn, yn, T, notch

    def box(s, th):
        r = 0.5 * (math.cos(th) + math.sin(th))
        return s.x0 + r, s.T - r, r, s.T - r

    def octagon(s, th):
        """N ⊕ Q0(th) as a CCW convex polygon (Q0 centred square)."""
        c, si = math.cos(th), math.sin(th)
        V = np.array([[x, y] for x in (s.x0, s.xn) for y in (s.yn, s.T)])
        # square vertices: ±e1/2 ± e2/2
        e1, e2 = np.array([c, si]), np.array([-si, c])
        sq = np.array([(a * e1 + b * e2) / 2 for a in (1, -1) for b in (1, -1)])
        pts = (V[:, None, :] + sq[None]).reshape(-1, 2)
        return hull(pts)

    def free(s, C, th, tol=1e-9):
        """C: (k,2) centres at angle th; True where the pose is feasible."""
        a, b, c, d = s.box(th)
        ok = (C[:, 0] >= a - tol) & (C[:, 0] <= b + tol) & (C[:, 1] >= c - tol) & (C[:, 1] <= d + tol)
        if s.notch:
            ok &= ~in_open_poly(C, s.octagon(th), tol)
        return ok

    def centres(s, th, h):
        """Centre sample for angle th: grid + box boundary + obstacle boundary + vertices."""
        a, b, c, d = s.box(th)
        if b < a or d < c:
            return np.zeros((0, 2))
        nx, ny = max(1, int(math.ceil((b - a) / h))), max(1, int(math.ceil((d - c) / h)))
        X, Y = np.meshgrid(np.linspace(a, b, nx + 1), np.linspace(c, d, ny + 1))
        G = [np.c_[X.ravel(), Y.ravel()]]
        if s.notch:
            O = s.octagon(th)
            for i in range(len(O)):
                p, q = O[i], O[(i + 1) % len(O)]
                L = np.linalg.norm(q - p)
                k = max(1, int(math.ceil(L / (h / 2))))
                t = np.linspace(0, 1, k + 1)[:, None]
                G.append(p + t * (q - p))
            # intersections of obstacle edges with box lines
            for i in range(len(O)):
                p, q = O[i], O[(i + 1) % len(O)]
                for xv in (a, b):
                    if (p[0] - xv) * (q[0] - xv) < 0:
                        t = (xv - p[0]) / (q[0] - p[0]); G.append((p + t * (q - p))[None])
                for yv in (c, d):
                    if (p[1] - yv) * (q[1] - yv) < 0:
                        t = (yv - p[1]) / (q[1] - p[1]); G.append((p + t * (q - p))[None])
        C = np.vstack(G)
        C[:, 0] = np.clip(C[:, 0], a, b); C[:, 1] = np.clip(C[:, 1], c, d)
        return C[s.free(C, th)]

    def inside_pts(s, Q, tol=1e-12):
        ok = (Q[:, 0] >= s.x0 - tol) & (Q[:, 0] <= s.T + tol) & (Q[:, 1] >= -tol) & (Q[:, 1] <= s.T + tol)
        if s.notch:
            ok &= ~((Q[:, 0] > s.x0 + tol) & (Q[:, 0] < s.xn - tol) & (Q[:, 1] > s.yn + tol) & (Q[:, 1] < s.T - tol))
        return ok


def hull(P):
    P = sorted(set(map(tuple, np.round(P, 15))))
    if len(P) < 3:
        return np.array(P)
    def cr(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lo, up = [], []
    for p in P:
        while len(lo) >= 2 and cr(lo[-2], lo[-1], p) <= 1e-15: lo.pop()
        lo.append(p)
    for p in reversed(P):
        while len(up) >= 2 and cr(up[-2], up[-1], p) <= 1e-15: up.pop()
        up.append(p)
    return np.array(lo[:-1] + up[:-1])


def in_open_poly(C, O, tol):
    ins = np.ones(len(C), bool)
    for i in range(len(O)):
        p, q = O[i], O[(i + 1) % len(O)]
        e = q - p; L = np.linalg.norm(e)
        crs = (e[0] * (C[:, 1] - p[1]) - e[1] * (C[:, 0] - p[0])) / L
        ins &= crs > tol
    return ins


def special_coords(lo, hi, anchors):
    ds = [0, 0.5, math.sqrt(2) / 2, 1, 1.5, 1 + math.sqrt(2) / 2, 2, 2.5, 2 + math.sqrt(2) / 2, 3,
          (1 + math.sqrt(2)) / 2, math.sqrt(2)]
    v = set()
    for a in anchors:
        for d in ds:
            for sgn in (1, -1):
                x = a + sgn * d
                if lo - 1e-12 <= x <= hi + 1e-12:
                    v.add(round(x, 12))
    return sorted(v)


def candidates(R, h):
    xs = set(np.round(np.arange(math.ceil(R.x0 / h) * h, R.T + 1e-9, h), 10)) | {R.x0, R.T}
    ys = set(np.round(np.arange(0, R.T + 1e-9, h), 10))
    sx = special_coords(R.x0, R.T, [R.x0, R.T] + ([R.xn] if R.notch else []))
    sy = special_coords(0, R.T, [0, R.T] + ([R.yn] if R.notch else []))
    xs = sorted(xs | set(sx)); ys = sorted(ys | set(sy))
    X, Y = np.meshgrid(xs, ys)
    Q = np.c_[X.ravel(), Y.ravel()]
    return Q[R.inside_pts(Q)]


METHOD = "highs"


def solve_lp(A, npts):
    res = linprog(np.ones(npts), A_ub=-A, b_ub=-np.ones(A.shape[0]), bounds=(0, None), method=METHOD,
                  options={"presolve": True})
    assert res.status == 0, res.message
    return res.fun, res.x, -res.ineqlin.marginals


def build_A(P, S, hs):
    """Sparse containment matrix rows=poses, cols=points."""
    rows, cols = [], []
    px, py = P[:, 0], P[:, 1]
    for k, (cx, cy, th) in enumerate(S):
        c, s = math.cos(th), math.sin(th)
        dx, dy = px - cx, py - cy
        m = (np.abs(dx * c + dy * s) <= hs) & (np.abs(-dx * s + dy * c) <= hs)
        j = np.nonzero(m)[0]
        rows.append(np.full(len(j), k)); cols.append(j)
    r = np.concatenate(rows); c = np.concatenate(cols)
    return sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(len(S), len(P)))


def pose_lattice(R, h, dth_deg, extra_axis=True):
    angs = list(np.arange(0, 90, dth_deg))
    if extra_axis:
        angs += [0.1, 0.25, 0.5, 1, 2, 89, 89.5, 89.75, 89.9]
    angs = sorted(set(round(a, 6) for a in angs if 0 <= a < 90))
    out = []
    for a in angs:
        th = math.radians(a)
        C = R.centres(th, h)
        out.append(np.c_[C, np.full(len(C), th)])
    return np.vstack(out)


def separate(R, P, w, hs, h, dth, rng, nlocal=4000, want=2000):
    """Return violated poses (coverage < 1 - 1e-7): lattice sweep + local refinement."""
    sup = w > 1e-12
    Ps, ws = P[sup], w[sup]
    S = pose_lattice(R, h, dth)
    cv = cov_poses(Ps, ws, S, hs)
    order = np.argsort(cv)
    worst = S[order[:200]]
    # local refinement: perturb the worst poses, project into the free set by rejection
    loc = []
    for scale in (h, h / 4, h / 16):
        k = nlocal
        base = worst[rng.integers(0, len(worst), k)]
        pert = base + rng.normal(0, 1, (k, 3)) * np.array([scale, scale, math.radians(dth) * scale / h])
        pert[:, 2] = np.mod(pert[:, 2], math.pi / 2)
        okm = np.zeros(k, bool)
        for i in range(k):
            okm[i] = R.free(pert[i:i + 1, :2], pert[i, 2])[0]
        # also snap to walls: push to box edges
        loc.append(pert[okm])
        snap = base.copy()
        for i in range(len(snap)):
            a, b, c, d = R.box(snap[i, 2])
            j = rng.integers(0, 4)
            if j == 0: snap[i, 0] = a
            elif j == 1: snap[i, 0] = b
            elif j == 2: snap[i, 1] = c
            else: snap[i, 1] = d
        okm = np.array([R.free(snap[i:i + 1, :2], snap[i, 2])[0] for i in range(len(snap))])
        loc.append(snap[okm])
    L = np.vstack(loc)
    cvl = cov_poses(Ps, ws, L, hs)
    allS = np.vstack([S, L]); allc = np.r_[cv, cvl]
    bad = np.nonzero(allc < 1 - 1e-7)[0]
    bad = bad[np.argsort(allc[bad])]
    # diversity: greedy pick, skip poses too close to picked ones
    pick = []
    for i in bad:
        p = allS[i]
        if pick:
            Pk = allS[pick]
            dd = np.abs(Pk[:, 0] - p[0]) + np.abs(Pk[:, 1] - p[1]) + 0.5 * np.abs(Pk[:, 2] - p[2])
            if dd.min() < 1e-9:
                continue
        pick.append(i)
        if len(pick) >= want:
            break
    return allS[pick], (allc.min() if len(allc) else 1.0), len(bad)


def arrangement_vertices(S):
    """Corners of the squares in S and all pairwise edge intersections (float)."""
    cx, cy, th = S[:, 0], S[:, 1], S[:, 2]
    c, s = np.cos(th), np.sin(th)
    e1 = np.c_[c, s] / 2; e2 = np.c_[-s, c] / 2
    Cn = np.c_[cx, cy]
    V = np.stack([Cn + e1 + e2, Cn - e1 + e2, Cn - e1 - e2, Cn + e1 - e2], 1)  # (m,4,2)
    E0 = V.reshape(-1, 2); E1 = np.roll(V, -1, axis=1).reshape(-1, 2)
    out = [E0]
    m = len(E0)
    for i in range(0, m, 400):
        a0, a1 = E0[i:i + 400, None, :], E1[i:i + 400, None, :]
        d1 = a1 - a0; d2 = (E1 - E0)[None]
        den = d1[..., 0] * d2[..., 1] - d1[..., 1] * d2[..., 0]
        r = E0[None] - a0
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (r[..., 0] * d2[..., 1] - r[..., 1] * d2[..., 0]) / den
            u = (r[..., 0] * d1[..., 1] - r[..., 1] * d1[..., 0]) / den
        ok = (np.abs(den) > 1e-12) & (t >= -1e-12) & (t <= 1 + 1e-12) & (u >= -1e-12) & (u <= 1 + 1e-12)
        ii, jj = np.nonzero(ok)
        out.append(a0[ii, 0] + t[ii, jj, None] * d1[ii, 0])
    return np.vstack(out)


def main():
    global LIB, NT, METHOD
    ap = argparse.ArgumentParser()
    ap.add_argument("tag"); ap.add_argument("x0", type=float); ap.add_argument("xn", type=float)
    ap.add_argument("yn", type=float)
    ap.add_argument("--T", type=float, default=4.0)
    ap.add_argument("--nonotch", action="store_true")
    ap.add_argument("--open", action="store_true", help="open semantics: boundary points do not count")
    ap.add_argument("--iters", type=int, default=40)
    ap.add_argument("--h0", type=float, default=0.1, help="initial point grid")
    ap.add_argument("--hc", type=float, default=0.01, help="candidate (pricing) point grid")
    ap.add_argument("--hs", type=float, default=0.02, help="pose centre lattice for separation")
    ap.add_argument("--dth", type=float, default=0.5, help="pose angle lattice (deg)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--addpts", type=int, default=1000)
    ap.add_argument("--method", default="highs")
    ap.add_argument("--resume", default=None, help="npz with P, S to start from")
    ap.add_argument("--keep_dc", type=float, default=0.7, help="prune points with dual coverage below this")
    ap.add_argument("--keep_cov", type=float, default=1.25, help="prune poses with primal coverage above this")
    a = ap.parse_args()
    NT = a.threads; METHOD = a.method
    LIB = kernel()
    R = Region(a.x0, a.xn, a.yn, a.T, notch=not a.nonotch)
    hsz = 0.5 - 1e-9 if a.open else 0.5 + 1e-9
    rng = np.random.default_rng(1)
    logf = open(os.path.join(OUT, a.tag + ".log"), "w")
    def log(*x):
        s = " ".join(str(y) for y in x); print(s, flush=True); logf.write(s + "\n"); logf.flush()
    log("args", vars(a))
    t0 = time.time()
    P = candidates(R, a.h0)
    Qc = candidates(R, a.hc)
    S = pose_lattice(R, 0.1, 5.0)
    if a.resume:
        st = np.load(a.resume); P, S = st["P"], st["S"]
    log(f"init points {len(P)}  candidates {len(Qc)}  poses {len(S)}")
    A = build_A(P, S, hsz)
    hist = []
    for it in range(a.iters):
        tl=time.time(); val, w, mu = solve_lp(A, len(P)); tl=time.time()-tl
        # pricing
        Smu = S[mu > 1e-12]; mus = mu[mu > 1e-12]
        Vx = arrangement_vertices(Smu); Vx = Vx[R.inside_pts(Vx, tol=1e-9)]
        Qi = np.vstack([Qc, Vx])
        dc = cov_points(Qi, Smu, mus, hsz)
        known = set(map(tuple, np.round(P, 10)))
        newp_idx = [i for i in np.argsort(-dc) if dc[i] > 1 + 1e-7 and tuple(np.round(Qi[i], 10)) not in known][:a.addpts]
        # separation
        ts=time.time(); newS, cmin, nbad = separate(R, P, w, hsz, a.hs, a.dth, rng); ts=time.time()-ts
        Lb = val / max(dc.max(), 1.0)
        log(f"it {it:2d}  LP {val:.5f}  Ldual {Lb:.5f}  pts {len(P)} (supp {int((w>1e-9).sum())})  poses {len(S)} "
            f"(dual supp {len(mus)})  maxdual {dc.max():.4f} +pts {len(newp_idx)}  mincov {cmin:.4f} "
            f"viol {nbad} +poses {len(newS)}  nnz {A.nnz} tlp {tl:.0f} tsep {ts:.0f} [{time.time()-t0:.0f}s]")
        hist.append(dict(it=it, lp=val, Ldual=Lb, maxdual=float(dc.max()), mincov=float(cmin)))
        if not newp_idx and len(newS) == 0:
            break
        # prune: drop points far from pricing in and poses far from binding
        dcP = cov_points(P, Smu, mus, hsz)
        kp = (w > 1e-9) | (dcP > a.keep_dc)
        cS = cov_poses(P[w > 1e-12], w[w > 1e-12], S, hsz)
        ks = (mu > 1e-12) | (cS < a.keep_cov)
        P, S = P[kp], S[ks]
        if newp_idx:
            P = np.vstack([P, Qi[newp_idx]])
        if len(newS):
            S = np.vstack([S, newS])
        A = build_A(P, S, hsz)
        np.savez(os.path.join(OUT, a.tag + "_state.npz"), P=P, S=S)
    # final: certified-ish numbers
    val, w, mu = solve_lp(A, len(P))
    log(f"final LP over {len(P)} points, {len(S)} poses: {val:.5f}")
    sup = w > 1e-9
    # (i) cover side: scale w by the min coverage over a finer lattice -> float upper estimate for this point set
    Sfin = pose_lattice(R, a.hs / 2, a.dth / 2)
    cf = cov_poses(P[sup], w[sup], Sfin, hsz)
    up = val / max(cf.min(), 1e-9)
    # (ii) dual side: mu / max dilated coverage on a fine grid  -> lower bound on any cover (float)
    hg = 0.0025
    xs = np.arange(R.x0, R.T + 1e-12, hg); ys = np.arange(0, R.T + 1e-12, hg)
    X, Y = np.meshgrid(xs, ys); G = np.c_[X.ravel(), Y.ravel()]
    G = G[R.inside_pts(G, tol=hg)]
    dil = 0.5 + hg * math.sqrt(2) / 2 + 1e-9
    Smu = S[mu > 1e-12]; mus = mu[mu > 1e-12]
    M = cov_points(G, Smu, mus, dil).max()
    Mc = max(cov_points(Qc, Smu, mus, hsz).max(), cov_points(G, Smu, mus, hsz).max())
    low = mus.sum() / M
    Vx = arrangement_vertices(Smu)
    Mv = cov_points(Vx, Smu, mus, 0.5 + 1e-9)
    Mx = max(Mv.max(), Mc)
    lowx = mus.sum() / Mx
    imax = np.argsort(-Mv)[:10]
    log(f"dual: {len(mus)} poses, arrangement vertices {len(Vx)}, max closed coverage {Mx:.5f} at "
        f"{np.round(Vx[imax[0]], 4)} -> L_exact = {lowx:.5f}  (lower bound on any closed point cover, float)")
    log(f"FINAL LP {val:.5f}  cover-scaled (finer lattice, min cov {cf.min():.4f}) {up:.5f}  "
        f"dual mass {mus.sum():.5f}  M_dilated {M:.4f}  -> L {low:.5f}   (M_closed_grid {Mc:.4f} -> {mus.sum()/Mc:.5f})")
    np.savez(os.path.join(OUT, a.tag + ".npz"), P=P, w=w, S=S, mu=mu)
    json.dump(dict(args=vars(a), lp=val, cover_scaled=up, dual_mass=float(mus.sum()), M=float(M), L=low,
                   L_exact=lowx, M_exact=float(Mx), hist=hist, secs=time.time() - t0), open(os.path.join(OUT, a.tag + ".json"), "w"), indent=1)


if __name__ == "__main__":
    main()
