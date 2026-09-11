#!/usr/bin/env python3
"""Uniform certificates by ILP with lazy placement cuts.

A uniform certificate (k, m, s) is a set of m points in [0,s]^2 such that every closed unit
square inside the container, at every angle, contains at least k of them.  With m <= 12k-1 it
proves s(12) >= s (twelve disjoint squares would need 12k points).  See UNIFORM.md.

Method
------
  * candidate points: the grid  i*s/n, i = 0..n  (n even, so the centre s/2 and the boundary
    are on the grid), reduced to orbits of a symmetry group G (d4 = full symmetry of the
    container, c4 = rotations, diag = reflection (x,y)->(y,x), none).
  * one binary variable per orbit; objective = number of points; constraint per placement
    (cx, cy, theta): the closed square of half-side HCUT (default 0.499 < sigma/2 of the
    verifier's N=2000 net) centred at (cx,cy) at angle theta contains >= k chosen points.
    The cuts are (slightly strengthened) NECESSARY conditions, so the ILP optimum is a lower
    bound on the true minimum m over this grid, and the ILP is exact once the verifier passes.
  * lazy generation: (a) a float scan over an (eta, dt) lattice of placements with the current
    solution, worst placement per block; (b) the exact verifier's witness placements
    (verify <cert> 12 N threads topk witnessfile).  Repeat until the verifier prints VERIFIED
    or the ILP proves m > mmax.
  * --near CERT --radius R restricts the candidate orbits to those within L_inf distance R of a
    point of CERT (a local refinement on a finer grid).

Usage
    python3 search/uniform/ilp_uniform.py K S [--fine 0.05] [--mmax 12K-1] [--sym d4]
        [--rounds 60] [--tl 300] [--tag TAG] [--N 2000] [--hcut 0.499] [--near CERT --radius R]
        [--threads T] [--dmult M]
Output: runs/uni_TAG.txt (certificate, only written when VERIFIED), log on stdout.
"""
import numpy as np, scipy.sparse as sp, math, sys, os, time, subprocess, argparse
from scipy.optimize import milp, LinearConstraint, Bounds
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
VERIFY = os.path.join(REPO, "verify", "target", "release", "verify")
sys.path.insert(0, os.path.join(REPO, "search"))
import lp_search as LS


def read_cert(path):
    t = open(path).read().split()
    sn, sd, D, W, m = map(int, t[:5])
    rows = [(int(t[5+3*i]), int(t[6+3*i]), int(t[7+3*i])) for i in range(m)]
    return Fraction(sn, sd), D, W, rows


class Model:
    def __init__(self, k, s, fine, sym, near=None, radius=0.0, log=print, lines=None):
        self.k = k; self.s = s; self.log = log; self.feas = False
        n = int(round(s / fine)); n += n % 2
        self.n = n; self.step = s / n
        g = np.arange(n + 1) * self.step
        maps = {
            'd4': lambda i, j: [(i, j), (n-i, j), (i, n-j), (n-i, n-j), (j, i), (n-j, i), (j, n-i), (n-j, n-i)],
            'c4': lambda i, j: [(i, j), (n-j, i), (n-i, n-j), (j, n-i)],
            'diag': lambda i, j: [(i, j), (j, i)],
            'none': lambda i, j: [(i, j)],
        }[sym]
        self.sym = sym
        seen = set(); orbits = []
        nearP = None
        if near is not None:
            _, D, _, rows = read_cert(near)
            nearP = np.array([[x / D, y / D] for x, y, _ in rows])
        lineidx = None
        if lines:
            # candidate points restricted to the union of the vertical and horizontal lines
            # x = a or y = a, a in `lines` (mirrored: a -> s-a), snapped to the grid
            lineidx = set()
            for a in lines:
                for v in (a, s - a):
                    lineidx.add(int(round(v / self.step)))
            log(f"  lines (grid index): {sorted(lineidx)} -> x = {[round(i*self.step,4) for i in sorted(lineidx)]}")
        for i in range(n + 1):
            for j in range(n + 1):
                orb = sorted(set(maps(i, j)))
                if orb[0] in seen: continue
                seen.add(orb[0])
                if lineidx is not None and not (i in lineidx or j in lineidx):
                    continue
                if nearP is not None:
                    a, b = orb[0]
                    if np.min(np.max(np.abs(nearP - np.array([g[a], g[b]])), axis=1)) > radius + 1e-9:
                        continue
                orbits.append(np.array([[g[a], g[b]] for a, b in orb]))
        self.orbits = orbits
        self.P = np.concatenate(orbits)
        self.own = np.concatenate([np.full(len(o), t) for t, o in enumerate(orbits)])
        self.size = np.array([len(o) for o in orbits], dtype=float)
        self.Own = sp.csr_matrix((np.ones(len(self.P)), (np.arange(len(self.P)), self.own)),
                                 shape=(len(self.P), len(orbits)))
        self.rows = []      # list of csr blocks
        self.keys = set(); self.nrows = 0
        log(f"  grid n={n} step={self.step:.5f} points={len(self.P)} orbits={len(orbits)} sym={sym}")

    def add_placements(self, cx, cy, th, h):
        """closed squares of half-side h; returns number of new rows"""
        cx = np.asarray(cx, float); cy = np.asarray(cy, float); th = np.asarray(th, float)
        key = np.stack([np.round(cx / 0.004), np.round(cy / 0.004), np.round(th / 0.0015)], 1).astype(np.int64)
        keep = []
        for i, kk in enumerate(map(tuple, key)):
            if kk in self.keys: continue
            self.keys.add(kk); keep.append(i)
        if not keep: return 0
        keep = np.array(keep); cx, cy, th = cx[keep], cy[keep], th[keep]
        added = 0
        for t in np.unique(th):
            sel = np.nonzero(th == t)[0]
            ct, st = math.cos(t), math.sin(t)
            q0 = self.P[:, 0] * ct + self.P[:, 1] * st; q1 = -self.P[:, 0] * st + self.P[:, 1] * ct
            u0 = cx[sel] * ct + cy[sel] * st; u1 = -cx[sel] * st + cy[sel] * ct
            for a in range(0, len(sel), 512):
                M = (np.abs(q0[None, :] - u0[a:a+512, None]) <= h) & (np.abs(q1[None, :] - u1[a:a+512, None]) <= h)
                C = sp.csr_matrix(M.astype(float)) @ self.Own
                self.rows.append(C.tocsr()); added += C.shape[0]
        self.nrows += added
        return added

    def solve(self, mmax, time_limit, threads):
        A = sp.vstack(self.rows).tocsr()
        cons = [LinearConstraint(A, lb=np.full(A.shape[0], self.k), ub=np.inf)]
        if mmax is not None:
            cons.append(LinearConstraint(sp.csr_matrix(self.size[None, :]), lb=0, ub=mmax))
        t0 = time.time()
        obj = np.zeros(len(self.size)) if self.feas else self.size
        res = milp(obj, constraints=cons, integrality=np.ones(len(self.size)),
                   bounds=Bounds(0, 1), options={'time_limit': time_limit, 'disp': False, 'mip_rel_gap': 0.0})
        dt = time.time() - t0
        if res.x is None:
            return None, res, dt
        x = np.round(res.x).astype(int)
        return x, res, dt

    def points(self, x):
        sel = np.nonzero(x[self.own])[0]
        return self.P[sel]

    def write(self, x, path, dmult=1):
        """exact grid coordinates: D = n*s_den*dmult, X = i*s_num*dmult"""
        fr = Fraction(self.s).limit_denominator(10**6)
        D = self.n * fr.denominator * dmult
        pts = self.points(x)
        I = np.round(pts / self.step).astype(np.int64)
        assert np.max(np.abs(I * self.step - pts)) < 1e-9
        X = I * fr.numerator * dmult
        with open(path, 'w') as f:
            f.write(f"{fr.numerator} {fr.denominator}\n{D}\n{self.k}\n{len(X)}\n")
            for a, b in X: f.write(f"{a} {b} 1\n")
        return len(X)


def scan_cuts(m, x, h, eta, dt, perang=150, B=8):
    """float scan of the current solution over an (eta,dt) lattice; returns worst placements"""
    P = m.points(x); w = np.ones(len(P))
    tmax = math.pi / 4 if m.sym == 'd4' else math.pi / 2
    out = []; nviol = 0; gmin = 1e9
    for tm in np.arange(dt / 2, tmax + dt, dt):
        r = LS.scan_angle(P, w, m.s, tm, h + 1e-9, eta)
        if r is None: continue
        vals, g0, g1, box = r
        G0, G1 = np.meshgrid(g0, g1, indexing='ij')
        ok = LS.admissible(G0.ravel(), G1.ravel(), box, 0.0)
        v = np.where(ok, vals.ravel(), 1e9)
        gmin = min(gmin, v.min())
        bad = v < m.k - 0.5
        nviol += int(bad.sum())
        if not bad.any(): continue
        VV = v.reshape(vals.shape); n0, n1 = VV.shape
        p0 = (-n0) % B; p1 = (-n1) % B
        VP = np.pad(VV, ((0, p0), (0, p1)), constant_values=1e9)
        bl = VP.reshape(VP.shape[0] // B, B, VP.shape[1] // B, B).transpose(0, 2, 1, 3).reshape(-1, B * B)
        mn = bl.min(axis=1); am = bl.argmin(axis=1)
        sb = np.nonzero(mn < m.k - 0.5)[0]
        if len(sb) > perang: sb = sb[np.argsort(mn[sb])[:perang]]
        nb1 = VP.shape[1] // B
        bi, bj = np.divmod(sb, nb1); ii, jj = np.divmod(am[sb], B)
        r0 = bi * B + ii; r1 = bj * B + jj
        okk = (r0 < n0) & (r1 < n1); r0 = r0[okk]; r1 = r1[okk]
        sel = r0 * n1 + r1
        lo, hi, ct, st = box
        U0 = G0.ravel()[sel]; U1 = G1.ravel()[sel]
        cx = U0 * ct - U1 * st; cy = U0 * st + U1 * ct
        out.append((cx, cy, np.full(len(cx), tm)))
    if out:
        cx = np.concatenate([o[0] for o in out]); cy = np.concatenate([o[1] for o in out]); th = np.concatenate([o[2] for o in out])
    else:
        cx = cy = th = np.zeros(0)
    return cx, cy, th, nviol, gmin


def run_verifier(cert, N, threads, topk, wit):
    r = subprocess.run([VERIFY, cert, "12", str(N), str(threads), str(topk), wit], capture_output=True, text=True)
    ok = "VERIFIED" in r.stdout and "NOT VERIFIED" not in r.stdout
    ln = [l for l in r.stdout.split("\n") if l.startswith("min covered")]
    mc = ln[0].split("=")[1].strip() if ln else "?"
    return ok, mc, r.stdout


def read_wit(path):
    pts = []
    for line in open(path):
        q = line.split()
        if len(q) == 4 and float(q[0]) < 1.0: pts.append((float(q[0]), float(q[1]), float(q[2]), float(q[3])))
    pts.sort()
    return pts


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int); ap.add_argument("s", type=float)
    ap.add_argument("--fine", type=float, default=0.05); ap.add_argument("--mmax", type=int, default=None)
    ap.add_argument("--sym", default="d4"); ap.add_argument("--rounds", type=int, default=60)
    ap.add_argument("--tl", type=float, default=300); ap.add_argument("--tag", default="x")
    ap.add_argument("--N", type=int, default=2000); ap.add_argument("--hcut", type=float, default=0.499)
    ap.add_argument("--near", default=None); ap.add_argument("--radius", type=float, default=0.06)
    ap.add_argument("--threads", type=int, default=8); ap.add_argument("--dmult", type=int, default=1)
    ap.add_argument("--eta", type=float, default=0.02); ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--topk", type=int, default=6); ap.add_argument("--maxwit", type=int, default=3000)
    ap.add_argument("--runs", default=os.path.join(REPO, "runs"))
    ap.add_argument("--lines", default=None, help="comma list of line coordinates a; candidates restricted to x=a or y=a (and s-a)")
    ap.add_argument("--feas", action="store_true", help="feasibility only (zero objective) with the m <= mmax constraint")
    a = ap.parse_args(argv)
    lines = [float(v) for v in a.lines.split(",")] if a.lines else None
    k = a.k; mmax = a.mmax if a.mmax is not None else 12 * k - 1
    os.makedirs(a.runs, exist_ok=True)
    out = os.path.join(a.runs, f"uni_{a.tag}.txt"); snap = os.path.join(a.runs, f"unisnap_{a.tag}.txt")
    wit = os.path.join(a.runs, f"uniwit_{a.tag}.txt")
    t0 = time.time()
    def log(msg): print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)
    log(f"k={k} s={a.s} mmax={mmax} fine={a.fine} sym={a.sym} hcut={a.hcut} near={a.near}")
    m = Model(k, a.s, a.fine, a.sym, a.near, a.radius, log, lines)
    m.feas = a.feas
    # initial cuts: axis-aligned squares on a lattice of centres, plus a few angles
    h = a.hcut
    tmax = math.pi / 4 if a.sym == 'd4' else math.pi / 2
    init = []
    for tm in np.linspace(0, tmax, 4):
        ct, st = math.cos(tm), math.sin(tm); wid = ct + abs(st)
        c = np.arange(wid / 2, a.s - wid / 2 + 1e-9, 0.08 if tm == 0 else 0.16)
        CX, CY = np.meshgrid(c, c, indexing='ij')
        init.append((CX.ravel(), CY.ravel(), np.full(CX.size, tm)))
    for cx, cy, th in init: m.add_placements(cx, cy, th, h)
    log(f"initial rows={m.nrows}")
    status = "UNKNOWN"
    for it in range(a.rounds):
        x, res, dts = m.solve(mmax, a.tl, a.threads)
        if x is None:
            log(f"it{it}: ILP {res.message.strip()} ({dts:.1f}s) rows={m.nrows}  -> {'INFEASIBLE' if res.status==2 else 'no solution'}")
            status = "INFEASIBLE" if res.status == 2 else "TIMEOUT"
            break
        mm = int(round(m.size @ x))
        cx, cy, th, nviol, gmin = scan_cuts(m, x, h, a.eta, a.dt)
        m.write(x, snap, a.dmult)
        ok, mc, so = run_verifier(snap, a.N, a.threads, a.topk, wit)
        log(f"it{it}: m={mm} ILP={res.fun:.2f} status={res.status} bound={getattr(res,'mip_dual_bound',float('nan')):.2f} ({dts:.1f}s) rows={m.nrows} scan viol={nviol} min={gmin:.0f} verifier={mc} {'VERIFIED' if ok else ''}")
        if ok:
            os.replace(snap, out); status = "VERIFIED"
            log(f"*** VERIFIED uniform certificate k={k} m={mm} s={a.s} -> {out}")
            break
        n1 = m.add_placements(cx, cy, th, h)
        pts = read_wit(wit) if os.path.exists(wit) else []
        pts = pts[:a.maxwit]
        n2 = m.add_placements([p[2] for p in pts], [p[3] for p in pts], [p[1] for p in pts], h) if pts else 0
        log(f"      +{n1} scan cuts, +{n2} verifier cuts")
        if n1 + n2 == 0:
            log("no new cuts but not verified?!"); break
    print(f"RESULT k={k} s={a.s} mmax={mmax} fine={a.fine} sym={a.sym} status={status} time={time.time()-t0:.0f}s")
    return status


if __name__ == "__main__":
    main()
