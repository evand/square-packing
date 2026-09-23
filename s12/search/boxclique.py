#!/usr/bin/env python3
"""boxclique.py -- certificates with BOX-CLIQUE columns (certificates/FORMAT.md, "Clique certificates").

The point set of a certificate is kept (D4 orbits, as in tighten.py) and the weights are
re-optimised by the same cutting-plane loop with the exact verifier as separation oracle, except
that the LP also has CLIQUE columns.  A clique is a union of pose boxes (bin k of the N-net and a
rectangle of rotated centres); it costs its weight once and covers every row whose pose lies in
one of its boxes.  To keep the LP D4-symmetric (its point columns are orbits) a clique column is
a clique ORBIT: the eight images of a clique under the symmetries of the container, sharing one
weight, cost 8 (coincident images are merged with their multiplicity).

Pricing.  The row duals y of the symmetric LP are only constrained through their D4-average
ybar = (1/8) sum_g g.y, which IS a fractional packing (mass <= 1 on the poses through every atom).
A clique K with ybar(K) > 1 is therefore a genuine non-Helly gain over the current atom set, and
the orbit column of K has reduced cost 8 (1 - ybar(K)) < 0.  The dual is taken from an
interior-point solve without crossover (a central dual; a vertex dual of this very degenerate LP
is a poor pricing signal), the primal from a simplex solve (a vertex: sparse weights).  Candidate
cliques: the max-mass clique of the closed-overlap graph of the support images (squares shrunk by
eps, so that boxes of half-size eps around the poses have meeting cores), then the best clique
through each of the heaviest images not yet in one.  Every image clique is boxed (bin, rotated
centre +- eps, outward-rounded over Q), grown to +- grow*eps if the cores still meet, and pruned
to a family that passes the exact pairwise core test with the verifier's own integer sigma_k --
so every clique written is one the verifier accepts.

    python3 search/boxclique.py CERT TAG [--Dp D --mul M] [--N 2000] [--threads 4] [--rounds 40]
                                [--colgen 30] [--per-round 8] [--eps 0.005] [--Q 1000000] [--n 12]
    -> runs/boxclique_TAG.txt  (points + cliques; verified at N by the Rust verifier; then
       re-check with  python3 xcheck.py runs/boxclique_TAG.txt N --all --n 12)

A clique certificate is tied to its N (the boxes are defined on the net), so the only check that
means anything is at that N; there is no 2N confirmation as in tighten.py.  With --colgen 0 the
loop is tighten.py's plain re-optimisation on the full [0,90) sweep (a control run).
"""
import sys, os, math, time, argparse, json
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from fractions import Fraction as F
import warnings

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import tighten as T
import clique_check as CC
import xcheck as XC          # exact core / pairwise test (independent implementation, exact sigma)
try:
    import highspy
except ImportError:
    highspy = None

HIGHS = {'random_seed': 0}


# ----------------------------------------------------------------------------- bins and boxes
def bin_rot(N, k):
    """(cos, sin) of theta_k, exact"""
    g = N * N + k * k
    return F(N * N - k * k, g), F(2 * k * N, g)

def sigma_rust(N, k):
    """the verifier's sigma_k: rounded DOWN to 1e-6, as an exact Fraction (same integer formula)"""
    c0, s0, g0 = N * N - k * k, 2 * k * N, N * N + k * k
    k2 = k + 1; c1, s1, g1 = N * N - k2 * k2, 2 * k2 * N, N * N + k2 * k2
    cd = c0 * c1 + s0 * s1; sd = c0 * s1 - s0 * c1
    return F((g0 * g1 * 10 ** 6) // (cd + sd), 10 ** 6)

def row_bin(th, N):
    """bin of a pose at angle th in [0, pi/2) (verifier witnesses sit exactly at theta_k)"""
    t = N * math.tan(th / 2)
    k = int(round(t))
    if abs(t - k) > 1e-6: k = int(math.floor(t))
    return min(max(k, 0), N - 1)

def core_rust(box, Q, N):
    """core of a box with the verifier's rounded sigma (a subset of xcheck's exact core)"""
    _, k, lo0, hi0, lo1, hi1 = box
    h = sigma_rust(N, k) / 2
    return (F(hi0, Q) - h, F(lo0, Q) + h, F(hi1, Q) - h, F(lo1, Q) + h), bin_rot(N, k)

def boxes_pairwise_ok(boxes, Q, N):
    """every box non-empty core and every pair meets -- with the verifier's sigma, exactly"""
    cores = [core_rust(b, Q, N) for b in boxes]
    for c, _ in cores:
        if c[0] > c[1] or c[2] > c[3]: return False
    for i in range(len(cores)):
        for j in range(i, len(cores)):
            if not XC.cores_meet(cores[i][0], cores[i][1], cores[j][0], cores[j][1]): return False
    return True

def d4_images(cx, cy, th, s):
    """the 8 images of a pose under the symmetries of [0,s]^2, angles normalised to [0, pi/2)"""
    out = []
    for g in range(8):
        x, y, a = cx, cy, th
        if g & 1: x, y, a = s - x, y, -a
        if g & 2: x, y, a = x, s - y, -a
        if g & 4: x, y, a = y, x, math.pi / 2 - a
        a = a % (math.pi / 2)
        if a > math.pi / 2 - 1e-9: a = 0.0
        out.append((x, y, a))
    return out

def box_of(cx, cy, th, N, Q, e):
    """the box (bin of th, rotated centre +- e, rounded outwards over Q) around a pose"""
    k = row_bin(th, N); c, s = bin_rot(N, k)
    u0 = float(c) * cx + float(s) * cy; u1 = -float(s) * cx + float(c) * cy
    return (0, k, int(math.floor((u0 - e) * Q)), int(math.ceil((u0 + e) * Q)), int(math.floor((u1 - e) * Q)), int(math.ceil((u1 + e) * Q)))


class Cliques:
    """clique-orbit columns.  A column is a list of (multiplicity, boxes) -- the distinct images of
    one clique -- and costs the sum of the multiplicities (8)."""
    def __init__(self, N, Q):
        self.N = N; self.Q = Q; self.cols = []; self.keys = set(); self._rot = {}

    def rot(self, k):
        if k not in self._rot:
            c, s = bin_rot(self.N, k); self._rot[k] = (float(c), float(s))
        return self._rot[k]

    def row_geom(self, rows):
        if not rows: return np.zeros(0, dtype=int), np.zeros(0), np.zeros(0)
        R = np.array(rows)
        K = np.array([row_bin(th, self.N) for th in R[:, 2]], dtype=int)
        cs = np.array([self.rot(int(k)) for k in K])
        U0 = R[:, 0] * cs[:, 0] + R[:, 1] * cs[:, 1]; U1 = -R[:, 0] * cs[:, 1] + R[:, 1] * cs[:, 0]
        return K, U0, U1

    def member_mask(self, boxes, K, U0, U1):
        """rows whose pose lies (strictly, by a hair) inside one of the boxes"""
        Q = self.Q; mask = np.zeros(len(K), dtype=bool)
        for _, k, lo0, hi0, lo1, hi1 in boxes:
            mask |= (K == k) & (U0 >= lo0 / Q + 1e-12) & (U0 <= hi0 / Q - 1e-12) & (U1 >= lo1 / Q + 1e-12) & (U1 <= hi1 / Q - 1e-12)
        return mask

    def cost(self):
        return np.array([sum(mu for mu, _ in col) for col in self.cols], dtype=float)

    def matrix(self, rows):
        K, U0, U1 = self.row_geom(rows)
        R = []; C = []; V = []
        for j, col in enumerate(self.cols):
            coef = np.zeros(len(K))
            for mu, bx in col: coef += mu * self.member_mask(bx, K, U0, U1)
            idx = np.nonzero(coef)[0]
            R.append(idx); C.append(np.full(len(idx), j)); V.append(coef[idx])
        if not R: return sp.csr_matrix((len(rows), 0))
        return sp.coo_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))), shape=(len(rows), len(self.cols))).tocsr()

    def add(self, images):
        """images: list of box lists (one per symmetry, already certified); merged by identity"""
        merged = {}
        for bx in images:
            key = tuple(sorted(tuple(b[1:]) for b in bx))
            if key in merged: merged[key][0] += 1
            else: merged[key] = [1, [tuple(b) for b in bx]]
        colkey = tuple(sorted(merged.keys()))
        if colkey in self.keys: return False
        self.keys.add(colkey); self.cols.append([(mu, bx) for mu, bx in merged.values()]); return True


def closed_overlap_graph(cx, cy, th, hs):
    """adjacency of the CLOSED squares (centre, angle, half-side hs): separating-axis test on the
    four edge normals, non-strict"""
    n = len(cx); C = np.c_[cx, cy]
    D = np.linalg.norm(C[:, None, :] - C[None, :, :], axis=2)
    cand = D <= (hs[:, None] + hs[None, :]) * math.sqrt(2) + 1e-9
    np.fill_diagonal(cand, False)
    adj = np.zeros((n, n), dtype=bool)
    for i in range(n):
        js = np.nonzero(cand[i])[0]
        if js.size == 0: continue
        ok = np.ones(js.size, dtype=bool)
        for src in (0, 1):
            for q in range(2):
                if src == 0:
                    a = th[i] + q * math.pi / 2; ax = np.repeat(np.array([[math.cos(a), math.sin(a)]]), js.size, 0)
                else:
                    a = th[js] + q * math.pi / 2; ax = np.c_[np.cos(a), np.sin(a)]
                phi = np.arctan2(ax[:, 1], ax[:, 0])
                ei = hs[i] * (np.abs(np.cos(th[i] - phi)) + np.abs(np.sin(th[i] - phi)))
                ej = hs[js] * (np.abs(np.cos(th[js] - phi)) + np.abs(np.sin(th[js] - phi)))
                ci = C[i, 0] * ax[:, 0] + C[i, 1] * ax[:, 1]; cj = C[js, 0] * ax[:, 0] + C[js, 1] * ax[:, 1]
                ok &= np.abs(ci - cj) <= ei + ej - 1e-9
        adj[i, js[ok]] = True
    return adj & adj.T


def price_cliques(rows, y, s, N, Q, eps, per_round=8, grow=3.0, time_limit=30, log=print):
    """clique orbits of negative reduced cost 8 (1 - ybar(K)) on the symmetrised dual ybar"""
    sup = np.nonzero(y > 1e-9)[0]
    if len(sup) == 0: return []
    # symmetrise: 8 images of every support row, mass y/8, coincident images merged
    acc = {}
    for r in sup:
        cx, cy, th, h = rows[r]
        for x, yy, a in d4_images(cx, cy, th, s):
            key = (round(x, 9), round(yy, 9), round(a, 9), round(h, 9))
            acc[key] = acc.get(key, 0.0) + y[r] / 8
    P = list(acc.keys()); ys = np.array([acc[k] for k in P])
    cx = np.array([p[0] for p in P]); cy = np.array([p[1] for p in P]); th = np.array([p[2] for p in P]); h = np.array([p[3] for p in P])
    n = len(P)
    adj = closed_overlap_graph(cx, cy, th, h - eps - 1e-9)
    def certify(poses, e):
        """boxes around the poses, pruned greedily (heaviest first) to a pairwise-certified family"""
        kept = []; bx = []
        for i in poses:
            b = box_of(cx[i], cy[i], th[i], N, Q, e)
            if boxes_pairwise_ok(bx + [b], Q, N): kept.append(i); bx.append(b)
        return kept, bx
    found = []; covered = np.zeros(n, dtype=bool)
    bw, bc, nodes, dt = CC.max_weight_clique(adj, ys, time_limit=time_limit)
    seeds = [[int(v) for v in bc]] if bc else []
    if bc: covered[bc] = True
    order = np.argsort(-ys)
    for i in order[:per_round * 4]:
        if len(seeds) >= per_round + 1: break
        if covered[i]: continue
        nb = np.nonzero(adj[i])[0]
        if nb.size == 0: covered[i] = True; continue
        sub = adj[np.ix_(nb, nb)]
        bw2, bc2, _, _ = CC.max_weight_clique(sub, ys[nb], time_limit=max(5, time_limit // 4))
        mem = [int(i)] + [int(nb[j]) for j in bc2]
        seeds.append(mem); covered[mem] = True
    for mem in seeds:
        mem = sorted(mem, key=lambda i: -ys[i])
        kept, bx = certify(mem, eps)
        if not kept: continue
        mass = float(ys[kept].sum()); rc = 8.0 * (1.0 - mass)
        if rc >= -1e-6: continue
        # the eight images of the kept clique, each boxed and certified on its own (grown if possible)
        images = []
        for g in range(8):
            poses = []
            for i in kept:
                x, yy, a = d4_images(cx[i], cy[i], th[i], s)[g]; poses.append((x, yy, a))
            bxg = []
            for e in (grow * eps, eps):
                cand = [box_of(x, yy, a, N, Q, e) for x, yy, a in poses]
                if boxes_pairwise_ok(cand, Q, N): bxg = cand; break
            if not bxg:                       # prune per image
                for x, yy, a in poses:
                    b = box_of(x, yy, a, N, Q, eps)
                    if boxes_pairwise_ok(bxg + [b], Q, N): bxg.append(b)
            if bxg: images.append(bxg)
        if images: found.append((rc, mass, len(kept), len(mem), images))
    found.sort(key=lambda t: t[0])
    if found: log(f"    priced {len(found)} clique orbits: best ybar(K)={found[0][1]:.4f} ({found[0][2]} of {found[0][3]} members kept); max-mass clique of ybar {bw:.4f} ({len(bc)} poses, {nodes} nodes, {dt:.0f}s)")
    else: log(f"    no clique with ybar(K) > 1 (max-mass clique of ybar {bw:.4f} over {len(bc)} poses)")
    return [f[4] for f in found]


def export(m, x, cl, z, path, WD=10 ** 7, factor=1.0, up=True, N=2000):
    s = F(m.K, m.D); lines = []; tot = 0
    for k, o in enumerate(m.orb_int):
        if x[k] <= 0: continue
        v = x[k] * factor * WD; w = int(math.ceil(v - 1e-9)) if up else int(math.floor(v))
        if w <= 0: continue
        for X, Y in o: lines.append(f"{X} {Y} {w}"); tot += w
    cls = []
    for j, col in enumerate(cl.cols):
        if z[j] <= 0: continue
        for mu, bx in col:
            v = z[j] * mu * factor * WD; w = int(math.ceil(v - 1e-9)) if up else int(math.floor(v))
            if w <= 0: continue
            cls.append(f"{w} {len(bx)}\n" + "\n".join(f"{b[1]} {b[2]} {b[3]} {b[4]} {b[5]}" for b in bx)); tot += w
    with open(path, 'w') as f:
        f.write(f"{s.numerator} {s.denominator}\n{m.D}\n{WD}\n{len(lines)}\n" + "\n".join(lines) + "\n")
        if cls: f.write(f"cliques {N} {cl.Q} {len(cls)}\n" + "\n".join(cls) + "\n")
    return tot / WD, len(lines), len(cls)


def solve(m, cl, margin, want_dual=True, force=0.0):
    """primal from a simplex solve (a vertex); dual from an interior-point solve without crossover.
    force > 0: the FIRST clique-orbit column is held at weight >= force (a format demonstration:
    the LP would otherwise use a clique only when it is cheaper than points)."""
    A = m.matrix(); B = cl.matrix(m.rows); nr = A.shape[0]
    M = sp.hstack([A, B]).tocsr() if B.shape[1] else A
    c = np.concatenate([m.sizes, cl.cost()])
    bounds = [(0, None)] * M.shape[1]
    if force > 0 and B.shape[1]: bounds[len(m.sizes)] = (force, None)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = linprog(c=c, A_ub=-M, b_ub=-np.full(nr, 1.0 + margin), bounds=bounds, method='highs', options=dict(HIGHS))
    if not res.success: return None
    y = np.maximum(-res.ineqlin.marginals[:nr], 0.0)
    if want_dual and highspy is not None:
        Mc = M.tocsc(); lp = highspy.HighsLp()
        lp.num_col_ = M.shape[1]; lp.num_row_ = nr
        lp.col_cost_ = c; lp.col_lower_ = np.array([b[0] for b in bounds], dtype=float); lp.col_upper_ = np.full(M.shape[1], highspy.kHighsInf)
        lp.row_lower_ = np.full(nr, 1.0 + margin); lp.row_upper_ = np.full(nr, highspy.kHighsInf)
        lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
        lp.a_matrix_.start_ = Mc.indptr; lp.a_matrix_.index_ = Mc.indices; lp.a_matrix_.value_ = Mc.data
        h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('solver', 'ipm'); h.setOptionValue('run_crossover', 'off')
        h.setOptionValue('random_seed', HIGHS['random_seed'])
        h.passModel(lp); h.run()
        if h.getModelStatus() == highspy.HighsModelStatus.kOptimal:
            y = np.maximum(np.array(h.getSolution().row_dual), 0.0)
    nx = len(m.sizes)
    return res.fun, res.x[:nx], res.x[nx:], y


def run(cert, tag, Dp=None, mul=1, N=2000, threads=4, topk=6, rounds=40, colgen=30, per_round=8, eps=0.005, Q=10 ** 6,
        margin=2e-6, n=12, log=print, warm=True, force=0.0):
    m, w0, s = T.build_model(cert, Dp, mul)
    cl = Cliques(N, Q); sf = float(s)
    log(f"[{tag}] {cert}: {len(m.orbits)} orbits / {len(m.P)} atoms, container {s} = {sf:.7f}, input total {float(m.sizes @ w0):.6f}; N={N} eps={eps} Q={Q}")
    t0 = time.time()
    if warm:
        lr = T.lattice_rows(m, w0[m.own], thr=1.05, nproc=threads); m.add_rows(lr)
        log(f"[{tag}] warm start: {len(m.rows)} lattice rows ({time.time()-t0:.0f}s)")
    tmp = f"runs/boxclique_{tag}_probe.txt"; sep = f"runs/boxclique_{tag}_sep.txt"
    pm = margin / 2; hist = []
    x = z = None
    for it in range(rounds):
        out = solve(m, cl, margin, want_dual=(it < colgen), force=force)
        if out is None: log(f"  it{it} LP failed"); return None
        val, x, z, y = out
        zc = float(cl.cost() @ z) if len(z) else 0.0; ncols = len(cl.cols)
        tw = float(m.sizes @ x) + zc
        export(m, x, cl, z, tmp, WD=10 ** 12, factor=1.0 / (1.0 + pm), up=False, N=N)
        mv, ok = T.run_verifier(tmp, N, topk=topk, sep=sep, threads=threads, n=n)
        wit = T.read_witnesses(sep, N) if os.path.exists(sep) else []
        if os.path.exists(sep): os.remove(sep)
        ncl = 0
        if it < colgen:
            for images in price_cliques(m.rows, y, sf, N, Q, eps, per_round=per_round, log=log):
                ncl += cl.add(images)
        added = m.add_rows(wit)
        nzc = int((z > 1e-9).sum()) if len(z) else 0
        log(f"  it{it} total={tw:.6f} (points {float(m.sizes @ x):.6f} + cliques {zc:.6f}, {nzc} active of {ncols}) obj={val:.6f} probe_min={float(mv):.7f} viol={len(wit)} new_rows={added} rows={len(m.rows)} new_cliques={ncl} t={time.time()-t0:.0f}s")
        hist.append(dict(it=it, total=tw, points=float(m.sizes @ x), cliques=zc, active=nzc, ncl=ncols, probe=float(mv), viol=len(wit)))
        if added == 0 and mv >= 1 and ncl == 0: break
    # finalize at W = 1e7, rounded up, checked at N (the clique N); feed back if needed
    out_path = f"runs/boxclique_{tag}.txt"
    for rnd in range(30):
        tw, npts, ncls = export(m, x, cl, z, out_path, WD=10 ** 7, up=True, N=N)
        mv1, ok1 = T.run_verifier(out_path, N, topk=topk, sep=sep, threads=threads, n=n)
        wit = T.read_witnesses(sep, N) if os.path.exists(sep) else []
        if os.path.exists(sep): os.remove(sep)
        log(f"  final[{rnd}] {out_path}: points={npts} cliques={ncls} total={tw:.7f} min@{N}={float(mv1):.7f} {'OK' if ok1 else 'FAIL'} (viol {len(wit)})")
        if ok1: break
        if tw >= n and mv1 >= 1: log("  covering holds but total >= n: not a certificate"); break
        m.add_rows(wit)
        out = solve(m, cl, margin, want_dual=False, force=force)
        if out is None: log("  LP failed in finalize"); return None
        val, x, z, y = out
    res = dict(cert=out_path, s=str(s), sf=sf, points=npts, cliques=ncls, total=tw, min=str(mv1), ok=ok1, hist=hist, t=time.time() - t0, N=N, Q=Q, eps=eps)
    json.dump(res, open(f"runs/boxclique_{tag}.json", 'w'), indent=1)
    log(f"[{tag}] RESULT s={s} ({sf:.7f}) points={npts} cliques={ncls} total={tw:.7f} min@{N}={mv1} {'VERIFIED' if ok1 else 'NOT VERIFIED'} t={time.time()-t0:.0f}s")
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cert'); ap.add_argument('tag')
    ap.add_argument('--Dp', type=int, default=None); ap.add_argument('--mul', type=int, default=1)
    ap.add_argument('--N', type=int, default=2000); ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--topk', type=int, default=6); ap.add_argument('--rounds', type=int, default=40)
    ap.add_argument('--colgen', type=int, default=30); ap.add_argument('--per-round', type=int, default=8)
    ap.add_argument('--eps', type=float, default=0.005); ap.add_argument('--Q', type=int, default=10 ** 6)
    ap.add_argument('--margin', type=float, default=2e-6); ap.add_argument('--n', type=int, default=12)
    ap.add_argument('--seed', type=int, default=0); ap.add_argument('--no-warm', action='store_true')
    ap.add_argument('--force', type=float, default=0.0, help='hold the first clique-orbit column at weight >= this (format demonstration)')
    a = ap.parse_args()
    HIGHS['random_seed'] = a.seed; T.HIGHS['random_seed'] = a.seed
    os.makedirs('runs', exist_ok=True)
    lf = open(f"runs/boxclique_{a.tag}.log", 'a')
    def log(msg):
        print(msg, flush=True); lf.write(msg + '\n'); lf.flush()
    log(f"boxclique.py {' '.join(sys.argv[1:])}")
    run(a.cert, a.tag, Dp=a.Dp, mul=a.mul, N=a.N, threads=a.threads, topk=a.topk, rounds=a.rounds, colgen=a.colgen,
        per_round=a.per_round, eps=a.eps, Q=a.Q, margin=a.margin, n=a.n, log=log, warm=not a.no_warm, force=a.force)


if __name__ == '__main__':
    main()
