"""Dual check for R5 strip covers: a fractional packing (pose weights y >= 0, sum Y) whose depth
max_p sum_q y_q [p in closed q] is D shows that every point cover has total >= Y / D.
Pose set: LP duals from a restricted cover LP (strip_lp atoms + cutting planes), then the depth is measured
over ALL points via the arrangement of the dual squares (vertices + offsets + a fine grid).
Usage: OMP_NUM_THREADS=1 python3 strip_dual.py H [L]"""
import os, sys, json, numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from multiprocessing import Pool
import ptcover as pc
from strip_lp import Lb, atoms_for

def corners(x, y, th):
    c, s = np.cos(th), np.sin(th)
    e1, e2 = np.array([c, s]) / 2, np.array([-s, c]) / 2
    m = np.array([x, y])
    return np.array([m + e1 + e2, m - e1 + e2, m - e1 - e2, m + e1 - e2])

def depth(poses, y, extra=None):
    """max over points p of sum y_q [p in q]; candidates: square vertices, pairwise edge intersections, offsets."""
    segs = []
    for (x, yy, th) in poses:
        C = corners(x, yy, th)
        for i in range(4): segs.append((C[i], C[(i + 1) % 4]))
    cand = [c for (x, yy, th) in poses for c in corners(x, yy, th)]
    S = np.array([[a[0], a[1], b[0], b[1]] for a, b in segs])
    n = len(S)
    for i in range(n):
        p, r = S[i, :2], S[i, 2:] - S[i, :2]
        q, sv = S[i + 1:, :2], S[i + 1:, 2:] - S[i + 1:, :2]
        den = r[0] * sv[:, 1] - r[1] * sv[:, 0]
        ok = np.abs(den) > 1e-14
        qp = q - p
        t = (qp[:, 0] * sv[:, 1] - qp[:, 1] * sv[:, 0]) / np.where(ok, den, 1)
        u = (qp[:, 0] * r[1] - qp[:, 1] * r[0]) / np.where(ok, den, 1)
        good = ok & (t >= -1e-12) & (t <= 1 + 1e-12) & (u >= -1e-12) & (u <= 1 + 1e-12)
        cand += list(p + t[good, None] * r)
    cand = np.array(cand)
    if extra is not None: cand = np.vstack([cand, extra])
    P = np.array(poses)
    best = 0
    for chunk in np.array_split(cand, max(1, len(cand) // 20000)):
        dx = chunk[:, None, 0] - P[None, :, 0]; dy = chunk[:, None, 1] - P[None, :, 1]
        c, s = np.cos(P[:, 2]), np.sin(P[:, 2])
        u = dx * c + dy * s; v = -dx * s + dy * c
        ins = (np.abs(u) <= 0.5 + 1e-9) & (np.abs(v) <= 0.5 + 1e-9)
        best = max(best, (ins.astype(float) @ y).max())
    return best

def run(L, H, rounds=14, procs=4):
    atoms = atoms_for(L, H)
    orb, n = pc.orbits(atoms, L, H)
    Pm = sp.csr_matrix((np.ones(len(orb)), (np.arange(len(orb)), orb)), shape=(len(orb), n))
    mult = np.bincount(orb, minlength=n).astype(float)
    poses = []
    for th in np.radians(np.arange(0, 90, 2.0)):
        k = (np.cos(th) + np.sin(th)) / 2
        if 2 * k > H: continue
        for x in np.linspace(k, L - k, 50):
            for yy in np.linspace(k, H - k, 9):
                poses.append((x, yy, th))
    # full (unsymmetrised) atoms for the dual: use orbit LP but keep pose rows unfolded
    G = sp.csr_matrix(pc.rows_of(atoms, poses)) @ Pm
    with Pool(procs) as pool:
        for rd in range(rounds + 1):
            r = linprog(mult, A_ub=-G, b_ub=-np.ones(G.shape[0]), bounds=(0, None), method='highs')
            w = r.x[orb]; sel = w > 1e-9
            res = pc.min_capture(atoms[sel], w[sel], L, H, nth=900, pool=pool)
            mc = min(v for v, _ in res)
            print(json.dumps(dict(L=round(L, 6), H=H, round=rd, lp=round(r.fun, 6), min_capture=round(mc, 6))),
                  flush=True)
            if mc >= 1 - 1e-9 or rd == rounds: break
            bad = [p for v, p in res if v < 1 - 1e-9]
            poses += bad
            G = sp.vstack([G, sp.csr_matrix(pc.rows_of(atoms, bad)) @ Pm]).tocsr()
    ydual = -r.ineqlin.marginals
    keep = ydual > 1e-9
    Pk = [poses[i] for i in np.where(keep)[0]]; yk = ydual[keep]
    # the orbit LP's dual weights live on the folded atom space; the true fractional packing is the
    # symmetrised family (each pose and its three mirrors, weight y/4 each) -- its depth is what counts
    sym = []
    ys = []
    for (x, yy, th), yv in zip(Pk, yk):
        for X, Y, T in ((x, yy, th), (L - x, yy, -th), (x, H - yy, -th), (L - x, H - yy, th)):
            sym.append((X, Y, T)); ys.append(yv / 4)
    ys = np.array(ys)
    D = depth(sym, ys)
    return r.fun, ys.sum(), D, sym, ys

if __name__ == '__main__' and sys.argv[1] != 'cg':
    H = float(sys.argv[1]); L = float(sys.argv[2]) if len(sys.argv) > 2 else Lb(H)
    val, Y, D, sym, ys = run(L, H)
    print(f"DUAL H={H} L={L:.6f} lp={val:.6f} sum_y={Y:.6f} depth={D:.6f} bound={Y / D:.6f} poses={len(sym)}")
    order = np.argsort(-ys)
    for i in order[:24]:
        x, yy, th = sym[i]
        print(f"   y={ys[i]:.4f}  ({x:.4f},{yy:.4f}) th={np.degrees(th) % 90:.2f}")
    np.savez(f"/home/evand/math/square-packing/public/runs/rectfront/r5pts/strip_dual_{H:.4f}_{L:.4f}.npz",
             poses=np.array(sym), y=ys)


def deep_points(poses, y, thr=1 + 1e-6, top=400):
    """points of maximal depth (candidates as in depth()), returned with their depth."""
    segs = []
    for (x, yy, th) in poses:
        C = corners(x, yy, th)
        for i in range(4): segs.append((C[i], C[(i + 1) % 4]))
    cand = [c for (x, yy, th) in poses for c in corners(x, yy, th)]
    S = np.array([[a[0], a[1], b[0], b[1]] for a, b in segs])
    for i in range(len(S)):
        p, r = S[i, :2], S[i, 2:] - S[i, :2]
        q, sv = S[i + 1:, :2], S[i + 1:, 2:] - S[i + 1:, :2]
        den = r[0] * sv[:, 1] - r[1] * sv[:, 0]
        ok = np.abs(den) > 1e-14
        qp = q - p
        t = (qp[:, 0] * sv[:, 1] - qp[:, 1] * sv[:, 0]) / np.where(ok, den, 1)
        u = (qp[:, 0] * r[1] - qp[:, 1] * r[0]) / np.where(ok, den, 1)
        good = ok & (t >= -1e-12) & (t <= 1 + 1e-12) & (u >= -1e-12) & (u <= 1 + 1e-12)
        cand += list(p + t[good, None] * r)
    cand = np.array(cand)
    # also centroids of small cells: average of candidate with nearby offsets
    P = np.array(poses)
    dx = cand[:, None, 0] - P[None, :, 0]; dy = cand[:, None, 1] - P[None, :, 1]
    c, s = np.cos(P[:, 2]), np.sin(P[:, 2])
    ins = (np.abs(dx * c + dy * s) <= 0.5 + 1e-9) & (np.abs(-dx * s + dy * c) <= 0.5 + 1e-9)
    d = ins.astype(float) @ y
    o = np.argsort(-d)
    o = o[d[o] > thr][:top]
    return cand[o], d[o]


def run2(L, H, iters=12, procs=4, nth=900):
    atoms = atoms_for(L, H)
    poses = []
    for th in np.radians(np.arange(0, 90, 2.0)):
        k = (np.cos(th) + np.sin(th)) / 2
        if 2 * k > H: continue
        for x in np.linspace(k, L / 2, 30):
            for yy in np.linspace(k, H / 2, 6):
                poses.append((x, yy, th))
    hist = []
    with Pool(procs) as pool:
        for it in range(iters):
            for inner in range(10):     # primal cutting planes
                orb, n = pc.orbits(atoms, L, H)
                Pm = sp.csr_matrix((np.ones(len(orb)), (np.arange(len(orb)), orb)), shape=(len(orb), n))
                mult = np.bincount(orb, minlength=n).astype(float)
                G = sp.csr_matrix(pc.rows_of(atoms, poses)) @ Pm
                r = linprog(mult, A_ub=-G, b_ub=-np.ones(G.shape[0]), bounds=(0, None), method='highs')
                w = r.x[orb]; sel = w > 1e-9
                res = pc.min_capture(atoms[sel], w[sel], L, H, nth=nth, pool=pool)
                mc = min(v for v, _ in res)
                if mc >= 1 - 1e-9: break
                poses += [p for v, p in res if v < 1 - 1e-9]
            ydual = -r.ineqlin.marginals
            keep = ydual > 1e-9
            sym, ys = [], []
            for (x, yy, th), yv in zip([poses[i] for i in np.where(keep)[0]], ydual[keep]):
                for X, Y, T in ((x, yy, th), (L - x, yy, -th), (x, H - yy, -th), (L - x, H - yy, th)):
                    sym.append((X, Y, T)); ys.append(yv / 4)
            ys = np.array(ys)
            pts, dep = deep_points(sym, ys)
            D = dep.max() if len(dep) else depth(sym, ys)
            hist.append(dict(it=it, lp=round(r.fun, 6), min_capture=round(mc, 6), depth=round(float(D), 6),
                             bound=round(ys.sum() / max(D, 1.0), 6), atoms=len(atoms), poses=len(poses)))
            print(json.dumps(hist[-1]), flush=True)
            if len(pts) == 0: break
            # fold deep points to the quadrant and add them (with tiny jitter-free copies)
            q = np.c_[np.minimum(pts[:, 0], L - pts[:, 0]), np.minimum(pts[:, 1], H - pts[:, 1])]
            atoms = np.unique(np.round(np.vstack([atoms, q, np.c_[L - q[:, 0], q[:, 1]],
                                                  np.c_[q[:, 0], H - q[:, 1]], np.c_[L - q[:, 0], H - q[:, 1]]]), 9),
                              axis=0)
    return r.fun, mc, D, ys.sum(), sym, ys, atoms[sel] if False else (atoms, w)

if __name__ == '__main__' and sys.argv[1] == 'cg':
    for hs in sys.argv[2:]:
        H = float(hs); L = Lb(H)
        val, mc, D, Y, sym, ys, (atoms, w) = run2(L, H)
        print(f"CG H={H} L={L:.6f} lp={val:.6f} min_capture={mc:.6f} dual_sum={Y:.6f} depth={D:.6f} "
              f"cover_cost_lower_bound={Y / max(D, 1):.6f}", flush=True)
        np.savez(f"/home/evand/math/square-packing/public/runs/rectfront/r5pts/strip_cg_{H:.4f}.npz",
                 poses=np.array(sym), y=ys, atoms=atoms, w=w)
