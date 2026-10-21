"""Weighted point covers near the R5 corner: LP over point atoms + arrangement-based min-capture check.

For fixed angle th, point p is in the closed square centred m iff m lies in the closed tilted unit square
Q_p(th) around p.  capture(m) = sum w_p [m in Q_p]; its minimum over the centre box [k, a-k] x [k, b-k] is
attained in an open cell of the arrangement of the Q_p's and the box, so we test points at small offsets
from every arrangement vertex (edge-edge, edge-box, box corners) in 16 directions.  Angles on a fine grid.
"""
import sys, json, numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from multiprocessing import Pool

EPS = 1e-7

def caps(P, w, M, th):
    c, s = np.cos(th), np.sin(th)
    dx = P[None, :, 0] - M[:, None, 0]; dy = P[None, :, 1] - M[:, None, 1]
    u = dx * c + dy * s; v = -dx * s + dy * c
    inside = (np.abs(u) <= 0.5 + 1e-12) & (np.abs(v) <= 0.5 + 1e-12)
    return inside.astype(float) @ w

def lines_at(P, th, a, b):
    """edge lines of every Q_p as (n, d): n.m = d, plus the centre-box lines."""
    c, s = np.cos(th), np.sin(th); k = (abs(c) + abs(s)) / 2
    L = []
    for n in (np.array([c, s]), np.array([-s, c])):
        for p in P:
            L.append((n, n @ p + 0.5)); L.append((n, n @ p - 0.5))
    for n, d in ((np.array([1., 0]), k), (np.array([1., 0]), a - k), (np.array([0, 1.]), k), (np.array([0, 1.]), b - k)):
        L.append((n, d))
    return L, k

def min_capture_angle(P, w, a, b, th):
    L, k = lines_at(P, th, a, b)
    if 2 * k > min(a, b): return np.inf, None
    N = np.array([l[0] for l in L]); D = np.array([l[1] for l in L])
    # pairwise intersections
    i, j = np.triu_indices(len(L), 1)
    det = N[i, 0] * N[j, 1] - N[i, 1] * N[j, 0]
    ok = np.abs(det) > 1e-12
    i, j, det = i[ok], j[ok], det[ok]
    x = (D[i] * N[j, 1] - D[j] * N[i, 1]) / det
    y = (N[i, 0] * D[j] - N[j, 0] * D[i]) / det
    V = np.stack([x, y], 1)
    keep = (V[:, 0] >= k - 1e-9) & (V[:, 0] <= a - k + 1e-9) & (V[:, 1] >= k - 1e-9) & (V[:, 1] <= b - k + 1e-9)
    V = V[keep]
    ang = np.arange(16) * np.pi / 8 + 0.1
    off = EPS * np.stack([np.cos(ang), np.sin(ang)], 1)
    M = (V[:, None, :] + off[None, :, :]).reshape(-1, 2)
    M = M[(M[:, 0] >= k) & (M[:, 0] <= a - k) & (M[:, 1] >= k) & (M[:, 1] <= b - k)]
    if len(M) == 0: return np.inf, None
    cv = caps(P, w, M, th)
    t = np.argmin(cv)
    return cv[t], (M[t, 0], M[t, 1], th)

def min_capture(P, w, a, b, nth=1800, pool=None):
    ths = np.linspace(0, np.pi / 2, nth, endpoint=False)
    args = [(P, w, a, b, th) for th in ths]
    res = pool.starmap(min_capture_angle, args) if pool else [min_capture_angle(*x) for x in args]
    return res

def rows_of(atoms, poses):
    out = np.zeros((len(poses), len(atoms)))
    for r, (x, y, th) in enumerate(poses):
        c, s = np.cos(th), np.sin(th)
        dx = atoms[:, 0] - x; dy = atoms[:, 1] - y
        out[r] = (np.abs(dx * c + dy * s) <= 0.5 + 1e-12) & (np.abs(-dx * s + dy * c) <= 0.5 + 1e-12)
    return out

def orbits(atoms, a, b):
    key = {}
    orb = []
    for x, y in atoms:
        kk = (round(min(x, a - x), 9), round(min(y, b - y), 9))
        orb.append(key.setdefault(kk, len(key)))
    return np.array(orb), len(key)

def run(a, b, row_pitch=0.02, grid_pitch=0.1, rounds=8, procs=8, nth=900):
    xs = np.unique(np.concatenate([np.arange(0, a + 1e-9, row_pitch), [1, a - 1], a - np.arange(0, a + 1e-9, row_pitch)]))
    atoms = [(x, y) for x in xs for y in (1.0, b - 1)]
    gx = np.linspace(0, a, round(a / grid_pitch) + 1); gy = np.linspace(0, b, round(b / grid_pitch) + 1)
    atoms += [(x, y) for x in gx for y in gy]
    atoms = np.unique(np.round(np.array(atoms), 10), axis=0)
    orb, n = orbits(atoms, a, b)
    Pm = sp.csr_matrix((np.ones(len(orb)), (np.arange(len(orb)), orb)), shape=(len(orb), n))
    mult = np.bincount(orb, minlength=n).astype(float)
    # initial poses: coarse grid in a quadrant (D2 symmetric weights)
    poses = []
    for th in np.radians(np.arange(0, 90, 2.0)):
        k = (np.cos(th) + np.sin(th)) / 2
        for x in np.linspace(k, a / 2, 14):
            for y in np.linspace(k, b / 2, 16):
                poses.append((x, y, th))
    G = sp.csr_matrix(rows_of(atoms, poses)) @ Pm
    hist = []
    with Pool(procs) as pool:
        for rd in range(rounds + 1):
            r = linprog(mult, A_ub=-G, b_ub=-np.ones(G.shape[0]), bounds=(0, None), method='highs')
            w = r.x[orb]
            sel = w > 1e-9
            res = min_capture(atoms[sel], w[sel], a, b, nth=nth, pool=pool)
            vals = np.array([v for v, _ in res])
            mc = vals.min()
            hist.append(dict(round=rd, lp=float(r.fun), min_capture=float(mc), rows=int(G.shape[0])))
            print(json.dumps(hist[-1]), flush=True)
            if mc >= 1 - 1e-9 or rd == rounds: break
            bad = [p for v, p in res if v < 1 - 1e-9]
            # fold to the quadrant (symmetry) -- weights are symmetric so any pose is fine as a row
            G = sp.vstack([G, sp.csr_matrix(rows_of(atoms, bad)) @ Pm]).tocsr()
    return atoms, w, hist

if __name__ == '__main__':
    a, b = float(sys.argv[1]), float(sys.argv[2])
    atoms, w, hist = run(a, b)
    sel = w > 1e-9
    print(f"total {w.sum():.6f}  ({sel.sum()} atoms)")
    for (x, y), ww in sorted(zip(atoms[sel], w[sel]), key=lambda t: (t[0][1], t[0][0])):
        print(f"  ({x:.4f}, {y:.4f})  w {ww:.4f}")
    np.savez(f"/home/evand/math/square-packing/public/runs/rectfront/r5pts/ptcover_{a:.4f}_{b:.4f}.npz",
             atoms=atoms, w=w, hist=json.dumps(hist))
