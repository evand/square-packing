"""Mixed cover LP for R5 near the corner: points + uniform segment densities, D2-symmetric.

capture(square) = sum_points w_p [p in closed square] + sum_segs rho_s * |s ∩ square|.
min total weight s.t. capture >= 1 on sampled poses; cutting planes from an escape optimiser.
Usage: python3 mixlp.py a b [pitch] [rounds]
"""
import sys, json, numpy as np, scipy.sparse as sp
from scipy.optimize import linprog, minimize
from multiprocessing import Pool

TOL = 1e-9

def build_atoms(a, b, pitch):
    N = max(2, round(a / pitch)); M = max(2, round(b / pitch))
    xs = np.linspace(0, a, N + 1); ys = np.linspace(0, b, M + 1)
    atoms = []  # (kind, x0, y0, x1, y1)
    for x in xs:
        for y in ys: atoms.append((0, x, y, x, y))
    for y in list(ys) + [1.0, b - 1]:
        for i in range(N): atoms.append((1, xs[i], y, xs[i + 1], y))
    for x in list(xs) + [1.0, a - 1]:
        for j in range(M): atoms.append((1, x, ys[j], x, ys[j + 1]))
    for p in [(1, 1), (a - 1, 1), (1, b - 1), (a - 1, b - 1)]:
        atoms.append((0, p[0], p[1], p[0], p[1]))
    A = np.array(atoms, float)
    # orbits under x -> a-x, y -> b-y
    def key(r):
        k, x0, y0, x1, y1 = r
        cands = []
        for mx in (0, 1):
            for my in (0, 1):
                X0, X1 = (a - x0, a - x1) if mx else (x0, x1)
                Y0, Y1 = (b - y0, b - y1) if my else (y0, y1)
                e = tuple(np.round(sorted([(X0, Y0), (X1, Y1)]), 9).ravel())
                cands.append((k,) + e)
        return min(cands)
    keys = [key(r) for r in A]
    uniq = {}
    orb = np.array([uniq.setdefault(k, len(uniq)) for k in keys])
    return A, orb, len(uniq)

def capture_rows(A, poses):
    """coefficient matrix rows (poses x atoms): 1 for contained points, length for segments."""
    K = A[:, 0]; x0, y0, x1, y1 = A[:, 1], A[:, 2], A[:, 3], A[:, 4]
    L = np.hypot(x1 - x0, y1 - y0)
    out = []
    for (cx, cy, th) in poses:
        c, s = np.cos(th), np.sin(th)
        u0 = (x0 - cx) * c + (y0 - cy) * s; v0 = -(x0 - cx) * s + (y0 - cy) * c
        du = (x1 - x0) * c + (y1 - y0) * s; dv = -(x1 - x0) * s + (y1 - y0) * c
        lo = np.zeros_like(u0); hi = np.ones_like(u0)
        for w0, dw in ((u0, du), (v0, dv)):
            with np.errstate(divide='ignore', invalid='ignore'):
                ta = (-0.5 - TOL - w0) / dw; tb = (0.5 + TOL - w0) / dw
            t1 = np.minimum(ta, tb); t2 = np.maximum(ta, tb)
            par = np.abs(dw) < 1e-15
            inside = np.abs(w0) <= 0.5 + TOL
            t1 = np.where(par, np.where(inside, -np.inf, np.inf), t1)
            t2 = np.where(par, np.where(inside, np.inf, -np.inf), t2)
            lo = np.maximum(lo, t1); hi = np.minimum(hi, t2)
        seg = np.maximum(0, hi - lo) * L
        pt = ((np.abs(u0) <= 0.5 + TOL) & (np.abs(v0) <= 0.5 + TOL)).astype(float)
        out.append(np.where(K == 0, pt, seg))
    return np.array(out)

def sample_poses(a, b, cpitch, dth):
    poses = []
    for th in np.radians(np.arange(0, 90, dth)):
        k = (np.cos(th) + np.sin(th)) / 2
        if 2 * k > min(a, b): continue
        nx = max(2, int(np.ceil((a / 2 - k) / cpitch)) + 1); ny = max(2, int(np.ceil((b / 2 - k) / cpitch)) + 1)
        for x in np.linspace(k, a / 2, nx):
            for y in np.linspace(k, b / 2, ny):
                poses.append((x, y, th))
    return poses

def orbit_matrix(R, orb, n):
    P = sp.csr_matrix((np.ones(len(orb)), (np.arange(len(orb)), orb)), shape=(len(orb), n))
    return sp.csr_matrix(R) @ P

def solve_lp(G, mult):
    n = G.shape[1]
    r = linprog(mult, A_ub=-G, b_ub=-np.ones(G.shape[0]), bounds=(0, None), method='highs')
    return r.x, r.fun

def esc_obj(p, a, b, A, w):
    x, y, th = p
    k = (abs(np.cos(th)) + abs(np.sin(th))) / 2
    X = np.clip(x, k, a - k); Y = np.clip(y, k, b - k)
    pen = 10 * (abs(x - X) + abs(y - Y))
    return capture_rows(A, [(X, Y, th)])[0] @ w + pen

def worst_poses(args):
    a, b, A, w, seeds = args
    out = []
    for s in seeds:
        r = minimize(esc_obj, s, args=(a, b, A, w), method='Nelder-Mead',
                     options=dict(xatol=1e-7, fatol=1e-9, maxiter=600))
        x, y, th = r.x; k = (abs(np.cos(th)) + abs(np.sin(th))) / 2
        out.append((r.fun, (float(np.clip(x, k, a - k)), float(np.clip(y, k, b - k)), float(th))))
    return out

def main(a, b, pitch=0.1, rounds=4, procs=8):
    A, orb, n = build_atoms(a, b, pitch)
    mult = np.bincount(orb, minlength=n).astype(float)  # weight per orbit var counts each member
    poses = sample_poses(a, b, 0.05, 2.0)
    with Pool(procs) as pool:
        chunks = np.array_split(np.arange(len(poses)), procs * 4)
        Rs = pool.starmap(capture_rows, [(A, [poses[i] for i in ch]) for ch in chunks])
        G = sp.vstack([orbit_matrix(R, orb, n) for R in Rs]).tocsr()
        hist = []
        for rd in range(rounds + 1):
            z, val = solve_lp(G, mult)
            w = z[orb]
            # escape: full-plane check via fine sample + Nelder–Mead polish
            fine = sample_poses(a, b, 0.02, 0.5)
            fine_c = np.concatenate(pool.starmap(capture_rows, [(A, [fine[i] for i in ch]) for ch in
                                                                np.array_split(np.arange(len(fine)), procs * 8)]))
            capf = fine_c @ w
            idx = np.argsort(capf)[:400]
            seeds = [fine[i] for i in idx]
            res = pool.map(worst_poses, [(a, b, A, w, seeds[i::procs]) for i in range(procs)])
            res = sorted([r for rr in res for r in rr])
            minc = min(res[0][0], capf.min())
            hist.append(dict(round=rd, lp=float(val), min_capture=float(minc), honest=float(val / minc) if minc > 0 else None,
                             rows=G.shape[0]))
            print(json.dumps(hist[-1]), flush=True)
            if rd == rounds: break
            new = [p for _, p in res[:200]] + [fine[i] for i in idx[:200]]
            Rn = capture_rows(A, new)
            G = sp.vstack([G, orbit_matrix(Rn, orb, n)]).tocsr()
    return A, w, hist

if __name__ == '__main__':
    a, b = float(sys.argv[1]), float(sys.argv[2])
    pitch = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
    rounds = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    A, w, hist = main(a, b, pitch, rounds)
    import os
    out = f"/home/evand/math/square-packing/public/runs/rectfront/r5pts/mixlp_{a:.4f}_{b:.4f}_{pitch}.npz"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out, A=A, w=w, hist=json.dumps(hist))
    # structure summary
    K = A[:, 0]; L = np.hypot(A[:, 3] - A[:, 1], A[:, 4] - A[:, 2])
    mass = np.where(K == 0, w, w * L)
    tot = mass.sum()
    print(f"total {tot:.4f}; points {mass[K==0].sum():.4f}, segments {mass[K==1].sum():.4f}")
    order = np.argsort(-mass)
    for i in order[:40]:
        if mass[i] < 1e-4: break
        k, x0, y0, x1, y1 = A[i]
        print(f"  {'pt ' if k==0 else 'seg'} ({x0:.3f},{y0:.3f})-({x1:.3f},{y1:.3f}) mass {mass[i]:.4f} dens {w[i]:.4f}")
