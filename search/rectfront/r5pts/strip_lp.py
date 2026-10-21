"""Weighted point-cover LP for R5's strip part (short side H < 2, long side L on the boundary curve).

Container [0, L] x [0, H] (long side horizontal).  Atoms: rows y in {H/2, 1, H-1, 1/2, H-1/2} at pitch `rp`
plus a grid at pitch `gp`; D2-symmetric weights (orbits under x -> L-x, y -> H-y).  Cutting planes from
ptcover.min_capture (arrangement-vertex check).  Usage: python3 strip_lp.py H [L] ...  (L default: curve).
"""
import sys, json, numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from multiprocessing import Pool
import ptcover as pc

R2 = np.sqrt(2); H0 = (6 * R2 - 1) / 4
import os
CONTACT = os.environ.get('CONTACT', '0') == '1'

def Lb(H):
    return 5.0 if H <= H0 else 4 + 6 * R2 - 4 * H

def atoms_for(L, H, rp=0.02, gp=0.1):
    xs = np.unique(np.round(np.concatenate([np.arange(0, L + 1e-9, rp), L - np.arange(0, L + 1e-9, rp),
                                            np.arange(0, 6)]), 9))
    xs = xs[(xs >= 0) & (xs <= L)]
    rows = sorted(set(np.round([H / 2, 1, H - 1, 0.5, H - 0.5], 9)))
    A = [(x, y) for x in xs for y in rows if 0 <= y <= H]
    gx = np.linspace(0, L, round(L / gp) + 1); gy = np.linspace(0, H, max(2, round(H / gp) + 1))
    A += [(x, y) for x in gx for y in gy]
    if CONTACT and H > 1.5:
        # zigzag contact lines: a-b along x + y = 2, b-c along y - x = H/2 - x*, x* = 2 + 2 sqrt2 - 1.5 H
        xs_ = 2 + 2 * R2 - 1.5 * H
        A += [(1, 1), (xs_, H / 2)]
        for t in np.arange(-0.6, 0.6001, 0.02):
            A += [(1 + t, 1 - t), (xs_ + t, H / 2 + t)]
    A = [(x, y) for x, y in A if 0 <= x <= L and 0 <= y <= H]
    A = A + [(L - x, y) for x, y in A] + [(x, H - y) for x, y in A] + [(L - x, H - y) for x, y in A]
    return np.unique(np.round(np.array(A), 9), axis=0)

def run(L, H, rounds=14, procs=int(__import__('os').environ.get('STRIP_PROCS', 3)), nth=900, rp=0.02, gp=0.1):
    atoms = atoms_for(L, H, rp, gp)
    orb, n = pc.orbits(atoms, L, H)
    Pm = sp.csr_matrix((np.ones(len(orb)), (np.arange(len(orb)), orb)), shape=(len(orb), n))
    mult = np.bincount(orb, minlength=n).astype(float)
    poses = []
    for th in np.radians(np.arange(0, 90, 2.0)):
        k = (np.cos(th) + np.sin(th)) / 2
        if 2 * k > H: continue
        for x in np.linspace(k, L / 2, 30):
            for y in np.linspace(k, H / 2, 6):
                poses.append((x, y, th))
    G = sp.csr_matrix(pc.rows_of(atoms, poses)) @ Pm
    with Pool(procs) as pool:
        for rd in range(rounds + 1):
            r = linprog(mult, A_ub=-G, b_ub=-np.ones(G.shape[0]), bounds=(0, None), method='highs')
            w = r.x[orb]; sel = w > 1e-9
            res = pc.min_capture(atoms[sel], w[sel], L, H, nth=nth, pool=pool)
            mc = min(v for v, _ in res)
            print(json.dumps(dict(L=round(L, 6), H=H, round=rd, lp=round(r.fun, 6), min_capture=round(mc, 6),
                                  rows=G.shape[0])), flush=True)
            if mc >= 1 - 1e-9 or rd == rounds: break
            bad = [p for v, p in res if v < 1 - 1e-9]
            G = sp.vstack([G, sp.csr_matrix(pc.rows_of(atoms, bad)) @ Pm]).tocsr()
    return atoms[sel], w[sel], r.fun, mc

if __name__ == '__main__':
    for s in sys.argv[1:]:
        if ',' in s: H, L = map(float, s.split(','))
        else: H = float(s); L = Lb(H)
        P, W, val, mc = run(L, H)
        print(f"RESULT H={H} L={L:.6f} total={W.sum():.6f} min_capture={mc:.6f}")
        for (x, y), w in sorted(zip(P, W), key=lambda t: (t[0][0], t[0][1])):
            if x <= L / 2 + 1e-9 and y <= H / 2 + 1e-9: print(f"   ({x:.4f},{y:.4f}) w {w:.4f}")
        np.savez(f"/home/evand/math/square-packing/public/runs/rectfront/r5pts/strip_lp_{H:.4f}_{L:.4f}.npz", P=P, W=W)
