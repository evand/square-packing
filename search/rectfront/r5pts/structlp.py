"""Structured weighted point cover LP on the R5 boundary: D2-symmetric atoms on rows y = 1, b-1 and a few
interior columns; cutting planes from the arrangement check (ptcover.min_capture)."""
import sys, json, numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from multiprocessing import Pool
import ptcover as pc

def atoms_for(a, b):
    xr = np.round(np.concatenate([np.arange(0.40, 1.0001, 0.02), [1.0]]), 6)
    A = [(x, y) for x in xr for y in (1.0, b - 1)]
    xi = [0.8, 0.9, 1.0, 1.1, a / 2]
    yi = np.round(np.arange(1.1, b / 2 + 1e-9, 0.05), 6)
    A += [(x, y) for x in xi for y in list(yi) + [b - t for t in yi]]
    A = np.array(A)
    full = np.concatenate([A, np.c_[a - A[:, 0], A[:, 1]]])
    full = np.unique(np.round(full, 9), axis=0)
    return full

def run(a, b, rounds=12, procs=4):
    atoms = atoms_for(a, b)
    orb, n = pc.orbits(atoms, a, b)
    Pm = sp.csr_matrix((np.ones(len(orb)), (np.arange(len(orb)), orb)), shape=(len(orb), n))
    mult = np.bincount(orb, minlength=n).astype(float)
    poses = []
    for th in np.radians(np.arange(0, 90, 3.0)):
        k = (np.cos(th) + np.sin(th)) / 2
        for x in np.linspace(k, a / 2, 12):
            for y in np.linspace(k, b / 2, 14):
                poses.append((x, y, th))
    G = sp.csr_matrix(pc.rows_of(atoms, poses)) @ Pm
    with Pool(procs) as pool:
        for rd in range(rounds + 1):
            r = linprog(mult, A_ub=-G, b_ub=-np.ones(G.shape[0]), bounds=(0, None), method='highs')
            w = r.x[orb]; sel = w > 1e-9
            res = pc.min_capture(atoms[sel], w[sel], a, b, nth=900, pool=pool)
            mc = min(v for v, _ in res)
            print(json.dumps(dict(a=a, b=b, round=rd, lp=round(r.fun, 6), min_capture=round(mc, 6))), flush=True)
            if mc >= 1 - 1e-9: break
            bad = [p for v, p in res if v < 1 - 1e-9]
            G = sp.vstack([G, sp.csr_matrix(pc.rows_of(atoms, bad)) @ Pm]).tocsr()
        res = pc.min_capture(atoms[sel], w[sel], a, b, nth=3600, pool=pool)
        mc = min(v for v, _ in res)
    return atoms[sel], w[sel], r.fun, mc

if __name__ == '__main__':
    for s in sys.argv[1:]:
        a, b = map(float, s.split(','))
        P, W, val, mc = run(a, b)
        print(f"RESULT a={a} b={b} total={W.sum():.6f} final_min_capture(3600 angles)={mc:.6f}")
        for (x, y), w in sorted(zip(P, W), key=lambda t: (t[0][1], t[0][0])):
            if x <= a / 2 + 1e-9 and y <= b / 2 + 1e-9: print(f"   ({x:.4f},{y:.4f}) w {w:.4f}")
        np.savez(f"/home/evand/math/square-packing/public/runs/rectfront/r5pts/structlp_{a:.4f}_{b:.4f}.npz", P=P, W=W)
