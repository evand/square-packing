"""Independent float check of a unit-square packing: containment (corner coords in [0,s]) and pairwise separating-axis
overlap depth (max over the 4 axes of the gap; >= 0 means disjoint interiors).  Prints worst values.
  satcheck_f64.py ID...   (10-23; negative control: a square moved past the edge onto a neighbour -> FAIL)"""
import os, sys, sqlite3, numpy as np
STORE = os.environ.get("PK_STORE", os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs/store.sqlite"))
db = sqlite3.connect(f'file:{STORE}?mode=ro', uri=True)
for idd in map(int, sys.argv[1:]):
    n, s, c = db.execute("SELECT n, s, coords FROM packing WHERE id=?", (idd,)).fetchone()
    P = np.frombuffer(c, '<f8').reshape(-1, 3); assert len(P) == n
    th = np.radians(P[:, 2]); u = np.stack([np.cos(th), np.sin(th)], 1); v = np.stack([-np.sin(th), np.cos(th)], 1)
    C = P[:, :2]
    corners = C[:, None, :] + 0.5 * (np.array([1, 1, -1, -1])[None, :, None] * u[:, None, :] + np.array([1, -1, -1, 1])[None, :, None] * v[:, None, :])
    cont = min(corners.min(), (s - corners).min())
    worst = np.inf
    for i in range(n):
        j = np.arange(i + 1, n)
        d = np.linalg.norm(C[j] - C[i], axis=1); j = j[d < 1.5]
        if not len(j): continue
        gaps = []
        for ax in (u[i], v[i]):
            pi = C[i] @ ax; pj = np.einsum('jkd,d->jk', corners[j], ax)
            gaps.append(np.maximum(pj.min(1) - (pi + 0.5), (pi - 0.5) - pj.max(1)))
        for A in (u[j], v[j]):
            pj = np.einsum('jd,jd->j', C[j], A); pi = np.einsum('kd,jd->jk', corners[i], A)
            gaps.append(np.maximum(pi.min(1) - (pj + 0.5), (pj - 0.5) - pi.max(1)))
        worst = min(worst, np.max(np.stack(gaps), 0).min())
    print(idd, n, f's {s:.10f} containment slack {cont:.2e}  worst pair separation {worst:.2e}', 'OK' if cont > -1e-9 and worst > -1e-9 else 'FAIL')
