import numpy as np, sys
from multiprocessing import Pool
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from escape import escape, rect_pts
r2 = np.sqrt(2)
def job(ab):
    a, b = ab
    e, x = escape(rect_pts(a, b), a, b, nx=50, ny=50, nth=90)
    return a, b, e, x
if __name__ == '__main__':
    rng = np.random.default_rng(1)
    jobs = []
    for a in np.linspace(4 - r2, 2 + 1/r2, 16):
        for t in np.linspace(0, 1, 8):
            jobs.append((a, a + t*(4 + r2 - 2*a)))
    for _ in range(80):
        a = rng.uniform(4 - r2, 2 + 1/r2); b = rng.uniform(a, 4 + r2 - a); jobs.append((a, b))
    with Pool(6) as p:
        res = p.map(job, jobs)
    w = max(res, key=lambda r: r[2])
    print(f"{len(res)} containers, max esc = {w[2]:.3e} at a={w[0]:.5f} b={w[1]:.5f}, pose {w[3]}")
    print("count esc > 1e-12:", sum(r[2] > 1e-12 for r in res))
