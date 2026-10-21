"""Scan cover-A-type families (weights be, 1-be on rows; ga interior) for validity along the R5 line."""
import sys, math, itertools, numpy as np, ptcover as pc
from multiprocessing import Pool
from cand3 import cover, line_b
AS = [1 + math.sqrt(2), 2.45, 2.48, 2.50, 2.52, 2.55, 4 - math.sqrt(2)]
def one(args):
    prm, a, nth = args
    P, W = cover(a, line_b(a), *prm)
    return min(pc.min_capture_angle(P, W, a, line_b(a), th)[0] for th in np.linspace(0, np.pi/2, nth, endpoint=False))
if __name__ == '__main__':
    be, ga = eval(sys.argv[1]), eval(sys.argv[2])
    c0s = [float(x) for x in sys.argv[3].split(',')]
    xis = [float(x) for x in sys.argv[4].split(',')]
    yis = [float(x) for x in sys.argv[5].split(',')]
    nth = int(sys.argv[6]) if len(sys.argv) > 6 else 360
    jobs = [((c0, be, xi, yi, ga), a, nth) for c0 in c0s for xi in xis for yi in yis for a in AS]
    with Pool(8) as pool:
        res = pool.map(one, jobs, chunksize=1)
    for i in range(0, len(jobs), len(AS)):
        r = res[i:i+len(AS)]
        print(jobs[i][0][:4], ' '.join(f"{v:.3f}" for v in r), 'OK' if min(r) > 1-1e-9 else '')
