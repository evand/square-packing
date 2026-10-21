"""Cover C as blocks: every closed unit square on the R5 line contains one full block (weight 1):
corner pairs {(c0,1),(1,1)} (x4 by mirrors), wall blocks {(c0,1),(1,b/2)}, {(c0,b-1),(1,b/2)} (x2 sides),
centre blocks {(1,1),(1,b/2),(a-1,b/2)} (x4).  Checked on all arrangement-vertex probes."""
import sys, math, numpy as np
from c10 import cover, line_b
from tightC import probes
from multiprocessing import Pool
import os
C0 = float(os.environ.get("C0", "0.9"))
def blocks(a, b):
    # indices into cover(): order r(c0,1), r(c0,b-1), r'(a-c0,1), r'(a-c0,b-1), s(1,1), s(1,b-1), s'(a-1,1), s'(a-1,b-1), m, m'
    r, rt, R, Rt, s, st, S, St, m, M = range(10)
    return [(r, s), (rt, st), (R, S), (Rt, St), (r, m), (rt, m), (R, M), (Rt, M), (s, m, M), (st, m, M), (S, m, M), (St, m, M)]
def run(args):
    a, th = args
    b = line_b(a); P, W = cover(a, b, C0); B = blocks(a, b)
    Mx, k = probes(P, a, b, th)
    c, s = np.cos(th), np.sin(th)
    dx = P[None, :, 0] - Mx[:, None, 0]; dy = P[None, :, 1] - Mx[:, None, 1]
    ins = (np.abs(dx*c + dy*s) <= .5 + 1e-12) & (np.abs(-dx*s + dy*c) <= .5 + 1e-12)
    ok = np.zeros(len(Mx), bool)
    for bl in B: ok |= ins[:, list(bl)].all(1)
    bad = np.where(~ok)[0]
    return (len(bad), a, th, Mx[bad[0]] if len(bad) else None, len(Mx))
if __name__ == '__main__':
    s, e = 1 + math.sqrt(2), 4 - math.sqrt(2)
    As = [s + (e - s) * i / 40 for i in range(41)]
    ths = np.linspace(0, np.pi / 2, 3600, endpoint=False)
    with Pool(8) as pool:
        res = pool.map(run, [(a, th) for a in As for th in ths], chunksize=200)
    nb = sum(r[0] for r in res); npb = sum(r[4] for r in res)
    print(f"probes {npb}, without a full block: {nb}")
    for r in res:
        if r[0]: print("  fail", r[1], math.degrees(r[2]), r[3]); break
