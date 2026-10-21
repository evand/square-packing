#!/usr/bin/env python3
"""Grid map of min_h(tilts, w0) - (4 + 2sqrt2 - 2 w0) over [0,90)^3 for n = 3 (tiltlp.py)."""
import itertools, json, math, os, sys
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tiltlp import min_h


def f(a):
    t, w0 = a
    return t, w0, min_h(t, w0)


if __name__ == '__main__':
    step, procs, out = float(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    w0s = [float(x) for x in sys.argv[4].split(',')]
    g = [i * step for i in range(int(round(90 / step)))]
    pts = [(t, w0) for w0 in w0s for t in itertools.product(g, g, g) if t[0] <= t[2]]   # A <-> C mirror
    with ProcessPoolExecutor(procs) as ex, open(out, 'w') as fo:
        for t, w0, h in ex.map(f, pts, chunksize=8):
            fo.write(json.dumps(dict(t=t, w0=w0, h=h, gap=h - (4 + 2 * math.sqrt(2) - 2 * w0))) + '\n')
