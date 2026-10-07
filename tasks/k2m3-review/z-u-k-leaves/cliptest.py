import sys, random
sys.path.insert(0, '/home/evand/math/square-packing/public/s12/search')
import zeromargin as zm
from fractions import Fraction as F
m = F(7); random.seed(5); bad = 0; n = 0; axis = 0; empty = 0
def best_centre(lo, hi):  # the centre coordinate in [lo,hi] closest to m/2 (most room)
    return min(max(m / 2, lo), hi)
def admissible_exists(box, u):
    c, s = zm.trig(u); w = c + s
    cx = best_centre(box[0], box[1]); cy = best_centre(box[2], box[3])
    return w / 2 <= cx <= m - w / 2 and w / 2 <= cy <= m - w / 2
for t in range(200000):
    k = random.randint(0, 14); wd = F(1, 10 * 2 ** k)
    def co():
        return random.choice([F(random.randint(0, 20), 40), F(1, 2) + F(random.randint(-200, 200), 2 ** random.randint(8, 20)), F(random.randint(0, 35), 10)])
    x0 = max(F(0), co()); y0 = max(F(0), co())
    x0 = min(x0, F(7, 2) - wd); y0 = min(y0, F(7, 2) - wd)
    ku = random.randint(3, 20); du = F(1, 2 ** ku)
    u0 = du * random.randint(0, int(F(1, 2) / du) - 1); u1 = u0 + du
    box = (x0, x0 + wd, y0, y0 + wd, u0, u1)
    cu1 = zm.clip_bin(box, m); nu1 = min(cu1, u1)
    B = zm.bin_data(u0, nu1)
    lo = B['wlo'] / 2
    is_empty = box[1] < lo or box[0] > m - lo or box[3] < lo or box[2] > m - lo
    # sample u in (nu1, u1] and, if EMPTY, in [u0, nu1]
    tests = [nu1 + (u1 - nu1) * F(j, 7) for j in range(1, 8)] if nu1 < u1 else []
    tests += [nu1 + (u1 - nu1) / 10 ** 9] if nu1 < u1 else []
    if is_empty:
        empty += 1; tests += [u0 + (nu1 - u0) * F(j, 7) for j in range(0, 8)]
    if nu1 == 0 and not is_empty: axis += 1; tests += [F(1, 10 ** 12)] if u1 > 0 else []
    for u in tests:
        n += 1
        if admissible_exists(box, u):
            bad += 1; print('BAD', box, cu1, u, is_empty)
print('tests', n, 'bad', bad, 'empty boxes', empty, 'axis boxes', axis)
