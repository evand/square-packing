"""Own re-proof of every LEB and CAP leaf (own code; Lemma U / Lemma K re-derived, crude chord range first).
LEB: all admissible Q of the box inside U = [a,b]^2  =>  mu(Q) >= area(Q) = 1.
CAP: mu(Q) >= 1 - c_y d_y - c_x d_x + rho_y c_y + rho_x c_x (width lemma; needs cy >= a when Q crosses y = a, etc.);
     d <= rho suffices.  rho = min density of the line y = a (x = a) over a range containing the chord.
     Chord range: crude [c0 - hmax, c1 + hmax]; if that fails, the triangle-cap range re-derived here:
     for d <= min(s, c) the part of Q below y = a is the triangle at the lowest vertex V_B = (cx - (c - s)/2,
     cy - (c + s)/2), with legs along (c, s) (to the right, slope s/c) and (-s, c) (to the left, slope -c/s);
     at depth d the chord is [x_B - d s / c, x_B + d c / s]."""
import json
from fractions import Fraction as F
from mass import PROF
a, b = F(14, 5), F(31, 5)
REC = '/home/evand/math/square-packing/public/s12/runs/qx2_k4x_k008/qxzm_full.jsonl'
SQ2H = F(70711, 100000)   # > 1/sqrt2


def trig(u):
    N = 1 + u * u
    return (1 - u * u) / N, 2 * u / N


def below45(u): return u * u + 2 * u - 1 <= 0


def hmax(u0, u1):
    if below45(u1): c, s = trig(u1); return (c + s) / 2
    if not below45(u0): c, s = trig(u0); return (c + s) / 2   # h decreasing past 45
    return SQ2H


def hmin(u0, u1):
    c0, s0 = trig(u0); c1, s1 = trig(u1)
    return min(c0 + s0, c1 + s1) / 2


def rho(key, t0, t1):
    P = PROF.get(key, ())
    # min density over [t0, t1]; 0 if any part of [t0,t1] is uncovered
    cur = t0; best = None
    for lo, hi, r in P:
        if hi <= t0 or lo >= t1: continue
        if lo > cur: return F(0)
        best = r if best is None else min(best, r)
        cur = max(cur, hi)
    if best is None or cur < t1: return F(0)
    return best


def leb(bx):
    x0, x1, y0, y1, u0, u1 = bx
    H = hmax(u0, u1)
    return x0 - H >= a and y0 - H >= a and x1 + H <= b and y1 + H <= b


def cap(bx, refined):
    x0, x1, y0, y1, u0, u1 = bx
    H = hmax(u0, u1)
    if x1 + H > b or y1 + H > b: return False
    for axis in (0, 1):
        lo, (o0, o1) = (y0, (x0, x1)) if axis == 0 else (x0, (y0, y1))
        d = a - lo + H
        if d <= 0: continue
        if lo < a: return False
        t0, t1 = o0 - H, o1 + H
        if refined and u0 > 0:   # c - s, c, s monotone on all of (0, 90deg): valid for bins past 45deg too
            # theta in [th0, th1] within (0, 45deg]: c in [c1, c0], s in [s0, s1]; d <= min(s0, c1) => triangle cap
            c0, s0 = trig(u0); c1, s1 = trig(u1)
            if d <= min(s0, c1):
                if axis == 0:   # line y = a, vertex B lowest: x_B = cx - (c - s)/2; chord [x_B - d s/c, x_B + d c/s]
                    t0 = max(t0, x0 - (c0 - s0) / 2 - d * s1 / c1)
                    t1 = min(t1, x1 - (c1 - s1) / 2 + d * c0 / s0)
                else:           # line x = a, leftmost vertex L = (cx - (c + s)/2, cy + (c - s)/2); legs to (s, -c)
                    # (downwards, slope -c/s: y drops c/s per unit x) and (c, s) (upwards, s/c per unit x)
                    t0 = max(t0, y0 + (c1 - s1) / 2 - d * c0 / s0)
                    t1 = min(t1, y1 + (c0 - s0) / 2 + d * s1 / c1)
        key = ('H', a) if axis == 0 else ('V', a)
        if not d <= rho(key, t0, t1): return False
    return True


fh = open(REC); fh.readline()
res = {'LEB': [0, 0], 'CAP': [0, 0, 0]}
fails = []
for ln in fh:
    d = json.loads(ln)
    for bb, k in d['leaves']:
        if k not in res: continue
        bx = tuple(F(v) for v in bb)
        if k == 'LEB':
            ok = leb(bx); res[k][0 if ok else 1] += 1
            if not ok: fails.append((k, bb))
        else:
            if cap(bx, False): res[k][0] += 1
            elif cap(bx, True): res[k][1] += 1
            else: res[k][2] += 1; fails.append((k, bb))
print('LEB [own-proved, not]:', res['LEB'], ' CAP [crude range, refined range, not proved]:', res['CAP'])
for f in fails[:20]: print('  not re-proved', f)
