#!/usr/bin/env python3
"""bandcut_scan_exact.py -- exact rational certificate that a band packing fits.

`delta > 0` measured in floating point proves nothing.  This script turns a float
configuration into an **exactly rational** one and evaluates the margin in `Fraction`
arithmetic with the pair/wall rows of `search/S6_SKELETON.md` §3.1:

 1. every angle is replaced by a nearby *rational* angle: pick `u in Q` close to
    `tan(theta/2)`, then `cos theta = (1-u^2)/(1+u^2)`, `sin theta = 2u/(1+u^2)` are rational
    and `|cos|+|sin|`, the wall half-widths and every pair width `1/2 + W(Dtheta)/2` stay in Q;
 2. the centres are re-optimised at those angles by the disjunctive LP (float), which recovers
    what the angle rounding cost;
 3. the centres are rounded to rationals with a fixed denominator;
 4. `bandcut_scan.exact_delta` evaluates `min(min pair gap, min wall slack)` exactly.

A positive exact value is a proof (modulo the rows, which are the separating-axis theorem)
that the `n` closed unit squares are pairwise disjoint inside `[0,T]^2`, i.e. `s(n) < T`.

    python3 search/bandcut_scan_exact.py --in runs/bandcut_scan_p132.json --T 12
"""
import argparse
import json
import math
import os
import sys
from fractions import Fraction

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from bandcut_scan import Q2, exact_delta, float_value, polish          # noqa: E402


def rationalise_angles(TH, den):
    """theta -> nearest angle whose half-tangent is a rational with denominator <= den."""
    us, cs, th2 = [], [], []
    for th in TH:
        t = math.fmod(th, math.pi / 2.0)
        if t < 0:
            t += math.pi / 2.0
        u = Fraction(math.tan(t / 2.0)).limit_denominator(den)
        c = (1 - u * u) / (1 + u * u)
        s = 2 * u / (1 + u * u)
        us.append(u); cs.append((c, s))
        th2.append(2.0 * math.atan(float(u)))
    return np.array(th2), cs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inp', required=True)
    ap.add_argument('--T', type=float, required=True)
    ap.add_argument('--aden', type=int, default=10 ** 7, help='denominator for tan(theta/2)')
    ap.add_argument('--cden', type=int, default=10 ** 9, help='denominator for the centres')
    ap.add_argument('--cut', type=float, default=3.2)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()

    d = json.load(open(a.inp))
    sq = d['sq'] if 'sq' in d else d['chain'][0]['sq']
    X = np.array([q[0] for q in sq]); Y = np.array([q[1] for q in sq])
    TH = np.array([q[2] for q in sq])
    s = d.get('s', a.T)
    lam = a.T / s
    X, Y = X * lam, Y * lam
    T = a.T
    n = len(X)
    print('n = %d   T = %g   (source side s = %.17g, dilation %.17g)' % (n, T, s, lam))
    print('float delta, published angles  : %+.9e' % float_value(X, Y, TH, T)[0])

    TH2, cs = rationalise_angles(TH, a.aden)
    print('max angle move by rationalising: %.3e rad' % float(np.abs(TH2 - (TH % (math.pi / 2))).max()))
    print('float delta, rational angles   : %+.9e' % float_value(X, Y, TH2, T)[0])
    v, X2, Y2 = polish(X, Y, TH2, T, cut=a.cut, iters=40, rounds=4)
    print('float delta, after LP re-solve : %+.9e' % v)

    cx = [Q2(Fraction(int(round(x * a.cden)), a.cden)) for x in X2]
    cy = [Q2(Fraction(int(round(y * a.cden)), a.cden)) for y in Y2]
    csq = [(Q2(c), Q2(sn)) for (c, sn) in cs]
    Tq = Fraction(int(round(T * 1)), 1) if abs(T - round(T)) < 1e-12 else Fraction(T)
    val, wit = exact_delta(cx, cy, csq, Tq, report=True)
    print()
    print('EXACT delta (Fraction)         : %s' % val.a)
    print('EXACT delta (float)            : %+.12e' % float(val))
    print('sign                           : %+d      binding row: %s' % (val.sign(), wit))
    print()
    if val.sign() > 0:
        print('=> %d closed unit squares at exactly rational angles and centres are pairwise'
              % n)
        print('   disjoint inside [0,%g]^2 with margin %s  ==>  s(%d) < %g, exactly.'
              % (T, val.a, n, T))
    else:
        print('=> NOT a certificate: the exact margin is not positive.')
    if a.out:
        json.dump({'T': T, 'n': n, 'delta_num': str(val.a.numerator),
                   'delta_den': str(val.a.denominator),
                   'cden': a.cden, 'aden': a.aden,
                   'sq': [[str(x.a), str(y.a), str(c[0].a), str(c[1].a)]
                          for x, y, c in zip(cx, cy, csq)]}, open(a.out, 'w'))
        print('wrote', a.out)


if __name__ == '__main__':
    main()
