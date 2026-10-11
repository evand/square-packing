#!/usr/bin/env python3
"""lifted.py m e d w q p [--timeout s] [--box x0 x1 y0 y1] [--srange s0 s1] [--tactic T]

Whole-problem exact decision on the division-free encoding, with if-then-else lifted to the Boolean level
(cofactor-term-ite), so NLSAT works in exactly the five pose variables (X, Y, c, s, lam).
Tilted case only (0 < s <= c); the axis-parallel slice is decided by nlsat.py --axis."""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import z3
from nagqe import counterexample_query_df

TACTICS = {
    'cofactor': lambda: z3.Then('simplify', 'cofactor-term-ite', 'simplify', 'qfnra-nlsat'),
    'elim':     lambda: z3.Then('simplify', 'elim-term-ite', 'simplify', 'qfnra-nlsat'),
}


def run(m, prm, timeout=600, box=None, srange=None, tactic='cofactor'):
    S, (X, Y, c, s, lam), sk, K = counterexample_query_df(m, **prm)
    R = lambda v: z3.RealVal(str(F(v)))
    fml = list(S.assertions())
    if box:
        fml += [X >= R(box[0]), X <= R(box[1]), Y >= R(box[2]), Y <= R(box[3])]
    if srange:
        fml += [s >= R(srange[0]), s <= R(srange[1])]
    g = z3.Goal(); g.add(*fml)
    t0 = time.time()
    try:
        r = list(z3.TryFor(TACTICS[tactic](), int(timeout * 1000))(g))
        if all(len(x) == 1 and z3.is_false(x[0]) for x in r):
            v = 'unsat'
        elif all(len(x) == 0 for x in r):
            v = 'sat'
        else:
            v = 'unknown'
    except z3.Z3Exception as ex:
        v = f'unknown ({ex})'
    return dict(m=m, prm={k: str(x) for k, x in prm.items()}, box=box and [str(b) for b in box],
                srange=srange and [str(x) for x in srange], tactic=tactic, verdict=v, secs=round(time.time() - t0, 2))


if __name__ == '__main__':
    a = argparse.ArgumentParser()
    a.add_argument('m', type=int)
    for k in 'edwqp':
        a.add_argument(k)
    a.add_argument('--timeout', type=float, default=600)
    a.add_argument('--box', nargs=4)
    a.add_argument('--srange', nargs=2)
    a.add_argument('--tactic', default='cofactor')
    o = a.parse_args()
    prm = {k: F(getattr(o, k)) for k in 'edwqp'}
    print(run(o.m, prm, o.timeout, o.box and [F(b) for b in o.box], o.srange and [F(x) for x in o.srange],
              o.tactic), flush=True)
