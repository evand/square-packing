#!/usr/bin/env python3
"""nlsat.py m e d w q p [--axis] [--timeout s] [--smt2 FILE] [--box X0 X1 Y0 Y1]

Same question as decide.py ("exists a fitting pose, 1 < lam <= 1.01, score <= 1"), but run through z3's complete
NLSAT procedure (CAD-based) after if-then-else elimination, instead of the default incremental solver.
--box restricts the centre to a rectangle (for splitting the pose space into windows)."""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import z3
from nagqe import counterexample_query


def nlsat(m, prm, tilted=True, timeout=600, box=None, smt2=None):
    S, (X, Y, c, s, lam), sc = counterexample_query(m, tilted=tilted, **prm)
    if box:
        x0, x1, y0, y1 = (z3.RealVal(str(F(b))) for b in box)
        S.add(X >= x0, X <= x1, Y >= y0, Y <= y1)
    if smt2:
        open(smt2, 'w').write(S.to_smt2())
    g = z3.Goal()
    g.add(*S.assertions())
    tac = z3.TryFor(z3.Then('simplify', 'elim-term-ite', 'purify-arith', 'simplify', 'qfnra-nlsat'),
                    int(timeout * 1000))
    t0 = time.time()
    try:
        r = tac(g)
        goals = list(r)
        if all(len(gg) == 1 and z3.is_false(gg[0]) for gg in goals):
            verdict = 'unsat'
        elif all(len(gg) == 0 for gg in goals):
            verdict = 'sat'
        else:
            verdict = 'unknown'
        model = None
        if verdict == 'sat':
            mc = r.as_expr  # noqa: model extraction below
            try:
                model = r[0].convert_model(z3.Model(r.ctx)) if False else None
            except Exception:
                model = None
    except z3.Z3Exception as ex:
        verdict, model = f'unknown ({ex})', None
    return dict(m=m, tilted=tilted, box=[str(b) for b in box] if box else None,
                verdict=verdict, secs=round(time.time() - t0, 2))


if __name__ == '__main__':
    a = argparse.ArgumentParser()
    a.add_argument('m', type=int)
    for k in 'edwqp':
        a.add_argument(k)
    a.add_argument('--axis', action='store_true')
    a.add_argument('--timeout', type=float, default=600)
    a.add_argument('--smt2')
    a.add_argument('--box', nargs=4)
    o = a.parse_args()
    prm = {k: F(getattr(o, k)) for k in 'edwqp'}
    print(nlsat(o.m, prm, not o.axis, o.timeout, o.box, o.smt2), flush=True)
