#!/usr/bin/env python3
"""decide.py m e d w q p [--axis] [--timeout s]: exact z3 verdict on "exists fitting pose, 1<lam<=1.01, score<=1".
sat => prints the witness and re-checks it with the Fraction evaluator (snapped to a rational angle)."""
import sys, os, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import z3
from nagqe import counterexample_query, score_exact, fits

def val(mdl, v):
    x = mdl.eval(v, model_completion=True)
    if z3.is_algebraic_value(x):
        x = x.approx(30)
    return F(x.as_fraction()) if hasattr(x, 'as_fraction') else F(str(x))

def run(m, prm, tilted, timeout):
    S, (X, Y, c, s, lam), sc = counterexample_query(m, tilted=tilted, **prm)
    S.set('timeout', int(timeout * 1000))
    t0 = time.time(); r = S.check(); dt = time.time() - t0
    out = dict(m=m, tilted=tilted, verdict=str(r), secs=round(dt, 2))
    if r == z3.sat:
        mdl = S.model()
        Xv, Yv, cv, sv, lv = (val(mdl, v) for v in (X, Y, c, s, lam))
        # snap the angle to a rational point on the circle via tan(a/2)
        tt = sv / (1 + cv) if tilted else F(0)
        tt = F(tt).limit_denominator(10 ** 12)
        cq, sq = (1 - tt * tt) / (1 + tt * tt), 2 * tt / (1 + tt * tt)
        ok = fits(m, Xv, Yv, cq, sq, lv)
        scq = score_exact(m, Xv, Yv, cq, sq, lv, **prm) if ok else None
        out.update(X=float(Xv), Y=float(Yv), deg=float(z3.RealVal(0).as_fraction()) if False else None,
                   s=float(sq), lam=float(lv), witness_fits=ok,
                   witness_score=(float(scq) if scq is not None else None),
                   witness_exact=(str(scq) if scq is not None else None))
    return out

if __name__ == '__main__':
    a = argparse.ArgumentParser()
    a.add_argument('m', type=int); a.add_argument('e'); a.add_argument('d'); a.add_argument('w'); a.add_argument('q'); a.add_argument('p')
    a.add_argument('--axis', action='store_true'); a.add_argument('--timeout', type=float, default=600)
    o = a.parse_args()
    prm = dict(e=F(o.e), d=F(o.d), w=F(o.w), q=F(o.q), p=F(o.p))
    print(run(o.m, prm, not o.axis, o.timeout), flush=True)
