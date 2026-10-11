#!/usr/bin/env python3
"""Box-wise branch pruning, then exact NLSAT on what is left.

For a box of poses (centre rectangle, optional s-range), every if-then-else condition in the score is decided
exactly by two small NLSAT calls (box & cond, box & not cond).  Decided conditions are substituted, so the final
CAD call sees only the branches that genuinely change inside the box.  Everything stays exact: a condition is
only replaced when NLSAT proves it constant on the box; undecided or timed-out conditions are left in place."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import z3
from nagqe import counterexample_query

TAC = lambda ms: z3.TryFor(z3.Then('simplify', 'elim-term-ite', 'purify-arith', 'simplify', 'qfnra-nlsat'), ms)


def decide(fmls, ms):
    g = z3.Goal()
    g.add(*fmls)
    try:
        r = list(TAC(ms)(g))
    except z3.Z3Exception:
        return 'unknown'
    if all(len(gg) == 1 and z3.is_false(gg[0]) for gg in r):
        return 'unsat'
    if all(len(gg) == 0 for gg in r):
        return 'sat'
    return 'unknown'


def ite_conds(e, seen=None, out=None):
    if seen is None:
        seen, out = set(), []
    if e.get_id() in seen:
        return out
    seen.add(e.get_id())
    if z3.is_app_of(e, z3.Z3_OP_ITE):
        out.append(e.arg(0))
    for ch in e.children():
        ite_conds(ch, seen, out)
    return out


def leaf(m, prm, box, srange=None, tilted=True, cond_ms=2000, final_ms=120000):
    """box = (x0, x1, y0, y1).  Returns dict(verdict, undecided, secs)."""
    t0 = time.time()
    S, (X, Y, c, s, lam), sc = counterexample_query(m, tilted=tilted, bound=None, **prm)
    R = lambda v: z3.RealVal(str(F(v)))
    base = list(S.assertions()) + [X >= R(box[0]), X <= R(box[1]), Y >= R(box[2]), Y <= R(box[3])]
    if srange:
        base += [s >= R(srange[0]), s <= R(srange[1])]
    if decide(base, cond_ms) == 'unsat':
        return dict(verdict='unsat', why='empty box', undecided=0, secs=round(time.time() - t0, 2))
    goal = sc <= 1
    undecided = 0
    for _ in range(4):                     # substitution can expose new conditions; iterate
        conds, subs = [], []
        for cnd in ite_conds(goal):
            if z3.is_true(cnd) or z3.is_false(cnd) or any(cnd.eq(x) for x in conds):
                continue
            conds.append(cnd)
            if decide(base + [cnd], cond_ms) == 'unsat':
                subs.append((cnd, z3.BoolVal(False)))
            elif decide(base + [z3.Not(cnd)], cond_ms) == 'unsat':
                subs.append((cnd, z3.BoolVal(True)))
        if not subs:
            undecided = len(conds)
            break
        goal = z3.simplify(z3.substitute(goal, *subs))
    else:
        undecided = len([c for c in ite_conds(goal) if not (z3.is_true(c) or z3.is_false(c))])
    v = decide(base + [goal], final_ms)
    return dict(verdict=v, undecided=undecided, secs=round(time.time() - t0, 2))


if __name__ == '__main__':
    import argparse
    a = argparse.ArgumentParser()
    a.add_argument('m', type=int)
    for k in 'edwqp':
        a.add_argument(k)
    a.add_argument('--box', nargs=4, required=True)
    a.add_argument('--srange', nargs=2)
    a.add_argument('--final', type=float, default=120)
    o = a.parse_args()
    prm = {k: F(getattr(o, k)) for k in 'edwqp'}
    print(leaf(o.m, prm, [F(b) for b in o.box], [F(x) for x in o.srange] if o.srange else None,
               final_ms=int(o.final * 1000)), flush=True)
