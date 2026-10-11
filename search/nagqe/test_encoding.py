#!/usr/bin/env python3
"""Cross-check the z3 score expression against the Fraction evaluator at random exact rational poses."""
import random, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import z3
from nagqe import counterexample_query, score_exact, fits, NAGAMOCHI

def check(m, prm, n, seed=1):
    rnd = random.Random(seed); bad = 0; done = 0
    while done < n:
        t = F(rnd.randint(0, 4142), 10000)            # tan(a/2), a in [0, 45deg]
        c, s = (1 - t * t) / (1 + t * t), 2 * t / (1 + t * t)
        lam = 1 + F(rnd.randint(1, 100), 10000)
        X = F(rnd.randint(0, 5000 * m), 10000); Y = F(rnd.randint(0, 5000 * m), 10000)
        if not fits(m, X, Y, c, s, lam):
            continue
        done += 1
        want = score_exact(m, X, Y, c, s, lam, **prm)
        S, (Xv, Yv, cv, sv, lv), sc = counterexample_query(m, tilted=(s > 0), bound=None, **prm)
        S.add(Xv == z3.RealVal(str(X)), Yv == z3.RealVal(str(Y)), cv == z3.RealVal(str(c)),
              sv == z3.RealVal(str(s)), lv == z3.RealVal(str(lam)))
        S.push(); S.add(sc != z3.RealVal(str(want))); r = S.check(); S.pop()
        if r != z3.unsat:
            bad += 1; print('MISMATCH', m, X, Y, c, s, lam, want, r)
    print(f'm={m}: {n} poses, {bad} mismatches')
    return bad

def check_df(m, prm, n, seed=1):
    from nagqe import counterexample_query_df
    rnd = random.Random(seed); bad = 0; done = 0
    while done < n:
        t = F(rnd.randint(1, 4142), 10000)
        c, s = (1 - t * t) / (1 + t * t), 2 * t / (1 + t * t)
        lam = 1 + F(rnd.randint(1, 100), 10000)
        X = F(rnd.randint(0, 5000 * m), 10000); Y = F(rnd.randint(0, 5000 * m), 10000)
        if not fits(m, X, Y, c, s, lam):
            continue
        done += 1
        want = score_exact(m, X, Y, c, s, lam, **prm) * (c * s) ** 2
        S, (Xv, Yv, cv, sv, lv), sc, K = counterexample_query_df(m, bound=None, **prm)
        S.add(Xv == z3.RealVal(str(X)), Yv == z3.RealVal(str(Y)), cv == z3.RealVal(str(c)),
              sv == z3.RealVal(str(s)), lv == z3.RealVal(str(lam)))
        S.push(); S.add(sc != z3.RealVal(str(want))); r = S.check(); S.pop()
        if r != z3.unsat:
            bad += 1; print('DF MISMATCH', m, X, Y, c, s, lam, want, r)
    print(f'df m={m}: {n} poses, {bad} mismatches')
    return bad



if __name__ == '__main__':
    tot = 0
    for m in (4, 5, 6):
        tot += check(m, NAGAMOCHI, 150, seed=m)
    tot += sum(check_df(m, NAGAMOCHI, 150, seed=10 + m) for m in (4, 5, 6))
    sys.exit(1 if tot else 0)
