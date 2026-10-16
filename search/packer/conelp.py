#!/usr/bin/env python3
"""Cone LP for weakcone.py, standalone (numpy + scipy only) so it can run under a newer HiGHS (ipenv.sh: scipy 1.18) when
the system scipy's HiGHS 1.8 returns model status 'Unknown' (postsolve leaves 3e-5 infeasibility on some s(110) points).
  max sum t  s.t.  J d >= t, 0 <= t <= 1, |d|_inf <= 1        python3 conelp.py IN.npz OUT.npy   (IN: JV without S column)
"""
import sys
import numpy as np
from scipy.optimize import linprog


def solve(JV):
    K, m = JV.shape
    c = np.r_[np.zeros(m), -np.ones(K)]
    for meth in ('highs', 'highs-ds', 'highs-ipm'):
        r = linprog(c, A_ub=np.hstack([-JV, np.eye(K)]), b_ub=np.zeros(K), bounds=[(-1, 1)] * m + [(0, 1)] * K, method=meth)
        if r.x is not None:
            return r.x[m:]
    return None


if __name__ == '__main__':
    t = solve(np.load(sys.argv[1])['JV'])
    np.save(sys.argv[2], t if t is not None else np.array([]))
