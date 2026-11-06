"""slack_separated / RELEASE_TOL: a pair is released only when another side line clears it by more than the contact
tolerance (11-05, s(303): a pair cleared by 4.9e-16 was released and the corner-corner MILP reported a spurious
first-order descent of -2.2e-4).

  python3 -m pytest search/exact/test_slack_release.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf
import exactsolve as E


def _pair(dy):
    """Two upright unit squares meeting corner to corner: B's lower-left corner on A's right side line (x = 3.5),
    and A's top line (y = 3.5) clearing B by dy."""
    mp.dps = 50
    P = E.Pack(2, mpf(10), [mpf(3), mpf(4)], [mpf(3), mpf(4) + mpf(dy)], [mpf(0), mpf(0)])
    C, Sn = P.trig()
    H = mpf(1) / 2
    # the incidence the KKT system would use: a corner of B at gap 0 on a side line of A
    eq = [('C', 1, a, 0, k) for k in range(4) for a in range(4)
          if abs(E.eval_contact(('C', 1, a, 0, k), P.X, P.Y, C, Sn, P.S, H, order=1)[0]) < mpf(10) ** -40
          and min(E.eval_contact(('C', 1, b, 0, k), P.X, P.Y, C, Sn, P.S, H, order=1)[0] for b in range(4)) > -mpf(10) ** -40]
    assert eq, 'no gap-0 incidence found'
    return P, eq[:1]


def test_tiny_clearance_not_released():
    P, eq = _pair('1e-15')
    assert E.slack_separated(P, eq, E.RELEASE_TOL) == {}
    # the old threshold (10^-(dps-15)) released it
    assert (0, 1) in E.slack_separated(P, eq, mpf(10) ** -65)


def test_real_clearance_released():
    P, eq = _pair('1e-8')
    rel = E.slack_separated(P, eq, E.RELEASE_TOL)
    assert (0, 1) in rel and abs(float(rel[(0, 1)][0]) - 1e-8) < 1e-12
