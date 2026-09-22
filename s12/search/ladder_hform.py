#!/usr/bin/env python3
"""ladder_hform: the THREE-level ladder H1, in closed form, checked square by square.

New file (task `tasks/counting-ladder`); nothing in `search/` is modified.  It imports
`t3_chain_hform` (the two-level Lemma H closed forms) and extends it by one rung.

Everything is for the disjunctive LP of `S6_SKELETON.md` sec 3.1 with the walls carrying delta
(`BANDCUT_K.md` sec 1.1); rows read `a.c - delta >= b` and a normalised non-negative combination
with `sum w a = 0` gives `delta <= -sum w b` (`t3-chain.md` sec 1.1).

Two ladders:

  UNIFORM (uladder / H1u).  Every square at tilt t.  Lemma H of `t3-chain.md` sec 1.3 --
  crossbar of T at weight 1, two legs of a and c links at weight tan t, six wall rows -- with
  `r in {0,1,2}` of the two closing MAIN-direction wall rows replaced by a RUNG: one further
  main-direction link at weight tan^2 t onto a fresh square, plus that square's two wall rows.
  r = 0 is Lemma H itself (`t3_chain_hform.H`).

  MIXED (mladder / H1m).  The `eps^3` cell of `T4_CYCLES.md` sec 3: a crossbar of T of which
  exactly two ADJACENT members are tilted by t and the rest are axis-parallel; A axial links
  before them, B after (A + B = T - 3, A >= 1); one leg of l1 links at the entry interface and
  one of l2 links at the exit interface (L = l1 + l2); one rung at the kink of the entry leg.

Usage:  python3 search/ladder_hform.py [check|cell89|pinwheel|tables|all]
        python3 search/ladder_hform.py t11  runs/bandcut_scan_exact110.json ...   # legs per kink
        python3 search/ladder_hform.py t11h runs/bandcut_scan_exact110.json ...   # best dual/chain
"""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t3_chain_hform as HF                                          # noqa: E402  (the old file)


def W(d):
    return abs(math.cos(d)) + abs(math.sin(d))


def m_of(d):
    return 0.5 + 0.5 * W(d)


def P_of(th):
    return 0.5 * W(th)


# ===================================================================== rows and their checking
class Sup(object):
    """a Farkas support: rows (name, {square: (ax, ay)}, b) with a weight each."""

    def __init__(self, T):
        self.T = T
        self.rows = []
        self.tilt = []

    def sq(self, tilt):
        self.tilt.append(tilt)
        return len(self.tilt) - 1

    def link(self, i, j, tau, w, name='link'):
        """separation of the pair (i, j) along the unit normal at angle tau, j on the + side."""
        n = (math.cos(tau), math.sin(tau))
        self.rows.append((name, {j: n, i: (-n[0], -n[1])},
                          m_of(self.tilt[j] - self.tilt[i]), w))

    def wall(self, i, ax, side, w, name=None):
        e = (1.0, 0.0) if ax == 'x' else (0.0, 1.0)
        if side == 'hi':
            e = (-e[0], -e[1])
        b = P_of(self.tilt[i]) - (self.T if side == 'hi' else 0.0)
        self.rows.append((name or ('%s-%s' % (side, ax)), {i: e}, b, w))

    def value(self, hard=False):
        """-> dict(resid, wmin, D, N, value).  `resid` is the max uncancelled centre coefficient;
        the certificate is valid iff resid = 0 and wmin >= 0."""
        ns = len(self.tilt)
        co = np.zeros((ns, 2))
        D = N = 0.0
        Dp = 0.0
        for (nm, d, b, w) in self.rows:
            for i, v in d.items():
                co[i] += w * np.array(v)
            D += w
            N += w * b
            if len(d) == 2:
                Dp += w
        den = Dp if hard else D
        return dict(resid=float(np.max(np.abs(co))), wmin=min(r[3] for r in self.rows),
                    D=D, Dp=Dp, N=N, value=-N / den)

    def lp(self, hard=False):
        """the BEST bound on this support (the weights re-optimised), for comparison."""
        from scipy.optimize import linprog
        ns = len(self.tilt)
        nr = len(self.rows)
        A = np.zeros((2 * ns + 1, nr))
        b = np.array([r[2] for r in self.rows])
        for q, (nm, d, bb, w) in enumerate(self.rows):
            for i, v in d.items():
                A[2 * i, q] += v[0]
                A[2 * i + 1, q] += v[1]
        A[-1, :] = [1.0 if (len(r[1]) == 2 or not hard) else 0.0 for r in self.rows]
        rhs = np.zeros(2 * ns + 1)
        rhs[-1] = 1.0
        res = linprog(-b, A_eq=A, b_eq=rhs, bounds=[(0, None)] * nr, method='highs')
        return float(res.fun) if res.success else math.inf


# ===================================================================== the uniform-tilt ladder
def uladder(T, a, c, t, r=2):
    """Lemma H (crossbar T, legs a and c) with r of its two closing main walls made into rungs."""
    C, S = math.cos(t), math.sin(t)
    nu = S * S / (C * C)                                   # rung weight = tan^2 t
    wc = S * S / C                                         # closing main-wall weight = S tan t
    s = Sup(T)
    ch = [s.sq(t) for _ in range(T)]
    for q in range(T - 1):
        s.link(ch[q], ch[q + 1], t, 1.0, 'cross%d' % q)
    s.wall(ch[0], 'x', 'lo', 1.0 / C)
    s.wall(ch[-1], 'x', 'hi', 1.0 / C)
    for (leg, at, side, sgn) in ((a, ch[0], 'lo', -1), (c, ch[-1], 'hi', +1)):
        prev = at
        for q in range(leg):
            i = s.sq(t)
            s.link(prev, i, t + sgn * math.pi / 2, S / C, 'leg')
            prev = i
        rung = r > 0
        r -= 1
        s.wall(prev, 'y', side, S + (nu * S if rung else 0.0))
        oside = 'hi' if side == 'lo' else 'lo'
        if rung:
            e = s.sq(t)
            if side == 'lo':
                s.link(prev, e, t, nu, 'rung')
            else:
                s.link(e, prev, t, nu, 'rung')
            s.wall(e, 'x', oside, nu * C)
            s.wall(e, 'y', oside, nu * S)
        else:
            s.wall(prev, 'x', oside, wc)
    return s


def H1u(T, L, t, r=2):
    """closed form of uladder.  r = 0 is t3_chain_hform.H(T, L, t)."""
    C, S = math.cos(t), math.sin(t)
    u = C + S
    P = u / 2.0
    Mv = HF.M_of(t)
    nu = S * S / (C * C)
    wc = S * S / C
    D = (T - 1) + L * S / C + 2 * Mv
    N = (T - 1) + L * S / C + (u - T) * Mv
    for _ in range(r):
        D += -wc + nu * (1.0 + C + S) + nu * S
        N += -wc * (P - T) + nu * 1.0 + nu * C * (P - T) + nu * S * (P - T) + nu * S * P
    return -N / D


# ===================================================================== the mixed-tilt ladder
def mladder(T, A, B, l1, l2, t, rung=True, legs=True):
    """crossbar of T with s_A, s_{A+1} at tilt t and the rest axis-parallel.

    link normals: 0..A-1 axial (owned by the axis-parallel square on the left), A and A+1 tilted
    (owned by a tilted square), A+2..T-2 axial again.  The normal changes at s_A (entry) and at
    s_{A+2} (exit), and a leg hangs at each; with legs=False those two residuals are taken by
    transverse WALL rows instead (mode `CW` of T4_CYCLES sec 3.2)."""
    assert A >= 1 and B >= 0 and A + B == T - 3
    C, S = math.cos(t), math.sin(t)
    s = Sup(T)
    ch = [s.sq(t if q in (A, A + 1) else 0.0) for q in range(T)]
    if not legs:      # CW: the two interface residuals taken by transverse WALL rows instead
        for q in range(T - 1):
            tau = t if q in (A, A + 1) else 0.0
            s.link(ch[q], ch[q + 1], tau, 1.0 / C if q in (A, A + 1) else 1.0, 'cross%d' % q)
        s.wall(ch[0], 'x', 'lo', 1.0)
        s.wall(ch[-1], 'x', 'hi', 1.0)
        s.wall(ch[A], 'y', 'lo', S / C)
        s.wall(ch[A + 2], 'y', 'hi', S / C)
        return s
    for q in range(T - 1):
        tau = t if q in (A, A + 1) else 0.0
        w = 1.0 if q < A else (C if q <= A + 1 else C * C)
        s.link(ch[q], ch[q + 1], tau, w, 'cross%d' % q)
    s.wall(ch[0], 'x', 'lo', 1.0)
    s.wall(ch[-1], 'x', 'hi', C * C)
    # entry leg: down from s_A to the lo-y wall.  first link normal tilted, the rest axial.
    prev = ch[A]
    kink = None
    for q in range(l1):
        i = s.sq(0.0)
        s.link(prev, i, (t if q == 0 else 0.0) - math.pi / 2, S if q == 0 else S * C, 'leg1')
        if q == 0:
            kink = i
        prev = i
    s.wall(prev, 'y', 'lo', S * C)
    if rung:
        e = s.sq(0.0)
        s.link(kink, e, 0.0, S * S, 'rung')
        s.wall(e, 'x', 'hi', S * S)
    else:
        s.wall(kink, 'x', 'hi', S * S)
    # exit leg: up from s_{A+2} to the hi-y wall, all links axial
    prev = ch[A + 2]
    for q in range(l2):
        i = s.sq(0.0)
        s.link(prev, i, math.pi / 2, C * S, 'leg2')
        prev = i
    s.wall(prev, 'y', 'hi', C * S)
    return s


def H1m(T, A, B, l1, l2, t, rung=True):
    C, S = math.cos(t), math.sin(t)
    u = C + S
    g = 0.5 * (1.0 + u)                                    # m of a 0 <-> t link
    rw = 2.0 if rung else 1.0
    D = (A + 1) + 2 * C + (B + 1) * C * C + S + (l1 + l2 + 1) * S * C + rw * S * S
    N = (0.5 + (A - 1) + g + C + C * g + B * C * C + C * C * (0.5 - T)
         + S * g + (l1 - 1) * S * C + 0.5 * S * C
         + (1.0 if rung else 0.0) * S * S + S * S * (0.5 - T)
         + l2 * C * S + C * S * (0.5 - T))
    return -N / D


def H1mw(T, l1, l2, t, rung=True):
    """the WALL-BLOCK mixed ladder: the two tilted crossbar members are s_0 and s_1, so the
    tilted square sits AGAINST the main wall (A = 0).  This is the T = 11 / T = 12 band shape.
    Series: N = (L + 2 - T) t + (3 - T) t^2 + (2T - 2L - 5/2) t^3/3 + ..., so even at L = T - 2
    it is +(T-3) t^2/(T+1) > 0 for T > 3."""
    C, S = math.cos(t), math.sin(t)
    u = C + S
    g = 0.5 * (1.0 + u)
    rw = 2.0 if rung else 1.0
    D = 1 + 2 * C + (T - 2) * C * C + S + (l1 + l2 + 1) * S * C + rw * S * S
    N = (u / 2 + C + C * g + (T - 3) * C * C + C * C * (0.5 - T) + S * g
         + (l1 - 1) * S * C + 0.5 * S * C
         + (1.0 if rung else 0.0) * S * S + S * S * (0.5 - T)
         + l2 * C * S + C * S * (0.5 - T))
    return -N / D


def mladderw(T, l1, l2, t, rung=True):
    """the wall-block ladder as an explicit support (for checking H1mw)."""
    C, S = math.cos(t), math.sin(t)
    s = Sup(T)
    ch = [s.sq(t if q in (0, 1) else 0.0) for q in range(T)]
    s.wall(ch[0], 'x', 'lo', 1.0)
    s.link(ch[0], ch[1], t, C, 'cross0')
    s.link(ch[1], ch[2], t, C, 'cross1')
    for q in range(2, T - 1):
        s.link(ch[q], ch[q + 1], 0.0, C * C, 'cross%d' % q)
    s.wall(ch[-1], 'x', 'hi', C * C)
    prev, kink = ch[0], None
    for q in range(l1):
        i = s.sq(0.0)
        s.link(prev, i, (t if q == 0 else 0.0) - math.pi / 2, S if q == 0 else S * C, 'leg1')
        if q == 0:
            kink = i
        prev = i
    s.wall(prev, 'y', 'lo', S * C)
    if rung:
        e = s.sq(0.0)
        s.link(kink, e, 0.0, S * S, 'rung')
        s.wall(e, 'x', 'hi', S * S)
    else:
        s.wall(kink, 'x', 'hi', S * S)
    prev = ch[2]
    for q in range(l2):
        i = s.sq(0.0)
        s.link(prev, i, math.pi / 2, C * S, 'leg2')
        prev = i
    s.wall(prev, 'y', 'hi', C * S)
    return s


def need_Lw(T, kink, cap=200):
    """least L for the WALL-BLOCK ladder at interface angle `kink`."""
    if kink <= 1e-12:
        return 0
    for L in range(2, cap):
        if H1mw(T, 1, L - 1, kink) <= 0:
            return L
    return None


def CWm(T, A, B, t):
    """the same crossbar with no legs at all (T4_CYCLES's `CW`)."""
    C, S = math.cos(t), math.sin(t)
    u = C + S
    g = 0.5 * (1.0 + u)
    D = (A + B + 2) + 2.0 / C + 2.0 * S / C
    N = (0.5 + (A - 1) + g + (1.0 / C) + (1.0 / C) * g + B + (0.5 - T)
         + (S / C) * (0.5 * u) + (S / C) * (0.5 - T))
    return -N / D


# ===================================================================== checks
def check(verbose=True):
    bad = [0]

    def rep(tag, ok, *a):
        if not ok:
            bad[0] += 1
            if verbose and bad[0] < 20:
                print('  FAIL', tag, *['%.12g' % x if isinstance(x, float) else x for x in a])

    for T in (3, 4, 5, 11):
        for a in (0, 1, 2):
            for c in (0, 1, 2):
                for td in (0.5, 5, 20, 45):
                    t = math.radians(td)
                    v = uladder(T, a, c, t, r=0).value()
                    cf = HF.H(T, a + c, t)
                    rep('LemmaH T=%d a=%d c=%d t=%g' % (T, a, c, td),
                        v['resid'] < 1e-12 and v['wmin'] > -1e-15
                        and abs(v['value'] - cf) < 1e-11 and abs(H1u(T, a + c, t, 0) - cf) < 1e-11,
                        v['resid'], v['wmin'], v['value'], cf)
                    for rr in (1, 2):
                        v = uladder(T, a, c, t, r=rr).value()
                        rep('H1u T=%d a=%d c=%d t=%g r=%d' % (T, a, c, td, rr),
                            v['resid'] < 1e-12 and v['wmin'] > -1e-15
                            and abs(v['value'] - H1u(T, a + c, t, rr)) < 1e-11,
                            v['resid'], v['wmin'], v['value'], H1u(T, a + c, t, rr))
    for T in (4, 5, 6, 11):
        for A in range(1, T - 2):
            B = T - 3 - A
            for (l1, l2) in ((1, 1), (2, 1), (1, 2), (max(1, T - 3), 1)):
                for td in (0.5, 1, 5, 20):
                    t = math.radians(td)
                    for rung in (True, False):
                        v = mladder(T, A, B, l1, l2, t, rung).value()
                        cf = H1m(T, A, B, l1, l2, t, rung)
                        rep('H1m T=%d A=%d B=%d l=%d,%d t=%g rung=%d'
                            % (T, A, B, l1, l2, td, rung),
                            v['resid'] < 1e-12 and v['wmin'] > -1e-15
                            and abs(v['value'] - cf) < 1e-11,
                            v['resid'], v['wmin'], v['value'], cf)
                for (l1, l2) in ((1, 1), (2, 3)):
                    for td in (1, 20):
                        t = math.radians(td)
                        for rung in (True, False):
                            v = mladderw(T, l1, l2, t, rung).value()
                            rep('H1mw T=%d l=%d,%d t=%g rung=%d' % (T, l1, l2, td, rung),
                                v['resid'] < 1e-12 and v['wmin'] > -1e-15
                                and abs(v['value'] - H1mw(T, l1, l2, t, rung)) < 1e-11,
                                v['resid'], v['wmin'], v['value'], H1mw(T, l1, l2, t, rung))
                v = mladder(T, A, B, 1, 1, math.radians(5), legs=False).value()
                rep('CWm T=%d A=%d' % (T, A),
                    v['resid'] < 1e-12 and v['wmin'] > -1e-15
                    and abs(v['value'] - CWm(T, A, B, math.radians(5))) < 1e-11,
                    v['resid'], v['wmin'], v['value'], CWm(T, A, B, math.radians(5)))
    print('check: %s   (Lemma H against t3_chain_hform; the 1- and 2-rung uniform ladder; the '
          'mixed ladder and its leg-free CW; T in {3,4,5,6,11}, all A, l1, l2, '
          't in {0.5,1,5,20,45} deg.  Every weight >= 0 and every centre coordinate cancelled '
          'to < 1e-12.)' % ('OK' if bad[0] == 0 else '%d FAILURES' % bad[0]))
    return bad[0]


# ===================================================================== the eps^3 cell
def cell89():
    meas = {1: dict(CW=+6.872476e-03, H=+5.955303e-05, H1=-5.285352e-07, d=-5.259969e-07),
            5: dict(CW=+3.226175e-02, H=+1.363023e-03, H1=-5.805332e-05, d=-6.188939e-05),
            10: dict(CW=+5.956405e-02, H=+4.902941e-03, H1=-4.061638e-04, d=-4.596117e-04)}
    print('\n# T = 4, the eps^3 cell (k = 8, 9 chain-free) of BANDCUT_K sec 3(b).')
    print('# ladder: crossbar of four, its two middle members tilted by eps, A = 1, B = 0,')
    print('#         one leg of one link at each interface (L = 2), one rung.')
    print('# closed form  vs  runs/t4_cycles_cell89.txt (LP duals on the same supports)')
    print('%7s | %12s %12s | %12s %12s | %12s %12s %12s'
          % ('eps deg', 'CW closed', 'CW meas', 'H closed', 'H meas', 'H1 closed', 'H1 meas',
             'delta* meas'))
    for d in sorted(meas):
        mm = meas[d]
        t = math.radians(d)
        print('%7g | %12.6e %12.6e | %12.6e %12.6e | %12.6e %12.6e %12.6e'
              % (d, CWm(4, 1, 0, t), mm['CW'], H1m(4, 1, 0, 1, 1, t, False), mm['H'],
                 H1m(4, 1, 0, 1, 1, t, True), mm['H1'], mm['d']))
    print('\n# divided by the power of eps (rad) the level structure predicts')
    print('%7s | %9s %9s | %9s %9s | %9s %9s %9s'
          % ('eps deg', 'CW/eps', 'meas', 'H/eps^2', 'meas', 'H1/eps^3', 'meas', 'delta*/eps^3'))
    for d in sorted(meas):
        mm = meas[d]
        t = math.radians(d)
        print('%7g | %9.4f %9.4f | %9.4f %9.4f | %9.4f %9.4f %9.4f'
              % (d, CWm(4, 1, 0, t) / t, mm['CW'] / t,
                 H1m(4, 1, 0, 1, 1, t, False) / t ** 2, mm['H'] / t ** 2,
                 H1m(4, 1, 0, 1, 1, t, True) / t ** 3, mm['H1'] / t ** 3, mm['d'] / t ** 3))
    print('\n# eps -> 0 limits of the closed form:')
    for (lab, f, p) in (('CW /eps  ', lambda t: CWm(4, 1, 0, t), 1),
                        ('H  /eps^2', lambda t: H1m(4, 1, 0, 1, 1, t, False), 2),
                        ('H1 /eps^3', lambda t: H1m(4, 1, 0, 1, 1, t, True), 3)):
        print('   %s -> %s   (eps = 1e-2, 1e-3, 1e-4, 1e-5 deg)'
              % (lab, ' '.join('%+.8f' % (f(math.radians(e)) / math.radians(e) ** p)
                               for e in (1e-2, 1e-3, 1e-4, 1e-5))))
    print('   predicted:   CW/eps -> (T-2)/(T+1) = %.8f ; H/eps^2 -> 1/(T+1) = %.8f ; '
          'H1/eps^3 -> -1/(2(T+1)) = %.8f' % (2.0 / 5, 1.0 / 5, -0.1))
    print('   (below eps ~ 1e-3 deg the H1 row is float cancellation, not the closed form: N is')
    print('    a difference of O(1) terms whose first three orders cancel.  The exact series is')
    print('    N = (L+2-T) e - B e^2 + (2T-2L-5/2) e^3/3 + (T/3-A/3-5/4) e^4, by sympy.)')


# ===================================================================== the pinwheel duals
def pinwheel():
    print('\n# S6_LOCAL sec 3 reads the two exact pinwheel duals as')
    print('#   n = 6,  T = 3: chain of three at 1/2 - t/4,  a transverse chain at t/2,  '
          'two links at t^2/2')
    print('#   n = 12, T = 4: chain of four  at 1/3 - 2t/9, a transverse chain at t/3, '
          'links at t^2/3')
    print('# The uniform ladder with L = T - 2 legs and r = 2 rungs, in HARD containment')
    print('# (sum of the PAIR weights = 1), has crossbar weight')
    print('#   w0 = 1 / [ (T-1) + L tan t + 2 tan^2 t ]  =  1/(T-1) - L t/(T-1)^2 + O(t^2)')
    for (T, a, c) in ((3, 1, 0), (4, 1, 1)):
        L = a + c
        print('\n  T = %d, L = %d = T-2, r = 2:   1/(T-1) = %.6f, -L/(T-1)^2 = %.6f per rad'
              % (T, L, 1.0 / (T - 1), -L / (T - 1.0) ** 2))
        print('   %6s %12s %12s %12s %12s %12s' % ('t deg', 'crossbar w0', 'predicted',
                                                   'leg/w0', 'rung/w0', 'value (hard)'))
        for td in (0.5, 1, 2, 4):
            t = math.radians(td)
            s = uladder(T, a, c, t, r=2)
            v = s.value(hard=True)
            w = {}
            for (nm, d, b, ww) in s.rows:
                w.setdefault(nm.split('-')[0].rstrip('0123456789'), []).append(ww)
            w0 = w['cross'][0]
            pred = 1.0 / ((T - 1) + L * math.tan(t) + 2 * math.tan(t) ** 2)
            print('   %6.2f %12.6f %12.6f %12.6f %12.6f %12.6e'
                  % (td, w0 / v['Dp'], pred, w['leg'][0] / w0, w['rung'][0] / w0,
                     v['value']))
        if T == 3:
            print('   exact n = 6 pinwheel value -3(u-1)^2/(u^2+3):  %s'
                  % ' '.join('%+.6e' % (-3 * (W(math.radians(td)) - 1) ** 2
                                        / (W(math.radians(td)) ** 2 + 3))
                             for td in (0.5, 1, 2, 4)))


# ===================================================================== tables
def tables():
    print('\n# H1m(T, L, eps), the MIXED ladder: N_1 = L + 2 - T is the first-order coefficient')
    print('# of the numerator, so the eps-ORDER is decided by L against T - 2 = L*(T, 0).')
    print('# (tilted block at the crossbar end: A = T-3, B = 0, which is what N_2 = -B needs)')
    print('%5s %4s %13s %13s %13s %8s  %s' % ('T', 'L', 'eps=0.5 deg', '1 deg', '5 deg',
                                              'order', 'sign'))
    for T in (4, 5, 6, 11):
        for L in range(max(2, T - 4), T + 1):
            A, B = T - 3, 0
            l1, l2 = 1, L - 1
            vs = [H1m(T, A, B, l1, l2, math.radians(d)) for d in (0.5, 1, 5)]
            o = math.log(abs(vs[1] / vs[0])) / math.log(2.0)
            print('%5d %4d %13.6e %13.6e %13.6e %8.2f  %s'
                  % (T, L, vs[0], vs[1], vs[2], o, 'KILLS' if vs[1] <= 0 else '-'))
    print('\n# the same without the rung (two-level H): order eps^2 at L = T - 2, wrong sign')
    print('%5s %4s %13s %13s %8s' % ('T', 'L', 'eps=0.5 deg', '1 deg', 'order'))
    for T in (4, 5, 6, 11):
        L = T - 2
        A, B = T - 3, 0
        vs = [H1m(T, A, B, 1, L - 1, math.radians(d), rung=False) for d in (0.5, 1)]
        print('%5d %4d %13.6e %13.6e %8.2f'
              % (T, L, vs[0], vs[1], math.log(abs(vs[1] / vs[0])) / math.log(2.0)))
    print('\n# the UNIFORM ladder H1u(T, L, t, r) against Lemma H (r = 0), and L*(T, t)')
    print('%5s %4s %8s %13s %13s %13s %9s' % ('T', 'L', 't deg', 'r=0 (Lemma H)', 'r=1', 'r=2',
                                              'L*(T,t)'))
    for T in (3, 4, 5, 11):
        for L in (T - 2, T - 1):
            for td in (1, 5, 20):
                t = math.radians(td)
                print('%5d %4d %8g %13.6e %13.6e %13.6e %9.3f'
                      % (T, L, td, H1u(T, L, t, 0), H1u(T, L, t, 1), H1u(T, L, t, 2),
                         HF.Lstar(T, t)))
    print('\n# the mixed ladder at a BAND tilt: how many transverse links a single interface')
    print('# would have to carry.  L_need = least L with H1m(T, 1, T-4, 1, L-1, t) <= 0.')
    print('%5s %8s %10s %10s %10s' % ('T', 't deg', 'L_need', 'T - 2', 'n - T'))
    for T in (4, 5, 11, 12):
        for td in (1, 5, 10, 20, 25, 29):
            print('%5d %8g %10s %10d %10d' % (T, td, need_L(T, math.radians(td)), T - 2,
                                              T * T - 2 * T))


# ===================================================================== T = 11 / T = 12
def tilt_of(th):
    t = math.degrees(th) % 90.0
    return math.radians(t if t <= 45.0 else t - 90.0)


def leg_graphs(T, delta, sq, tr, tol=1e-9):
    """For each distinct tilt class th of the packing, the directed graph of pairs separated
    along the SINGLE transverse normal n(th) (+tr direction), n an edge normal of one of the two
    squares -- i.e. the links a transverse chain with a COMMON normal may use.

    -> dict th -> adjacency list (i -> list of j)."""
    n = len(sq)
    ths = sorted(set(round(tilt_of(s[2]), 9) for s in sq))
    out = {}
    for th in ths:
        g = [[] for _ in range(n)]
        # transverse normal at angle th from the transverse axis
        if tr == 'y':
            nv = (-math.sin(th), math.cos(th))
        else:
            nv = (math.cos(th), -math.sin(th))
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                if (abs(tilt_of(sq[i][2]) - th) > tol and abs(tilt_of(sq[j][2]) - th) > tol):
                    continue                      # not an edge normal of either square
                dx, dy = sq[j][0] - sq[i][0], sq[j][1] - sq[i][1]
                m = 0.5 + 0.5 * W(sq[j][2] - sq[i][2])
                if nv[0] * dx + nv[1] * dy - m >= delta - tol:
                    g[i].append(j)
        out[th] = g
    return out


def longest_from(g, v, avoid, order):
    """longest path in g starting at v, avoiding `avoid`.  For a FIXED normal n every edge
    strictly increases n.c, so g is acyclic and `order` (the vertices sorted by n.c) is a
    topological order: this is a linear-time DAG longest path, no simple-path search."""
    memo = [0] * len(g)
    for u in reversed(order):
        if u in avoid and u != v:
            continue
        best = 0
        for j in g[u]:
            if j in avoid:
                continue
            if 1 + memo[j] > best:
                best = 1 + memo[j]
        memo[u] = best
    return memo[v]


def t11(path, cap=1e-6, top=400, out=None):
    """the ladders available on the tight wall-to-wall chains of a verified positive packing."""
    import t11_chains as TC11
    T, delta, sq = TC11.load(path)
    n = len(sq)
    Tn = int(round(T))
    print('\n=== %s :  T = %g  n = %d  delta = %+.9e  (need L = T-2 = %d transverse links) ==='
          % (os.path.basename(path), T, n, delta, Tn - 2))
    res = []
    for axis in 'xy':
        tr = 'y' if axis == 'x' else 'x'
        edges = TC11.build(T, delta, sq, axis, 1e-9)
        lo, hi = TC11.walls(T, delta, sq, axis)
        gg = leg_graphs(T, delta, sq, tr)
        ordr = {}
        for th in gg:
            nv = ((-math.sin(th), math.cos(th)) if tr == 'y'
                  else (math.cos(th), -math.sin(th)))
            ordr[th] = sorted(range(n), key=lambda i: nv[0] * sq[i][0] + nv[1] * sq[i][1])
        rg = {}
        for th, g in gg.items():
            r = [[] for _ in range(n)]
            for i in range(n):
                for j in g[i]:
                    r[j].append(i)
            rg[th] = r
        paths = TC11.k_best_paths(n, edges, lo, hi, top, cap)
        seen = set()
        for (g0, path0, links) in paths:
            d = TC11.describe(T, delta, sq, axis, g0, path0, links)
            key = tuple(i for i, t in zip(path0, d['tilts']) if abs(t) > 1.0)
            if key in seen:
                continue
            seen.add(key)
            taus = [0.0] + [math.radians(x) for x in d['taus']] + [0.0]
            changes = [s for s in range(len(path0)) if abs(taus[s + 1] - taus[s]) > 1e-9]
            legs = []
            for s in changes:
                v = path0[s]
                thv = round(tilt_of(sq[v][2]), 9)
                av = set(path0)
                # the leg's first link must use the transverse normal OWNED BY the interface
                # square (that is what cancels the kink); the rest must share it.
                up = longest_from(gg[thv], v, av, ordr[thv]) if thv in gg else 0
                dn = longest_from(rg[thv], v, av, ordr[thv][::-1]) if thv in rg else 0
                # and, for comparison, the best over any single normal class
                anyup = max(longest_from(gg[th], v, av, ordr[th]) for th in gg)
                anydn = max(longest_from(rg[th], v, av, ordr[th][::-1]) for th in rg)
                legs.append((s, math.degrees(tilt_of(sq[v][2])), up, dn, anyup, anydn))
            res.append(dict(axis=axis, slack=g0, k=len(path0), path=list(path0),
                            tilts=d['tilts'], taus=d['taus'], Lbudget=d['L'], C=d['C'],
                            bound=d['bound'], rise=d['rise'], changes=changes, legs=legs,
                            Lown=sum(max(l[2], l[3]) for l in legs),
                            Lany=sum(max(l[4], l[5]) for l in legs),
                            kink=max((abs(taus[s + 1] - taus[s]) for s in changes), default=0.0)))
    print('%3s %3s %8s %8s %5s %6s %6s %7s %7s %9s  %s'
          % ('ax', 'k', 'slack', 'maxkink', 'nchg', 'L own', 'L any', 'L need', 'H1mw',
             'rise', 'legs: at, tilt, links up/down (own normal | best single normal)'))
    for r in res:
        need = need_Lw(Tn, r['kink'])
        r['Lneed'] = need
        r['H1mw'] = H1mw(Tn, 1, max(r['Lown'] - 1, 1), r['kink'])
        print('%3s %3d %8.1e %8.1f %5d %6d %6d %7s %7.4f %9.2f  %s'
              % (r['axis'], r['k'], r['slack'], math.degrees(r['kink']), len(r['changes']),
                 r['Lown'], r['Lany'], need, r['H1mw'], r['rise'],
                 '; '.join('s%d t=%+.1f %d/%d (%d/%d)' % (l[0], l[1], l[2], l[3], l[4], l[5])
                           for l in r['legs'])))
    print('   L need = least L making the wall-block ladder H1mw(T, 1, L-1, maxkink) <= 0;')
    print('   H1mw = that ladder evaluated at the L the packing actually offers (generously: the')
    print('   longer of the two directions at each kink, and no check that the directions are the')
    print('   ones the kink signs demand).  A positive H1mw means no such certificate closes.')
    if out:
        json.dump(res, open(out, 'w'), indent=1)
        print('wrote %s' % out)
    return res


def t11h(path, cap=1e-6, top=400, tol=1e-9, out=None):
    """The DECISIVE test: the best H-shaped / ladder dual on each tight chain of a verified
    positive packing, using the same restricted-dual instrument as `t3_chain.py` at T = 3.

    Modes (t3_chain.h_support): C = the chain's links + its two end walls; CW = + every wall row;
    H = + EVERY transverse-type pair row (so it is an upper bound on what any leg structure of
    any length can achieve); HP = + the main-direction rows between chain squares; H1 = H + the
    single best further main-direction pair row (the rung)."""
    import t11_chains as TC11
    import t3_chain as TC
    T, delta, sq = TC11.load(path)
    n = len(sq)
    z = np.zeros(1 + 3 * n)
    for i, (x, y, th) in enumerate(sq):
        z[1 + 3 * i], z[2 + 3 * i], z[3 + 3 * i] = x, y, th
    TH = np.array(z[3::3])
    A, b, tags = TC.all_rows(z, n, T, delta, tol=tol)
    print('\n=== %s : best ladder dual per tight chain (%d rows in the system) ==='
          % (os.path.basename(path), len(b)))
    print('%3s %3s %8s  %12s %12s %12s %12s   %s'
          % ('ax', 'k', 'slack', 'C', 'CW', 'H', 'H1', 'tilts along the chain'))
    res = []
    for axis in 'xy':
        edges = TC11.build(T, delta, sq, axis, 1e-9)
        lo, hi = TC11.walls(T, delta, sq, axis)
        seen = set()
        for (g0, path0, links) in TC11.k_best_paths(n, edges, lo, hi, top, cap):
            d = TC11.describe(T, delta, sq, axis, g0, path0, links)
            key = tuple(i for i, t in zip(path0, d['tilts']) if abs(t) > 1.0)
            if key in seen:
                continue
            seen.add(key)
            vals = {}
            for mode in ('C', 'CW', 'H'):
                keep = TC.h_support(tags, TH, axis, list(path0), mode)
                vals[mode] = TC.dual_bound(A, b, keep)[0]
            # H1 = H + the best single further main-direction pair row (the rung).  The full
            # dual is delta* > 0 here, so every mode must be > 0; the point is HOW FAR.
            base = TC.h_support(tags, TH, axis, list(path0), 'H')
            extra = [q for q, tg in enumerate(tags)
                     if tg[0] == 'pair' and q not in base and TC.row_axis(tg, TH) == axis
                     and (tg[1] in path0 or tg[2] in path0)]
            best = math.inf
            for q in extra[:24]:
                v = TC.dual_bound(A, b, base | {q})[0]
                best = min(best, v)
            vals['H1'] = best
            vals['nextra'] = len(extra)
            print('%3s %3d %8.1e  %12.5e %12.5e %12.5e %12.5e   %s'
                  % (axis, len(path0), g0, vals['C'], vals['CW'], vals['H'], vals['H1'],
                     ' '.join('%+.0f' % t for t in d['tilts'])))
            res.append(dict(axis=axis, k=len(path0), slack=g0, path=list(path0),
                            tilts=d['tilts'], **vals))
    print('   (delta = %+.9e; a mode is a valid upper bound on delta, so a value > delta is '
          'no contradiction and a value > 0 means that shape of certificate does not exist here)'
          % delta)
    if out:
        json.dump(res, open(out, 'w'), indent=1)
        print('wrote %s' % out)
    return res


def need_L(T, kink):
    """least L = l1 + l2 for which the BEST-PLACED mixed ladder (tilted block at the crossbar's
    far end: A = T-3, B = 0) with interface angle `kink` is <= 0."""
    if kink <= 1e-12:
        return 0
    for L in range(2, 8 * T):
        if H1m(T, T - 3, 0, 1, L - 1, kink) <= 0:
            return L
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', nargs='?', default='all')
    ap.add_argument('arg', nargs='*')
    ap.add_argument('--cap', type=float, default=1e-6)
    ap.add_argument('--top', type=int, default=400)
    a = ap.parse_args()
    if a.cmd in ('check', 'all'):
        check()
    if a.cmd in ('cell89', 'all'):
        cell89()
    if a.cmd in ('pinwheel', 'all'):
        pinwheel()
    if a.cmd in ('tables', 'all'):
        tables()
    if a.cmd in ('t11', 't11h'):
        for p in a.arg:
            if a.cmd == 't11h':
                t11h(p, cap=a.cap, top=a.top,
                     out='runs/ladder_h_%s.json' % os.path.basename(p).split('.')[0])
                continue
            t11(p, cap=a.cap, top=a.top,
                out='runs/ladder_%s.json' % os.path.basename(p).split('.')[0])


if __name__ == '__main__':
    main()
