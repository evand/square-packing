#!/usr/bin/env python3
"""t3_chord_ks: Kearney-Shiu 2002 (EJC 9 #R14) section 3 rewritten in the repo's row
language, with chord rows.  Nothing in search/ is modified; s6skel and t3_chord_geom
are imported.

Source read for this script: the paper itself (doi 10.37236/1631), section 3 verbatim.

Commands
    points      the 13 points and their A/B/C classification
    unavoid     is the 7-point set (6) unavoidable for [0,3]^2?          [measured]
    lemma1      any unit square covering C covers a B-point              [measured]
    ineq7       Lemma 2's two intercepts and the three inequalities (7)  [measured/proved]
    lemma3      Lemma 3 and the form the proof actually uses             [measured]
    certb       configuration (b): the Farkas certificate, residual and value
    all         everything
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import t3_chord_geom as G                                            # noqa: E402

R2 = math.sqrt(2.0)
GREEN = [(R2 - 0.5, 1.0), (1.5, 1.0), (3.5 - R2, 1.0), (1.5, 1.5),
         (R2 - 0.5, 2.0), (1.5, 2.0), (3.5 - R2, 2.0)]
RED = [(3.0 - q[1], q[0]) for q in GREEN]          # 90 deg rotation about (3/2, 3/2)
CPT = (1.5, 1.5)
BPTS = [(1.5, 1.0), (1.5, 2.0), (1.0, 1.5), (2.0, 1.5)]
APTS = [p for p in set(GREEN) | set(RED) if p != CPT and p not in BPTS]


# ---------------------------------------------------------------- containment functional
def inside(c, th, p):
    """>= 0 iff the point p lies in the closed unit square (centre c, angle th)."""
    C, S = math.cos(th), math.sin(th)
    dx, dy = p[0] - c[0], p[1] - c[1]
    return 0.5 - max(abs(dx * C + dy * S), abs(-dx * S + dy * C))


def inside_v(cx, cy, th, p):
    C, S = np.cos(th), np.sin(th)
    dx, dy = p[0] - cx, p[1] - cy
    return 0.5 - np.maximum(np.abs(dx * C + dy * S), np.abs(-dx * S + dy * C))


def cmd_points():
    print('C-point  ', CPT)
    print('B-points ', sorted(BPTS))
    print('A-points ', sorted(APTS), '  (%d)' % len(APTS))
    print('green    ', sorted(GREEN))
    print('red      ', sorted(RED))
    d = sorted({round(math.hypot(p[0] - CPT[0], p[1] - CPT[1]), 9)
                for p in set(GREEN) | set(RED)})
    print('distances from C:', d)
    print('|green u red| =', len(set(GREEN) | set(RED)))


# ---------------------------------------------------------------- unavoidability (measured)
def cmd_unavoid(nth=721, ngrid=361):
    """min over unit squares inside [0,3]^2 of max_p inside(c, th, p) over the 7 green
    points.  >= 0 says the set is unavoidable.  Grid + local polish."""
    worst = (1e9, None)
    for k in range(nth):
        th = k * (math.pi / 2) / nth
        p = 0.5 * (abs(math.cos(th)) + abs(math.sin(th)))
        xs = np.linspace(p, 3 - p, ngrid)
        X, Y = np.meshgrid(xs, xs, indexing='ij')
        best = None
        for q in GREEN:
            v = inside_v(X, Y, th, q)
            best = v if best is None else np.maximum(best, v)
        i = int(np.argmin(best))
        if best.flat[i] < worst[0]:
            worst = (float(best.flat[i]), (float(X.flat[i]), float(Y.flat[i]), th))
    print('unavoidable set (6): min over squares of max_p inside = %+.6f' % worst[0])
    print('  argmin  c = (%.5f, %.5f)  th = %.4f deg' %
          (worst[1][0], worst[1][1], math.degrees(worst[1][2])))
    print('  verdict:', 'UNAVOIDABLE (>= 0)' if worst[0] >= -1e-9 else 'AVOIDABLE (< 0)')
    return worst[0]


# ---------------------------------------------------------------- Lemma 1 (measured)
def cmd_lemma1(nth=901, ngrid=361):
    """min over unit squares covering C of max_B inside(B).  >= 0 is Lemma 1."""
    worst = (1e9, None)
    for k in range(nth):
        th = k * (math.pi / 2) / nth
        r = 0.7072
        xs = np.linspace(CPT[0] - r, CPT[0] + r, ngrid)
        ys = np.linspace(CPT[1] - r, CPT[1] + r, ngrid)
        X, Y = np.meshgrid(xs, ys, indexing='ij')
        ok = inside_v(X, Y, th, CPT) >= 0.0
        best = None
        for q in BPTS:
            v = inside_v(X, Y, th, q)
            best = v if best is None else np.maximum(best, v)
        best = np.where(ok, best, 1e9)
        i = int(np.argmin(best))
        if best.flat[i] < worst[0]:
            worst = (float(best.flat[i]), (float(X.flat[i]), float(Y.flat[i]), th))
    print('Lemma 1: min over squares covering C of max_B inside(B) = %+.6f' % worst[0])
    print('  argmin  c = (%.5f, %.5f)  th = %.4f deg' %
          (worst[1][0], worst[1][1], math.degrees(worst[1][2])))
    print('  verdict:', 'HOLDS' if worst[0] >= -1e-9 else 'FAILS')
    return worst[0]


# ---------------------------------------------------------------- Lemma 2 and (7)
def lemma2_square(t):
    """The unit square of Lemma 2: a corner on the x-axis, edge angle th = 2 arctan t,
    the point (0,1) on the opposite edge.  -> (c, th)."""
    th = 2.0 * math.atan(t)
    C, S = math.cos(th), math.sin(th)
    # corner on the x-axis at (x0, 0); the square occupies the wedge above it.  Its
    # lowest vertex is the corner, so c = (x0, 0) + (p_x, p_y) with the half-diagonal
    # pointing up: c = corner + ((C - S)/2 + S*?, ...).  Solve instead from the two
    # stated incidences, which determine the square: put c = (cx, cy) with
    #   lowest vertex on y = 0  ->  cy = (C + S)/2
    #   (0,1) on an edge        ->  one of the two slab equalities is tight.
    cy = 0.5 * (C + S)
    # (0,1) on the LEFT/UPPER edge: the slab value -dx*S + dy*C = +1/2 with
    # dx = -cx, dy = 1 - cy.
    cx = (0.5 - (1.0 - cy) * C) / S if S > 1e-12 else 0.5
    return (cx, cy), th


def cmd_ineq7(n=200001):
    """Lemma 2's two coordinates, checked against the exact chord geometry, and (7)."""
    w1 = w2 = 0.0
    m1 = m2 = m3 = 1e9
    arg = [None, None, None]
    for k in range(n):
        t = k / (n - 1.0)
        c, th = lemma2_square(t)
        u1 = (1 + t * t) / (1 + t)             # claimed: on the line y = 1
        u2 = (1 + 2 * t - t * t) / 2.0         # claimed: on the line x = 1
        w1 = max(w1, abs(inside(c, th, (u1, 1.0))))
        w2 = max(w2, abs(inside(c, th, (1.0, u2))))
        if u1 < m1:
            m1, arg[0] = u1, t
        if u2 < m2:
            m2, arg[1] = u2, t
        if u1 + u2 < m3:
            m3, arg[2] = u1 + u2, t
    print('Lemma 2: the two stated points lie ON the square (|inside| = 0)')
    print('   max |inside((1+t^2)/(1+t), 1)|      = %.3e' % w1)
    print('   max |inside(1, (1+2t-t^2)/2)|       = %.3e' % w2)
    print('inequality (7), 0 <= t <= 1:')
    print('   min (1+t^2)/(1+t)        = %.9f  at t = %.5f   (2 sqrt2 - 2 = %.9f)'
          % (m1, arg[0], 2 * R2 - 2))
    print('   min (1+2t-t^2)/2         = %.9f  at t = %.5f   (1/2)' % (m2, arg[1]))
    print('   min of the sum           = %.9f  at t = %.5f   (3/2)' % (m3, arg[2]))
    # the first inequality IS the repo's wall-strip chord lemma at a = 1
    best = (1e9, None)
    for k in range(200001):
        th = k * (math.pi / 2) / 200000
        C, S = math.cos(th), math.sin(th)
        if C * S < 1e-12:
            v = 1.0 if False else (C + S - 1.0) / max(C * S, 1e-300)
            continue
        v = (C + S - 1.0) / (C * S)            # Stromquist f(th, 1), = chord at d = 1 - p
        if v < best[0]:
            best = (v, th)
    print('cross-check: min_th (C + S - 1)/(CS) = %.9f at %.4f deg  == 2 sqrt2 - 2'
          % (best[0], math.degrees(best[1])))
    return m1, m2, m3


# ---------------------------------------------------------------- Lemma 3 (measured)
def cmd_lemma3(nth=1801, ngrid=501):
    """Lemma 3 as the proof uses it: V covers C = (3/2,3/2) and none of
    (1,3/2), (1,2), (2,3/2), (2,2); how high can the BOTTOM of V's chord on x = 1 be,
    and must V meet x = 1 at all?"""
    blocked = [(1.0, 1.5), (1.0, 2.0), (2.0, 1.5), (2.0, 2.0)]
    best = (-1e9, None)
    nmeet = 0
    tot = 0
    for k in range(nth):
        th = k * (math.pi / 2) / nth
        p = 0.5 * (abs(math.cos(th)) + abs(math.sin(th)))
        xs = np.linspace(max(p, CPT[0] - 0.71), min(3 - p, CPT[0] + 0.71), ngrid)
        ys = np.linspace(max(p, CPT[1] - 0.71), min(3 - p, CPT[1] + 0.71), ngrid)
        X, Y = np.meshgrid(xs, ys, indexing='ij')
        ok = inside_v(X, Y, th, CPT) >= 0.0
        ok &= inside_v(X, Y, th, (1.5, 2.0)) >= 0.0      # V covers the green B-point in region 2
        for q in blocked:
            ok &= inside_v(X, Y, th, q) < 0.0
        if not ok.any():
            continue
        # bottom of the chord on x = 1, for those that meet it
        meets = ok & (np.abs(X - 1.0) <= p)
        tot += int(ok.sum())
        nmeet += int(meets.sum())
        if meets.any():
            Xi, Yi = X[meets], Y[meets]
            C, S = math.cos(th), math.sin(th)
            e = 1.0 - Xi
            b0 = (-0.5 - e * C) / S if S > 1e-12 else np.full(Xi.shape, -np.inf)
            b1 = (e * S - 0.5) / C if C > 1e-12 else np.full(Xi.shape, -np.inf)
            lo = Yi + np.maximum(b0, b1)
            j = int(np.argmax(lo))
            if lo[j] > best[0]:
                best = (float(lo[j]), (float(Xi[j]), float(Yi[j]), th))
    print('Lemma 3 (as used): V covers C and (3/2,2), avoids (1,3/2),(1,2),(2,3/2),(2,2)')
    print('   admissible grid points: %d;  of them meeting x = 1: %d  (%.1f %%)'
          % (tot, nmeet, 100.0 * nmeet / max(tot, 1)))
    print('   max bottom of V-chord on x = 1 = %.6f   (Lemma 3 claims <= 5/3 = %.6f)'
          % (best[0], 5.0 / 3.0))
    print('   argmax  c = (%.5f, %.5f)  th = %.4f deg'
          % (best[1][0], best[1][1], math.degrees(best[1][2])))
    return best[0], nmeet, tot


# ---------------------------------------------------------------- the certificate
def cmd_certb():
    """Configuration (b), final step, as a Farkas certificate in the centre variables.

    Squares:  W covers the red A-point (1, sqrt2 - 1/2)   [region 7]
              X covers the green A-point (sqrt2 - 1/2, 1) [region 7] and the
                red B-point (1, 3/2)                      [region 4]
              V covers C and the green B-point (3/2, 2)   [region 2]
    Line L: x = 1.  alpha_i, beta_i are the ends of the chord of square i on L; on a
    leaf that fixes the active slab branch each is affine in (x_i, y_i).

      r1  alpha_X - beta_W          >= delta            chord disjointness, X above W
      r2  beta_W                    >= sqrt2 - 1/2      W covers (1, sqrt2 - 1/2)
      r3  alpha_V - beta_X          >= delta            chord disjointness, V above X
      r4  -alpha_V                  >= -5/3             Lemma 3
      r5  beta_X - alpha_X          >= 2 sqrt2 - 2      Lemma 2 + (7) first inequality

    weights (1,1,1,1,1): every centre variable cancels identically, for ANY branch
    assignment, and
      0 >= 2 delta + (sqrt2 - 1/2) - 5/3 + (2 sqrt2 - 2),
      delta <= (13/6 + 2 - 3 sqrt2)/2 = -0.037986...
    """
    # explicit representatives, only to exhibit the affine forms; the cancellation is
    # independent of them.
    reps = {'W': (1.30, 0.70, math.radians(20.0)),
            'X': (1.30, 1.62, math.radians(35.0)),
            'V': (1.62, 2.10, math.radians(10.0))}
    idx = {'W': 0, 'X': 1, 'V': 2}
    nv = 6
    forms = {}
    for nm, (x, y, th) in reps.items():
        al, be, bl, bh = G.chord(x, y, th, 1.0, 'v')
        lo, hi = G._affine_endpoints(th, 1.0, 'v', bl, bh)
        for tag, f, val in (('alpha', lo, al), ('beta', hi, be)):
            a = np.zeros(nv)
            a[2 * idx[nm]] = f[0]
            a[2 * idx[nm] + 1] = f[1]
            forms[(nm, tag)] = (a, f[2])
            chk = f[0] * x + f[1] * y + f[2]
            assert abs(chk - val) < 1e-12, (nm, tag, chk, val)
    rows = []          # (name, coefficient vector, rhs constant, delta coefficient)

    def add(name, expr, rhs, dcoef):
        a = np.zeros(nv); k = 0.0
        for (nm, tag, sgn) in expr:
            f, c = forms[(nm, tag)]
            a = a + sgn * f
            k += sgn * c
        rows.append((name, a, rhs - k, dcoef))

    add('r1 alpha_X - beta_W >= delta', [('X', 'alpha', +1), ('W', 'beta', -1)], 0.0, 1.0)
    add('r2 beta_W >= sqrt2 - 1/2', [('W', 'beta', +1)], R2 - 0.5, 0.0)
    add('r3 alpha_V - beta_X >= delta', [('V', 'alpha', +1), ('X', 'beta', -1)], 0.0, 1.0)
    add('r4 -alpha_V >= -5/3', [('V', 'alpha', -1)], -5.0 / 3.0, 0.0)
    add('r5 beta_X - alpha_X >= 2 sqrt2 - 2', [('X', 'beta', +1), ('X', 'alpha', -1)],
        2 * R2 - 2, 0.0)
    w = np.ones(5)
    A = np.array([r[1] for r in rows])
    b = np.array([r[2] for r in rows])
    dc = np.array([r[3] for r in rows])
    resid = float(np.max(np.abs(w @ A)))
    dtot = float(w @ dc)
    val = -float(w @ b) / dtot     # rows are a.c >= b + dcoef*delta
    print('rows (coefficients on (x_W,y_W,x_X,y_X,x_V,y_V), rhs, delta-coef):')
    for (nm, a, bb, d) in rows:
        print('   %-32s  %s   b = %+.9f   dcoef = %g'
              % (nm, np.array2string(a, precision=4, suppress_small=True,
                                     max_line_width=200), bb, d))
    print('weights          w = %s' % w)
    print('residual  max|w.A| = %.3e        (exact cancellation)' % resid)
    print('sum of delta-coefficients = %g' % dtot)
    print('bound      delta <= %.12f' % val)
    exact = (13.0 / 6.0 + 2.0 - 3.0 * R2) / 2.0
    print('closed form  (13/6 + 2 - 3 sqrt2)/2 = %.12f   |diff| = %.2e'
          % (exact, abs(exact - val)))
    print("K-S's own numbers: 13/6 - sqrt2 = %.9f  <  2 sqrt2 - 2 = %.9f  (gap %.9f)"
          % (13 / 6 - R2, 2 * R2 - 2, (2 * R2 - 2) - (13 / 6 - R2)))
    return resid, val


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'all'
    order = ['points', 'unavoid', 'lemma1', 'ineq7', 'lemma3', 'certb']
    todo = order if cmd == 'all' else [cmd]
    for c in todo:
        print('\n===== %s =====' % c)
        globals()['cmd_' + c]()


if __name__ == '__main__':
    main()
