#!/usr/bin/env python3
"""t3_chord_scan: the chord-row pigeonhole, measured.

For a line L (axis 'v' = x = c, 'h' = y = c) and a set K of k squares whose INTERIORS
meet L, the open chords are pairwise disjoint subintervals of L n [delta, T-delta]^2, so

    sum_{i in K} h_i(L)  <=  T - (k+1) delta            (the chord-chain row, C3)

and hence  delta <= (T - sum h) / (k+1).  `best_line_bound` minimises the right-hand
side over lines, over k, and over the k longest chords on the line; it is the strongest
bound any chord-chain certificate can give at the configuration, because it uses the
EXACT chords (any provable lower bound H_i <= h_i can only be weaker).

F(c) = sum over all squares whose interior meets the line; it is piecewise linear in c
with breakpoints
only at c = x_i +- p_i and c = x_i +- (C_i - S_i)/2, so the maximum over c is attained
at one of 4n candidates per axis and the scan is exact, not a sample.

Commands
    packings [files...]   the 749 margin-0 packings: does a chord certificate exist?
    cone                  the coherent cone / uniform tilt: where obstruction X lives
    unif                  the uniform-tilt family, order of the chord gap in t
    col3                  the 1-D sub-problem: can three chords be forced onto one line?
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s6skel                                                        # noqa: E402

T = 3.0


def geom(z, n):
    X = np.array(z[1::3]); Y = np.array(z[2::3]); TH = np.array(z[3::3])
    C = np.abs(np.cos(TH)); S = np.abs(np.sin(TH))
    return X, Y, TH, C, S, 0.5 * (C + S)


def chords_on(z, n, a, axis, tol=1e-12):
    """(h_i), counting ONLY the squares whose INTERIOR meets the line.

    Closed semantics allows two squares to share an edge, and then the two closed
    chords on that edge's line coincide; what is pairwise disjoint is the family of
    OPEN chords {L n int(Q_i)}, which is empty unless d_i < p_i strictly.  This is the
    correction that makes the chord-chain row (C3) valid under the repo's semantics."""
    X, Y, TH, C, S, P = geom(z, n)
    ctr = X if axis == 'v' else Y
    out = np.zeros(n)
    for i in range(n):
        d = abs(ctr[i] - a)
        if d < P[i] - tol:
            out[i] = s6skel.chord_len(TH[i], d)
    return out


def candidates(z, n, axis):
    X, Y, TH, C, S, P = geom(z, n)
    ctr = X if axis == 'v' else Y
    band = 0.5 * np.abs(C - S)
    cs = set()
    for i in range(n):
        for d in (0.0, P[i], -P[i], band[i], -band[i]):
            for eps in (0.0, 1e-9, -1e-9):
                c = ctr[i] + d + eps
                if -1e-9 <= c <= T + 1e-9:
                    cs.add(round(min(max(c, 0.0), T), 12))
    return sorted(cs)


def max_F(z, n):
    """(max_c F(c), c, axis, k) -- exact, over the breakpoints."""
    best = (-1.0, None, None, 0)
    for axis in ('v', 'h'):
        for c in candidates(z, n, axis):
            h = chords_on(z, n, c, axis)
            tot = float(h.sum()); k = int((h > 1e-12).sum())
            if tot > best[0]:
                best = (tot, c, axis, k)
    return best


def best_line_bound(z, n):
    """min over (axis, c, k) of (T - sum of the k longest chords)/(k+1).

    -> (bound, c, axis, k, sum_h)
    """
    best = (1e9, None, None, 0, 0.0)
    for axis in ('v', 'h'):
        for c in candidates(z, n, axis):
            h = np.sort(chords_on(z, n, c, axis))[::-1]
            run = 0.0
            for k in range(1, n + 1):
                if h[k - 1] <= 1e-12:
                    break
                run += h[k - 1]
                v = (T - run) / (k + 1)
                if v < best[0]:
                    best = (v, c, axis, k, run)
    return best


# ------------------------------------------------------------------ the recorded packings
def load(paths):
    out = []
    for p in paths:
        if not os.path.exists(p):
            continue
        for line in open(p):
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            if 'z' not in r:
                continue
            out.append(r)
    return out


def cmd_packings(argv):
    paths = argv or ['runs/t3_exist_row1.jsonl', 'runs/t3_exist_hunt1.jsonl',
                     'runs/t3_exist_scope1.jsonl', 'runs/t3_chain_scan1.jsonl']
    rows = load(paths)
    n = 6
    pk = [r for r in rows if r.get('delta', -1) >= -1e-9]
    print('rows read %d;  in scope (delta >= -1e-9): %d' % (len(rows), len(pk)))
    fired = 0
    worst = (1e9, None)
    tab = {}
    for r in pk:
        z = r['z']
        b, c, axis, k, sh = best_line_bound(z, n)
        mf = max_F(z, n)
        nz = sum(1 for th in z[3::3]
                 if abs(((th + math.pi / 4) % (math.pi / 2)) - math.pi / 4) < 1e-9)
        key = nz
        tab.setdefault(key, [0, 0, 1e9, -1e9])
        tab[key][0] += 1
        if b <= 1e-9:
            fired += 1
            tab[key][1] += 1
        tab[key][2] = min(tab[key][2], b)
        tab[key][3] = max(tab[key][3], b)
        if b < worst[0]:
            worst = (b, (r.get('fam'), r.get('label'), c, axis, k, sh, mf[0]))
    print('chord certificate (bound <= 0) exists at %d / %d packings' % (fired, len(pk)))
    print('%-6s %6s %8s %12s %12s' % ('#axpar', 'count', 'fires', 'min bound', 'max bound'))
    for k in sorted(tab):
        v = tab[k]
        print('%-6d %6d %8d %12.6f %12.6f' % (k, v[0], v[1], v[2], v[3]))
    print('best (most negative) bound: %.9f  at %s' % (worst[0], worst[1]))
    # and the max chord sum distribution
    ms = sorted(max_F(r['z'], n)[0] for r in pk)
    print('max_c F(c) over the packings: min %.6f  median %.6f  max %.6f'
          % (ms[0], ms[len(ms) // 2], ms[-1]))
    print('  (a chord-chain certificate needs max_c F >= T = 3 at delta = 0)')


# ------------------------------------------------------------------ the coherent cone
def cmd_cone(argv):
    rows = load(['runs/t3_chain_scan1.jsonl'])
    n = 6
    fams = ['unif1', 'unif5', 'unif10', 'unif20', 'unif30', 'unif40', 'unif45',
            'coh0.5', 'coh1', 'coh2', 'coh5', 'coh10', 'coh20',
            'near0.5', 'near1', 'near2', 'near5', 'near10', 'near20',
            'Z3', 'Z4', 'Z5', 'Z6']
    print('NOTE: the chord-chain row (C3) is valid only at delta >= 0.  Every family below')
    print('      except Z3..Z6 has delta* < 0, so its "chordbound" column is NOT a valid')
    print('      bound -- it is the size of the (C3) violation, amplified by ~1/(C S).')
    print('      The meaningful column at delta* < 0 is max_c F versus T = 3.')
    print('%-9s %-26s %12s %12s %12s %6s %5s' %
          ('fam', 'label', 'delta*', 'max_c F', 'chordbound', 'k', 'axis'))
    for f in fams:
        sel = [r for r in rows if r['fam'] == f]
        if not sel:
            continue
        r = max(sel, key=lambda q: q['delta'])      # the best (largest-margin) optimum
        z = r['z']
        mf = max_F(z, n)
        b = best_line_bound(z, n)
        print('%-9s %-26s %12.3e %12.6f %12.6f %6d %5s' %
              (f, r['label'][:26], r['delta'], mf[0], b[0], b[3], b[2]))


# ------------------------------------------------------------------ the uniform family
def cmd_unif(argv):
    """delta*(t) against the chord gap T - max_c F(t), along the uniform-tilt line,
    recomputed from scratch with s6skel's instrument."""
    rows = load(['runs/t3_chain_scan1.jsonl'])
    print('NOTE: delta* < 0 throughout, so "chordbound" is not a valid bound (see `cone`).')
    print('%8s %14s %14s %14s %10s' %
          ('t (deg)', 'delta*', 'T - max_c F', 'gap / t', 'chordbound'))
    for lab, tdeg in (('uniform 1.0 deg', 1.0), ('uniform 5.0 deg', 5.0),
                      ('uniform 10.0 deg', 10.0), ('uniform 20.0 deg', 20.0),
                      ('uniform 30.0 deg', 30.0), ('uniform 40.0 deg', 40.0),
                      ('uniform 45.0 deg', 45.0)):
        sel = [r for r in rows if r['label'] == lab]
        if not sel:
            continue
        r = max(sel, key=lambda q: q['delta'])
        z = r['z']
        mf = max_F(z, 6)
        b = best_line_bound(z, 6)
        t = math.radians(tdeg)
        print('%8.1f %14.6e %14.6e %14.4f %10.6f'
              % (tdeg, r['delta'], T - mf[0], (T - mf[0]) / t, b[0]))


# ------------------------------------------------------------------ the 1-D sub-problem
def cmd_col3(argv):
    """The pigeonhole in its strongest form: six trapezoids h_i(c) of area 1, height
    sec t, plateau half-width (C-S)/2, support half-width p, with centres free in
    [p, T-p].  What is  min over centre placements of  max_c sum_i h_i(c)?
    (This is an upper bound on what ANY chord-chain certificate can extract from the
    x-positions alone -- it ignores every other constraint, so if it is < T there is no
    chord-chain certificate that pigeonhole alone can produce.)"""
    print('%8s %10s %10s %12s %12s %12s' %
          ('t(deg)', 'p', 'sec t', 'min-max F', '3 sec t', 'area/T'))
    for tdeg in (0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 45.0):
        t = math.radians(tdeg)
        C, S = math.cos(t), math.sin(t)
        p = 0.5 * (C + S)
        band = 0.5 * (C - S)
        hmax = 1.0 / C

        def F(xs, c):
            s_ = 0.0
            for x in xs:
                d = abs(x - c)
                if d <= band:
                    s_ += hmax
                elif d < p:
                    s_ += (p - d) / (C * S) if C * S > 1e-12 else 0.0
            return s_

        def maxF(xs):
            cand = set()
            for x in xs:
                for dd in (0.0, band, -band, p, -p):
                    cand.add(x + dd)
            return max(F(xs, c) for c in cand)

        # 2 + 2 + 2 columns is the obvious minimiser; also try a random search
        lo, hi = p, T - p
        best = (1e9, None)
        rng = np.random.default_rng(11)
        for trial in range(4000):
            if trial == 0:
                xs = [lo, lo, 0.5 * (lo + hi), 0.5 * (lo + hi), hi, hi]
            else:
                xs = list(rng.uniform(lo, hi, 6))
            cur = maxF(xs)
            for _ in range(300):
                j = int(rng.integers(6))
                old = xs[j]
                xs[j] = min(max(old + rng.normal(0, 0.06), lo), hi)
                v = maxF(xs)
                if v <= cur:
                    cur = v
                else:
                    xs[j] = old
            if cur < best[0]:
                best = (cur, sorted(round(x, 4) for x in xs))
        print('%8.1f %10.6f %10.6f %12.6f %12.6f %12.6f'
              % (tdeg, p, 1.0 / C, best[0], 3.0 / C, 6.0 / T))
        print('         placement %s' % (best[1],))


def cmd_band(argv):
    """The leaf hypothesis the chord certificate needs, against what pigeonhole gives.

    A chord-chain certificate on a vertical line needs THREE squares whose centres are
    within (C-S)/2 of the line, i.e. three x-coordinates inside a window of length
    C - S (then each chord is exactly sec t and the three sum to 3 sec t > T = 3).
    Six centres lie in an interval of length T - u, so pigeonhole delivers three inside
    a window of length (T - u)/2 and no less.  The certificate therefore needs
        (T - u)/2  <=  C - S     <=>     3 cos t - sin t >= 3,
    which holds only at t = 0.  Columns:
      needed      C - S
      pigeonhole  (T - u)/2                        (sharp: the 2-2-2 placement)
      gap         (T-u)/2 - (C-S)                  = t/2 + O(t^2) at T = 3
      actual      the smallest 3-window actually realised at the recorded optimum,
                  minus (C - S)  -- the amount by which the TRUTH fails
    """
    rows = load(['runs/t3_chain_scan1.jsonl'])
    print('%8s %10s %10s %12s %12s %12s %12s' %
          ('t(deg)', 'C-S', '(T-u)/2', 'pigeon gap', 'actual w3', 'actual gap', 'delta*'))
    for tdeg in (0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 45.0):
        lab = 'uniform %.1f deg' % tdeg
        sel = [r for r in rows if r['label'] == lab]
        t = math.radians(tdeg)
        C, S = math.cos(t), math.sin(t)
        need = C - S
        pig = (T - C - S) / 2.0
        if not sel:
            print('%8.1f %10.6f %10.6f %12.6f %12s %12s %12s'
                  % (tdeg, need, pig, pig - need, '-', '-', '-'))
            continue
        r = max(sel, key=lambda q: q['delta'])
        z = r['z']
        best = 1e9
        for ctr in (z[1::3], z[2::3]):
            v = sorted(ctr)
            for i in range(len(v) - 2):
                best = min(best, v[i + 2] - v[i])
        print('%8.1f %10.6f %10.6f %12.6f %12.6f %12.3e %12.3e'
              % (tdeg, need, pig, pig - need, best, best - need, r['delta']))


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'packings'
    globals()['cmd_' + cmd](sys.argv[2:])


if __name__ == '__main__':
    main()
