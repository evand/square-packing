#!/usr/bin/env python3
"""linecuts.py -- the general-line chord inequality, as a family of packing-LP rows.

`notes/chord-lemma.md` proves, for the BOTTOM WALL of the container, that at most 3 pairwise
CLOSED-disjoint unit squares inside `[0,t]^2`, `t <= 4`, can have their centre within 1 of the wall.
The proof goes through the horizontal line `y = 9/10`: every such square meets it in a chord of
length `>= 1`, the chords are pairwise disjoint closed intervals, and four closed intervals of
length `>= 1` do not fit in a closed segment of length `4`.

Nothing in that argument uses that the line is near a wall, or that it is horizontal.  For ANY
line `L`, with `Lam` = the length of `L` intersected with the container,

    mu( { S : |chord of S on L| >= 1 } )  <=  ceil(Lam) - 1                       (the LINE CUT)

is valid for every packing measure.  This module builds those rows.

The predicate (`notes/chord-lemma.md` Lemma 1 + Corollary 2, rotated into the line's frame): write
`phi` for the direction angle of `L`, `psi = theta - phi` for the square's angle relative to it and
`d` for the distance from the square's centre to `L`.  Then

    |chord| >= 1   <=>   d <= D(psi),      D(psi) = (|cos psi| + |sin psi|)/2 - |cos psi sin psi| .

`D` decreases from `1/2` at `psi = 0` (mod 90 deg) to `(sqrt 2 - 1)/2 = 0.20711` at `psi = 45 deg`:
a tilted square can hide between two lines, an axis-parallel one cannot hide from a family of
horizontal lines of pitch `< 1`.

BOUNDARY CONVENTION.  `mu(B) <= N` is valid for every `B` that is a SUBSET of the true
`{chord >= 1}` set, and NOT for a superset (`mu` is monotone).  So the safe direction is to
EXCLUDE poses on the boundary `d = D(psi)`, which is what `eps < 0` in `member()` does; `eps = 0`
is the exact closed criterion (mathematically valid: a chord of length exactly 1 is a chord of
length `>= 1`) and `eps > 0` is a superset and is NOT valid.  The default is `eps = -1e-12`.

Usage:
    python3 search/linecuts.py --grid 4x4                    # the sanity check of LINECUTS.md 4
    python3 search/linecuts.py --measure runs/tl_B40K_measure_pure_nochord.txt --t 4 --top 20
"""
import argparse
import math

import numpy as np


# ====================================================================== the family
def clip_length(nx, ny, c0, t):
    """length of {x : n.x = c0} intersected with [0,t]^2 (0 if it misses)."""
    dx, dy = -ny, nx                       # direction of the line
    px, py = c0 * nx, c0 * ny              # the point of L closest to the origin
    lo, hi = -1e18, 1e18
    for p, d in ((px, dx), (py, dy)):
        if abs(d) < 1e-15:
            if p < -1e-12 or p > t + 1e-12:
                return 0.0
        else:
            a, b = (0.0 - p) / d, (t - p) / d
            if a > b:
                a, b = b, a
            lo, hi = max(lo, a), min(hi, b)
    return max(hi - lo, 0.0)


def cap_of(Lam, tol=1e-9):
    """the exact bound: N < Lam for integer N, i.e. N <= ceil(Lam) - 1."""
    r = round(Lam)
    if abs(Lam - r) < tol:
        return int(r) - 1
    return int(math.floor(Lam))


class LineFamily:
    """arrays `nx, ny, c0, phi, cap` -- one entry per line."""

    def __init__(self, t):
        self.t = float(t)
        self.nx, self.ny, self.c0, self.phi, self.cap, self.lam = [], [], [], [], [], []
        self.tag = []

    def add(self, nx, ny, c0, tag='', maxcap=3):
        n = math.hypot(nx, ny)
        nx, ny, c0 = nx / n, ny / n, c0 / n
        Lam = clip_length(nx, ny, c0, self.t)
        if Lam <= 1.0 + 1e-12:
            return False                    # cap 0: no unit chord fits.  Valid but vacuous here
        k = cap_of(Lam)
        if k > maxcap or k < 0:
            return False
        self.nx.append(nx); self.ny.append(ny); self.c0.append(c0)
        self.phi.append(math.atan2(nx, -ny))          # direction (-ny, nx)
        self.cap.append(float(k)); self.lam.append(Lam); self.tag.append(tag)
        return True

    def finish(self):
        for k in ('nx', 'ny', 'c0', 'phi', 'cap', 'lam'):
            setattr(self, k, np.asarray(getattr(self, k), dtype=float))
        return self

    def __len__(self):
        return len(self.cap)

    def describe(self, i):
        return (f'{self.tag[i]:>14s}  n=({self.nx[i]:+.4f},{self.ny[i]:+.4f}) c={self.c0[i]:+.6f} '
                f'Lam={self.lam[i]:.4f} cap={self.cap[i]:.0f}')


def build(t, pitch=0.05, exact=(1.0, 2.0, 3.0), diag_pitch=0.0, maxcap=3, hv=True):
    """horizontal and vertical lines on a grid of the given pitch (plus the exact heights
    `exact`), and optionally the two diagonal directions."""
    F = LineFamily(t)
    cs = []
    if hv and pitch > 0:
        n = int(round(t / pitch))
        cs = [i * pitch for i in range(n + 1)]
    cs = sorted(set([round(c, 12) for c in cs] + [float(e) for e in exact]))
    for c in cs:
        F.add(0.0, 1.0, c, f'y={c:g}', maxcap)
        F.add(1.0, 0.0, c, f'x={c:g}', maxcap)
    if diag_pitch > 0:
        n = int(round(2 * t / diag_pitch))
        for i in range(n + 1):
            c = i * diag_pitch                      # x + y = c
            F.add(1.0, 1.0, c, f'x+y={c:g}', maxcap)
            d = -t + i * diag_pitch                 # x - y = d
            F.add(1.0, -1.0, d, f'x-y={d:g}', maxcap)
    return F.finish()


# ====================================================================== the predicate
def Dfun(psi):
    c, s = np.abs(np.cos(psi)), np.abs(np.sin(psi))
    return (c + s) / 2.0 - c * s


def chord_length(F, P):
    """(nlines x nposes) exact chord length of pose `P[j]` on line `F[i]` (0 if they miss),
    `notes/chord-lemma.md` Lemma 1 in the line's frame."""
    P = np.asarray(P, dtype=float).reshape(-1, 3)
    d = np.abs(np.outer(F.nx, P[:, 0]) + np.outer(F.ny, P[:, 1]) - F.c0[:, None])
    psi = P[None, :, 2] - F.phi[:, None]
    C, S = np.abs(np.cos(psi)), np.abs(np.sin(psi))
    p = (C + S) / 2.0
    big = 1e300
    a = np.where(C > 1e-14, 1.0 / np.maximum(C, 1e-300), big)
    b = np.where(S > 1e-14, 1.0 / np.maximum(S, 1e-300), big)
    c = np.where(C * S > 1e-14, (p - d) / np.maximum(C * S, 1e-300), big)
    L = np.minimum(np.minimum(a, b), c)
    L = np.where(L >= big / 2, 1.0, L)                 # axis-parallel: the chord is the side
    return np.where(d <= p, np.maximum(L, 0.0), 0.0)


def cap_tau(Lam, tau, tol=1e-9):
    """`N` chords of length >= tau, pairwise disjoint and closed, inside a segment of length `Lam`
    satisfy `N < Lam/tau` (the proof of Theorem L with `1` replaced by `tau`)."""
    q = Lam / tau
    r = round(q)
    return (int(r) - 1) if abs(q - r) < tol else int(math.floor(q))


def member(F, P, eps=-1e-12):
    """boolean (nlines x nposes): does pose `P[j]` have a chord of length >= 1 on line `F[i]`?

    `P` is an (n, 3) array of `(cx, cy, theta)`.  `eps < 0` excludes poses exactly on the
    boundary of the set (a strict subset of it, hence unconditionally valid); `eps = 0` is the
    exact closed criterion; `eps > 0` is a SUPERSET and is not valid."""
    P = np.asarray(P, dtype=float).reshape(-1, 3)
    if not len(F) or not len(P):
        return np.zeros((len(F), len(P)), dtype=bool)
    d = np.abs(np.outer(F.nx, P[:, 0]) + np.outer(F.ny, P[:, 1]) - F.c0[:, None])
    D = Dfun(P[None, :, 2] - F.phi[:, None])
    return d <= D + eps


def mass(F, P, mu, eps=-1e-12):
    """mu(B_L) for every line of the family"""
    M = member(F, P, eps)
    return M @ np.asarray(mu, dtype=float)


# ====================================================================== CLI
def read_measure(path):
    P, mu = [], []
    t = None
    for line in open(path):
        s = line.strip()
        if s.startswith('#'):
            for tok in s.replace(',', ' ').split():
                if tok.startswith('t='):
                    try:
                        t = float(tok[2:])
                    except ValueError:
                        pass
            continue
        q = s.split()
        if len(q) >= 5 and q[0] == 'pose':
            P.append((float(q[1]), float(q[2]), math.radians(float(q[3]))))
            mu.append(float(q[4]))
    return t, np.array(P, dtype=float).reshape(-1, 3), np.array(mu, dtype=float)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--measure', action='append', default=[])
    ap.add_argument('--grid', default=None, help='"4x4" or "4x4-4" (grid minus the four centres)')
    ap.add_argument('--t', type=float, default=4.0)
    ap.add_argument('--pitch', type=float, default=0.05)
    ap.add_argument('--diag-pitch', type=float, default=0.0)
    ap.add_argument('--eps', type=float, default=-1e-12)
    ap.add_argument('--top', type=int, default=15)
    ap.add_argument('--check', action='store_true', help='the y = 0.5 / 1 / 1.5 / 2 sanity table')
    ap.add_argument('--tau', default=None,
                    help='comma-separated thresholds: also screen mu({chord >= tau}) <= cap_tau')
    a = ap.parse_args()

    if a.check or a.grid:
        t = a.t
        P = [(0.5 + i, 0.5 + j, 0.0) for i in range(int(t)) for j in range(int(t))]
        mu = np.ones(len(P))
        if a.grid and a.grid.endswith('-4'):
            keep = [k for k, (x, y, _) in enumerate(P)
                    if not (1.0 < x < 3.0 and 1.0 < y < 3.0)]
            P = [P[k] for k in keep]; mu = np.ones(len(P))
        P = np.array(P, dtype=float)
        print(f'# {len(P)} axis-parallel unit squares at the {int(t)}x{int(t)} grid centres, mass 1 each')
        F = LineFamily(t)
        for c in (0.5, 0.9, 1.0, 1.0000001, 1.5, 2.0, 2.5, 3.0):
            F.add(0.0, 1.0, c, f'y={c:g}', maxcap=99)
        F.finish()
        for e, name in ((0.0, 'closed  d<=D'), (-1e-12, 'strict  d<D'), (1e-9, 'slack   d<=D+1e-9')):
            m = mass(F, P, mu, e)
            print(f'  {name:22s}  ' + '  '.join(f'{F.tag[i]}:{m[i]:.0f}' for i in range(len(F))))
        print('  cap for every one of these lines: Lam = 4 -> cap 3')
        return

    F = build(a.t, a.pitch, diag_pitch=a.diag_pitch)
    print(f'# {len(F)} lines, t = {a.t}, pitch {a.pitch}, diag pitch {a.diag_pitch}')
    for path in a.measure:
        t, P, mu = read_measure(path)
        t = t or a.t
        m = mass(F, P, mu, a.eps)
        v = m - F.cap
        o = np.argsort(-v)
        print(f'\n{path}\n  mass {mu.sum():.6f}, {len(P)} poses; '
              f'{int((v > 1e-9).sum())} of {len(F)} lines violated, total excess {v[v>0].sum():.6f}, '
              f'max excess {v[o[0]]:+.6f}')
        for i in o[:a.top]:
            print(f'    {F.describe(i)}  mu(B_L) = {m[i]:.6f}   excess {v[i]:+.6f}')
        if a.tau:
            Lch = chord_length(F, P)
            for tau in [float(s) for s in a.tau.split(',')]:
                mm = (Lch >= tau - 1e-12).astype(float) @ mu
                cp = np.array([cap_tau(F.lam[i], tau) for i in range(len(F))], dtype=float)
                vv = mm - cp
                i = int(np.argmax(vv))
                print(f'    tau={tau:.4f}: {int((vv > 1e-7).sum())} of {len(F)} violated; worst '
                      f'{F.tag[i]} (Lam={F.lam[i]:.3f} cap={cp[i]:.0f}) mu={mm[i]:.6f} '
                      f'excess {vv[i]:+.6f}')


if __name__ == '__main__':
    main()
