#!/usr/bin/env python3
"""jlevy/squares' packing colours, ported to Python, plus a renderer for our regularized records.

Colour = angle class (0 deg -> family 0, 45 deg -> family 1, the rest from family 2 by descending class size);
shade = number of full-side contacts, walls included for axis-parallel squares (4 -> darkest, 0 -> lightest).
Ported by reading jlevy/squares packing/src/sqpack/render/{color,style}.py (not executed).  Check: every fill of
jlevy's atlas/known-best/rendering/n-029.svg is in PAL.  Simplification: angle classes are tolerance-seeded only
(jlevy also merges classes proven aligned by full-side contacts; within 1e-6 rad that changes nothing here).

  python3 search/regularize/jstyle.py N [out.png]     draw record N from lists/certs (regularized view)

Ported code: Copyright (c) 2026 Joshua Levy, MIT License:
  Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
  documentation files (the "Software"), to deal in the Software without restriction, including without limitation
  the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to the following conditions:
  The above copyright notice and this permission notice shall be included in all copies or substantial portions of
  the Software.
  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
  WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
  COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
  OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import math
from colorsys import rgb_to_hls

HUES = ("#1faa8e","#c3c45f","#aa5585","#166eac","#b3543b","#c9a13a","#23b4e8","#158655","#a7539d","#75951c",
        "#c8691e","#1990a2","#714fad","#67b45c","#e26e82","#8286e9","#ce871b","#147e7c","#4571c9","#a86cc6")
RAMPS = {0: (0.50, 0.70, 0.080), 1: (0.68, 0.862, 0.150)}
QT = math.pi / 2

def s2l(c): return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
def l2s(c): return 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055

def hex_oklch(h):
    r, g, b = (s2l(int(h[i:i+2], 16) / 255) for i in (1, 3, 5))
    L = (0.4122214708*r + 0.5363325363*g + 0.0514459929*b) ** (1/3)
    M = (0.2119034982*r + 0.6806995451*g + 0.1073969566*b) ** (1/3)
    S = (0.0883024619*r + 0.2817188376*g + 0.6299787005*b) ** (1/3)
    a = 1.9779984951*L - 2.4285922050*M + 0.4505937099*S
    bb = 0.0259040371*L + 0.7827717662*M - 0.8086757660*S
    return (0.2104542553*L + 0.7936177850*M - 0.0040720468*S, math.hypot(a, bb), math.degrees(math.atan2(bb, a)) % 360)

def oklch_lin(l, c, h):
    a, b = c*math.cos(math.radians(h)), c*math.sin(math.radians(h))
    L = (l + 0.3963377774*a + 0.2158037573*b) ** 3
    M = (l - 0.1055613458*a - 0.0638541728*b) ** 3
    S = (l - 0.0894841775*a - 1.2914855480*b) ** 3
    return (4.0767416621*L - 3.3077115913*M + 0.2309699292*S,
            -1.2684380046*L + 2.6097574011*M - 0.3413193965*S,
            -0.0041960863*L - 0.7034186147*M + 1.7076147010*S)

def oklch_hex(l, c, h):
    return '#' + ''.join(f'{min(255, max(0, round(l2s(min(1.0, max(0.0, v))) * 255))):02x}' for v in oklch_lin(l, c, h))

def maxc(l, h):
    lo, hi = 0.0, 0.45
    for _ in range(48):
        m = (lo + hi) / 2
        if all(-1e-4 <= v <= 1 + 1e-4 for v in oklch_lin(l, m, h)): lo = m
        else: hi = m
    return lo

def spread(ls): return [x + 0.015 * max(0, i - 1) for i, x in enumerate(ls)]

def hex_hsl(h):
    r, g, b = (int(h[i:i+2], 16) / 255 for i in (1, 3, 5))
    mx, mn = max(r, g, b), min(r, g, b); L = (mx + mn) / 2; d = mx - mn
    if d == 0: return 0.0, 0.0, L
    s = d / (1 - abs(2*L - 1))
    if mx == r: hs = (g - b) / d; hs += 6 if hs < 0 else 0
    elif mx == g: hs = (b - r) / d + 2
    else: hs = (r - g) / d + 4
    return hs, s, L

def hsl_hex(hs, s, L):
    sec = int(hs) % 6; C = (1 - abs(2*L - 1)) * s; X = C * (1 - abs((hs % 2) - 1))
    r, g, b = [(C, X, 0), (X, C, 0), (0, C, X), (0, X, C), (X, 0, C), (C, 0, X)][sec]
    m = L - C / 2
    ch = lambda v: int(math.floor(min(1, max(0, v)) * 255 + 0.5))
    return f'#{ch(r+m):02x}{ch(g+m):02x}{ch(b+m):02x}'

def palette(k=5, span=0.2):
    fams = []
    for i, hx in enumerate(HUES):
        hs, s, L = hex_hsl(hx)
        if i in RAMPS:
            dk, lt, cr = RAMPS[i]; hue = hex_oklch(hsl_hex(hs, s, L))[2]
            step = (lt - dk) / (k - 1)
            fams.append([oklch_hex(l, min(cr, 0.95 * maxc(l, hue)), hue) for l in spread([dk + step*j for j in range(k)])])
            continue
        cen = 0.52 + (L - 0.52) * 0.85; step = span / (k - 1); lo = cen - span / 4
        ls = spread([lo + step*j for j in range(k)])
        sat = min(max(s, 0.5), 0.85)
        ss = [min(1, max(0.05, sat * (1 + 0.12 * j / (k - 1)))) for j in range(k)]
        fams.append([hsl_hex(hs, a, b) for a, b in zip(ss, ls)])
    return fams

PAL = palette()

def corners(x, y, t):
    c, s = math.cos(t), math.sin(t)
    return [(x + c*u - s*v, y + s*u + c*v) for u, v in ((-.5, -.5), (.5, -.5), (.5, .5), (-.5, .5))]

def odist(a, b):
    d = abs(a - b) % QT
    return min(d, QT - d)

def colours(sq, S, atol=1e-6, ctol=2e-6):
    """sq: list of (x, y, theta_radians). Returns (fills, hue_idx, contact_counts, class_of)."""
    n = len(sq); ori = [t % QT for _, _, t in sq]
    E = [[(P[i], P[(i+1) % 4]) for i in range(4)] for P in (corners(*q) for q in sq)]
    cand = [[[] for _ in range(4)] for _ in range(n)]
    for i in range(n):
        if odist(ori[i], 0) <= atol:
            for e, ((x0, y0), (x1, y1)) in enumerate(E[i]):
                for nm, ax, bd in (('wl', 0, 0), ('wr', 0, S), ('wb', 1, 0), ('wt', 1, S)):
                    if max(abs((x0, y0)[ax] - bd), abs((x1, y1)[ax] - bd)) <= ctol: cand[i][e].append(nm)
    for i in range(n):
        for j in range(i+1, n):
            if odist(ori[i], ori[j]) > atol: continue
            for a, ea in enumerate(E[i]):
                for b, eb in enumerate(E[j]):
                    (p0, p1), (q0, q1) = ea, eb
                    d1 = max(abs(p0[0]-q0[0]), abs(p0[1]-q0[1]), abs(p1[0]-q1[0]), abs(p1[1]-q1[1]))
                    d2 = max(abs(p0[0]-q1[0]), abs(p0[1]-q1[1]), abs(p1[0]-q0[0]), abs(p1[1]-q0[1]))
                    if min(d1, d2) <= ctol: cand[i][a].append(j); cand[j][b].append(i)
    cnt = [sum(1 for e in c if e) for c in cand]
    # angle classes, tolerance seeded (full-side merges only join classes already within atol here)
    reps, cls = [], []
    for i in range(n):
        m = next((k for k, r in enumerate(reps) if odist(ori[i], r) <= atol), None)
        if m is None: reps.append(ori[i]); cls.append([i])
        else: cls[m].append(i)
    hue = [0] * len(cls); unp = []
    for k, r in enumerate(reps):
        if odist(r, 0) <= atol: hue[k] = 0
        elif odist(r, QT/2) <= atol: hue[k] = 1
        else: unp.append(k)
    unp.sort(key=lambda k: (-len(cls[k]), min(cls[k])))
    for p, k in enumerate(unp): hue[k] = 2 + p % 18
    cof = {i: k for k, m in enumerate(cls) for i in m}
    fills = [PAL[hue[cof[i]]][4 - cnt[i]] for i in range(n)]
    return fills, [hue[cof[i]] for i in range(n)], cnt, cof


CERTS = __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'lists', 'certs')

def load(n):
    """Regularized record n from lists/certs: (S, [(x, y, theta_radians)]), floats (certificate rationals)."""
    import gzip
    from fractions import Fraction as F
    with gzip.open(f'{CERTS}/n-{n:03d}.cert.gz', 'rt') as f:
        L = [l for l in f.read().split('\n') if l.strip() and not l.startswith('#')]
    S = float(F(L[0].split()[1])); sq = []
    for l in L[1:]:
        x, y, t = map(F, l.split())
        c, s = (1 - t*t) / (1 + t*t), 2*t / (1 + t*t)
        sq.append((float(x), float(y), math.atan2(float(s), float(c))))
    return S, sq

def draw(ax, S, sq, lw=1.0, labels=None, label_fs=9, fills=None, alpha=1.0):
    from matplotlib.patches import Polygon, Rectangle
    f, hue, cnt, cof = colours(sq, S)
    fills = fills or f
    for i, q in enumerate(sq):
        ax.add_patch(Polygon(corners(*q), closed=True, fc=fills[i], ec='#000000', lw=lw, joinstyle='round', alpha=alpha))
        if labels:
            ax.text(q[0], q[1], labels[i], ha='center', va='center', fontsize=label_fs, color='#17202a', fontweight='bold')
    ax.add_patch(Rectangle((0, 0), S, S, fill=False, ec='#000000', lw=lw * 1.6))
    ax.set_xlim(-0.04*S, 1.04*S); ax.set_ylim(-0.04*S, 1.04*S); ax.set_aspect('equal'); ax.axis('off')
    return f, hue, cnt, cof

if __name__ == '__main__':
    import sys
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    n = int(sys.argv[1]); out = sys.argv[2] if len(sys.argv) > 2 else f'n{n}.png'
    S, sq = load(n)
    fig, ax = plt.subplots(figsize=(8, 8), dpi=200)
    draw(ax, S, sq, lw=1.3)
    fig.savefig(out, bbox_inches='tight', pad_inches=0.15, facecolor='white')
    print(out)
