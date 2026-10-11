import numpy as np, html
from rectpack import half

def corners(x, y, t):
    c, s = np.cos(t), np.sin(t)
    d = [(-.5,-.5),(.5,-.5),(.5,.5),(-.5,.5)]
    return [(x + c*u - s*v, y + s*u + c*v) for u, v in d]

def tiltdeg(t):
    a = np.degrees(t) % 90
    return min(a, 90-a)

def panel(squares, box, extra='', title='', caption='', scale=58, pad=0.25):
    """squares: list of (x,y,t,role); box: (W,H) container drawn from origin"""
    W, H = box
    xs = [p[0] for sq in squares for p in corners(*sq[:3])] + [0, W]
    ys = [p[1] for sq in squares for p in corners(*sq[:3])] + [0, H]
    x0, x1, y0, y1 = min(xs)-pad, max(xs)+pad, min(ys)-pad, max(ys)+pad
    vw, vh = (x1-x0)*scale, (y1-y0)*scale
    T = lambda p: f"{(p[0]-x0)*scale:.2f},{(y1-p[1])*scale:.2f}"
    out = [f'<svg viewBox="0 0 {vw:.1f} {vh:.1f}" width="{vw:.0f}" role="img" aria-label="{html.escape(title)}">']
    out.append(f'<polygon class="box" points="{T((0,0))} {T((W,0))} {T((W,H))} {T((0,H))}"/>')
    out.append(extra.replace('@T', '') if not callable(extra) else extra(T))
    for (x, y, t, role) in squares:
        pts = ' '.join(T(p) for p in corners(x, y, t))
        out.append(f'<polygon class="sq {role}" points="{pts}"/>')
        d = tiltdeg(t)
        if d > 1.5 and role != 'ghost':
            cx, cy = T((x, y)).split(',')
            out.append(f'<text class="ang" x="{cx}" y="{float(cy)+4:.1f}" text-anchor="middle">{d:.0f}°</text>')
    out.append('</svg>')
    return f'<figure><div class="fig">{"".join(out)}</div><figcaption><b>{title}</b>{caption}</figcaption></figure>'

def load(f, n, roles=None):
    v = np.load(f); x, y, t = v[:n], v[n:2*n], v[2*n:]
    roles = roles or ['strip']*n
    return [(x[k], y[k], t[k], roles[k]) for k in range(n)]

def width(sq):
    return max(p[0] for s in sq for p in corners(*s[:3]))

def polyT(T, pts, cls):
    return f'<polygon class="{cls}" points="{" ".join(T(p) for p in pts)}"/>'

secs = []

# A. strip shape classes
H = 3.99
A1 = load('cfg_n4_H3.99.npy', 4); W1 = width(A1)
q = np.load('cfg_n4_H3.99_intr.npy')
A1i = A1 + [(q[0], q[1], q[2], 'ghost')]
# WLOG convention: flat end square flush right at the TOP (flip y)
flip = lambda L, HH: [(x, HH - y, -t, r) for (x, y, t, r) in L]
A1 = flip(A1, H); A1i = flip(A1i, H)
A2 = load('stag_W1.95_H3.99.npy', 4); A3 = load('stag_W1.99_H3.99.npy', 4)
A4 = load('cfg_n4_H3.95.npy', 4)
A5 = load('wall_W1.95_H3.99.npy', 4)
secs.append(('Four-square strips at height just under 4',
 'Minimum-width and near-minimum packings of 4 unit squares in a w × H box, H = 3.99 unless noted. Three shape classes show up so far; this is not a census. The first two put flat end squares in opposite corners; the wall-column class does not.',
 [panel(A1i, (W1, H), title='Tilted-middle class, w ≈ 1.914', caption='Two flat end squares, two tilted middle squares (the Stromquist and Stenlund shape), drawn with the flush-right end at the top. Dashed: the deepest single intruding square, reaching c ≈ 0.89 into the box, nearly to x = 1.'),
  panel(A4, (width(A4), 3.95), title='Tilted-middle class, H = 3.95, w ≈ 1.939', caption='Same class at a lower height; the middle pair tilts about 24°. Our search is loose: Stromquist fits 1.9 × 3.9475.'),
  panel(A2, (1.95, H), title='Staggered class, w = 1.95', caption='Three flat squares alternating sides, one tilted 18°. Near-axis; slots beside the flush-left squares are almost a unit high.'),
  panel(A3, (1.99, H), title='Staggered class, w = 1.99', caption='Tilt down to 7.5°. As w → 2 the tilt goes to 0 and this becomes two tiling columns: the margin-zero core.'),
  panel(A5, (1.95, 3.99), title='Wall-column class, w = 1.95', caption='Three flat squares stacked against the wall, one square tilted about 26° beside them, dipping into the gap. No end square is flush right, so the top-unit block does not apply. Found inside a margin-zero hit (next row, right).')]))

# B. margin-zero leaf hits
def leafp(f, w, title, cap):
    sq = load(f, 12, ['strip']*4 + ['other']*8)
    return panel(sq, (4, 4), extra=lambda T: polyT(T, [(0,0),(w,0),(w,4),(0,4)], 'band'), title=title, caption=cap, scale=40)
secs.append(('The strip leaf at margin zero',
 '12 squares with 4 of them (amber) wholly inside a wall strip of width w < 2 (shaded), minimising the side. The search lands on tilings minus 4 at side 4.000, with touching squares.',
 [leafp('leaf_w1.9_s1.npy', 1.9, 'Slid tile, strip width 1.87', 'Three tiles of column 0 plus a column-1 tile slid partway into column 0’s empty cell. Side 4.002.'),
  leafp('leaf_w1.95_s1.npy', 1.95, 'Two-column stagger, strip width 1.93', 'Two tiles from each of columns 0 and 1, interleaved into each other’s empty cells. Side 4.000.'),
  leafp('leaf_w1.99_s2.npy', 1.99, 'Stagger with tilts, strip width 1.99', 'A looser variant with two tilted squares, still at side 4.000.')]))

# C. 8-square
C1 = load('cfg_n8_H3.99.npy', 8, ['other']*8); C2 = load('cfg_n8_H3.95.npy', 8, ['other']*8)
secs.append(('Eight squares at height just under 4',
 'Best widths found for 8 squares in a w × H box (upper bounds). Flat squares in every corner, a tilted block of four in the middle.',
 [panel(C1, (width(C1), 3.99), title='H = 3.99, w ≈ 2.86', caption='Middle block all at 57°. The top and bottom rows leave a gap of about 0.86 between the corner squares, which the block pokes into.', scale=48),
  panel(C2, (width(C2), 3.95), title='H = 3.95, w ≈ 2.94', caption='Same structure, block at 42–47°.', scale=48)]))

# D. notched region
def notchp(f, c, title, cap):
    W, H = 2.99, 3.99
    sq = load(f, 8, ['other']*8) if f else []
    ex = lambda T: polyT(T, [(0,H-1),(c,H-1),(c,H),(0,H)], 'notch')
    if not sq: sq = []
    return panel(sq, (W, H), extra=ex, title=title, caption=cap, scale=48)
secs.append(('Notched region: the 8-square side of the stitch',
 'The region left for the 8 squares under the generous seam: width H − 1 = 2.99 below the top unit, with the top unit blocked out to a notch of width n (hatched). Proved n ≥ 0.69; the true strip gives n ≈ 0.9.',
 [notchp('notchfit_H3.99_c0.3_W2.99.npy', 0.3, 'n = 0.3: fits', 'The top-left corner square moves to just right of the notch.'),
  notchp('notchfit_H3.99_c0.5_W2.99.npy', 0.5, 'n = 0.5: fits', 'Still fits, with the middle block tilted about 16–20°.'),
  notchp(None, 0.69, 'n = 0.6, 0.69, 0.8, 0.9: no fit found', '300 seeded starts each (float search, not a proof). The threshold is between 0.5 and 0.6.')]))

# E. stitch
H = 3.99; a = W1 - 1
def stitch(T):
    X = lambda x: T((x, 0)).split(",")[0]
    Y = lambda y: T((0, y)).split(",")[1]
    out = polyT(T, [(1,0),(1+a,0),(1+a,H-1),(1,H-1)], 'overlap')
    out += polyT(T, [(0,H-1),(1+a,H-1),(1+a,H),(0,H)], 'border') + polyT(T, [(0,0),(1,0),(1,H-1),(0,H-1)], 'border')
    out += f'<line class="cut" x1="{X(1+a)}" y1="{Y(-0.05)}" x2="{X(1+a)}" y2="{Y(H+0.05)}"/>'
    out += f'<line class="cut" x1="{X(1)}" y1="{Y(-0.05)}" x2="{X(1)}" y2="{Y(H+0.05)}"/>'
    # dimension: a and c along the bottom, below the box
    yb = float(Y(0)) + 14
    out += f'<line class="dimline" x1="{X(1)}" y1="{yb}" x2="{X(1+a)}" y2="{yb}"/>'
    out += f'<text class="dim" x="{(float(X(1))+float(X(1+a)))/2:.1f}" y="{yb+14:.1f}" text-anchor="middle">c = a = {a:.3f}</text>'
    out += f'<text class="dim" x="{X(1)}" y="{float(Y(H))-8:.1f}" text-anchor="end">x = 1 </text>'
    out += f'<text class="dim" x="{X(1+a)}" y="{float(Y(H))-8:.1f}" text-anchor="start"> x = 1 + a = {1+a:.3f}</text>'
    out += f'<text class="dim" x="{float(X(1+a))+8:.1f}" y="{float(Y(H-0.5))+4:.1f}">top unit: blocked to 1 + a</text>'
    out += f'<text class="dim ov" x="{(float(X(1))+float(X(1+a)))/2:.1f}" y="{Y(1.75)}" text-anchor="middle" transform="rotate(-90 {(float(X(1))+float(X(1+a)))/2:.1f} {Y(1.75)})">overlap allowed</text>'
    return out
secs.append(('The stitched picture',
 f'Evan’s optimistic seam on the tilted-middle strip from the first row, with the flush-right end square at the top. Values: H = 3.99, a = {a:.3f} (strip width {1+a:.3f}), c = a, so the lower seam sits at x = 1 and the notch width is n = a. Hatched: the border B the 8 squares must avoid. Red tint: the strip box below the top unit, where the 8 squares may overlap the strip’s own squares.',
 [panel(A1, (H, H), extra=stitch, title='Strip plus border, H = 3.99', caption='Everything right of the hatching is available to the 8 squares, including the red overlap zone. That is the notched region of the previous row, placed at x = 1. The 8 squares never see the strip squares themselves.', scale=60, pad=0.45)]))

body = []
for (h, intro, figs) in secs:
    body.append(f'<section><h2>{h}</h2><p class="intro">{intro}</p><div class="grid">{"".join(figs)}</div></section>')

page = f'''<title>Strip Shape Classes</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;600&family=IBM+Plex+Mono:wght@400;500&family=Fraunces:opsz,wght@9..144,600&display=swap">
<style>
/* layout: one column of labelled figure rows, drafting-sheet feel */
:root {{
  --bg: #f6f7f4; --fg: #1d2421; --muted: #5b6660; --line: #c9cfc9;
  --strip: #d08a2c; --strip-f: #f2d3a6; --other: #3c6e8f; --other-f: #c9dcea;
  --ghost: #b8323a; --band: #e9e2d4; --notch: #9aa39d; --card: #ffffff;
  --display: "Fraunces", Georgia, serif; --body: "IBM Plex Sans", system-ui, sans-serif; --mono: "IBM Plex Mono", ui-monospace, monospace;
}}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{
  --bg: #151a18; --fg: #e3e8e4; --muted: #9aa6a0; --line: #3a443f;
  --strip: #e3a14a; --strip-f: #5a4223; --other: #7fb2d4; --other-f: #26394a;
  --ghost: #ff6b72; --band: #2a2a22; --notch: #6c7570; --card: #1c2220; color-scheme: dark }} }}
:root[data-theme="dark"] {{
  --bg: #151a18; --fg: #e3e8e4; --muted: #9aa6a0; --line: #3a443f;
  --strip: #e3a14a; --strip-f: #5a4223; --other: #7fb2d4; --other-f: #26394a;
  --ghost: #ff6b72; --band: #2a2a22; --notch: #6c7570; --card: #1c2220; color-scheme: dark }}
body {{ background: var(--bg); color: var(--fg); font-family: var(--body); font-size: 15px; line-height: 1.5; }}
main {{ max-width: 1180px; margin: 0 auto; padding-inline: 20px; padding-block: 28px 56px; display: grid; gap: 40px; }}
header h1 {{ font-family: var(--display); font-size: 2rem; margin: 0 0 6px; text-wrap: balance; }}
header p {{ margin: 0; max-width: 70ch; color: var(--muted); }}
.key {{ display: flex; flex-wrap: wrap; gap: 8px 20px; margin-top: 14px; font-size: 13px; color: var(--muted); }}
.key span {{ display: inline-flex; align-items: center; gap: 6px; }}
.sw {{ width: 14px; height: 14px; display: inline-block; border: 1.5px solid; }}
h2 {{ font-family: var(--display); font-size: 1.35rem; margin: 0 0 4px; text-wrap: balance; }}
.intro {{ margin: 0 0 16px; max-width: 75ch; color: var(--muted); }}
.grid {{ display: flex; flex-wrap: wrap; gap: 18px; align-items: flex-start; }}
figure {{ margin: 0; background: var(--card); border: 1px solid var(--line); padding: 12px; display: grid; gap: 8px; max-width: 100%; min-width: 0; flex: 0 1 auto; }}
.fig {{ overflow-x: auto; }}
figure svg {{ display: block; max-width: 100%; height: auto; }}
figcaption {{ font-size: 13px; max-width: 34ch; color: var(--muted); }}
figcaption b {{ display: block; color: var(--fg); font-weight: 600; font-family: var(--mono); font-size: 12.5px; }}
.box {{ fill: none; stroke: var(--fg); stroke-width: 1.5; }}
.sq {{ stroke-width: 1.5; stroke-linejoin: round; }}
.strip {{ fill: var(--strip-f); stroke: var(--strip); }}
.other {{ fill: var(--other-f); stroke: var(--other); }}
.ghost {{ fill: none; stroke: var(--ghost); stroke-dasharray: 5 3; stroke-width: 2; }}
.band {{ fill: var(--band); stroke: none; }}
.notch, .border {{ fill: url(#hatch); stroke: var(--notch); stroke-width: 1; }}
.overlap {{ fill: var(--ghost); fill-opacity: 0.16; stroke: var(--ghost); stroke-width: 1; stroke-dasharray: 3 3; }}
.dimline {{ stroke: var(--fg); stroke-width: 1.2; }}
.dim {{ font-family: var(--mono); font-size: 11.5px; fill: var(--fg); }}
.dim.ov {{ fill: var(--ghost); font-weight: 500; }}
.cut {{ stroke: var(--muted); stroke-dasharray: 4 4; stroke-width: 1.2; }}
.ang {{ font-family: var(--mono); font-size: 11px; fill: var(--fg); }}
</style>
<svg width="0" height="0" style="position:absolute" aria-hidden="true"><defs><pattern id="hatch" width="7" height="7" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><rect width="7" height="7" fill="var(--bg)"/><line x1="0" y1="0" x2="0" y2="7" stroke="var(--notch)" stroke-width="2"/></pattern></defs></svg>
<main>
<header>
<h1>Strip shape classes for s(12)</h1>
<p>Float search results from the 4 + 8 strip-split leaf (notes/strip-split.md, 2026-10-11). Every packing here is a float upper bound, not a certificate. Angles inside squares are tilts from the axes, shown when above 1.5°.</p>
<div class="key"><span><i class="sw" style="background:var(--strip-f);border-color:var(--strip)"></i>strip squares</span><span><i class="sw" style="background:var(--other-f);border-color:var(--other)"></i>other squares</span><span><i class="sw" style="border-color:var(--ghost);border-style:dashed"></i>deepest intruder</span><span><i class="sw" style="background:repeating-linear-gradient(45deg,var(--notch) 0 2px,var(--bg) 2px 5px);border-color:var(--notch)"></i>blocked border or notch</span><span><i class="sw" style="background:var(--band);border-color:var(--band)"></i>strip band</span><span><i class="sw" style="background:color-mix(in srgb, var(--ghost) 16%, transparent);border-color:var(--ghost);border-style:dashed"></i>overlap allowed</span></div>
</header>
{"".join(body)}
</main>'''
open('strip-shapes.html', 'w').write(page)
print(len(page))
