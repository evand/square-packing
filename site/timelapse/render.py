#!/usr/bin/env python3
"""Render the records time-lapse: out/history.json + out/views/ -> out/timelapse.mp4 (+ final still, 1080p copy).

Layout: jlevy's Atlas triangle.  Row k = ceil(sqrt n) holds n = (k-1)^2 + 1 .. k^2, right-aligned (n = k^2, the full
grid, on the right edge; column c = k^2 - n counts from the right), with a small gap before each row's grid suffix
(the n whose record is the k x k grid minus squares, as of the last frame), as in jlevy's poster.
Colours: jlevy's (search/regularize/jstyle.py, a port of jlevy/squares packing/src/sqpack/render/color.py): hue = angle
class (0 deg, 45 deg, then the other classes of that packing by size), shade = number of full-side contacts.

Time: piecewise-linear "regularised" clock (SEGMENTS: decades compressed, Sept-Oct 2026 stretched), then every step
(the states sharing one timestamp) at least MIN_GAP seconds after the previous one.  A changed cell cross-fades and
gets a crimson frame that fades out; the caption lists the latest records.

  render.py [--fps 30] [--still-only] [--frames-dir DIR] [--strip N]     (--strip N: one n's states side by side)
"""
import argparse, json, math, os, subprocess, sys
from datetime import datetime, timezone
import numpy as np
from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
OUT = os.path.join(HERE, 'out')
sys.path.insert(0, os.path.join(ROOT, 'search', 'regularize'))
import jstyle  # noqa: E402

W, H = 3840, 2160
NMAX, ROWS = 324, 18
FONT = '/nix/store/bz36mpwdhrcc0nyyl24zh86hfpc9md6a-dejavu-fonts-2.37/share/fonts/truetype/DejaVuSans.ttf'
FONTB = FONT.replace('DejaVuSans.ttf', 'DejaVuSans-Bold.ttf')
BG = (255, 255, 255)
INK = (23, 32, 42)
MUTED = (110, 116, 124)
ACCENT = (178, 24, 60)          # jlevy's poster marks recent results in crimson
SS = 4                          # tile supersampling

# real time -> video seconds (piecewise linear between knots)
KNOTS = [('1979-01-01', 3.0), ('2024-01-01', 16.0), ('2026-09-01', 36.0), ('2026-10-05', 52.0), ('2026-10-15', 100.0)]
INTRO, OUTRO = 3.0, 6.0
MIN_GAP = 0.27
FADE_IN, HILITE = 0.3, 2.0


def font(sz, bold=False):
    try:
        return ImageFont.truetype(FONTB if bold else FONT, sz)
    except OSError:
        return ImageFont.load_default()


def ts(s):
    return datetime.fromisoformat(s.rstrip('Z')).replace(tzinfo=timezone.utc).timestamp()


def clock(t):
    K = [(ts(d + 'T00:00'), v) for d, v in KNOTS]
    x = ts(t)
    if x <= K[0][0]:
        return K[0][1]
    for (a, va), (b, vb) in zip(K, K[1:]):
        if x <= b:
            return va + (vb - va) * (x - a) / (b - a)
    return K[-1][1]


# ------------------------------------------------------------------------------------------------ layout

def row_of(n):
    return math.ceil(math.sqrt(n) - 1e-12)


def layout(final_side):
    """cell top-left corners: {n: (x, y)}, cell size, label font size"""
    margin_r, margin_b, top = 70, 40, 34
    gapf = 0.45
    pitch = (W - 2 * margin_r) / (35 + gapf)
    cell = int(pitch * 0.86)
    rpitch = (H - top - margin_b) / ROWS
    pos = {}
    for k in range(1, ROWS + 1):
        ns = list(range((k - 1) ** 2 + 1, k * k + 1))
        g0 = next((m for m in ns if all(final_side[q] >= k - 1e-12 for q in ns if q >= m)), k * k)
        for n in ns:
            c = k * k - n
            gap = gapf * pitch if n < g0 else 0.0
            x = W - margin_r - cell - c * pitch - gap
            y = top + (k - 1) * rpitch
            pos[n] = (int(round(x)), int(round(y)))
    return pos, cell, int(rpitch - cell - 4)


# ------------------------------------------------------------------------------------------------ tiles

def tile(view, cell):
    """RGB tile of one view at cell x cell, jlevy colours"""
    S = view['s']
    sq = [(x, y, math.radians(t)) for x, y, t in view['squares']]
    fills, _, _, _ = jstyle.colours(sq, S)
    N = cell * SS
    im = Image.new('RGB', (N, N), BG)
    d = ImageDraw.Draw(im)
    u = (N - 2 * SS) / S                     # px per unit (1 px border room)
    o = SS
    lw = max(1, int(round(SS * min(1.4, max(0.55, 0.045 * u / SS)))))
    for (x, y, t), f in zip(sq, fills):
        P = [(o + px * u, N - o - py * u) for px, py in jstyle.corners(x, y, t)]
        d.polygon(P, fill=f)
    for (x, y, t), f in zip(sq, fills):
        P = [(o + px * u, N - o - py * u) for px, py in jstyle.corners(x, y, t)]
        d.line(P + [P[0]], fill=(0, 0, 0), width=lw, joint='curve')
    bw = max(SS, int(round(SS * 1.5)))
    d.rectangle([o - bw // 2, o - bw // 2, N - o + bw // 2 - 1, N - o + bw // 2 - 1], outline=(0, 0, 0), width=bw)
    return im.resize((cell, cell), Image.LANCZOS)


# ------------------------------------------------------------------------------------------------ timeline

def steps(H, V):
    """[(t, [(n, state index)])] in time order, plus the initial state index per n"""
    by_t = {}
    for n, sts in H.items():
        for i, st in enumerate(sts):
            if st['t'] is not None:
                by_t.setdefault(st['t'], []).append((int(n), i))
    return sorted(by_t.items())


def schedule(S):
    out, last = [], INTRO - MIN_GAP
    for t, ch in S:
        v = max(INTRO + clock(t) - KNOTS[0][1], last + MIN_GAP)
        out.append((v, t, ch)); last = v
    return out, last + OUTRO


def fmt_date(t):
    d = datetime.fromisoformat(t.rstrip('Z'))
    return d


def short(s, k=34):
    return s if len(s) <= k else s[:k - 1] + '…'


# ------------------------------------------------------------------------------------------------ main render

class Renderer:
    def __init__(self):
        D = json.load(open(os.path.join(OUT, 'history.json')))
        self.H = {int(n): v for n, v in D['hist'].items()}
        self.V = {int(n): v for n, v in json.load(open(os.path.join(OUT, 'views.json'))).items()}
        final = {n: v[-1]['s'] for n, v in self.H.items()}
        self.pos, self.cell, self.lab = layout(final)
        self.tiles = {}
        self.fn = font(max(12, int(self.lab * 0.95)))
        self.grid_tile = {}

    def view(self, n, i):
        k = self.V[n][i]
        if k is None:                        # no drawing: grid of the right size as a placeholder
            return None
        if k not in self.tiles:
            v = json.load(open(os.path.join(OUT, 'views', k + '.json')))
            self.tiles[k] = tile(v, self.cell)
        return self.tiles[k]

    def base(self, cur):
        im = Image.new('RGB', (W, H), BG)
        d = ImageDraw.Draw(im)
        for n, (x, y) in self.pos.items():
            t = self.view(n, cur[n])
            if t is not None:
                im.paste(t, (x, y))
            d.text((x + self.cell / 2, y + self.cell + 2), str(n), fill=MUTED, font=self.fn, anchor='mt')
        return im

    def caption(self, im, t, recent, frac):
        d = ImageDraw.Draw(im)
        x0, y0 = 70, 60
        d.text((x0, y0), 'Best known packings of n unit squares in a square', fill=INK, font=font(60, True))
        d.text((x0, y0 + 84), 'n = 1 … 324, record by record, 1979 – October 2026', fill=MUTED, font=font(38))
        if t is None:
            big = 'before 1979'
        else:
            dt = fmt_date(t)
            big = dt.strftime('%Y') if dt.year < 2024 else dt.strftime('%B %Y') if dt.year < 2026 or dt.month < 9 \
                else dt.strftime('%-d %B %Y')
        d.text((x0, y0 + 170), big, fill=INK, font=font(120, True))
        # timeline bar
        by = y0 + 340
        bw = 1500
        d.line([(x0, by), (x0 + bw, by)], fill=(200, 200, 200), width=6)
        span = KNOTS[-1][1] - KNOTS[0][1]
        for dd, v in KNOTS:
            xx = x0 + bw * (v - KNOTS[0][1]) / span
            d.line([(xx, by - 14), (xx, by + 14)], fill=(160, 160, 160), width=3)
            lab = dd[:4] if dd[5:] == '01-01' else dd[:7] if dd[8:] == '01' else dd
            d.text((xx, by + 22), lab, fill=MUTED, font=font(26), anchor='mt')
        d.ellipse([x0 + bw * frac - 14, by - 14, x0 + bw * frac + 14, by + 14], fill=ACCENT)
        # legend
        ly = by + 90
        sw = 34
        for k, (lab, fam) in enumerate((('0°', 0), ('45°', 1), ('other tilts', 2), ('', 3), ('', 4))):
            d.rectangle([x0 + k * (sw + 6), ly, x0 + k * (sw + 6) + sw, ly + sw], fill=jstyle.PAL[fam][1], outline=(0, 0, 0), width=2)
        d.text((x0 + 5 * (sw + 6) + 14, ly + sw / 2), 'colour = tilt angle (0°, 45°, others)', fill=INK, font=font(30), anchor='lm')
        for k in range(5):
            d.rectangle([x0 + k * (sw + 6), ly + 50, x0 + k * (sw + 6) + sw, ly + 50 + sw], fill=jstyle.PAL[0][4 - k],
                        outline=(0, 0, 0), width=2)
        d.text((x0 + 5 * (sw + 6) + 14, ly + 50 + sw / 2), 'shade = full-side contacts (0 … 4)', fill=INK, font=font(30), anchor='lm')
        # recent records
        ry = ly + 140
        d.text((x0, ry), 'Latest records', fill=INK, font=font(36, True))
        ry += 56
        for line, age in recent[:9]:
            a = max(0.0, min(1.0, 1 - age / 6.0))
            col = tuple(int(ACCENT[c] * a + MUTED[c] * (1 - a)) for c in range(3))
            d.text((x0, ry), line, fill=col, font=font(30))
            ry += 42
        ry += 20
        d.text((x0, ry), 'Data: Friedman / Ellsworth catalogue (kingbird), jlevy/squares register and issues, '
               'regularised views', fill=MUTED, font=font(24))
        d.text((x0, ry + 34), 'Colours after jlevy\'s Atlas.  evand/square-packing site/timelapse', fill=MUTED, font=font(24))

    def run(self, fps, frames_dir=None, still_only=False, out_mp4=None):
        S = steps(self.H, self.V)
        sched, total = schedule(S)
        cur = {n: 0 for n in self.H}
        nfr = int(math.ceil(total * fps))
        canvas = self.base(cur)
        prev_tiles = {}
        changed_at = {}                      # n -> (video time, old tile)
        recent = []
        si = 0
        t_label = None
        enc = None
        if not still_only:
            out_mp4 = out_mp4 or os.path.join(OUT, 'timelapse.mp4')
            enc = subprocess.Popen(['ffmpeg', '-y', '-loglevel', 'error', '-f', 'rawvideo', '-pix_fmt', 'rgb24', '-s', f'{W}x{H}',
                                    '-r', str(fps), '-i', '-', '-c:v', 'libx264', '-preset', 'medium', '-crf', '20', '-threads', '3',
                                    '-pix_fmt', 'yuv420p', '-movflags', '+faststart', out_mp4], stdin=subprocess.PIPE)
        frames = range(nfr) if not still_only else [nfr - 1]
        if still_only:
            for v, t, ch in sched:
                for n, i in ch:
                    cur[n] = i
            canvas = self.base(cur)
            t_label = sched[-1][1]
            fr = canvas.copy()
            self.caption(fr, t_label, self.recent_lines(sched, len(sched) - 1, sched[-1][0] + 10), 1.0)
            fr.save(os.path.join(OUT, 'final.png'))
            return
        for f in frames:
            T = f / fps
            while si < len(sched) and sched[si][0] <= T:
                v, t, ch = sched[si]
                for n, i in ch:
                    old = self.view(n, cur[n])
                    cur[n] = i
                    new = self.view(n, i)
                    if new is not None:
                        canvas.paste(new, self.pos[n])
                    changed_at[n] = (v, old)
                t_label = t
                si += 1
            fr = canvas.copy()
            d = ImageDraw.Draw(fr)
            for n, (v, old) in list(changed_at.items()):
                age = T - v
                if age > HILITE:
                    del changed_at[n]; continue
                x, y = self.pos[n]
                if age < FADE_IN and old is not None:
                    new = self.view(n, cur[n])
                    fr.paste(Image.blend(old, new, age / FADE_IN), (x, y))
                a = 1 - age / HILITE
                w = max(2, int(round(9 * a)))
                col = tuple(int(ACCENT[c] * a + 255 * (1 - a)) for c in range(3))
                d.rectangle([x - 6, y - 6, x + self.cell + 5, y + self.cell + 5], outline=col, width=w)
            frac = (T - INTRO) / (KNOTS[-1][1] - KNOTS[0][1])
            frac = max(0.0, min(1.0, (min(T, sched[-1][0]) - INTRO) / max(1e-9, sched[-1][0] - INTRO)))
            self.caption(fr, t_label, self.recent_lines(sched, si - 1, T), self.bar_pos(t_label))
            enc.stdin.write(fr.tobytes())
            if f % 300 == 0:
                print(f'frame {f}/{nfr}  t={t_label}  tiles={len(self.tiles)}', flush=True)
        enc.stdin.close(); enc.wait()
        canvas.save(os.path.join(OUT, 'final-canvas.png'))
        fr.save(os.path.join(OUT, 'final.png'))

    def bar_pos(self, t):
        if t is None:
            return 0.0
        return (clock(t) - KNOTS[0][1]) / (KNOTS[-1][1] - KNOTS[0][1])

    def recent_lines(self, sched, upto, T):
        out = []
        for j in range(upto, -1, -1):
            v, t, ch = sched[j]
            own = [(n, i) for n, i in ch if self.H[n][i]['src'] != 'derived']
            for n, i in sorted(own):
                st = self.H[n][i]
                prev = self.H[n][i - 1]['s'] if i else None
                d = datetime.fromisoformat(t.rstrip('Z'))
                ds = d.strftime('%Y') if d.year < 2024 else d.strftime('%Y-%m') if d.year < 2026 or d.month < 9 else d.strftime('%m-%d')
                who = short(st['who'] or '?', 30)
                via = st['via'] if st['via'].startswith('#') else ''
                out.append((f'{ds}  s({n}) = {st["s"]:.6f}  {who} {via}'.rstrip(), T - v))
            if len(out) >= 9:
                break
        return out


def strip(n, cell=360):
    """debug: every state of n, raw source (top) and view (bottom)"""
    sys.path.insert(0, HERE)
    import regview
    H = json.load(open(os.path.join(OUT, 'history.json')))['hist'][str(n)]
    V = json.load(open(os.path.join(OUT, 'views.json')))[str(n)]
    im = Image.new('RGB', (len(H) * (cell + 20), 2 * cell + 80), BG)
    d = ImageDraw.Draw(im)
    for i, (st, k) in enumerate(zip(H, V)):
        S, A = regview.load(st['draw'], n)
        im.paste(tile(dict(s=S, squares=A.tolist()), cell), (i * (cell + 20), 0))
        if k:
            im.paste(tile(json.load(open(os.path.join(OUT, 'views', k + '.json'))), cell), (i * (cell + 20), cell + 10))
        d.text((i * (cell + 20), 2 * cell + 20), f'{(st["t"] or "-")[:10]} {st["s"]:.6f} {short(st["who"], 24)}', fill=INK, font=font(16))
    p = os.path.join(OUT, f'strip-{n}.png'); im.save(p); print(p)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fps', type=int, default=30)
    ap.add_argument('--still-only', action='store_true')
    ap.add_argument('--strip', type=int, nargs='*')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    if a.strip:
        for n in a.strip:
            strip(n)
    else:
        Renderer().run(a.fps, still_only=a.still_only, out_mp4=a.out)
