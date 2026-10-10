"""Packing: n unit squares in [0, s]^2, centres (x, y) lower-left origin, angle in degrees mod 90.

One type for the whole toolchain; readers for every format we meet, writers for ours and the common third-party ones,
an f64 feasibility check, a grid test (full lines / staggered axis chains), and a D4-canonical hash for dedupe.
"""
from __future__ import annotations
import gzip, hashlib, json, math, os, re
from dataclasses import dataclass, field
from fractions import Fraction
import numpy as np


def side_val(x) -> float:
    x = str(x).strip()
    return float(Fraction(x)) if '/' in x else float(x)


@dataclass
class Packing:
    s: float
    sq: np.ndarray                      # (n, 3): x, y, deg (mod 90)
    s_str: str | None = None            # side exactly as reported by the source (decimal or n/d)
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.sq = np.asarray(self.sq, float).reshape(-1, 3).copy()
        self.sq[:, 2] %= 90.0

    @property
    def n(self) -> int:
        return len(self.sq)

    @property
    def k(self) -> int:
        """trivial bound: ceil(sqrt n)."""
        return math.ceil(math.sqrt(self.n) - 1e-12)

    def tuples(self):
        return [tuple(map(float, r)) for r in self.sq]

    # ---------- checks ----------
    def max_pen(self) -> float:
        """max over walls and pairs of the SAT penetration depth (0 = feasible), f64."""
        P = self.sq; th = np.radians(P[:, 2]); c, si = np.cos(th), np.sin(th)
        off = np.array([[.5, .5], [-.5, .5], [-.5, -.5], [.5, -.5]])
        C = P[:, None, :2] + np.stack([off[:, 0] * c[:, None] - off[:, 1] * si[:, None],
                                       off[:, 0] * si[:, None] + off[:, 1] * c[:, None]], -1)
        worst = max(0.0, -C.min(), C.max() - self.s)
        d = np.linalg.norm(P[:, None, :2] - P[None, :, :2], axis=-1)
        I, J = np.nonzero(np.triu(d < math.sqrt(2) + 1e-9, 1))
        for i, j in zip(I, J):
            pen = math.inf
            for kk in (i, j):
                for ax in ((c[kk], si[kk]), (-si[kk], c[kk])):
                    a, b = C[i] @ ax, C[j] @ ax
                    pen = min(pen, min(a.max() - b.min(), b.max() - a.min()))
                    if pen <= 0: break
                if pen <= 0: break
            worst = max(worst, pen)
        return float(worst)

    def ntilted(self, thr=1.0) -> int:
        a = self.sq[:, 2]
        return int(np.sum(np.minimum(a, 90 - a) > thr))

    def grid_lines(self, k=None) -> int:
        """> 0 if a full line or staggered chain of >= k axis squares forces side >= k (layout.full_lines)."""
        from layout import _straight_lines, axis_chain
        k = self.k if k is None else k
        ax = [(x, y) for x, y, a in self.sq if min(a, 90 - a) < 1.0]
        c = _straight_lines(ax, k)
        if c == 0 and max(axis_chain(ax, 0), axis_chain(ax, 1)) >= k:
            return 1
        return c

    # ---------- identity ----------
    def canon_hash(self, pos=1e-6, ang=1e-4, sd=1e-9) -> str:
        """Hash invariant under the 8 symmetries of the box and square relabelling (rounded coordinates)."""
        s = self.s; X, Y, A = self.sq[:, 0], self.sq[:, 1], self.sq[:, 2]
        best = None
        for sw in (False, True):
            for fx in (False, True):
                for fy in (False, True):
                    x, y = (Y, X) if sw else (X, Y)
                    a = (90 - A) if sw else A        # swap = reflection: angle -> -angle
                    if fx: x = s - x; a = 90 - a
                    if fy: y = s - y; a = 90 - a
                    a = np.round(a / ang) * ang % 90
                    a = np.where(np.abs(a - 90) < ang / 2, 0.0, a)
                    rows = sorted(zip(np.round(x / pos).astype(np.int64), np.round(y / pos).astype(np.int64),
                                      np.round(a / ang).astype(np.int64)))
                    h = hashlib.sha1(json.dumps([self.n, round(s / sd)] + rows, default=int).encode()).hexdigest()[:20]
                    best = h if best is None or h < best else best
        return best

    # ---------- writers ----------
    def text(self, fmt='ours') -> str:
        s_out = self.s_str if self.s_str and fmt == 'ours' else repr(float(self.s))
        if fmt == 'ours':
            return f'{self.n} {s_out}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(a)!r}\n' for x, y, a in self.sq)
        if fmt == 'couzo':
            return f'# n = {self.n}\n# s = {float(self.s)!r}\n# x y theta(rad)\n' + ''.join(
                f'{x:.17e} {y:.17e} {math.radians(a):.17e}\n' for x, y, a in self.sq)
        if fmt == 'ellsworth':
            h = float(self.s) / 2
            return f's: {float(self.s)!r}\n' + ''.join(
                f'Square {i + 1}: x={float(x) - h!r}, y={float(y) - h!r}, deg={float(a)!r}\n'
                for i, (x, y, a) in enumerate(self.tuples()))
        if fmt == 'json':
            return json.dumps(dict(n=self.n, s=float(self.s), squares=self.sq.tolist()))
        raise ValueError(fmt)

    def write(self, path, fmt='ours'):
        tmp = f'{path}.tmp'
        with open(tmp, 'w') as f:
            f.write(self.text(fmt))
        os.replace(tmp, path)


# ---------- readers ----------
def _open(p):
    return gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p)


def parse(txt: str) -> Packing | None:
    """Any known format -> Packing, or None if not recognised."""
    t0 = txt.lstrip()
    if t0.startswith('{'):
        d = json.loads(txt)
        if 's_exact' in d and 'squares' in d:                   # SQUISH cert.json: rational x, y, t = tan(theta/2)
            sq = []
            for q in d['squares']:
                x, y, t = (Fraction(str(v)) for v in (q[:3] if isinstance(q, list) else (q['x'], q['y'], q['t'])))
                sq.append((float(x), float(y), math.degrees(2 * math.atan(float(t)))))
            return Packing(side_val(d['s_exact']), sq, str(d['s_exact']))
        if 'squares' in d and 's' in d:
            return Packing(float(d['s']), d['squares'])
        return None
    if 'Square 1:' in txt and re.search(r'^(?:Final )?s:\s', txt, re.M):
        s_str = re.search(r'^(?:Final )?s:\s*([0-9.eE+-]+)', txt, re.M).group(1); s = float(s_str)
        sq = [(float(a) + s / 2, float(b) + s / 2, float(c)) for a, b, c in
              re.findall(r'Square \d+:\s*x=([^,\s]+),\s*y=([^,\s]+),\s*deg=([^,\s]+)', txt)]
        return Packing(s, sq, s_str)
    if t0.startswith('# n ='):
        s_str = re.search(r'# s =\s*(\S+)', txt).group(1)
        rows = [l.split() for l in txt.splitlines() if l.strip() and not l.startswith('#')]
        return Packing(float(s_str), [(float(x), float(y), math.degrees(float(t))) for x, y, t, *_ in rows], s_str)
    if 'witness:' in txt and 'representation' in txt:
        import yaml
        w = yaml.safe_load(txt)['witness']
        if w['coordinates'].get('origin') != 'lower-left':
            return None
        if w.get('representation') == 'center-basis':           # basis = direction (cos, sin), unnormalised
            sq = [(side_val(q['center'][0]), side_val(q['center'][1]),
                   math.degrees(math.atan2(side_val(q['basis'][1]), side_val(q['basis'][0])))) for q in w['squares']]
            return Packing(side_val(w['side']), sq, str(w['side']))
        if w.get('representation') == 'corners':                # 4 corners in order: centre = mean, angle of edge 0->1
            sq = []
            for q in w['squares']:
                C = [(side_val(a), side_val(b)) for a, b in q['corners']]
                sq.append((sum(c[0] for c in C) / 4, sum(c[1] for c in C) / 4,
                           math.degrees(math.atan2(C[1][1] - C[0][1], C[1][0] - C[0][0]))))
            return Packing(side_val(w['side']), sq, str(w['side']))
        if w.get('representation') != 'center-angle':
            return None
        f = 1.0 if w['coordinates']['angle_unit'] == 'degrees' else 180 / math.pi
        sq = [(float(q['center'][0]), float(q['center'][1]), f * float(q['angle'])) for q in w['squares']]
        return Packing(side_val(w['side']), sq, str(w['side']))
    lines = [l for l in txt.splitlines() if l.strip() and not l.startswith('#')]
    if not lines:
        return None
    first = lines[0].split()
    if len(first) == 2 and first[0].isdigit():
        n = int(first[0]); rows = [l.split() for l in lines[1:n + 1]]
        if len(rows) == n and all(len(r) >= 3 for r in rows):
            if any('/' in v for r in rows for v in r[:3]):       # exact cert: rational x y t, t = tan(theta/2)
                sq = [(float(Fraction(a)), float(Fraction(b)), math.degrees(2 * math.atan(float(Fraction(t)))))
                      for a, b, t, *_ in rows]
            else:
                sq = [(float(a), float(b), float(c)) for a, b, c, *_ in rows]
            return Packing(side_val(first[1]), sq, first[1])
    return None


def read(path) -> Packing | None:
    with _open(path) as f:
        p = parse(f.read())
    if p is not None:
        p.meta.setdefault('path', str(path))
    return p
