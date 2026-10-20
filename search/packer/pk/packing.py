"""Packing: n unit squares in [0, s + dx] x [0, s] (dx = 0: the square box; dx > 0: rectangle pieces), centres (x, y) lower-left origin, angle in degrees mod 90.

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
    dx: float = 0.0                     # container [0, s + dx] x [0, s] (rectangle pieces, fq --rect-dx)

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
        worst = max(0.0, -C.min(), C[..., 0].max() - self.s - self.dx, C[..., 1].max() - self.s)
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
        """Hash invariant under the symmetries of the box (8 for a square, 4 for a rectangle) and square relabelling
        (rounded coordinates).  dx = 0 hashes exactly as before rectangles existed."""
        s = self.s; X, Y, A = self.sq[:, 0], self.sq[:, 1], self.sq[:, 2]
        best = None
        for sw in ((False, True) if self.dx == 0 else (False,)):
            for fx in (False, True):
                for fy in (False, True):
                    x, y = (Y, X) if sw else (X, Y)
                    a = (90 - A) if sw else A        # swap = reflection: angle -> -angle
                    if fx: x = s + self.dx - x; a = 90 - a
                    if fy: y = s - y; a = 90 - a
                    a = np.round(a / ang) * ang % 90
                    a = np.where(np.abs(a - 90) < ang / 2, 0.0, a)
                    rows = sorted(zip(np.round(x / pos).astype(np.int64), np.round(y / pos).astype(np.int64),
                                      np.round(a / ang).astype(np.int64)))
                    h = hashlib.sha1(json.dumps([self.n, round(s / sd)] + ([['dx', round(self.dx / sd)]] if self.dx else []) + rows, default=int).encode()).hexdigest()[:20]
                    best = h if best is None or h < best else best
        return best

    # ---------- writers ----------
    def text(self, fmt='ours') -> str:
        if self.dx and fmt not in ('ours', 'json'):
            raise ValueError(f'{fmt}: square boxes only (dx = {self.dx})')
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
            return json.dumps(dict(n=self.n, s=float(self.s), squares=self.sq.tolist(), **({'dx': self.dx} if self.dx else {})))
        raise ValueError(fmt)

    def write(self, path, fmt='ours'):
        tmp = f'{path}.tmp'
        with open(tmp, 'w') as f:
            f.write(self.text(fmt))
        os.replace(tmp, path)


# ---------- readers ----------
def _open(p):
    return gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p)


MAX_N = 20000                 # sanity bound on n for strict parsing (untrusted input)
MAX_READ = 64 << 20           # read(): refuse files (decompressed) larger than this
MAX_TOKEN = 4000              # longest number token accepted by the strict parsers (digits of an exact rational)
_RAT = re.compile(r'[+-]?\d+(?:/\d+)?')
_DEC = re.compile(r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?')


def _num(tok: str, exact: bool) -> Fraction:
    """One token of an untrusted file -> Fraction.  exact: integer or p/q only; else also a plain decimal.  Raises."""
    if len(tok) > MAX_TOKEN:
        raise ValueError(f'number token too long ({len(tok)} chars)')
    if _RAT.fullmatch(tok):
        f = Fraction(tok)                                     # ZeroDivisionError on q = 0
    elif not exact and _DEC.fullmatch(tok):
        f = Fraction(tok)
    else:
        raise ValueError(f'not a {"rational" if exact else "number"}: {tok[:40]!r}')
    return f


def parse_sq(txt: str) -> Packing | None:
    """dostashev's formats (github.com/dostashev/square-packing README, 'File formats'), parsed strictly:
      .sqcert  '# sqcert v1' header comments, then `n <n>`, `s <p/q>`, then exactly n lines `x y t` (exact rationals:
               centre (x, y), t = tan(theta/2)).
      .sqpk    '# sqpk v1' ..., `n <n>`, `s <decimal>`, n lines `cx cy theta_deg` (decimals).
    Anything else in the file (extra rows, short rows, junk tokens, missing header) raises ValueError.  The side is kept
    exactly in s_str (the rational for .sqcert)."""
    lines = [l.strip() for l in txt.splitlines()]
    body = [l for l in lines if l and not l.startswith('#')]
    head = next((l for l in lines if l), '')
    exact = head.startswith('# sqcert')
    if not (exact or head.startswith('# sqpk')):
        return None
    if len(body) < 2 or not re.fullmatch(r'n\s+\d{1,6}', body[0]) or not re.fullmatch(r's\s+\S+', body[1]):
        raise ValueError('sqcert/sqpk: expected `n <n>` then `s <side>`')
    n = int(body[0].split()[1]); s_str = body[1].split()[1]
    if not 1 <= n <= MAX_N:
        raise ValueError(f'sqcert/sqpk: n = {n} out of range')
    rows = body[2:]
    if len(rows) != n:
        raise ValueError(f'sqcert/sqpk: {len(rows)} square lines, header says n = {n}')
    S = _num(s_str, exact)
    sq = []
    for r in rows:
        tok = r.split()
        if len(tok) != 3:
            raise ValueError(f'sqcert/sqpk: square line with {len(tok)} fields')
        x, y, a = (_num(v, exact) for v in tok)
        sq.append((float(x), float(y), math.degrees(2 * math.atan(float(a))) if exact else float(a)))
    return Packing(float(S), sq, s_str)


def validate(p: Packing, strict=True) -> Packing:
    """Sanity checks for packings read from untrusted files; raises ValueError.  Always: 1 <= n, every number finite,
    s > 0.  strict: n <= MAX_N, s >= sqrt(n) (area) - 1e-9 and s <= n + 1, every centre within one unit of the box
    (a mis-read format, not a near-feasible packing).  Overlap is not checked here (the store does that: max_pen)."""
    if p.n < 1 or not np.all(np.isfinite(p.sq)) or not math.isfinite(p.s) or p.s <= 0 or not math.isfinite(p.dx):
        raise ValueError(f'invalid packing (n={p.n}, s={p.s})')
    if strict:
        if p.n > MAX_N:
            raise ValueError(f'n = {p.n} > {MAX_N}')
        if not math.sqrt(p.n) - 1e-9 <= p.s <= p.n + 1:
            raise ValueError(f'side {p.s} impossible for n = {p.n}')
        X, Y = p.sq[:, 0], p.sq[:, 1]
        if X.min() < -1 or Y.min() < -1 or X.max() > p.s + p.dx + 1 or Y.max() > p.s + 1:
            raise ValueError('centres far outside the box')
    return p


def parse(txt: str, strict=False) -> Packing | None:
    """Any known format -> Packing, or None if not recognised.  strict (untrusted input): validate() the result, and
    raise on a recognised but malformed file instead of silently truncating."""
    p = _parse(txt, strict)
    return None if p is None else validate(p, strict) if strict else p


def _parse(txt: str, strict=False) -> Packing | None:
    t0 = txt.lstrip()
    if t0.startswith('# sqcert') or t0.startswith('# sqpk'):
        return parse_sq(txt)
    if t0.startswith('{'):
        d = json.loads(txt)
        if not isinstance(d, dict):
            return None
        if strict and len(d.get('squares') or ()) > MAX_N:
            raise ValueError('too many squares')
        if 's_exact' in d and 'squares' in d:                   # SQUISH cert.json: rational x, y, t = tan(theta/2)
            sq = []
            for q in d['squares']:
                x, y, t = (Fraction(str(v)) for v in (q[:3] if isinstance(q, list) else (q['x'], q['y'], q['t'])))
                sq.append((float(x), float(y), math.degrees(2 * math.atan(float(t)))))
            return Packing(side_val(d['s_exact']), sq, str(d['s_exact']))
        if 'squares' in d and 's' in d:
            return Packing(float(d['s']), d['squares'], dx=float(d.get('dx', 0.0)))
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
        if strict and len(w.get('squares') or ()) > MAX_N:
            raise ValueError('too many squares')
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
        if strict and (len(lines) != n + 1 or n > MAX_N or any(len(r) != 3 for r in rows)
                       or any(len(v) > MAX_TOKEN for r in rows for v in r)):
            raise ValueError(f'"n s" format: header n = {n}, {len(lines) - 1} rows (or a row without 3 fields)')
        if len(rows) == n and all(len(r) >= 3 for r in rows):
            if any('/' in v for r in rows for v in r[:3]):       # exact cert: rational x y t, t = tan(theta/2)
                sq = [(float(Fraction(a)), float(Fraction(b)), math.degrees(2 * math.atan(float(Fraction(t)))))
                      for a, b, t, *_ in rows]
            else:
                sq = [(float(a), float(b), float(c)) for a, b, c, *_ in rows]
            return Packing(side_val(first[1]), sq, first[1])
    return None


def read(path, strict=False, max_bytes=MAX_READ) -> Packing | None:
    """parse() a file (.gz transparently); a file larger than max_bytes (after decompression) raises ValueError.
    An explorer archive file removed by pk/archive's dedupe is written back from the store first."""
    if not os.path.exists(path):
        from .archive import ensure
        ensure(path)
    with _open(path) as f:
        txt = f.read(max_bytes + 1)
    if len(txt) > max_bytes:
        raise ValueError(f'{path}: larger than {max_bytes} bytes')
    p = parse(txt, strict)
    if p is not None:
        p.meta.setdefault('path', str(path))
    return p
