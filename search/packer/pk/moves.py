"""Moves with explicit parameters.

A move is f(p: Packing, rng: random.Random, **params) -> (sq ndarray (n, 3), info dict).  Each declares its parameter
space; unspecified parameters are drawn from it, and every value actually used is returned in info['params'], so the
outcome of any proposal can be regressed on its parameters afterwards (noise level, region size, pressure ...).

An *arm* = (move name, fixed params): `ARMS` reproduces the 10-08 explorer's move kinds as presets; new arms are just
new dicts (e.g. kick at fixed sigma bins for a noise sweep).

Param spec: ('log', lo, hi) log-uniform, ('uni', lo, hi), ('int', lo, hi) inclusive, ('choice', [..]), or a constant.
"""
from __future__ import annotations
import math, os, random, subprocess, tempfile
import numpy as np
from .packing import Packing

H = 0.5
HERE = os.path.dirname(os.path.abspath(__file__))
ANNEAL = os.path.join(HERE, '..', 'target-dev', 'release', 'anneal')
if not os.path.exists(ANNEAL):
    ANNEAL = os.path.join(HERE, '..', 'target', 'release', 'anneal')


def draw(spec, rng):
    if not isinstance(spec, tuple):
        return spec
    kind = spec[0]
    if kind == 'log':
        return math.exp(rng.uniform(math.log(spec[1]), math.log(spec[2])))
    if kind == 'uni':
        return rng.uniform(spec[1], spec[2])
    if kind == 'int':
        return rng.randint(spec[1], spec[2])
    if kind == 'choice':
        return rng.choice(spec[1])
    raise ValueError(spec)


MOVES = {}


def move(name, **space):
    """Register f(p, rng, **params) under name with its parameter space."""
    def deco(f):
        def run(p: Packing, rng, **fixed):
            prm = {k: (fixed[k] if k in fixed else draw(v, rng)) for k, v in space.items()}
            prm.update({k: v for k, v in fixed.items() if k not in space})
            sq, info = f(p, rng, **prm)
            info = dict(info or {}); info['params'] = {k: (round(v, 6) if isinstance(v, float) else v) for k, v in prm.items()}
            return np.asarray(sq, float), info
        run.space, run.name, run.__doc__ = space, name, f.__doc__
        MOVES[name] = run
        return run
    return deco


def _region(p, rng, r):
    """indices within r of a random point (r <= 0: all)."""
    if r is None or r <= 0:
        return np.arange(p.n)
    c = np.array([rng.uniform(0, p.s), rng.uniform(0, p.s)])
    return np.nonzero(np.hypot(*(p.sq[:, :2] - c).T) < r)[0]


# ---------- noise ----------
@move('kick', sigma=('log', 0.003, 0.05), ang='x20', r=0.0, bg=0.0)
def kick(p, rng, sigma, ang, r, bg):
    """Gaussian kick of positions (sigma) and angles (ang: 'x20' = 20 sigma degrees, 'sym' = equal corner displacement,
    or a number of degrees) of the squares within r of a random point (r = 0: all); others get sigma bg."""
    sa = 20 * sigma if ang == 'x20' else math.degrees(sigma / (H * math.sqrt(2))) if ang == 'sym' else float(ang)
    sq = p.sq.copy()
    idx = set(_region(p, rng, r).tolist())
    for i in range(p.n):
        s_, a_ = (sigma, sa) if i in idx else (bg, 20 * bg)
        if s_ > 0:
            sq[i, 0] += rng.gauss(0, s_); sq[i, 1] += rng.gauss(0, s_); sq[i, 2] += rng.gauss(0, a_)
    return sq, dict(moved=len(idx))


# ---------- structured ----------
@move('crot', r=('uni', 1.0, 2.5), phi=('choice', [None]))
def crot(p, rng, r, phi):
    """Rotate the cluster within r of a random point by phi degrees (None: 90 w.p. 0.15, else +-U(2, 45))."""
    if phi is None:
        phi = 90.0 if rng.random() < 0.15 else rng.choice([-1, 1]) * rng.uniform(2, 45)
    c0 = np.array([rng.uniform(0, p.s), rng.uniform(0, p.s)])
    idx = np.nonzero(np.hypot(*(p.sq[:, :2] - c0).T) < r)[0]
    c, sn = math.cos(math.radians(phi)), math.sin(math.radians(phi))
    sq = p.sq.copy()
    d = sq[idx, :2] - c0
    sq[idx, 0] = c0[0] + c * d[:, 0] - sn * d[:, 1]; sq[idx, 1] = c0[1] + sn * d[:, 0] + c * d[:, 1]; sq[idx, 2] += phi
    return sq, dict(moved=len(idx), phi=round(phi, 2))


@move('mirror', kind=('choice', ['h', 'v', 'd']), r=('uni', 1.2, 3.0))
def mirror(p, rng, kind, r):
    """Reflect the cluster within r of a random point across a horizontal / vertical / diagonal line through it."""
    c0 = (rng.uniform(0, p.s), rng.uniform(0, p.s))
    sq = p.sq.copy(); m = 0
    for i, (x, y, t) in enumerate(p.sq):
        if math.hypot(x - c0[0], y - c0[1]) < r:
            m += 1; dx, dy = x - c0[0], y - c0[1]
            sq[i] = (x, c0[1] - dy, -t) if kind == 'h' else (c0[0] - dx, y, -t) if kind == 'v' else (c0[0] + dy, c0[1] + dx, 90 - t)
    return sq, dict(moved=m)


def _legacy(fn_name, mod='explore'):
    """Legacy move (draws its own parameters from rng); kept for parity until ported."""
    def f(p, rng, **_):
        import importlib
        m = importlib.import_module(mod)
        sq, desc = getattr(m, fn_name)(p.s, p.tuples(), rng, _A)
        return sq, dict(desc=desc)
    return f


for _name, _fn, _mod in [('aswap', 'mv_aswap', 'mcmin'), ('band', 'mv_band', 'mcmin'), ('reinsert', 'mv_reinsert', 'mcmin'),
                         ('rowslide', 'mv_rowslide', 'explore'), ('chainshift', 'mv_chainshift', 'explore')]:
    move(_name)(_legacy(_fn, _mod))


class _A:                                    # attribute bag the legacy mcmin moves read
    alpha = 2.0; reinsert = 'clear'; remove = 'uniform'; tau = 0.0


# ---------- annealing ----------
@move('melt', rad=('uni', 1.5, 3.5), rmax=('uni', 0.05, 0.2), sweeps=3000, bp=2000.0)
def melt(p, rng, rad, rmax, sweeps, bp):
    """Regional melt (anneal melt): squares within rad of a random square get rounded corners up to rmax, hard-particle
    NPT MC at pressure bp for `sweeps` sweeps (ramp up / hold / down profile), squared up again."""
    i = rng.randrange(p.n)
    with tempfile.TemporaryDirectory() as tmp:
        p.write(f'{tmp}/in.txt')
        r = subprocess.run([ANNEAL, 'melt', '--in', f'{tmp}/in.txt', '--out', f'{tmp}/out.txt', '--cx', repr(float(p.sq[i, 0])),
                            '--cy', repr(float(p.sq[i, 1])), '--rad', f'{rad:.4f}', '--rmax', f'{rmax:.4f}',
                            '--sweeps', str(int(sweeps)), '--bp', repr(float(bp)), '--seed', str(rng.randrange(1 << 30))],
                           capture_output=True, text=True, timeout=600)
        if r.returncode != 0 or not os.path.exists(f'{tmp}/out.txt'):
            return p.sq.copy(), dict(failed=True)
        from .packing import read
        q = read(f'{tmp}/out.txt')
    return q.sq, dict(s_new=q.s)


# ---------- arms ----------
# 10-08 explorer kinds as presets (weights = explore.WEIGHTS; melt off by default there)
ARMS = {
    'kick': ('kick', {}), 'kicksym': ('kick', dict(ang='sym')),
    'lkick': ('kick', dict(r=('uni', 1, 3), sigma=('uni', 0.05, 0.25))),
    'bigkick': ('kick', dict(sigma=('log', 0.05, 0.2))),
    'crot': ('crot', {}), 'aswap': ('aswap', {}), 'band': ('band', {}), 'reinsert': ('reinsert', {}),
    'rowslide': ('rowslide', {}), 'chainshift': ('chainshift', {}), 'mirror': ('mirror', {}), 'melt': ('melt', {}),
}
WEIGHTS = dict(kick=2, kicksym=1, lkick=2, bigkick=2, crot=1, aswap=2, band=1, reinsert=0, rowslide=2, chainshift=0,
               mirror=2)


def parse_arm(spec: str):
    """'kick:sigma=log:0.001:0.3,ang=sym,r=2' -> ('kick', {...}); a bare ARMS name -> ARMS[name]."""
    if spec in ARMS:
        return ARMS[spec]
    mv, _, rest = spec.partition(':')
    prm = {}
    for kv in filter(None, rest.split(',')):
        k, v = kv.split('=', 1)
        parts = v.split(':')
        if parts[0] in ('log', 'uni', 'int') and len(parts) == 3:
            f = int if parts[0] == 'int' else float
            prm[k] = (parts[0], f(parts[1]), f(parts[2]))
        else:
            try:
                prm[k] = float(v) if any(c in v for c in '.e') else int(v)
            except ValueError:
                prm[k] = v
    if mv not in MOVES:
        raise KeyError(f'unknown move {mv}')
    return mv, prm


def apply(p: Packing, arm, seed: int, **override):
    """Proposal from parent p under arm (name in ARMS, or (move, params)) with a reproducible seed."""
    mv, prm = parse_arm(arm) if isinstance(arm, str) else arm
    rng = random.Random(seed)
    fixed = {k: (draw(v, rng) if isinstance(v, tuple) else v) for k, v in {**prm, **override}.items()}
    sq, info = MOVES[mv](p, rng, **fixed)
    info['move'] = mv
    info['arm'] = arm if isinstance(arm, str) else mv
    return Packing(info.get('s_new', p.s), sq, meta=dict(parent=p.meta.get('id'))), info
