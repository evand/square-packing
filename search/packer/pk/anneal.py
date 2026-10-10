"""Schedule-driven annealing (hard rounded squares, NPT MC; `anneal sched` in src/anneal.rs).

A schedule is a path through the (shape, pressure) plane plus rotation share and region, as curves over progress t in
[0, 1]: bp(t) (pressure, log-interpolated), r(t) (corner radius: 0.5 disk ... 0 square), prot(t) (share of rotation
moves).  `Schedule.from_path` builds the named path families with a few parameters; any curve can also be given directly
("t:v,t:v,...").

  sch = Schedule.from_path('disk-first', sweeps=20000, bp1=3000)
  fin, info, traj = run(sch, n=110, seed=1, snaps=20)          # from scratch
  fin, info, traj = run(sch.but(region=(cx, cy, 2.5)), start=p)  # from a packing (shapes start square)
"""
from __future__ import annotations
import json, math, os, subprocess, tempfile
from dataclasses import dataclass, replace, asdict
import numpy as np
from .packing import Packing, read

HERE = os.path.dirname(os.path.abspath(__file__))
ANNEAL = os.environ.get('PK_ANNEAL', os.path.join(HERE, '..', 'target', 'release', 'anneal'))


def curve(knots):
    return ','.join(f'{t:g}:{v:g}' for t, v in knots)


@dataclass(frozen=True)
class Schedule:
    bp: str = '0:1,1:3000'
    r: str = '0:0.5,0.2:0.5,0.8:0,1:0'
    prot: str = '0.5'
    sweeps: int = 20000
    box_every: int = 1
    phi0: float = 0.25
    region: tuple | None = None          # (cx, cy, rad)
    lag_tol: float = 1e9                  # expand the box while shapes lag their radius by > this (off by default: the final uniform
    #                                       growth preserves the arrangement; 0.002 = shapes on schedule, bench lag110)
    bias: str = 'lin'                     # Q4 bias: lin (U = lam n Q4) or harm (U = lam n (Q4 - q0)^2 / 2)
    bias_lam: str = '0'                   # curve
    bias_q0: str = '0'                    # curve
    name: str = ''

    def but(self, **kw):
        return replace(self, **kw)

    def args(self):
        a = ['--bp', self.bp, '--r', self.r, '--prot', self.prot, '--sweeps', str(self.sweeps),
             '--box-every', str(self.box_every), '--phi0', repr(self.phi0), '--lag-tol', repr(self.lag_tol)]
        if self.bias_lam != '0':
            a += ['--bias', self.bias, '--bias-lam', self.bias_lam, '--bias-q0', self.bias_q0]
        if self.region:
            a += ['--cx', repr(float(self.region[0])), '--cy', repr(float(self.region[1])), '--rad', repr(float(self.region[2]))]
        return a

    @staticmethod
    def from_path(path, sweeps=20000, bp0=1.0, bp1=3000.0, rmax=0.5, t0=0.2, t1=0.8, prot=0.5, resoft=0.1):
        """Named paths through the (shape, pressure) plane (pressure always ramps bp0 -> bp1 log-linearly over [0, 1]):
          disk-first   r = rmax until t0, -> 0 by t1       (compress as disks, square up late; the 10-08 constructor)
          square-first r: rmax -> 0 by t0, then squares    (square up while dilute, compress as squares)
          diagonal     r: rmax -> 0 linearly over [0, t1]  (shape and pressure together)
          late-square  r = rmax until t1, -> 0 by 0.95     (square up only near jamming)
          resoft       diagonal, then r back up to `resoft` around (t1 + 1)/2 and down again by 0.97
          hold         r = rmax throughout until 0.95 (melt-like hold for rounding rmax < 0.5)"""
        R = {'disk-first': [(0, rmax), (t0, rmax), (t1, 0), (1, 0)],
             'square-first': [(0, rmax), (t0, 0), (1, 0)],
             'diagonal': [(0, rmax), (t1, 0), (1, 0)],
             'late-square': [(0, rmax), (t1, rmax), (0.95, 0), (1, 0)],
             'resoft': [(0, rmax), (t1, 0), ((t1 + 1) / 2, resoft), (0.97, 0), (1, 0)],
             'hold': [(0, rmax), (0.9, rmax), (0.95, 0), (1, 0)]}[path]
        return Schedule(bp=curve([(0, bp0), (1, bp1)]), r=curve(R), prot=f'{prot:g}', sweeps=int(sweeps),
                        name=f'{path}')


def parse_traj(path):
    """-> list of dict(t, s, bp, r, p: Packing, radii: ndarray)."""
    out = []
    if not os.path.exists(path):
        return out
    lines = open(path).read().splitlines()
    i = 0
    while i < len(lines):
        h = lines[i].split()
        meta = dict(t=float(h[2]), s=float(h[4]), bp=float(h[6]), r=float(h[8]))
        n, s = lines[i + 1].split()
        rows = np.array([list(map(float, l.split())) for l in lines[i + 2:i + 2 + int(n)]])
        out.append(dict(meta, p=Packing(float(s), rows[:, :3]), radii=rows[:, 3]))
        i += 2 + int(n)
    return out


def run(sch: Schedule, start: Packing | None = None, n: int | None = None, seed: int = 1, snaps: int = 0, timeout=1800):
    """Anneal from `start` (shapes start square) or from scratch (n squares, random placement at phi0, shapes r(0)).
    Returns (final squared-up Packing, summary dict, trajectory list)."""
    with tempfile.TemporaryDirectory() as tmp:
        cmd = [ANNEAL, 'sched', '--out', f'{tmp}/out.txt', '--seed', str(seed), *sch.args()]
        if start is not None:
            start.write(f'{tmp}/in.txt'); cmd += ['--in', f'{tmp}/in.txt']
        else:
            cmd += ['--n', str(n)]
        if snaps:
            cmd += ['--snaps', str(snaps)]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if r.returncode != 0:
            return None, dict(error=r.stderr[-300:]), []
        info = json.loads(r.stdout.strip().splitlines()[-1])
        fin = read(f'{tmp}/out.txt')
        traj = parse_traj(f'{tmp}/out.txt.traj') if snaps else []
    return fin, info, traj


def line_occupancy(p: Packing, axis_tol=1.0):
    """Max number of near-axis squares crossed by one horizontal or vertical line (>= k: a full line, side >= k)."""
    a = p.sq[:, 2]
    ax = p.sq[np.minimum(a, 90 - a) < axis_tol][:, :2]
    if len(ax) == 0:
        return 0
    best = 0
    for c in (0, 1):
        v = np.sort(ax[:, c])
        # a line at y0 crosses squares with |v - y0| < 0.5: max over y0 of points in an open window of width 1
        j = 0
        for i in range(len(v)):
            while v[i] - v[j] >= 1.0 - 1e-9:
                j += 1
            best = max(best, i - j + 1)
    return int(best)


def line_load(p: Packing, step=0.005):
    """Exact obstruction diagnostic: max over horizontal and vertical lines of the summed chord lengths of the squares
    the line crosses (chords are disjoint, so load <= s; a line with load = s is a grid lock; load >= k forces s >= k).
    Evaluated on a grid of line positions (step), so it slightly underestimates the max."""
    th = np.radians(p.sq[:, 2])
    c, s_ = np.cos(th), np.sin(th)
    off = np.array([[.5, .5], [-.5, .5], [-.5, -.5], [.5, -.5]])
    V = p.sq[:, None, :2] + np.stack([off[:, 0] * c[:, None] - off[:, 1] * s_[:, None],
                                      off[:, 0] * s_[:, None] + off[:, 1] * c[:, None]], -1)       # (n, 4, 2)
    best = 0.0
    for ax in (1, 0):                      # ax = coordinate held fixed by the line (1: horizontal lines y = y0)
        o = 1 - ax
        A, B = V, np.roll(V, -1, axis=1)   # edges A -> B
        ys = np.arange(step / 2, p.s, step)
        a, b = A[..., ax][None], B[..., ax][None]                  # (1, n, 4)
        y0 = ys[:, None, None]
        cross = (np.minimum(a, b) <= y0) & (np.maximum(a, b) > y0)
        u = np.where(cross, (y0 - a) / np.where(b - a == 0, 1, b - a), np.nan)
        xo = A[..., o][None] + u * (B[..., o] - A[..., o])[None]   # intersection coordinate along the line
        chord = np.nanmax(xo, axis=2) - np.nanmin(xo, axis=2)       # (ny, n), nan if not crossed
        load = np.nansum(chord, axis=1)
        best = max(best, float(load.max()))
    return best
