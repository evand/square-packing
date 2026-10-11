#!/usr/bin/env python3
"""Ring DP conditioned on interior squares (float demonstrator).

For fixed poses of some interior squares (obstacles), the wall ring of a closed-cell state is still a path problem:
a later square needs only the least position of its predecessor (gaps are up-sets, see wallchain.py), so the state
stays F(q, t) = least ALLOWED position, where 'allowed' now also excludes the obstacles.  Obstacles make a square's
allowed positions a non-interval set, handled by a fine position grid ("first allowed p >= L").

Relaxation (sound for exclusion once made rigorous): drops ring pairs that are not consecutive (cross-corner pairs,
i/i+2 pairs, the pair across a vacant ring cell) and the interior squares not conditioned on.
FLOAT: sampled grids; 'infeasible' is evidence, not proof.

Usage: ringcond.py MASK --cond interior-S [--cond interior-W] [--step 0.04] [--tstep 3] [--jobs 8]
Reports, per conditioned square (or jointly), how many sampled obstacle poses leave the ring feasible."""
import sys, os, json, time, argparse, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wallchain
from wallchain import U, RING, CORNER, FRAME, rect_in_frame, gap, wall_floor, load_cells, ring_path
from subsets import load as load_polys

INV = {'S': lambda p, q: (p, q), 'E': lambda p, q: (U - q, p), 'N': lambda p, q: (U - p, U - q),
       'W': lambda p, q: (q, U - p)}


def disjoint(x1, y1, t1, x2, y2, t2):
    """closed unit squares disjoint-or-touching (separating axis); broadcasting."""
    dx, dy = x2 - x1, y2 - y1
    ok = np.zeros(np.broadcast(x1, y1, t1, x2, y2, t2).shape, bool)
    for phi in (t1, t1 + np.pi / 2, t2, t2 + np.pi / 2):
        nx, ny = np.cos(phi), np.sin(phi)
        H = 0.5 * wallchain.SIGMA * (np.abs(np.cos(phi - t1)) + np.abs(np.sin(phi - t1)) +
                                     np.abs(np.cos(phi - t2)) + np.abs(np.sin(phi - t2)))
        ok |= np.abs(nx * dx + ny * dy) >= H - 1e-12
    return ok


class Ring:
    def __init__(self, chain, cells, dq=0.02, dp=0.004, tstep=2.0):
        self.chain, self.cells = chain, cells
        self.t = np.deg2rad(np.arange(0, 90, tstep))
        self.dq, self.dp = dq, dp
        self.sq = []       # per square: dict(wall, pgrid, qgrid, base_allowed[p,q,t], (frame2 data for corners))
        wall = None
        for k, name in enumerate(chain):
            if name.startswith('corner'):
                win, wout = CORNER[name]
            else:
                win = wout = name.split('-')[1][0]
            if k == 0:
                wall = wout if not name.startswith('corner') else wout
                # a chain starting at a corner starts in the corner's outgoing frame
            entry = dict(name=name, wall=wall)
            p0, p1, q0, q1 = rect_in_frame(cells[name], wall)
            entry['p'] = np.arange(p0, p1 + 1e-12, dp)
            entry['q'] = np.linspace(q0, q1, max(2, int(round((q1 - q0) / dq)) + 1))
            entry['base'] = self._base(wall, entry['p'], entry['q'])
            if name.startswith('corner') and k > 0:
                P0, P1, Q0, Q1 = rect_in_frame(cells[name], wout)
                entry['wall2'] = wout
                entry['q2'] = np.linspace(Q0, Q1, max(2, int(round((Q1 - Q0) / dq)) + 1))
                wall = wout
            self.sq.append(entry)

    def _xy(self, wall, p, q):
        return INV[wall](p, q)

    def _base(self, wall, pg, qg):
        P, Q, T = np.meshgrid(pg, qg, self.t, indexing='ij')
        x, y = self._xy(wall, P, Q)
        h = wall_floor(T)
        return (x >= h - 1e-12) & (x <= U - h + 1e-12) & (y >= h - 1e-12) & (y <= U - h + 1e-12)

    def allowed(self, entry, obstacles):
        A = entry['base'].copy()
        if obstacles:
            P, Q, T = np.meshgrid(entry['p'], entry['q'], self.t, indexing='ij')
            x, y = self._xy(entry['wall'], P, Q)
            for (ox, oy, ot) in obstacles:
                A &= disjoint(x, y, T, ox, oy, ot)
        return A

    @staticmethod
    def first_allowed(A, pg, L):
        """A[p, q, t] bool, L[q, t] lower bounds -> least allowed p >= L (inf if none)."""
        idx0 = np.searchsorted(pg, L - 1e-12)                       # first grid index >= L
        npg = len(pg)
        # next allowed index at or after i: suffix scan
        nxt = np.full(A.shape, npg, dtype=np.int32)
        cur = np.full(A.shape[1:], npg, dtype=np.int32)
        for i in range(npg - 1, -1, -1):
            cur = np.where(A[i], i, cur)
            nxt[i] = cur
        idx0c = np.minimum(idx0, npg - 1)
        j = np.take_along_axis(nxt, idx0c[None, :, :], axis=0)[0]
        j = np.where(idx0 >= npg, npg, j)
        out = np.where(j < npg, pg[np.minimum(j, npg - 1)], np.inf)
        return out

    def run(self, obstacles=()):
        t = self.t
        F = q = None
        for k, e in enumerate(self.sq):
            A = self.allowed(e, obstacles)
            if k == 0:
                L = np.full((len(e['q']), len(t)), -np.inf)
                F = self.first_allowed(A, e['p'], L)
                q = e['q']
            else:
                dq = e['q'][:, None] - q[None, :]                     # (nB, nA)
                L = np.full((len(e['q']), len(t)), np.inf)
                fin = np.isfinite(F)
                for jA in range(len(t)):
                    if not fin[:, jA].any():
                        continue
                    g = gap(dq[:, :, None], t[jA], t[None, None, :])      # (nB, nA, nt)
                    L = np.minimum(L, (F[:, jA][None, :, None] + g).min(axis=1))
                F = self.first_allowed(A, e['p'], L)
                q = e['q']
            if 'wall2' in e:
                # corner turn: p2 = q1, q2 = U - p1; G(q2, t) = least q1 with p1 = U - q2 allowed and >= F(q1, t)
                q2 = e['q2']
                G = np.full((len(q2), len(t)), np.inf)
                pg = e['p']
                # feasible corner poses: allowed and p1 >= F(q1, t); a depth grid point q2 on the new wall accepts any
                # p1 within half a depth step of U - q2 (errs towards feasibility, i.e. towards the relaxation)
                feas = A & (pg[:, None, None] >= F[None, :, :] - 1e-12)          # (np, nq1, nt)
                h = self.dq / 2
                for i2, qq in enumerate(q2):
                    lo = np.searchsorted(pg, U - qq - h - 1e-12)
                    hi = np.searchsorted(pg, U - qq + h + 1e-12)
                    if hi <= lo:
                        continue
                    ok = feas[lo:hi].any(axis=0)                                   # (nq1, nt)
                    G[i2] = np.where(ok.any(axis=0), np.where(ok, q[:, None], np.inf).min(axis=0), np.inf)
                F, q = G, q2
            if not np.isfinite(F).any():
                return False, k
        return True, len(self.sq)


def interior_poses(name, polys, step, tstep):
    A, B, V = polys[name]
    xs = np.arange(V[:, 0].min(), V[:, 0].max() + 1e-12, step)
    ys = np.arange(V[:, 1].min(), V[:, 1].max() + 1e-12, step)
    pts = [(x, y) for x in xs for y in ys if np.all(A @ np.array([x, y]) <= B + 1e-12)]
    ts = np.deg2rad(np.arange(0, 90, tstep))
    out = []
    for (x, y) in pts:
        for tt in ts:
            h = wall_floor(tt)
            if h <= x <= U - h and h <= y <= U - h:
                out.append((x, y, tt))
    return out


def chains_of(mask, names):
    path = ring_path(mask, names)
    chains, cur = [], []
    for c in path:
        if c is None:
            if cur:
                chains.append(cur)
            cur = []
        else:
            cur.append(c)
    if cur:
        chains.append(cur)
    return chains


_RINGS = None


def _work(args):
    pose = args
    ok = True
    for R in _RINGS:
        f, _ = R.run([pose])
        ok &= f
        if not ok:
            break
    return ok


def _init(mask, dq, dp, tstep, sigma=1.0):
    global _RINGS
    wallchain.SIGMA = sigma
    cells, names = load_cells()
    _RINGS = [Ring(ch, cells, dq, dp, tstep) for ch in chains_of(mask, names)]


if __name__ == '__main__':
    import multiprocessing as mp
    a = argparse.ArgumentParser()
    a.add_argument('mask', type=int)
    a.add_argument('--cond', action='append', default=[])
    a.add_argument('--step', type=float, default=0.04)
    a.add_argument('--tstep', type=float, default=3.0)
    a.add_argument('--dq', type=float, default=0.02)
    a.add_argument('--dp', type=float, default=0.004)
    a.add_argument('--ring-tstep', type=float, default=2.0)
    a.add_argument('--jobs', type=int, default=8)
    a.add_argument('--shrink', type=float, default=0.0, help='use squares of side 1 - shrink')
    a.add_argument('--control-family', action='store_true', help='obstacles = the family packing (site data)')
    a.add_argument('--out')
    o = a.parse_args()
    cells, names = load_cells()
    polys, _ = load_polys()
    t0 = time.time()
    _init(o.mask, o.dq, o.dp, o.ring_tstep, 1 - o.shrink)
    f, k = (True, None)
    res = {}
    for R in _RINGS:
        ff, kk = R.run([])
        print('unconditioned chain', [e['name'] for e in R.sq][0], '..', [e['name'] for e in R.sq][-1],
              'feasible' if ff else f'INFEASIBLE at {kk}', flush=True)
    for c in o.cond:
        poses = interior_poses(c, polys, o.step, o.tstep)
        with mp.Pool(o.jobs, initializer=_init, initargs=(o.mask, o.dq, o.dp, o.ring_tstep, 1 - o.shrink)) as pool:
            oks = pool.map(_work, poses, chunksize=8)
        good = [p for p, ok in zip(poses, oks) if ok]
        res[c] = dict(poses=len(poses), ring_feasible=len(good),
                      sample=[[round(x, 3), round(y, 3), round(float(np.rad2deg(tt)), 1)] for x, y, tt in good[:20]])
        print(f'{c}: {len(good)} of {len(poses)} sampled poses leave the ring feasible  ({time.time() - t0:.0f} s)',
              flush=True)
    if o.out:
        json.dump(dict(mask=o.mask, args=vars(o), result=res), open(o.out, 'w'), indent=1)
