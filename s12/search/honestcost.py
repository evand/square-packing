#!/usr/bin/env python3
"""honest cost of a `t = 4` QSTAB + odd-polygon dual read as a weighted cover (task honest-cost).

Doc: `search/HONEST.md`.  Reads a `cliquelever.py` checkpoint state (`cl_TAG_poses.txt`,
`_rows.txt`, `_cliques.txt`, `_pgons.json`), re-solves its LP on the FIXED pose set to get the
dual, and then treats that dual as a cover of `[0,4]^2`:

    Theta(w) = sum_p y_p  +  sum_K z_K  +  sum_Pi (k-1)/2 * zeta_Pi
    capture_w(S) = sum_{p in S} y_p  +  sum_K z_K [S credited K]  +  sum_Pi zeta_Pi [S credited Pi]
    honest(w) = Theta(w) / min over ALL admissible poses S of capture_w(S)

Semantics: `t = 4`, closed unit squares, closed containment, a point on the boundary of `S` counts
(`search/ZEROMARGIN.md` §1).  A pose is admissible iff its closed square lies in `[0,4]^2`, i.e.
`cx, cy in [w(theta)/2, 4 - w(theta)/2]` with `w(theta) = |cos| + |sin|`; only `theta in [0, 90)`
is needed because a unit square is invariant under a quarter turn.  D4 symmetry of the CONTAINER
is not assumed anywhere (these runs pin nothing, `sym 1` in the pose files), so the whole domain
is scanned.

**Crediting rules for the clique rows** (`--rule`):

* `meet` (default; the rule `cliquelever.price`'s `clique_cost` uses): `S` is credited `z_K` iff
  `S` closed-meets every member of `K`.  This is the *generous* rule: on the loaded poses it
  agrees with the LP row (the rows are maximal over the loaded set), but off them it credits any
  square that could be added to the row.  It is generous in a way worth stating plainly: the set
  `{S : S meets every member of K}` need NOT be a clique -- two disjoint squares can each meet
  every member of a non-Helly clique -- so a cover credited this way is not by itself a valid
  weighting.  `min capture` under `meet` is therefore an UPPER bound on the min capture of any
  sound rule, and `honest` under `meet` is a LOWER bound on the honest cost.
* `core`: `S` is credited `z_K` iff `S` contains the *core* of `K`, the common intersection
  `C_K = intersection of the member squares`; rows with `C_K = {}` credit nobody.  This is the
  rule a `clique_of_cores` verifier can check (`search/BOXCLIQUE.md`: cores pairwise meet, and
  here they coincide, so the credited set really is a clique) and it is SOUND.  It is the
  box/anchor-clique rule of the brief specialised to rows whose "boxes" are single poses.

The polygon rows are credited by their own rule in both cases: `Pi` with anchors `A_0..A_{k-1}`
counts iff `S` contains both endpoints of some side `[A_i, A_{i+1}]` (`search/rankfamily.py`).
That rule is sound over the continuum -- the lemma is proved for arbitrary admissible squares.

Everything here is FLOAT (the brief: this is a packing-side measurement of a cover, nothing is
certified).  The point atoms and the containment tests use `packing_dual.capture`'s `tol = 1e-8`,
the clique/polygon tests a `--geo-tol` slack; both slacks make `capture` larger, hence `honest`
smaller, i.e. they err on the side of the cover looking good.

Commands
--------
    dual   TAG  re-solve the state's LP, dump y/z/zeta to runs/hc_TAG_dual.npz, and self-test
                (Theta = LP, and capture >= 1 on every pose of P)
    scan   TAG  minimise capture over the admissible poses: the 0.04/2.5deg lattice, the
                `family_rows.py` families, the poses of P, the row points as centres, a pattern
                search, a knife-edge (+-eps) refinement and a 1.69M-pose dense confirmation
    fine   TAG  the knife-edge stage alone, from a finished scan's worst poses
    sound  TAG  is the `meet` rule a valid cover rule?  looks for two DISJOINT admissible squares
                that both credit one clique row
    repair TAG  add the worst poses found by `scan` as columns of the packing LP (= rows of the
                cover), re-solve, then re-scan (step 4 of the brief)
"""
import argparse
import json
import math
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import leaf_ceiling as lc          # noqa: E402
import packing_dual as pd          # noqa: E402
import cliquelever as cl           # noqa: E402
import rankfamily as rf            # noqa: E402

RUNS = cl.RUNS
T = 4.0


# ====================================================================== the state and its dual
def run_args(**kw):
    """the `cliquelever run` argument namespace this module needs (defaults = the pure instance:
    no pinned regions, no chord rows)"""
    a = argparse.Namespace(
        cmd='run', TAG='hc', t='4', r='1', bnd_delta=1e-6, clique_time=30.0, cq_top=30,
        anatomy=0, procs=2, exact=[], load_poses=[], load_rows=[], resume_rows=[],
        resume_cliques=[], resume_pgons=[], Q=10 ** 5, Dc=10 ** 6, corners='....',
        patterns='........', chord=False, extend=True, threads=4, lp_tlim=0.0, iters=1,
        row_cap=3000, rows0=0, row_pitch=0.0, cq_want=0, cq_age=0, ktol=1e-6, time=1e9,
        ckpt=0, price=0, price_pitch=0.04, price_dth=2.5, price_tol=1e-7, cg_want=400,
        finalize=False, pent=[], pent_want=0, pent_cands=240, pent_restarts=60,
        pent_time=90.0, pent_seed=20260911, no_warm=True, master=False, master_add=2000,
        master_passes=0, master_trim=0, lattice_every=0)
    for k, v in kw.items():
        setattr(a, k, v)
    return a


def load_state(tag, base=RUNS, log=print, **kw):
    """the pose set, coverage rows, clique rows and polygon rows of a `cl_TAG_*` checkpoint"""
    a = run_args(TAG=tag, **kw)
    ps = cl.Poses(Fr(a.t), Fr(a.r), a.bnd_delta)
    n = ps.add_exact_file(os.path.join(base, f'cl_{tag}_poses.txt'))
    ps.finish()
    log(f'   {n} exact poses, {ps.ncol} columns')
    lever = cl.Lever(ps, a, log)
    nr = lever.resume_rows(os.path.join(base, f'cl_{tag}_rows.txt'))
    nc = lever.resume_cliques(os.path.join(base, f'cl_{tag}_cliques.txt'))
    pgpath = os.path.join(base, f'cl_{tag}_pgons.json')
    npg, nbad = rf.resume(lever, pgpath) if os.path.exists(pgpath) else (0, 0)
    log(f'   {nr} coverage rows, {nc} clique rows, {npg} polygon rows ({nbad} dropped)')
    assert lever.target is None, lever.target
    return ps, lever


def solve_state(lever, log=print, master=False):
    """the LP on the loaded column set, solved to the exact FULL-LP optimum: the dual has to be
    feasible for every loaded column or the self-test below fails by construction.  `master`
    routes through `Lever.solve_master` with `force_full_master` set, which prices every loaded
    column after every pass and only stops when none is improving -- the same optimum, reached by
    a sequence of small ipm solves instead of one huge one (search/CLMASTER.md)."""
    t0 = time.time()
    if master:
        lever.master = True
        lever.force_full_master = True
        lever.a.master_passes = 0
        sol = lever.solve_master()
        if sol is not None:
            log(f'   master closed: {lever.master_n} of {lever.ps.ncol} columns, '
                f'{lever.master_passes_used} passes, {lever.rc_npos} improving left')
    else:
        sol = lever.solve_full()
    if sol is None:
        raise SystemExit('LP failed')
    x, obj = sol
    lp = float(obj)
    log(f'   LP = {lp:.9f} in {time.time() - t0:.0f}s')
    rc = lever.price_columns()
    log(f'   dual feasibility on the {lever.ps.ncol} loaded columns: max rc = {rc.max():+.3e} '
        f'({int((rc > 1e-7).sum())} columns with rc > 1e-7)')
    return x, lp, rc


def dump_dual(lever, x, lp, path):
    """write y (with its row points), z (with the clique member lists) and zeta (with the exact
    polygon anchors) to an npz"""
    mem = np.concatenate([np.array(c['members'], dtype=np.int64) for c in lever.cliques]) \
        if lever.cliques else np.zeros(0, np.int64)
    off = np.cumsum([0] + [len(c['members']) for c in lever.cliques]).astype(np.int64)
    np.savez_compressed(
        path, lp=np.array(lp), x=x,
        rx=np.array([float(u) / float(w) for u, _, w in lever.rows]),
        ry=np.array([float(v) / float(w) for _, v, w in lever.rows]),
        y=lever.y, z=lever.z, pgz=lever.pgz,
        cmem=mem, coff=off,
        pgrhs=np.array([g['rhs'] for g in lever.pgons]),
        pganch=np.array(json.dumps([[[str(u) for u in p] for p in g['anchors']]
                                    for g in lever.pgons])),
        pf=lever.ps.F)
    return path


# ====================================================================== the cover
def poly_clip(P, nx, ny, c, tol):
    """Sutherland-Hodgman: the part of the convex polygon P with nx*x + ny*y <= c + tol"""
    out = []
    m = len(P)
    for i in range(m):
        x0, y0 = P[i]
        x1, y1 = P[(i + 1) % m]
        d0 = nx * x0 + ny * y0 - c - tol
        d1 = nx * x1 + ny * y1 - c - tol
        if d0 <= 0:
            out.append((x0, y0))
        if (d0 < 0 < d1) or (d1 < 0 < d0):
            s = d0 / (d0 - d1)
            out.append((x0 + s * (x1 - x0), y0 + s * (y1 - y0)))
    return out


def square_polygon(cx, cy, co, si):
    h = 0.5
    return [(cx + co * a - si * b, cy + si * a + co * b)
            for a, b in ((h, h), (-h, h), (-h, -h), (h, -h))]


def clique_core(F, members, tol=1e-9):
    """the common intersection of the member squares, as a convex polygon (float).  Each square
    contributes its four half planes in its own frame.

    `tol` matters: the binding cliques here are *point cliques* (all loaded squares through one
    point `p`), whose exact core is the single point `p`, and an exact clipper loses it to
    rounding after a few hundred clips.  With `tol = 1e-9` a point core survives as a blob of
    diameter `~1e-9`, which is what the credit test wants (`S` contains `p`)."""
    i0 = members[0]
    P = square_polygon(F[i0, 0], F[i0, 1], F[i0, 2], F[i0, 3])
    for i in members:
        cx, cy, co, si = F[i]
        u = cx * co + cy * si
        v = -cx * si + cy * co
        for (nx, ny, c) in ((co, si, u + 0.5), (-co, -si, -(u - 0.5)),
                            (-si, co, v + 0.5), (si, -co, -(v - 0.5))):
            P = poly_clip(P, nx, ny, c, tol)
            if not P:
                return []
    return P


class Cover:
    """a dual read as a cover: the point atoms, the clique rows (with their cores) and the
    polygon rows, with a batched `capture` over arbitrary poses"""

    def __init__(self, d, geo_tol=1e-9, ytol=1e-9, log=print):
        self.geo_tol = geo_tol
        self.log = log
        y = d['y']
        sel = np.nonzero(y > ytol)[0]
        o = np.argsort(d['rx'][sel], kind='stable')          # pd.capture wants ax ascending
        sel = sel[o]
        self.ax = np.ascontiguousarray(d['rx'][sel])
        self.ay = np.ascontiguousarray(d['ry'][sel])
        self.aw = np.ascontiguousarray(y[sel])
        self.n_atoms = len(sel)
        self.point_total = float(y.sum())
        self.F = d['pf']
        z = d['z']
        cmem, coff = d['cmem'], d['coff']
        self.cz = []
        self.cmem = []
        for i in range(len(z)):
            if z[i] > ytol:
                self.cz.append(float(z[i]))
                self.cmem.append(cmem[coff[i]:coff[i + 1]])
        self.cz = np.array(self.cz)
        self.clique_total = float(z.sum())
        self.n_cliques = len(self.cz)
        # cores, for the sound (box/anchor) rule
        self.ccore = [clique_core(self.F, m) for m in self.cmem]
        self.n_core_empty = sum(1 for p in self.ccore if not p)
        self.core_empty_mass = float(sum(zz for zz, p in zip(self.cz, self.ccore) if not p))
        self.core_diam = [0.0 if not p else
                          float(max(abs(u[0] - v[0]) + abs(u[1] - v[1]) for u in p for v in p))
                          for p in self.ccore]
        self.n_core_point = sum(1 for p, dd in zip(self.ccore, self.core_diam)
                                if p and dd <= 1e-6)
        pgz = d['pgz']
        anch = json.loads(str(d['pganch']))
        rhs = d['pgrhs']
        self.pz, self.ppts, self.prhs = [], [], []
        for i in range(len(pgz)):
            if pgz[i] > ytol:
                self.pz.append(float(pgz[i]))
                self.prhs.append(float(rhs[i]))
                self.ppts.append(np.array([[float(Fr(u)) for u in p] for p in anch[i]]))
        self.pz = np.array(self.pz)
        self.n_pgons = len(self.pz)
        self.pgon_total = float((pgz * rhs).sum())
        self.theta = self.point_total + self.clique_total + self.pgon_total
        self.lp = float(d['lp'])

    def describe(self):
        return (f'cover: Theta = {self.theta:.9f} (points {self.point_total:.6f} + cliques '
                f'{self.clique_total:.6f} + polygons {self.pgon_total:.6f}); LP = {self.lp:.9f}; '
                f'{self.n_atoms} atoms, {self.n_cliques} clique rows with dual '
                f'({self.n_core_empty} of them with an EMPTY core, carrying '
                f'{self.core_empty_mass:.6f}; {self.n_core_point} with a POINT core), '
                f'{self.n_pgons} polygon rows with dual')

    # ------------------------------------------------------------------ the three terms
    def point_capture(self, poses, threads=4):
        if not hasattr(self, '_lib'):
            self._lib = pd.build_lib()
        return pd.capture(self._lib, self.ax, self.ay, self.aw, poses, threads)

    def clique_capture(self, poses, rule='meet'):
        """sum of z_K over the credited clique rows.  `meet`: S closed-meets every member (the
        pricer's rule).  `core`: S contains the whole core of K (sound)."""
        P = np.ascontiguousarray(np.asarray(poses, dtype=float).reshape(-1, 3))
        out = np.zeros(len(P))
        if not len(P) or not self.n_cliques:
            return out
        cx, cy = P[:, 0], P[:, 1]
        co, si = np.cos(P[:, 2]), np.sin(P[:, 2])
        tol = self.geo_tol
        F = self.F
        for ci in range(self.n_cliques):
            mem = self.cmem[ci]
            if rule == 'core':
                V = self.ccore[ci]
                if not V:
                    continue
                Vx = np.array([p[0] for p in V])
                Vy = np.array([p[1] for p in V])
                sel = np.nonzero((cx >= Vx.min() - 0.7072) & (cx <= Vx.max() + 0.7072)
                                 & (cy >= Vy.min() - 0.7072) & (cy <= Vy.max() + 0.7072))[0]
                for k in range(len(Vx)):
                    if not len(sel):
                        break
                    dx = Vx[k] - cx[sel]
                    dy = Vy[k] - cy[sel]
                    u = np.abs(dx * co[sel] + dy * si[sel])
                    v = np.abs(-dx * si[sel] + dy * co[sel])
                    sel = sel[np.maximum(u, v) <= 0.5 + tol]
                out[sel] += self.cz[ci]
                continue
            # ---- `meet`: the pricer's arithmetic (cliquelever.price.clique_cost), uncapped
            bx0, bx1 = F[mem, 0].min() - 1.4143, F[mem, 0].max() + 1.4143
            by0, by1 = F[mem, 1].min() - 1.4143, F[mem, 1].max() + 1.4143
            sel = np.nonzero((cx >= bx0) & (cx <= bx1) & (cy >= by0) & (cy <= by1))[0]
            for k in mem:
                if not len(sel):
                    break
                Fk = F[k]
                g0, g1, g2, g3 = cx[sel], cy[sel], co[sel], si[sel]
                cr = Fk[2] * g2 + Fk[3] * g3
                sr = Fk[2] * g3 - Fk[3] * g2
                wj = 0.5 * (np.abs(cr) + np.abs(sr))
                ddx = g0 - Fk[0]
                ddy = g1 - Fk[1]
                m1 = 0.5 + wj - np.abs(ddx * Fk[2] + ddy * Fk[3])
                m2 = 0.5 + wj - np.abs(-ddx * Fk[3] + ddy * Fk[2])
                m3 = 0.5 + wj - np.abs(ddx * g2 + ddy * g3)
                m4 = 0.5 + wj - np.abs(-ddx * g3 + ddy * g2)
                sel = sel[np.minimum(np.minimum(m1, m2), np.minimum(m3, m4)) >= -tol]
            out[sel] += self.cz[ci]
        return out

    def pgon_capture(self, poses):
        """sum of zeta_Pi over the polygon rows S is credited: S contains both endpoints of some
        side (rankfamily.pgon_cost's rule, in floats)"""
        P = np.ascontiguousarray(np.asarray(poses, dtype=float).reshape(-1, 3))
        out = np.zeros(len(P))
        if not len(P) or not self.n_pgons:
            return out
        cx, cy = P[:, 0], P[:, 1]
        co, si = np.cos(P[:, 2]), np.sin(P[:, 2])
        tol = self.geo_tol
        for gi in range(self.n_pgons):
            A = self.ppts[gi]
            sel = np.nonzero((cx >= A[:, 0].min() - 0.7072) & (cx <= A[:, 0].max() + 0.7072)
                             & (cy >= A[:, 1].min() - 0.7072) & (cy <= A[:, 1].max() + 0.7072))[0]
            if not len(sel):
                continue
            cs, ss, xs, ys = co[sel], si[sel], cx[sel], cy[sel]
            cont = []
            for (px, py) in A:
                dx, dy = px - xs, py - ys
                u = np.abs(dx * cs + dy * ss)
                v = np.abs(-dx * ss + dy * cs)
                cont.append(np.maximum(u, v) <= 0.5 + tol)
            k = len(A)
            inrow = np.zeros(len(sel), bool)
            for i in range(k):
                inrow |= cont[i] & cont[(i + 1) % k]
            out[sel[inrow]] += self.pz[gi]
        return out

    def capture(self, poses, rule='meet', threads=4, parts=False):
        P = np.asarray(poses, dtype=float).reshape(-1, 3)
        pt = self.point_capture(P, threads)
        cq = self.clique_capture(P, rule)
        pg = self.pgon_capture(P)
        if parts:
            return pt + cq + pg, pt, cq, pg
        return pt + cq + pg


# ====================================================================== the pose domain
def admissible_clamp(P, t=T, margin=0.0):
    """clamp the centres into the admissible box of their own angle"""
    P = np.array(P, dtype=float).reshape(-1, 3)
    w2 = 0.5 * (np.abs(np.cos(P[:, 2])) + np.abs(np.sin(P[:, 2])))
    lo, hi = w2 + margin, t - w2 - margin
    P[:, 0] = np.clip(P[:, 0], lo, hi)
    P[:, 1] = np.clip(P[:, 1], lo, hi)
    return P


def lattice(pitch, dth_deg, t=T):
    """the full-domain lattice: theta in [0, 90) (a unit square is quarter-turn invariant), and
    for each theta an even grid of admissible centres"""
    out = []
    for th in np.arange(0.0, math.pi / 2 - 1e-12, math.radians(dth_deg)):
        w2 = 0.5 * (abs(math.cos(th)) + abs(math.sin(th)))
        lo, hi = w2, t - w2
        k = max(int(math.floor((hi - lo) / pitch)), 1)
        g = np.linspace(lo, hi, k + 1)
        A, B = np.meshgrid(g, g, indexing='ij')
        out.append(np.column_stack([A.ravel(), B.ravel(), np.full(A.size, th)]))
    return np.concatenate(out)


def pattern_search(cover, seeds, rule, threads, steps=None, rounds=3, log=print):
    """batched coordinate/pattern search on `capture` (piecewise constant in the centre, so a
    fixed pattern of steps with a geometric ladder, not a gradient method)"""
    if steps is None:
        steps = [(0.04, 2.5), (0.02, 1.0), (0.01, 0.5), (4e-3, 0.25), (1e-3, 0.1),
                 (2.5e-4, 0.04), (6e-5, 0.01), (1.5e-5, 2.5e-3)]
    cur = admissible_clamp(seeds)
    val = cover.capture(cur, rule, threads)
    for (h, dd) in steps:
        dth = math.radians(dd)
        moves = np.array([[h, 0, 0], [-h, 0, 0], [0, h, 0], [0, -h, 0], [0, 0, dth], [0, 0, -dth],
                          [h, h, 0], [-h, -h, 0], [h, -h, 0], [-h, h, 0],
                          [h, 0, dth], [-h, 0, -dth], [0, h, dth], [0, -h, -dth],
                          [h, h, dth], [-h, -h, -dth], [h, -h, dth], [-h, h, -dth]])
        for _ in range(rounds):
            cand = admissible_clamp((cur[:, None, :] + moves[None, :, :]).reshape(-1, 3))
            v = cover.capture(cand, rule, threads).reshape(len(cur), len(moves))
            j = np.argmin(v, axis=1)
            best = v[np.arange(len(cur)), j]
            better = best < val - 1e-13
            if not better.any():
                break
            cand = cand.reshape(len(cur), len(moves), 3)
            cur[better] = cand[np.arange(len(cur)), j][better]
            val = np.minimum(val, best)
        log(f'      pattern h={h:g} dth={dd:g}: min capture {val.min():.9f}')
    return cur, val


def eps_refine(cover, seeds, rule, threads, scales=(1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9),
               rounds=2, log=print):
    """knife-edge refinement.  `capture` is DISCONTINUOUS at the heavy grid lines: a square whose
    edge sits exactly on `x = 1` captures every atom on that line, and one `1e-7` off it captures
    none of them (`packing_dual.capture`'s tolerance is `1e-8`).  That is `ZEROMARGIN.md` §1 /
    `RUNG2.md` §4.3's knife edge, and the infimum of `capture` over the admissible poses is
    attained just OFF the lines, at a scale no pattern ladder with a `1e-5` floor can see.  So
    probe every sign pattern of `+-eps` in `(cx, cy, theta)` at each scale."""
    cur = admissible_clamp(seeds)
    val = cover.capture(cur, rule, threads)
    for e in scales:
        moves = np.array([[a, b, c] for a in (0.0, e, -e) for b in (0.0, e, -e)
                          for c in (0.0, e, -e) if a or b or c])
        for _ in range(rounds):
            cand = admissible_clamp((cur[:, None, :] + moves[None, :, :]).reshape(-1, 3))
            v = cover.capture(cand, rule, threads).reshape(len(cur), len(moves))
            j = np.argmin(v, axis=1)
            best = v[np.arange(len(cur)), j]
            better = best < val - 1e-15
            if not better.any():
                break
            cand = cand.reshape(len(cur), len(moves), 3)
            cur[better] = cand[np.arange(len(cur)), j][better]
            val = np.minimum(val, best)
        log(f'      eps {e:g}: min capture {val.min():.9f}')
    return cur, val


def cmd_fine(a):
    """the knife-edge stage on its own: re-refine the worst poses a `scan` already found (its
    `hc_TAG_RULE_worst.txt`) plus the item-3 families, at `1e-4 .. 1e-9`"""
    t0 = time.time()

    def log(m):
        print(f'[{time.time() - t0:7.0f}s] {m}', flush=True)

    d = np.load(os.path.join(RUNS, f'hc_{a.TAG}_dual.npz'))
    cov = Cover(d, a.geo_tol, log=log)
    log(f'# honestcost fine {a.TAG} rule={a.rule}')
    log(cov.describe())
    seeds = []
    wpath = os.path.join(RUNS, f'hc_{a.TAG}_{a.rule}_worst.txt')
    if os.path.exists(wpath):
        seeds.append(np.loadtxt(wpath).reshape(-1, 3))
    if a.families and os.path.exists(a.families):
        FP = np.array([[float(v) for v in q[:3]] for q in
                       (line.split() for line in open(a.families))
                       if len(q) >= 3 and not q[0].startswith('#')])
        FP = admissible_clamp(FP)
        FP[:, 2] = np.mod(FP[:, 2], math.pi / 2)
        FP = admissible_clamp(FP)
        VF = cov.capture(FP, a.rule, a.threads)
        log(f'   families: {len(FP)} poses, min capture {VF.min():.9f}')
        seeds.append(FP[np.argsort(VF)[:a.nseed]])
    S = np.concatenate(seeds)
    V0 = cov.capture(S, a.rule, a.threads)
    o = np.argsort(V0)[:a.nseed]
    log(f'   refining {len(o)} seeds (best {V0[o].min():.9f})')
    P, V = eps_refine(cov, S[o], a.rule, a.threads, rounds=a.rounds, log=log)
    i = int(np.argmin(V))
    mn = float(V[i])
    h = cov.theta / mn if mn > 0 else float('inf')
    log(f'   MIN capture {mn:.9f} at {P[i].tolist()} -> HONEST = {h:.6f}')
    out = dict(tag=a.TAG, rule=a.rule, theta=cov.theta, min_capture=mn, honest=h,
               minimisers=[anatomy(cov, p, a.rule)
                           for p in P[np.argsort(V)[:a.nmin]]])
    for m in out['minimisers']:
        log(f"   min at ({m['pose'][0]:.9f}, {m['pose'][1]:.9f}, {m['theta_deg']:.6f} deg): "
            f"capture {m['capture']:.6f} = points {m['points']:.6f} ({m['n_points']}) + cliques "
            f"{m['cliques']:.6f} ({m['n_cliques']}) + polygons {m['pgons']:.6f} ({m['n_pgons']})")
    o = np.argsort(V)[:a.nworst]
    np.savetxt(os.path.join(RUNS, f'hc_{a.TAG}_{a.rule}_fine_worst.txt'), P[o],
               header='cx cy theta (the worst poses after the knife-edge refinement)')
    json.dump(out, open(os.path.join(RUNS, f'hc_{a.TAG}_{a.rule}_fine.json'), 'w'), indent=1)
    log(f'   wrote runs/hc_{a.TAG}_{a.rule}_fine.json')


# ====================================================================== soundness of `meet`
def pair_meets(p, Q, tol=0.0):
    """float closed-intersection margin of the pose `p = (cx, cy, th)` against the poses `Q`
    (the 4-axis separating-axis test of `cliquelever.Poses.meet_margin`): >= -tol iff they meet"""
    Q = np.asarray(Q, dtype=float).reshape(-1, 3)
    ci, si = math.cos(p[2]), math.sin(p[2])
    cj, sj = np.cos(Q[:, 2]), np.sin(Q[:, 2])
    ddx = Q[:, 0] - p[0]
    ddy = Q[:, 1] - p[1]
    cr = ci * cj + si * sj
    sr = ci * sj - si * cj
    w = 0.5 * (np.abs(cr) + np.abs(sr))
    m1 = 0.5 + w - np.abs(ddx * ci + ddy * si)
    m2 = 0.5 + w - np.abs(-ddx * si + ddy * ci)
    m3 = 0.5 + w - np.abs(ddx * cj + ddy * sj)
    m4 = 0.5 + w - np.abs(-ddx * sj + ddy * cj)
    return np.minimum(np.minimum(m1, m2), np.minimum(m3, m4)) >= -tol


def one_row_cover(cover, ci):
    sub = Cover.__new__(Cover)
    sub.__dict__.update(cover.__dict__)
    sub.cz = cover.cz[ci:ci + 1]
    sub.cmem = cover.cmem[ci:ci + 1]
    sub.ccore = cover.ccore[ci:ci + 1]
    sub.n_cliques = 1
    sub.pz = np.zeros(0)
    sub.ppts = []
    sub.n_pgons = 0
    return sub


def cmd_sound(a):
    """is the `meet` rule a valid COVER rule?  For each clique row with dual, the credited set is
    `{S : S closed-meets every member of K}`; the accounting `sum_i capture(S_i) <= Theta` over a
    packing needs that set to be a CLIQUE.  Search it for two disjoint admissible squares."""
    t0 = time.time()

    def log(m):
        print(f'[{time.time() - t0:7.0f}s] {m}', flush=True)

    d = np.load(os.path.join(RUNS, f'hc_{a.TAG}_dual.npz'))
    cov = Cover(d, a.geo_tol, log=log)
    log(f'# honestcost sound {a.TAG}')
    log(cov.describe())
    L = lattice(a.pitch, a.dth)
    rng = np.random.default_rng(20260912)
    bad, badmass, tested = 0, 0.0, 0
    worst = []
    order = np.argsort(-cov.cz)
    for ci in order[:a.ncliques]:
        sub = one_row_cover(cov, int(ci))
        v = sub.clique_capture(L, 'meet')
        cred = L[v > 0]
        if len(cred) > a.sample:
            cred = cred[rng.choice(len(cred), a.sample, replace=False)]
        tested += 1
        hit = None
        for i in range(len(cred)):
            m = pair_meets(cred[i], cred)
            j = np.nonzero(~m)[0]
            if len(j):
                hit = (cred[i], cred[int(j[0])])
                break
        if hit is not None:
            bad += 1
            badmass += float(cov.cz[ci])
            if len(worst) < 6:
                worst.append(dict(clique=int(ci), z=float(cov.cz[ci]),
                                  size=int(len(cov.cmem[ci])),
                                  credited_lattice=int((v > 0).sum()),
                                  disjoint_pair=[[float(u) for u in hit[0]],
                                                 [float(u) for u in hit[1]]]))
        log(f'   clique {int(ci)} z={cov.cz[ci]:.6f} size={len(cov.cmem[ci])} '
            f'credited {int((v > 0).sum())} lattice poses -> '
            f'{"DISJOINT PAIR FOUND (rule unsound for this row)" if hit is not None else "clique"}')
    out = dict(tag=a.TAG, tested=tested, unsound=bad, unsound_dual=badmass,
               clique_total=cov.clique_total, theta=cov.theta, examples=worst)
    log(f'   {bad} of {tested} clique rows tested are NOT cliques under the `meet` rule, '
        f'carrying dual {badmass:.6f} of {cov.clique_total:.6f}')
    json.dump(out, open(os.path.join(RUNS, f'hc_{a.TAG}_sound.json'), 'w'), indent=1)


# ====================================================================== report
def anatomy(cover, pose, rule, threads=1):
    """what one pose captures, per column type"""
    P = np.array(pose, dtype=float).reshape(1, 3)
    tot, pt, cq, pg = cover.capture(P, rule, threads, parts=True)
    inp = []
    for i in range(cover.n_atoms):
        dx, dy = cover.ax[i] - P[0, 0], cover.ay[i] - P[0, 1]
        co, si = math.cos(P[0, 2]), math.sin(P[0, 2])
        if max(abs(dx * co + dy * si), abs(-dx * si + dy * co)) <= 0.5 + 1e-8:
            inp.append((float(cover.ax[i]), float(cover.ay[i]), float(cover.aw[i])))
    inp.sort(key=lambda r: -r[2])
    cqs = []
    for ci in range(cover.n_cliques):
        if one_row_cover(cover, ci).clique_capture(P, rule)[0] > 0:
            cqs.append((int(ci), float(cover.cz[ci]), int(len(cover.cmem[ci])),
                        bool(cover.ccore[ci])))
    pgs = []
    for gi in range(cover.n_pgons):
        sub = Cover.__new__(Cover)
        sub.__dict__.update(cover.__dict__)
        sub.pz = cover.pz[gi:gi + 1]
        sub.ppts = cover.ppts[gi:gi + 1]
        sub.n_pgons = 1
        if sub.pgon_capture(P)[0] > 0:
            pgs.append((int(gi), float(cover.pz[gi]), int(len(cover.ppts[gi]))))
    return dict(pose=[float(v) for v in P[0]], theta_deg=math.degrees(P[0, 2]) % 90.0,
                capture=float(tot[0]), points=float(pt[0]), cliques=float(cq[0]),
                pgons=float(pg[0]), n_points=len(inp), top_points=inp[:8],
                n_cliques=len(cqs), clique_rows=cqs[:10], n_pgons=len(pgs), pgon_rows=pgs[:10])


# ====================================================================== commands
def cmd_dual(a):
    def log(m):
        print(m, flush=True)
    log(f'# honestcost dual {a.TAG}')
    ps, lever = load_state(a.TAG, a.base, log, threads=a.threads, lp_tlim=0.0)
    x, lp, rc = solve_state(lever, log, master=a.master)
    tag = a.out_tag or a.TAG
    p = os.path.join(RUNS, f'hc_{tag}_dual.npz')
    dump_dual(lever, x, lp, p)
    theta = float(lever.y.sum() + lever.z.sum()
                  + sum(g['rhs'] * lever.pgz[i] for i, g in enumerate(lever.pgons)))
    log(f'   Theta = {theta:.9f}, LP = {lp:.9f}, gap = {theta - lp:+.3e}')
    log(f'   wrote {p}')
    # ---- self-test: every loaded pose must be captured >= 1 under the `meet` rule
    d = np.load(p, allow_pickle=True)
    cov = Cover(d, a.geo_tol, log=log)
    log('   ' + cov.describe())
    P = np.column_stack([ps.F[:, 0], ps.F[:, 1], np.arctan2(ps.F[:, 3], ps.F[:, 2])])
    for rule in ('meet', 'core'):
        c = cov.capture(P, rule, a.threads)
        log(f'   SELFTEST {rule}: min capture over the {len(P)} poses of P = {c.min():.9f} '
            f'({int((c < 1 - 1e-9).sum())} poses under 1 - 1e-9)')
    json.dump(dict(tag=tag, lp=lp, theta=theta, rc_max=float(rc.max()),
                   atoms=cov.n_atoms, cliques=cov.n_cliques, pgons=cov.n_pgons,
                   point_total=cov.point_total, clique_total=cov.clique_total,
                   pgon_total=cov.pgon_total, core_empty=cov.n_core_empty,
                   core_empty_mass=cov.core_empty_mass),
              open(os.path.join(RUNS, f'hc_{tag}_dual.json'), 'w'), indent=1)


def cmd_scan(a):
    t0 = time.time()

    def log(m):
        print(f'[{time.time() - t0:7.0f}s] {m}', flush=True)

    d = np.load(os.path.join(RUNS, f'hc_{a.TAG}_dual.npz'), allow_pickle=True)
    cov = Cover(d, a.geo_tol, log=log)
    log(f'# honestcost scan {a.TAG} rule={a.rule}')
    log(cov.describe())
    out = dict(tag=a.TAG, rule=a.rule, theta=cov.theta, lp=cov.lp, geo_tol=a.geo_tol,
               atoms=cov.n_atoms, cliques=cov.n_cliques, pgons=cov.n_pgons)

    def honest_of(v):
        return cov.theta / v if v > 0 else float('inf')

    def report(name, P, V):
        i = int(np.argmin(V))
        out[name] = dict(n=int(len(P)), min=float(V[i]),
                         honest=honest_of(float(V[i])), pose=[float(v) for v in P[i]])
        log(f'   {name}: {len(P)} poses, min capture {V[i]:.9f} -> honest '
            f'{honest_of(float(V[i])):.6f}')
        return i

    best_P, best_V = [], []

    # ---- (a) the 0.04 / 2.5 deg lattice
    L = lattice(a.pitch, a.dth)
    VL = cov.capture(L, a.rule, a.threads)
    report('lattice', L, VL)
    o = np.argsort(VL)[:a.nseed]
    best_P.append(L[o])
    best_V.append(VL[o])

    # ---- (a') family_rows.py item-3 families at pitch 0.004, --tilted
    if a.families and os.path.exists(a.families):
        FP = []
        for line in open(a.families):
            q = line.split()
            if len(q) >= 3 and not q[0].startswith('#'):
                FP.append((float(q[0]), float(q[1]), float(q[2])))
        FP = admissible_clamp(np.array(FP))
        FP[:, 2] = np.mod(FP[:, 2], math.pi / 2)
        FP = admissible_clamp(FP)
        VF = cov.capture(FP, a.rule, a.threads)
        report('families', FP, VF)
        o = np.argsort(VF)[:a.nseed]
        best_P.append(FP[o])
        best_V.append(VF[o])

    # ---- (a'') the poses of P itself and the arrangement/row points as centres
    PF = d['pf']
    PP = np.column_stack([PF[:, 0], PF[:, 1], np.arctan2(PF[:, 3], PF[:, 2])])
    VP = cov.capture(PP, a.rule, a.threads)
    report('poses_of_P', PP, VP)
    best_P.append(PP)
    best_V.append(VP)

    nth = max(int(round(90.0 / a.dth)), 1)
    RC = []
    rx, ry = d['rx'], d['ry']
    yv = d['y']
    hv = np.argsort(-yv)[:a.nrowseed]
    for i in hv:
        for j in range(nth):
            RC.append((rx[i], ry[i], math.radians(90.0 * j / nth)))
    RC = admissible_clamp(np.array(RC))
    VR = cov.capture(RC, a.rule, a.threads)
    report('row_centres', RC, VR)
    o = np.argsort(VR)[:a.nseed]
    best_P.append(RC[o])
    best_V.append(VR[o])

    # ---- (b) local pattern search from the worst of everything above
    seeds = np.concatenate(best_P)
    sv = np.concatenate(best_V)
    o = np.argsort(sv)[:a.ndescend]
    seeds = seeds[o]
    log(f'   descent from {len(seeds)} seeds (worst capture {sv[o].min():.9f})')
    DP, DV = pattern_search(cov, seeds, a.rule, a.threads, log=log)
    report('descent', DP, DV)

    # ---- (b') the knife-edge stage: capture is discontinuous at the heavy grid lines
    o = np.argsort(DV)[:a.nfine]
    FP2, FV2 = eps_refine(cov, DP[o], a.rule, a.threads, log=log)
    report('knife_edge', FP2, FV2)
    DP = np.concatenate([DP, FP2])
    DV = np.concatenate([DV, FV2])

    # ---- (c) a dense confirmation scan
    if a.dense > 0:
        DL = lattice(a.dense_pitch, a.dense_dth)
        log(f'   dense scan: {len(DL)} poses (pitch {a.dense_pitch}, dtheta {a.dense_dth})')
        VD = np.empty(len(DL))
        CH = 400000
        for lo in range(0, len(DL), CH):
            VD[lo:lo + CH] = cov.capture(DL[lo:lo + CH], a.rule, a.threads)
        report('dense', DL, VD)
        o2 = np.argsort(VD)[:a.ndescend]
        DP2, DV2 = pattern_search(cov, DL[o2], a.rule, a.threads, log=log)
        DP2, DV2 = eps_refine(cov, DP2[np.argsort(DV2)[:a.nfine]], a.rule, a.threads, log=log)
        report('dense_descent', DP2, DV2)
        DP = np.concatenate([DP, DP2])
        DV = np.concatenate([DV, DV2])

    # ---- the verdict
    i = int(np.argmin(DV))
    mn = float(DV[i])
    out['min_capture'] = mn
    out['honest'] = honest_of(mn)
    log(f'   MIN capture {mn:.9f} at {DP[i].tolist()} -> HONEST = {honest_of(mn):.6f}')
    # the distinct minimisers, with their anatomy
    o = np.argsort(DV)
    keep = []
    for j in o:
        if DV[j] > mn + a.near:
            break
        p = DP[j]
        if all(abs(p[0] - q[0]) + abs(p[1] - q[1]) > 0.05
               or abs(math.degrees(p[2] - q[2])) > 2.0 for q in keep):
            keep.append(p)
        if len(keep) >= a.nmin:
            break
    out['minimisers'] = [anatomy(cov, p, a.rule) for p in keep]
    for m in out['minimisers']:
        log(f"   min at ({m['pose'][0]:.5f}, {m['pose'][1]:.5f}, {m['theta_deg']:.3f} deg): "
            f"capture {m['capture']:.6f} = points {m['points']:.6f} ({m['n_points']}) + cliques "
            f"{m['cliques']:.6f} ({m['n_cliques']}) + polygons {m['pgons']:.6f} ({m['n_pgons']})")
    # the worst poses, for --repair
    o = np.argsort(DV)[:a.nworst]
    np.savetxt(os.path.join(RUNS, f'hc_{a.TAG}_{a.rule}_worst.txt'), DP[o],
               header='cx cy theta (the worst poses found: capture ascending)')
    json.dump(out, open(os.path.join(RUNS, f'hc_{a.TAG}_{a.rule}_scan.json'), 'w'), indent=1)
    log(f'   wrote runs/hc_{a.TAG}_{a.rule}_scan.json')


def cmd_repair(a):
    """step 4: add the worst poses found as coverage rows (their centres) and as columns, re-solve,
    re-scan.  One round only."""
    t0 = time.time()

    def log(m):
        print(f'[{time.time() - t0:7.0f}s] {m}', flush=True)

    ps, lever = load_state(a.TAG, a.base, log, threads=a.threads, lp_tlim=0.0)
    W = np.loadtxt(os.path.join(RUNS, f'hc_{a.TAG}_{a.rule}_worst.txt')).reshape(-1, 3)[:a.nadd]
    log(f'   {len(W)} worst poses read')
    idx = ps.add_float_poses([tuple(w) for w in W], lever.a.Q, lever.a.Dc)
    ps.finish()
    log(f'   +{len(idx)} new poses -> {ps.n} poses, {ps.ncol} columns')
    A2, _ = ps.contains_rows(lever.rows)
    lever.A = A2
    lever._Acsc = None
    K = []
    nmem = 0
    for c in lever.cliques:
        mem = c['members']
        add = ps.joins_all(idx, mem)
        kept = []
        for j in add:
            j = int(j)
            if not kept or ps.meets_all(j, np.array(kept)).all():
                kept.append(j)
        c['members'] = mem + kept
        nmem += len(kept)
        cols = [cc for i in c['members'] for cc in ps.cols_of[i]]
        K.append(sp.csr_matrix((np.ones(len(cols)), (np.zeros(len(cols), dtype=int), cols)),
                               shape=(1, ps.ncol)))
    lever.ckeys = set(frozenset(c['members']) for c in lever.cliques)
    lever.K = sp.vstack(K, format='csr') if K else sp.csr_matrix((0, ps.ncol))
    lever._Kcsc = None
    log(f'   clique rows extended by {nmem} memberships')
    if lever.pgons:
        log(f'   polygon rows re-derived, +{rf.regrow(lever)} memberships')
    lever.hi.h = None
    x, lp, rc = solve_state(lever, log)
    p = os.path.join(RUNS, f'hc_{a.TAG}rep_dual.npz')
    dump_dual(lever, x, lp, p)
    log(f'   wrote {p} (LP {lp:.9f}); now run: scan {a.TAG}rep')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)

    c = sub.add_parser('dual')
    c.add_argument('TAG')
    c.add_argument('--base', default=RUNS)
    c.add_argument('--threads', type=int, default=4)
    c.add_argument('--geo-tol', type=float, default=1e-9)
    c.add_argument('--master', action='store_true',
                   help='close the restricted master exactly instead of one full-LP solve '
                        '(same optimum, far cheaper on the 14k-column states)')
    c.add_argument('--out-tag', default=None, help='write hc_OUTTAG_dual.npz instead')

    c = sub.add_parser('scan')
    c.add_argument('TAG')
    c.add_argument('--rule', default='meet', choices=('meet', 'core'))
    c.add_argument('--threads', type=int, default=4)
    c.add_argument('--geo-tol', type=float, default=1e-9)
    c.add_argument('--pitch', type=float, default=0.04)
    c.add_argument('--dth', type=float, default=2.5)
    c.add_argument('--families', default=os.path.join(RUNS, 'hc_families.txt'))
    c.add_argument('--nseed', type=int, default=400)
    c.add_argument('--nrowseed', type=int, default=400)
    c.add_argument('--ndescend', type=int, default=400)
    c.add_argument('--nfine', type=int, default=100)
    c.add_argument('--dense', type=int, default=1)
    c.add_argument('--dense-pitch', type=float, default=0.02)
    c.add_argument('--dense-dth', type=float, default=1.0)
    c.add_argument('--near', type=float, default=0.02)
    c.add_argument('--nmin', type=int, default=8)
    c.add_argument('--nworst', type=int, default=400)

    c = sub.add_parser('fine')
    c.add_argument('TAG')
    c.add_argument('--rule', default='meet', choices=('meet', 'core'))
    c.add_argument('--threads', type=int, default=4)
    c.add_argument('--geo-tol', type=float, default=1e-9)
    c.add_argument('--families', default=os.path.join(RUNS, 'hc_families.txt'))
    c.add_argument('--nseed', type=int, default=100)
    c.add_argument('--rounds', type=int, default=2)
    c.add_argument('--nmin', type=int, default=6)
    c.add_argument('--nworst', type=int, default=400)

    c = sub.add_parser('sound')
    c.add_argument('TAG')
    c.add_argument('--threads', type=int, default=4)
    c.add_argument('--geo-tol', type=float, default=1e-9)
    c.add_argument('--pitch', type=float, default=0.02)
    c.add_argument('--dth', type=float, default=2.5)
    c.add_argument('--ncliques', type=int, default=40)
    c.add_argument('--sample', type=int, default=400)

    c = sub.add_parser('repair')
    c.add_argument('TAG')
    c.add_argument('--rule', default='meet', choices=('meet', 'core'))
    c.add_argument('--base', default=RUNS)
    c.add_argument('--threads', type=int, default=4)
    c.add_argument('--nadd', type=int, default=400)

    a = ap.parse_args()
    sys.set_int_max_str_digits(0)
    {'dual': cmd_dual, 'scan': cmd_scan, 'fine': cmd_fine, 'sound': cmd_sound,
     'repair': cmd_repair}[a.cmd](a)


if __name__ == '__main__':
    main()
