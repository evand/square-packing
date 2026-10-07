#!/usr/bin/env python3
"""farfield_strong.py -- the near-axis mass cap `mu(tilt < eps) <= K` in the STRONGEST certifiable
`t = T` family: points (every arrangement vertex) + odd polygons + region equalities (corner
`k = 4`) + the chord rows.  Task `tasks/farfield-strong/README.md`; doc `search/FARFIELD_STRONG.md`.

It **imports and never modifies** `search/cliquelever.py` (the family's LP: `Poses`, `Lever`, the
exact row generation and the exact finalisation), `search/leaf_ceiling.py` (exact arrangement,
regions), `search/rankfamily.py` (the odd-polygon rows) and `search/arch_farfield.py` (the exact
near-axis test).  The only new object is one extra LP row,

    sum_{poses of tilt < eps} mu  <=  K ,

plus the guards that keep the exact finalisation inside it.

DIRECTION, the whole discipline (search/ARCH_FARFIELD.md 1).  Every number this file prints is
**packing-side**: an exactly certified measure that satisfies every row of the family AND the cap
is feasible for the continuum problem, so its mass is a LOWER bound on the capped value.  A cell
`>= n` therefore refutes the cap at that `eps`; a cell `< n` proves NOTHING (it may be the pose
pool -- read the uncapped value of the same pool first, ARCH_FARFIELD.md 5).  No cover-side
(upper) bound is produced here: that needs a dual feasible for every pose of the continuum.

The near-axis test is exact and errs towards a STRONGER cap: `sin(tilt) = min(|cos|,|sin|)` in the
rational parameterisation `theta = 2 arctan(p/q)`, compared against `sin(eps)` rounded UP, so the
set called near-axis is a superset of the true one (`arch_farfield.near_axis_exact`).

Modes
    pool  OUT --t 4 --exact F... --float F... --rot A,B,... [--rot-src F...]
          build a sym-1 exact pose pool (`pose p q cx cy mass`) out of exact measure files
          (sym 1 or sym 8, the 8 dihedral images expanded), float support files, and rigid
          rotations of either by `alpha` degrees about the container centre (a rotation preserves
          coverage exactly and moves every angle by alpha: the far-field columns of
          `arch_farfield.py rotseed`).
    cap   TAG --pool P --t 4 --eps E --cap K [--uncapped] ...
          the capped LP of the strong family on that pool, row-generated on the EXACT arrangement
          of its own support until nothing is violated, then finalised exactly (masses rounded
          down; regions topped up only on exact coverage / polygon / chord / cap slack).
    check FILE --t 4 --eps E --cap K [--pgons J] [--corners 1111] [--chord]
          re-derive mass, max coverage M, near-axis mass, region and chord masses (and the polygon
          rows if a pgon file is given) from the measure file alone, in integers.
    table TAG...
          collect `runs/farfield_strong_TAG.json` into one table.
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

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import leaf_ceiling as lc                                            # noqa: E402
import cliquelever as cl                                             # noqa: E402
import rankfamily                                                    # noqa: E402
import arch_farfield as af                                           # noqa: E402

RUNS = os.path.join(REPO, 'runs')
LABS = cl.LABS


# ====================================================================== pool building
def _d4(lst, t):
    out = []
    for (p, q, cx, cy, m) in lst:
        for (x, y, pp, qq) in lc.d4_images(cx, cy, p, q, t):
            out.append((pp, qq, x, y, m / 8))
    return out


def _read_any(path, t):
    """-> [(p, q, cx, cy, mass)] exact, D4 images expanded when the file says sym 8, and float
    support files snapped."""
    with open(path) as f:
        head = f.read(4000)
    out = []
    if 'pose' in head and len(
            [l for l in head.splitlines() if l.startswith('pose')] or ['']) and \
            any(len(l.split()) == 6 for l in head.splitlines() if l.startswith('pose')):
        tt, sym, poses = lc.read_measure(path, t_override=t)
        for (p, q, cx, cy, m) in poses:
            if sym == 8:
                for (x, y, pp, qq) in lc.d4_images(cx, cy, p, q, t):
                    out.append((pp, qq, x, y, m / 8))
            else:
                out.append((p, q, cx, cy, m))
    else:
        _, src = lc.read_float_poses(path)
        for (cx, cy, th, mu) in src:
            p, q, x, y = lc.snap_pose(cx, cy, th, t, 10 ** 6, 10 ** 8)
            out.append((p, q, x, y, Fr(0)))
    return out


def _rotate(poses, t, alpha_deg):
    """rigid rotation of the whole container about its centre: coverage is preserved exactly and
    every angle moves by alpha.  Images that leave [0,t]^2 are dropped."""
    c = float(t) / 2
    a = math.radians(alpha_deg)
    ca, sa = math.cos(a), math.sin(a)
    out = []
    for (p, q, cx, cy, m) in poses:
        x, y = float(cx) - c, float(cy) - c
        nx, ny = c + ca * x - sa * y, c + sa * x + ca * y
        th = 2 * math.atan2(p, q) + a
        w = (abs(math.cos(th)) + abs(math.sin(th))) / 2
        if nx < w - 1e-12 or nx > float(t) - w + 1e-12 or ny < w - 1e-12 or ny > float(t) - w + 1e-12:
            continue
        pp, qq, X, Y = lc.snap_pose(nx, ny, th, t, 10 ** 6, 10 ** 8)
        out.append((pp, qq, X, Y, Fr(0)))
    return out


def cmd_pool(a):
    t = Fr(a.t)
    acc, seen = [], set()

    def push(lst):
        n = 0
        for (p, q, cx, cy, m) in lst:
            k = (p, q, cx, cy)
            if k in seen:
                continue
            seen.add(k)
            acc.append((p, q, cx, cy, m))
            n += 1
        return n

    for f in a.exact + a.float:
        n = push(_read_any(f, t))
        print(f'   +{n} poses from {os.path.basename(f)}  (pool {len(acc)})')
    for f in a.float8:
        n = push(_d4(_read_any(f, t), t))
        print(f'   +{n} poses from {os.path.basename(f)} (D4-expanded)  (pool {len(acc)})')
    rot_src = []
    for f in (a.rot_src or (a.exact + a.float)):
        rot_src += _read_any(f, t)
    for f in a.rot_src8:
        rot_src += _d4(_read_any(f, t), t)
    if a.rot:
        for al in [float(v) for v in a.rot.split(',')]:
            n = push(_rotate(rot_src, t, al))
            print(f'   +{n} poses rotated by {al:+g} deg  (pool {len(acc)})')
    if a.tilt:
        # band-edge columns: every source pose KEPT WHERE IT IS but re-angled to tilt exactly A.
        # The cap is evaded by tilting just past eps, and this is that move applied pose by pose
        # rather than to the container as a whole (which drops poses out of it).
        for al in [float(v) for v in a.tilt.split(',')]:
            out = []
            for (p, q, cx, cy, m) in rot_src:
                th = math.radians(al)
                w = (abs(math.cos(th)) + abs(math.sin(th))) / 2
                x = min(max(float(cx), w), float(t) - w)
                y = min(max(float(cy), w), float(t) - w)
                pp, qq, X, Y = lc.snap_pose(x, y, th, t, 10 ** 6, 10 ** 8)
                out.append((pp, qq, X, Y, Fr(0)))
            n = push(out)
            print(f'   +{n} poses re-angled to {al:+g} deg  (pool {len(acc)})')
    with open(a.OUT, 'w') as f:
        f.write(f'# farfield_strong pool; t = {t}; sym 1; {len(acc)} poses\n')
        f.write('# pose p q cx cy mass : theta = 2 arctan(p/q); mass = the source mass (seed only)\n')
        for (p, q, cx, cy, m) in acc:
            f.write(f'pose {p} {q} {cx} {cy} {Fr(m)}\n')
    print(f'{len(acc)} poses -> {a.OUT}')


# ====================================================================== the capped LP
class CapLever(cl.Lever):
    """cliquelever.Lever + the row `sum_{tilt < eps} mu <= K`.  Nothing else changes."""

    def __init__(self, ps, args, log):
        super().__init__(ps, args, log)
        self.eps = 0.0
        self.capK = None
        self.se = af.sin_eps_rational(1.0, up=True)
        self.lam_cap = 0.0
        self.chord_rhs = float(args.chord_rhs)
        self.cell = ''
        self._near_pose = None
        self._near_col = None

    def set_cell(self, eps, cap, cell):
        """switch to the cell (eps, cap); the pose set and the coverage / polygon / region / chord
        rows are shared across cells (rows are only ever added, and every row is valid for every
        cell), so only the cap row and the HiGHS model change."""
        self.eps = float(eps)
        self.capK = None if cap is None else float(cap)
        self.se = af.sin_eps_rational(self.eps, up=True)
        self.cell = cell
        self._near_pose = None
        self._near_col = None
        self.hi.h = None                 # force a rebuild with the new cap row

    # -------------------------------------------------------------- the near-axis set (exact)
    def near_pose(self):
        if self._near_pose is None or len(self._near_pose) != self.ps.n:
            self._near_pose = np.array(
                [af.near_axis_exact(p, q, self.se) for (p, q, _, _) in self.ps.P], dtype=bool)
            self._near_col = None
        return self._near_pose

    def near_col(self):
        if self._near_col is None or len(self._near_col) != self.ps.ncol:
            self._near_col = self.near_pose()[self.ps.col_pose_a]
        return self._near_col

    # -------------------------------------------------------------- model
    def rebuild(self, cols=None):
        super().rebuild(cols)
        if self.chord and abs(self.chord_rhs - 3.0) > 1e-12:
            self.hi.bounds([('D', i) for i in range(4)], -self.hi.INF, self.chord_rhs)
        if self.capK is not None:
            v = self.near_col().astype(float)
            if cols is not None:
                v = v[cols]
            self.hi.append([('X', 0)], [-self.hi.INF], [self.capK],
                           sp.csr_matrix(v.reshape(1, -1)))

    def read_duals(self, d):
        super().read_duals(d)
        self.lam_cap = 0.0
        for mi, k in enumerate(self.hi.rk):
            if k[0] == 'X':
                self.lam_cap = max(-d[mi], 0.0)

    def price_columns(self):
        rc = super().price_columns()
        if self.capK is not None and self.lam_cap:
            rc = rc - self.lam_cap * self.near_col().astype(float)
        return rc

    # -------------------------------------------------------------- exact finalisation
    def finalize(self, x, mu, procs):
        """round the LP masses DOWN to 1e-9, scale down if the exact maximum coverage exceeds 1 or
        the cap is exceeded, then top the pinned regions back up ONLY on exact slack of every row
        of the family (coverage, polygons, chord) and of the cap.  Every comparison below is in
        Python integers / Fractions."""
        ps, DM = self.ps, self.DM
        nearp = self.near_pose()
        xc = [int(math.floor(v * DM)) for v in x]

        def pose_int():
            m = [0] * ps.n
            for c, v in enumerate(xc):
                m[ps.col_pose[c]] += v
            return m

        def arrangement(mi):
            sup = [i for i in range(ps.n) if mi[i] > 0]
            squares = [ps.sq[i][:8] + (mi[i],) for i in sup]
            V = lc.enumerate_vertices(squares, self.t, verbose=False, procs=procs)
            Inc = lc.incidences(squares, V, verbose=False, procs=procs)
            w = np.array([s[8] for s in squares], dtype=np.int64)
            return sup, squares, V, Inc, w, Inc.cov(w)

        mi = pose_int()
        sup, squares, V, Inc, w, cov = arrangement(mi)
        M = Fr(int(cov.max()), DM) if len(cov) else Fr(0)
        self.log(f'   final: rounded down: mass {sum(mi) / DM:.9f}, exact M = {float(M):.12f}')
        if M > 1:
            self.log('   final: M > 1 after rounding down; scaling by 1/M exactly')
            xc = [(v * M.denominator) // M.numerator for v in xc]
            mi = pose_int()
            sup, squares, V, Inc, w, cov = arrangement(mi)
            M = Fr(int(cov.max()), DM) if len(cov) else Fr(0)
        if self.capK is not None:
            NEAR = sum(mi[i] for i in range(ps.n) if nearp[i])
            capI = int(round(self.capK * DM))
            if NEAR > capI:
                self.log(f'   final: near-axis {NEAR / DM:.9f} > cap; scaling the near-axis '
                         f'columns down exactly')
                for c in range(ps.ncol):
                    if nearp[ps.col_pose[c]]:
                        xc[c] = (xc[c] * capI) // NEAR
                mi = pose_int()
                sup, squares, V, Inc, w, cov = arrangement(mi)
                M = Fr(int(cov.max()), DM) if len(cov) else Fr(0)

        # ---------------- region top-up on exact slack (coverage, polygons, chord, cap)
        pg = rankfamily.TopUpGuard(self, mi)
        pos = {i: k for k, i in enumerate(sup)}
        capI = None if self.capK is None else int(round(self.capK * DM))
        nearmass = sum(mi[i] for i in range(ps.n) if nearp[i])
        chordM = (self.chord_matrix_cached().toarray() > 0) if self.chord else None
        chordmass = ([int(sum(xc[c] for c in range(ps.ncol) if chordM[dd, c])) for dd in range(4)]
                     if self.chord else None)
        chordI = int(round(self.chord_rhs * DM))
        resid = {}
        if self.target:
            for lab, kv in sorted(self.target.items()):
                li = LABS.index(lab)
                cols = [c for c in range(ps.ncol) if ps.col_reg[c] == li]
                want = int(kv * DM)
                d = want - sum(xc[c] for c in cols)
                tries = 0
                while d > 0 and tries < 60:
                    tries += 1
                    bestc, bests = None, 0
                    for c in cols:
                        if xc[c] <= 0:
                            continue
                        i = ps.col_pose[c]
                        if i not in pos:
                            continue
                        k = pos[i]
                        s = int(DM - cov[Inc.cols[k]].max()) if len(Inc.cols[k]) else DM
                        s = min(s, pg.slack(i))
                        if capI is not None and nearp[i]:
                            s = min(s, capI - nearmass)
                        if self.chord:
                            for dd in range(4):
                                if chordM[dd, c]:
                                    s = min(s, chordI - chordmass[dd])
                        if s > bests:
                            bestc, bests = c, s
                        if bests >= d:
                            break
                    if bestc is None or bests <= 0:
                        break
                    add = min(d, bests)
                    xc[bestc] += add
                    d -= add
                    i = ps.col_pose[bestc]
                    pg.credit(i, add)
                    if nearp[i]:
                        nearmass += add
                    if self.chord:
                        for dd in range(4):
                            if chordM[dd, bestc]:
                                chordmass[dd] += add
                    k = pos[i]
                    w[k] += add
                    cov[Inc.cols[k]] += add
                resid[lab] = d / DM
                self.log(f'   final: region {lab} -> {Fr(sum(xc[c] for c in cols), DM)} '
                         f'(target {kv}, residual {d / DM:+.3e})')

        mi = pose_int()
        sup, squares, V, Inc, w, cov = arrangement(mi)
        M = Fr(int(cov.max()), DM) if len(cov) else Fr(0)
        mass = Fr(sum(mi), DM)
        nearm = Fr(sum(mi[i] for i in range(ps.n) if nearp[i]), DM)
        meta = [(Fr(s[0], s[2]), Fr(s[1], s[2]), s[8]) for s in squares]
        ok, assigned, amb, note = lc.region_report(meta, ps.boxes, self.target, DM)
        pgok, pgworst, pgrep = rankfamily.check_final(self, mi)
        ch = [Fr(0)] * 4
        if self.chord:
            t, r = ps.t, ps.r
            for i in range(ps.n):
                if mi[i] <= 0:
                    continue
                p, q, cx, cy = ps.P[i]
                for dd, cond in enumerate((cy <= r, cy >= t - r, cx <= r, cx >= t - r)):
                    if cond:
                        ch[dd] += Fr(mi[i], DM)
        capok = (self.capK is None or nearm <= Fr(capI, DM))
        chok = (not self.chord or all(v <= Fr(chordI, DM) for v in ch))
        self.log(f'   FINAL EXACT: mass = {mass} = {float(mass):.9f}; M = {float(M):.12f} '
                 f'{"OK" if M <= 1 else "FAIL"}; near-axis(eps={self.eps}) = {float(nearm):.9f} '
                 f'{"OK" if capok else "FAIL"} (cap {self.capK}); regions '
                 f'{"OK" if ok else "FAIL"} ({note}); chord {[float(v) for v in ch]} '
                 f'{"OK" if chok else "FAIL"}; {len(self.pgons)} polygon rows '
                 f'{"OK" if pgok else "FAIL"} (worst {pgworst:+.3e}); {len(sup)} poses')
        res = dict(mass=str(mass), mass_float=float(mass), M=str(M), M_float=float(M),
                   near=str(nearm), near_float=float(nearm), cap=self.capK, eps=self.eps,
                   sin_eps=str(self.se), regions_ok=bool(ok), region_resid=resid,
                   regions={lab: str(Fr(assigned.get(lab, 0), DM)) for lab in LABS},
                   chord=[str(v) for v in ch], chord_ok=bool(chok), chord_rhs=self.chord_rhs,
                   pgons=len(self.pgons), pgons_ok=bool(pgok), pgons_worst=pgworst,
                   cap_ok=bool(capok), poses=len(sup), vertices=len(V),
                   certified=bool(M <= 1 and ok and pgok and capok and chok))
        path = os.path.join(HERE, f'farfield_strong_{self.a.TAG}{self.cell}_exact.txt')
        self.write_measure(path, mi, note=(
            f'farfield-strong: corner k=4 + polygons + regions + chord, cap mu(tilt < {self.eps} '
            f'deg) <= {self.capK}; M = {M}, near-axis = {nearm}'))
        res['file'] = path
        # the heaviest poses, for the note
        order = sorted(range(ps.n), key=lambda i: -mi[i])[:20]
        res['heaviest'] = [dict(cx=float(ps.P[i][2]), cy=float(ps.P[i][3]),
                                theta=math.degrees(2 * math.atan2(ps.P[i][0], ps.P[i][1])),
                                mass=mi[i] / DM, near=bool(nearp[i]))
                           for i in order if mi[i] > 0]
        return res


def cmd_cap(a):
    os.makedirs(RUNS, exist_ok=True)
    logf = open(os.path.join(RUNS, f'farfield_strong_{a.TAG}.log'), 'a')

    def log(m):
        print(m, flush=True)
        logf.write(m + '\n')
        logf.flush()

    T0 = time.time()
    a.corners = a.corners if a.corners != '' else None
    a.patterns = a.patterns if a.patterns != '' else None
    a.pent = [int(v) for v in a.pent.split(',') if v.strip()] if a.pent else []
    ps = cl.Poses(Fr(a.t), Fr(a.r), a.bnd_delta)
    for f in a.pool:
        n = ps.add_exact_file(f)
        log(f'   +{n} exact poses from {os.path.basename(f)}')
    ps.finish()
    lever = CapLever(ps, a, log)
    log(f'# farfield_strong {a.TAG}: t={ps.t} corners={a.corners} patterns={a.patterns} '
        f'chord={a.chord}(rhs {a.chord_rhs}): {ps.n} poses, {ps.ncol} columns; args={vars(a)}')
    if a.row_pitch > 0:
        log(f'   +{lever.add_grid_rows(a.row_pitch)} grid rows (pitch {a.row_pitch})')
    for f in a.resume_rows:
        log(f'   +{lever.resume_rows(f)} exact rows resumed from {os.path.basename(f)}')
    for f in a.resume_pgons:
        nadd, nbad = rankfamily.resume(lever, f)
        log(f'   +{nadd} polygon rows from {os.path.basename(f)} (membership re-derived over this '
            f'pose set; {nbad} dropped)')
    for f in a.seed_rows:
        sp_ = _read_any(f, Fr(a.t))
        DM = 10 ** 9
        sq = [lc.make_square(cx, cy, p, q, Fr(a.t), max(int(m * DM), 1), tag=i)
              for i, (p, q, cx, cy, m) in enumerate(sp_)]
        V = lc.enumerate_vertices(sq, Fr(a.t), verbose=False, procs=a.procs)
        Inc = lc.incidences(sq, V, verbose=False, procs=a.procs)
        cov = Inc.cov(np.array([s[8] for s in sq], dtype=np.int64))
        o = np.argsort(-cov)[:a.rows0]
        log(f'   +{lever.add_rows([V[i] for i in o])} rows at the heaviest arrangement vertices '
            f'of {os.path.basename(f)} ({len(sq)} squares, {len(V)} vertices)')
    log(f'   {len(lever.rows)} rows, {time.time() - T0:.0f}s')

    cells = []
    if a.uncapped:
        cells.append((0.0, None))
    for e in a.eps:
        for k in a.cap:
            cells.append((e, k))
    if a.recheck:
        cells.append((0.0, None))
    out = []
    jpath = os.path.join(RUNS, f'farfield_strong_{a.TAG}.json')
    seen_cells = set()
    for (eps, cap) in cells:
        name = ('_uncapped' if cap is None
                else f'_e{str(eps).replace(".", "p")}_c{str(cap).replace(".", "p")}')
        while name in seen_cells:
            name += 'b'
        seen_cells.add(name)
        lever.set_cell(eps, cap, name)
        nn = int(lever.near_pose().sum())
        log(f'=== cell{name}: eps = {eps} deg (sin(eps) <= {lever.se}; the near-axis set is a '
            f'SUPERSET, so the cap is if anything too strong), cap = {cap}; {nn} of {ps.n} poses '
            f'near-axis')
        hist = []
        x = mu = None
        Tc = time.time()
        for it in range(a.iters):
            sol = lever.solve()
            if sol is None:
                log('   LP failed'); break
            x, val = sol
            mu = lever.pose_mass(x)
            sup, squares, V, Inc, M, cov, w = lever.certify(mu, a.procs)
            bad = np.nonzero(cov > lever.DM * (1.0 + 1e-9))[0]
            nrow = 0
            if len(bad):
                o = bad[np.argsort(-cov[bad])[:a.row_cap]]
                nrow = lever.add_rows([V[i] for i in o])
            npg = 0
            if a.pent:
                npg, pbest, pg_secs, _nf = rankfamily.separate(
                    lever, sup, squares, V, Inc, w, mu, it)
            R = lever.region_split(x)
            near = float(x[lever.near_col()].sum())
            hist.append(dict(it=it, LP=float(val), M=M, rows=len(lever.rows), bad=int(len(bad)),
                             near=near, lam_cap=lever.lam_cap, interior=float(R[12]),
                             corners=R[:4].tolist(), chord_dual=lever.chord_dual.tolist(),
                             secs=time.time() - Tc))
            log(f'   [{it}] LP={val:.9f} M={M:.9f} rows={len(lever.rows)}(+{nrow}) '
                f'pg={len(lever.pgons)}(+{npg}) sup={len(sup)} verts={len(V)} bad={len(bad)} '
                f'near={near:.6f} lam_cap={lever.lam_cap:.6f} int={R[12]:.6f} '
                f'corners={np.round(R[:4], 6).tolist()} ({time.time() - Tc:.0f}s)')
            if nrow == 0 and npg == 0:
                log('   converged: every violated arrangement vertex of the support is already '
                    'a row (the exact optimum of this pose set)')
                break
            if time.time() - Tc > a.time:
                log('   time limit for this cell'); break
        # the raw LP support, for the polish below and for a second `cap` pass restricted to it
        spath = os.path.join(RUNS, f'farfield_strong_{a.TAG}{name}_lpsupport.txt')
        lever.write_measure(spath, [int(round(v * lever.DM)) for v in mu],
                            note=f'LP support of cell{name} (float masses rounded to 1e-9)')
        # ---- POLISH (`LEAF_CEILING.md` §5.4).  On a large pool the row generation sits in a
        # degenerate tail: the LP moves the mass onto new poses whose arrangement vertices are new
        # rows, so `M` oscillates just above 1 and never closes, and `finalize` would then have to
        # scale the measure down by 1/M.  Restricted to the poses of the support, the same exact-row
        # LP closes in a few rounds, and its value is a value those poses really attain.
        lev = lever
        if a.polish_iters > 0 and (mu > 1e-12).any():
            ps2 = cl.Poses(Fr(a.t), Fr(a.r), a.bnd_delta)
            for i in np.nonzero(mu > 1e-12)[0]:
                p, q, cx, cy = ps.P[int(i)]
                ps2.add(p, q, cx, cy, float(mu[int(i)]))
            ps2.finish()
            lev2 = CapLever(ps2, a, log)
            # carry EVERY polygon row across (resumed and separated alike), with its membership
            # re-derived from its exact rational anchors over the support pose set
            for c in lever.pgons:
                mem = rankfamily.row_members(ps2, c['pts'])
                if mem and not rankfamily.verify_row(ps2, c['pts'], mem):
                    lev2.add_pgon(mem, c['rhs'], c, c['pts'])
            lev2.add_rows(list(lever.rows))
            lev2.set_cell(eps, cap, name)
            log(f'   polish: {ps2.n} support poses, {len(lev2.rows)} rows, '
                f'{len(lev2.pgons)} polygon rows')
            nrows0 = len(lever.rows)
            for it in range(a.polish_iters):
                sol = lev2.solve()
                if sol is None:
                    log('   polish: LP failed'); break
                x2, val2 = sol
                mu2 = lev2.pose_mass(x2)
                s2, sq2, V2, I2, M2, cov2, w2 = lev2.certify(mu2, a.procs)
                bad2 = np.nonzero(cov2 > lev2.DM * (1.0 + 1e-9))[0]
                n2 = 0
                if len(bad2):
                    o2 = bad2[np.argsort(-cov2[bad2])[:a.row_cap]]
                    n2 = lev2.add_rows([V2[i] for i in o2])
                log(f'   polish[{it}] LP={val2:.9f} M={M2:.9f} rows={len(lev2.rows)}(+{n2}) '
                    f'sup={len(s2)} bad={len(bad2)} lam_cap={lev2.lam_cap:.6f} '
                    f'({time.time() - Tc:.0f}s)')
                if n2 == 0:
                    log('   polish: converged (exact optimum of the support pose set)')
                    break
            x, mu, lev = x2, mu2, lev2
            lever.add_rows(lev2.rows[nrows0:])      # the new rows help the next cell too
        fin = lev.finalize(x, mu, a.procs)
        fin.update(cell=name, iters=len(hist), LP=hist[-1]['LP'] if hist else None,
                   lam_cap=hist[-1]['lam_cap'] if hist else None, secs=time.time() - Tc)
        out.append(fin)
        json.dump(dict(tag=a.TAG, args=vars(a), poses=ps.n, cells=out,
                       seconds=time.time() - T0), open(jpath, 'w'), indent=1)
        if lever.pgons:
            json.dump([dict(anchors=c['anchors']) for c in lever.pgons],
                      open(os.path.join(RUNS, f'farfield_strong_{a.TAG}_pgons.json'), 'w'))
        log(f'RESULT tag={a.TAG}{name} eps={eps} cap={cap} LP={fin["LP"]:.9f} '
            f'exact={fin["mass_float"]:.9f} near={fin["near_float"]:.9f} '
            f'M={fin["M_float"]:.12f} certified={fin["certified"]} '
            f'secs={time.time() - Tc:.0f}')


# ====================================================================== independent check
def cmd_check(a):
    t = Fr(a.t)
    tt, sym, poses = lc.read_measure(a.FILE, t_override=t)
    assert sym == 1, 'this checker wants a sym-1 measure file'
    DM = 1
    for (_, _, _, _, m) in poses:
        DM = DM * m.denominator // math.gcd(DM, m.denominator)
    squares = []
    for i, (p, q, cx, cy, m) in enumerate(poses):
        squares.append(lc.make_square(cx, cy, p, q, t, int(m * DM), tag=i))
    V = lc.enumerate_vertices(squares, t, verbose=False, procs=a.procs)
    Inc = lc.incidences(squares, V, verbose=False, procs=a.procs)
    w = np.array([s[8] for s in squares], dtype=np.int64)
    cov = Inc.cov(w)
    M = Fr(int(cov.max()), DM)
    mass = Fr(sum(int(s[8]) for s in squares), DM)
    se = af.sin_eps_rational(a.eps, up=True)
    near = Fr(sum(int(s[8]) for s, (p, q, _, _, _) in zip(squares, poses)
                  if af.near_axis_exact(p, q, se)), DM)
    boxes = lc.region_boxes(t, Fr(a.r))
    meta = [(cx, cy, int(m * DM)) for (_, _, cx, cy, m) in poses]

    class _A:
        corners = a.corners
        patterns = a.patterns
        quad = None
    tgt = lc.region_targets(t, _A) if (a.corners or a.patterns) else None
    ok, assigned, amb, note = lc.region_report(meta, boxes, tgt, DM)
    ch = [Fr(0)] * 4
    r = Fr(a.r)
    for (p, q, cx, cy, m) in poses:
        for dd, cond in enumerate((cy <= r, cy >= t - r, cx <= r, cx >= t - r)):
            if cond:
                ch[dd] += m
    print(f'{a.FILE}')
    print(f'  {len(poses)} poses, {len(V)} arrangement vertices (whole container)')
    print(f'  mass       = {mass} = {float(mass):.9f}')
    print(f'  M          = {M} = {float(M):.12f}   {"OK <= 1" if M <= 1 else "FAIL > 1"}')
    print(f'  near-axis  = {near} = {float(near):.9f}   (tilt < {a.eps} deg, sin(eps) <= {se})'
          + (f'   {"OK" if near <= Fr(a.cap) else "FAIL"} <= cap {a.cap}' if a.cap is not None else ''))
    print(f'  regions    = {"OK" if ok else "FAIL"}  ({note})')
    print(f'  corner masses = {[str(Fr(assigned.get(f"C{i}", 0), DM)) for i in range(4)]}')
    print(f'  chord (wall strips) = {[float(v) for v in ch]}  <= {a.chord_rhs}: '
          f'{"OK" if all(v <= Fr(a.chord_rhs) for v in ch) else "FAIL"}')
    if a.pgons:
        ps = cl.Poses(t, Fr(a.r), 1e-6)
        idxs = [ps.add(p, q, cx, cy, float(m)) for (p, q, cx, cy, m) in poses]
        ps.finish()
        sa = argparse.Namespace(chord=False, threads=1, lp_tlim=0.0, corners=None, patterns=None,
                                quad=None, pent_seed=0, master=False, TAG='check')
        L = cl.Lever(ps, sa, lambda m: None)
        L.DM = DM
        nadd, nbad = rankfamily.resume(L, a.pgons)
        mi = [0] * ps.n
        for j, (_, _, _, _, m) in enumerate(poses):
            mi[idxs[j]] += int(m * DM)
        pgok, worst, rep = rankfamily.check_final(L, mi)
        print(f'  polygon rows: {nadd} re-derived from {os.path.basename(a.pgons)} '
              f'({nbad} dropped), {"OK" if pgok else "FAIL"}, worst violation {worst:+.9f}')


def cmd_table(a):
    rows = []
    for tag in a.TAG:
        p = os.path.join(RUNS, f'farfield_strong_{tag}.json')
        if not os.path.exists(p):
            print(f'   [{tag}: no json]'); continue
        d = json.load(open(p))
        for e in d.get('cells', []):
            rows.append((tag + e['cell'], e['eps'], e['cap'], e['mass_float'], e['near_float'],
                         e['M_float'], e['certified'], e['regions_ok'], e['pgons_ok'],
                         e['chord_ok']))
    print(f'{"tag":36s} {"eps":>6s} {"cap":>6s} {"exact mass":>14s} {"near":>10s} '
          f'{"M":>14s}  cert reg pg chord')
    for r in sorted(rows, key=lambda r: (r[1], -1 if r[2] is None else r[2])):
        print(f'{r[0]:36s} {r[1]:6.1f} {str(r[2]):>6s} {r[3]:14.9f} {r[4]:10.6f} '
              f'{r[5]:14.12f}  {int(r[6])}    {int(r[7])}   {int(r[8])}  {int(r[9])}')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)

    c = sub.add_parser('pool')
    c.add_argument('OUT')
    c.add_argument('--t', default='4')
    c.add_argument('--exact', action='append', default=[])
    c.add_argument('--float', action='append', default=[])
    c.add_argument('--float8', action='append', default=[],
                   help='float support file that is D4-symmetrised: expand the 8 images')
    c.add_argument('--rot', default='')
    c.add_argument('--rot-src', action='append', default=[])
    c.add_argument('--rot-src8', action='append', default=[])
    c.add_argument('--tilt', default='', help='band-edge columns: re-angle every source pose to these tilts')

    c = sub.add_parser('cap')
    c.add_argument('TAG')
    c.add_argument('--pool', action='append', required=True)
    c.add_argument('--t', default='4')
    c.add_argument('--r', default='1')
    c.add_argument('--eps', type=float, nargs='+', default=[])
    c.add_argument('--cap', type=float, nargs='+', default=[])
    c.add_argument('--uncapped', action='store_true')
    c.add_argument('--recheck', action='store_true',
                   help='re-run the uncapped cell LAST, on the final (larger) row set')
    c.add_argument('--corners', default='1111')
    c.add_argument('--patterns', default='')
    c.add_argument('--chord', action='store_true')
    c.add_argument('--chord-rhs', type=float, default=3.0)
    c.add_argument('--resume-rows', action='append', default=[])
    c.add_argument('--resume-pgons', action='append', default=[])
    c.add_argument('--seed-rows', action='append', default=[],
                   help='sym-1 exact measure file whose arrangement vertices seed the row set')
    c.add_argument('--row-pitch', type=float, default=0.0)
    c.add_argument('--rows0', type=int, default=30000)
    c.add_argument('--row-cap', type=int, default=20000)
    c.add_argument('--iters', type=int, default=60)
    c.add_argument('--polish-iters', type=int, default=25)
    c.add_argument('--time', type=float, default=7200)
    c.add_argument('--threads', type=int, default=2)
    c.add_argument('--procs', type=int, default=2)
    c.add_argument('--lp-tlim', type=float, default=0.0)
    c.add_argument('--bnd-delta', type=float, default=1e-6)
    c.add_argument('--clique-time', type=float, default=1.0)
    c.add_argument('--cq-want', type=int, default=0)
    c.add_argument('--cq-age', type=int, default=0)
    c.add_argument('--ktol', type=float, default=1e-6,
                   help='violation tolerance of the odd-polygon separation (rankfamily)')
    c.add_argument('--cq-top', type=int, default=0)
    c.add_argument('--pent', default='', help='comma-separated odd k (e.g. 5,7,9); empty = off')
    c.add_argument('--pent-want', type=int, default=6)
    c.add_argument('--pent-cands', type=int, default=240)
    c.add_argument('--pent-restarts', type=int, default=400)
    c.add_argument('--pent-time', type=float, default=90.0)
    c.add_argument('--anatomy', type=int, default=0)
    c.set_defaults(master=False, quad=None, pent_seed=20260921)

    c = sub.add_parser('check')
    c.add_argument('FILE')
    c.add_argument('--t', default='4')
    c.add_argument('--r', default='1')
    c.add_argument('--eps', type=float, required=True)
    c.add_argument('--cap', type=float, default=None)
    c.add_argument('--corners', default='1111')
    c.add_argument('--patterns', default=None)
    c.add_argument('--chord-rhs', type=float, default=3.0)
    c.add_argument('--pgons', default=None)
    c.add_argument('--procs', type=int, default=2)

    c = sub.add_parser('table')
    c.add_argument('TAG', nargs='+')

    a = ap.parse_args()
    dict(pool=cmd_pool, cap=cmd_cap, check=cmd_check, table=cmd_table)[a.cmd](a)


if __name__ == '__main__':
    main()
