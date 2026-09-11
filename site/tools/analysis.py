#!/usr/bin/env python3
"""Geometric analysis of a packing of unit squares in a square container.

Input: container side s and squares (cx, cy, theta_deg) in a y-up frame, as produced by parse_svg.
Everything here is float64 (inputs carry 30 digits, so doubles are good to ~1e-14 absolute);
contacts that are borderline are re-checked in mpmath at 50 digits.

Computed:
  angles        theta mod 90 per square, tilt in (-45,45], exact-equal groups
  contacts      square-square and square-wall contacts with type (corner-edge, edge-edge, corner-corner)
  rigidity      first-order: which squares can move at all (collective motion, LP) and alone (others fixed)
  free regions  sampled translation region (theta fixed) for squares free alone, plus pure-rotation range
  near misses   smallest positive gaps (the slivers), refined in mpmath when tiny
  symmetry      D4 elements preserving the packing as drawn
  category      trivial / diagonal (0 and 45 only) / tilted
"""
import numpy as np, itertools, math
from scipy.optimize import linprog, nnls
from scipy.spatial import cKDTree
from mpmath import mp, mpf, cos as mcos, sin as msin, radians as mrad
mp.dps = 50

CONTACT_TOL = 1e-8       # gap below this is a contact
EXACT_TOL = mpf('1e-20')  # mp gap below this is an exact (analytic) contact
NEAR_MAX = 0.25           # report near misses up to this size (the site lists the smallest few)
SLIDE_TOL = 1e-6          # a slide shorter than this is numerical noise, not a motion

def corners(cx, cy, th):
    c, s = math.cos(math.radians(th)), math.sin(math.radians(th))
    return np.array([[cx + c*dx - s*dy, cy + s*dx + c*dy] for dx, dy in ((-.5,-.5),(.5,-.5),(.5,.5),(-.5,.5))])

def normals(P):
    """Outward unit normals of a CCW convex polygon's edges (edge k from P[k] to P[k+1])."""
    E = np.roll(P, -1, axis=0) - P
    N = np.stack([E[:,1], -E[:,0]], axis=1)
    return N / np.linalg.norm(N, axis=1)[:, None]

def sat(A, B):
    """Return (gap, axis, owner) — largest separating gap over the edge normals of A (owner 0) and B (owner 1).
    Axis points from owner's polygon toward the other."""
    best = (-np.inf, None, None)
    for owner, P in ((0, A), (1, B)):
        Q = B if owner == 0 else A
        for n in normals(P):
            g = np.min(Q @ n) - np.max(P @ n)
            if g > best[0]: best = (g, n, owner)
    return best

def mp_gap(sq_i, sq_j):
    def crn(q):
        cx, cy, th = (mpf(v) for v in q); c, s = mcos(mrad(th)), msin(mrad(th)); h = mpf('0.5')
        return [(cx + c*dx - s*dy, cy + s*dx + c*dy) for dx, dy in ((-h,-h),(h,-h),(h,h),(-h,h))]
    A, B = crn(sq_i), crn(sq_j); best = mpf('-inf')
    for P, Q in ((A, B), (B, A)):
        for k in range(4):
            x1, y1 = P[k]; x2, y2 = P[(k+1) % 4]; nx, ny = y2 - y1, x1 - x2
            L = (nx*nx + ny*ny) ** mpf('0.5'); nx, ny = nx/L, ny/L
            g = min(nx*x + ny*y for x, y in Q) - max(nx*x + ny*y for x, y in P)
            if g > best: best = g
    return best

# ---------------------------------------------------------------------------
# Nonlinear verification of a proposed (first-order) motion.
#
# The LP below works with the linearised non-penetration conditions, so any motion that is
# *tangent* to a constraint -- a cluster turning about a corner contact, a square sliding while
# its neighbour rotates -- penetrates by O(eps^2) when the whole configuration is simply stepped
# along the LP direction, and stepping alone would reject nearly every real mechanism.  Instead we
# step and then project back onto the true feasible set: find a genuinely feasible configuration
# near q0 + eps*x using the exact (SAT) gaps, correcting by the smallest displacement that repairs
# them.  Only the squares the direction moves are variables; only pairs close enough to touch get
# constraints.  If a feasible configuration survives with some square displaced by >= MIN_DISP the
# motion is real; a direction blocked at second order cannot be repaired and is rejected.
# ---------------------------------------------------------------------------
_OFF = np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]])
PAIR_R = math.sqrt(2) + 0.1     # centre distance beyond which two unit squares cannot touch
MIN_DISP = 1e-4                 # a motion counts as finite if some square moves this much (units / rad)
FEAS_TOL = 1e-10                # accepted penetration of a verified configuration
CLUSTER_CAP = 800               # largest contact-graph neighbourhood searched for one square's motion
CLUSTER_BUDGET = 120000          # ... shrunk when many squares need one, to bound the total work

def near_pairs(C, R):
    """Index arrays (i, j), i < j, of centres within R of each other."""
    if len(C) < 2: return np.zeros(0, int), np.zeros(0, int)
    pr = cKDTree(C).query_pairs(R, output_type='ndarray')
    if not len(pr): return np.zeros(0, int), np.zeros(0, int)
    return pr[:, 0], pr[:, 1]

def _corners_all(state):
    """state (n,3) = (cx, cy, theta in RADIANS) -> corners (n,4,2), CCW."""
    c, s_ = np.cos(state[:, 2]), np.sin(state[:, 2])
    x = state[:, 0:1] + np.outer(c, _OFF[:, 0]) - np.outer(s_, _OFF[:, 1])
    y = state[:, 1:2] + np.outer(s_, _OFF[:, 0]) + np.outer(c, _OFF[:, 1])
    return np.stack([x, y], -1)

def _edge_normals_all(Cn):
    E = np.roll(Cn, -1, axis=1) - Cn
    N = np.stack([E[:, :, 1], -E[:, :, 0]], -1)
    return N / np.linalg.norm(N, axis=2)[:, :, None]

def _sat_batch(Cn, N, I, J):
    """SAT over pairs (I[k], J[k]).  Returns gap (K,), contact vertex q (K,2), axis a (K,2)
    (outward normal of the face owner, pointing towards the other body), and the owner /
    other body indices (K,).  gap = (q - p).a for the maximising axis."""
    K = len(I)
    g = np.empty((K, 8)); qi = np.empty((K, 8), int)
    for t, (Ow, Ot) in enumerate(((I, J), (J, I))):
        pa = np.einsum('kad,kcd->kac', N[Ow], Cn[Ow])
        pb = np.einsum('kad,kcd->kac', N[Ow], Cn[Ot])
        g[:, 4*t:4*t+4] = pb.min(2) - pa.max(2)
        qi[:, 4*t:4*t+4] = pb.argmin(2)
    k = g.argmax(1); r = np.arange(K)
    own = k >= 4; e = k % 4
    A = np.where(own, J, I); B = np.where(own, I, J)
    return g[r, k], Cn[B, qi[r, k]], N[A, e], A, B

_WALL_AX = np.array([0, 0, 1, 1]); _WALL_SG = np.array([1., -1., 1., -1.])

def _project_motion(state0, s, M, xdir, eps, pairs, walls, steps=1, iters=12):
    """Follow the motion xdir out to a displacement of eps by predictor-corrector continuation:
    advance eps/steps along xdir, then pull back onto the true feasible set with the smallest
    correction that repairs the contacts it has driven into (min ||u|| s.t. g + J u >= 0 over the
    active ones).  A motion merely tangent to its contacts is violated by O(step^2) and one
    correction fixes it, and the trajectory may bend away from xdir as it goes; a motion blocked at
    second order stalls at the first step it cannot repair.  Returns the furthest feasible
    displacement reached (m,3) -- all zeros if even the first step is blocked."""
    m = len(M); nsq = len(state0)
    I, J = pairs; wi, wc, ww = walls
    posarr = np.full(nsq, -1); posarr[M] = np.arange(m)
    cols3 = np.arange(3); nw = len(wi)
    wsg = _WALL_SG[ww]; wax = _WALL_AX[ww]; wpos = posarr[wi]
    wcon = np.where(wsg < 0, s, 0.0)
    wr = np.arange(nw)
    cr = lambda v, w: v[:, 0]*w[:, 1] - v[:, 1]*w[:, 0]
    margin = max(0.02, 4 * eps / steps)   # a corner further back than this cannot become binding
    def gaps_jac(d):
        """Signed gaps and their exact gradients at state0 + d.

        Each near pair is linearised on its currently separating axis, as one row per corner of the
        other body that could bind on it -- NOT as the single SAT gap.  The SAT gap is a min over
        those corners, so linearising it through the current argmin over-estimates it, and a
        correction can pivot the far corner of an edge-edge contact straight into the neighbour
        while the model says the gap is growing.  The owner's own extent along its face normal is
        exactly centre.a + 1/2, so each row is smooth: (corner - centre_A).a - 1/2 >= 0."""
        st = state0.copy(); st[M] = state0[M] + d
        Cn = _corners_all(st); N = _edge_normals_all(Cn)
        blocks = []
        if len(I):
            g, q, a, A, B = _sat_batch(Cn, N, I, J)
            off = np.einsum('kcd,kd->kc', Cn[B], a) - (np.einsum('kd,kd->k', st[A, :2], a) + 0.5)[:, None]
            kk, cc = np.where(off < off.min(1)[:, None] + margin)
            qv = Cn[B[kk], cc]; av = a[kk]; Av = A[kk]; Bv = B[kk]
            Jp = np.zeros((len(kk), 3*m)); rows = np.arange(len(kk))
            gA = np.stack([-av[:, 0], -av[:, 1], -cr(qv - st[Av, :2], av)], 1)
            gB = np.stack([av[:, 0], av[:, 1], cr(qv - st[Bv, :2], av)], 1)
            for bod, gg in ((Av, gA), (Bv, gB)):
                p = posarr[bod]; mk = p >= 0
                if mk.any(): Jp[rows[mk][:, None], 3*p[mk][:, None] + cols3] = gg[mk]
            blocks.append((off[kk, cc], Jp))
        if nw:
            co = Cn[wi, wc]
            gw = np.zeros((nw, 3))
            gw[wr, wax] = wsg
            gw[:, 2] = wsg * np.where(wax == 0, -(co[:, 1] - st[wi, 1]), co[:, 0] - st[wi, 0])
            Jw = np.zeros((nw, 3*m))
            Jw[wr[:, None], 3*wpos[:, None] + cols3] = gw
            blocks.append((wcon + wsg * co[wr, wax], Jw))
        if not blocks: return np.zeros(0), np.zeros((0, 3*m))
        return np.concatenate([b[0] for b in blocks]), np.vstack([b[1] for b in blocks])
    inc = eps / steps
    d = np.zeros((m, 3)); best = np.zeros((m, 3))
    for _ in range(steps):
        d = d + inc * xdir; prev = np.inf
        for _ in range(iters):
            val, Jc = gaps_jac(d)
            if not np.isfinite(val).all(): return best
            viol = -val.min()
            if viol <= 1e-14: break
            if viol > 0.9 * prev: break                   # correction has stalled: really blocked
            prev = viol
            # Active set: the contacts the motion rides on (gap ~ 0) plus the violated ones; the rest
            # have slack to spare over a correction this small.  Smallest u with G u >= -val over
            # them -- an inequality, so the repair may open a riding contact as well as close one.
            # Solved through its dual: u = G' lam with lam >= 0, a non-negative least squares fit of
            # G' lam to the equality (minimum-norm) solution.
            act = val < 0.05 * inc
            G = Jc[act]; b = -val[act]
            try:
                u = G.T @ np.linalg.lstsq(G @ G.T, b, rcond=None)[0]
                u = G.T @ nnls(G.T, u, maxiter=8*len(b))[0]
            except Exception:
                break
            big = np.abs(u).max()
            if big > inc: u *= inc / big                  # keep the correction local
            d = d + u.reshape(m, 3)
        if gaps_jac(d)[0].min() < -1e-14: break           # cannot get past this point
        best = d.copy()
    return best

def _still_ok(state0, s, M, d, pairs, walls, tol):
    """Exact check of the moved part: every pair involving a moved square, and its walls."""
    st = state0.copy(); st[M] = state0[M] + d
    Cn = _corners_all(st)
    if Cn[M][:, :, 0].min() < -tol or Cn[M][:, :, 1].min() < -tol: return False
    if Cn[M][:, :, 0].max() > s + tol or Cn[M][:, :, 1].max() > s + tol: return False
    I, J = pairs
    if not len(I): return True
    return bool(_sat_batch(Cn, _edge_normals_all(Cn), I, J)[0].min() >= -tol)

def verify_motion(sq, s, x, allpairs=None, state0=None, target=None):
    """Is the first-order motion x (n,3; theta component in radians) realisable as a finite motion?
    Returns the per-square displacement (n,3) of a genuinely feasible configuration, or None.
    With a `target` square the step sizes are measured on that square (an LP direction is often
    dominated by other squares, so scaling by the largest component leaves the one we are asking
    about barely moving) and the answer only counts if the target itself ends up displaced."""
    n = len(sq)
    amp = np.abs(x).max() if target is None else np.abs(x[target]).max()
    if amp <= 0: return None
    if target is not None and np.abs(x).max() > 500 * amp: return None   # others move 500x as far
    xn = x / amp
    moved = np.abs(xn).max(axis=1) > 1e-7
    if not moved.any(): return None
    M = np.where(moved)[0]
    if state0 is None: state0 = np.column_stack([sq[:, 0], sq[:, 1], np.radians(sq[:, 2])])
    ii, jj = near_pairs(state0[:, :2], PAIR_R) if allpairs is None else allpairs
    keep = moved[ii] | moved[jj]
    ii, jj = ii[keep], jj[keep]                      # every pair the motion can possibly disturb
    Cn0 = _corners_all(state0)
    Cm = Cn0[M]                                      # containment: only walls a mover is near
    D = np.stack([Cm[:, :, 0].min(1), s - Cm[:, :, 0].max(1),
                  Cm[:, :, 1].min(1), s - Cm[:, :, 1].max(1)], 1)
    mi, mw = np.where(D < 0.25)
    walls = (np.repeat(M[mi], 4), np.tile(np.arange(4), len(mi)), np.repeat(mw, 4))
    # constraints for the optimiser: only pairs close enough to interact over the step we take
    if len(ii):
        near = _sat_batch(Cn0, _edge_normals_all(Cn0), ii, jj)[0] < 0.3
        cpairs = (ii[near], jj[near])
    else:
        cpairs = (ii, jj)
    # work on the sub-configuration that the motion and its constraints touch
    sub = np.unique(np.concatenate([M, cpairs[0], cpairs[1], walls[0]]).astype(int))
    ren = np.full(n, -1); ren[sub] = np.arange(len(sub))
    sub0 = state0[sub]
    Ms, cps, wls = ren[M], (ren[cpairs[0]], ren[cpairs[1]]), (ren[walls[0]], walls[1], walls[2])
    tpos = None if target is None else int(np.where(M == target)[0][0])
    for eps in (1e-2, 1e-3, 1e-4):
        d = _project_motion(sub0, s, Ms, xn[M], eps, cps, wls)
        if np.abs(d if tpos is None else d[tpos]).max() < MIN_DISP: continue
        if _still_ok(state0, s, M, d, (ii, jj), walls, FEAS_TOL):
            dx = np.zeros((n, 3)); dx[M] = d
            return dx
    return None

def analyze(s, squares_str, want_regions=True):
    sq = np.array([[float(a), float(b), float(c)] for a, b, c in squares_str])
    n = len(sq)
    state0 = np.column_stack([sq[:, 0], sq[:, 1], np.radians(sq[:, 2])])
    P = _corners_all(state0)
    th = sq[:, 2] % 90.0
    tilt = np.where(th > 45, th - 90, th)

    # ---- angle groups (exact within 1e-9) ----
    order = np.argsort(th); groups = []
    for i in order:
        if groups and abs(th[i] - groups[-1]['theta']) < 1e-9: groups[-1]['members'].append(int(i))
        else: groups.append({'theta': float(th[i]), 'members': [int(i)]})
    for g in groups: g['count'] = len(g['members'])

    # ---- contacts ----
    contacts = []   # dict(i, j or None, wall, type, points [[x,y]..], normal [nx,ny]) normal: direction i moves to separate
    near = []
    # walls: 0 left x=0, 1 right x=s, 2 bottom y=0, 3 top y=s
    wallN = [np.array([1.,0.]), np.array([-1.,0.]), np.array([0.,1.]), np.array([0.,-1.])]
    for w, (ax, val) in enumerate(((0, 0.), (0, s), (1, 0.), (1, s))):
        D = P[:, :, ax] - val if val == 0 else val - P[:, :, ax]      # (n,4)
        G = D.min(1)
        for i in np.where(G < NEAR_MAX)[0]:
            i = int(i); g = float(G[i])
            if g < CONTACT_TOL:
                idx = np.where(D[i] < CONTACT_TOL)[0]
                contacts.append({'i': i, 'j': None, 'wall': w, 'type': 'edge-wall' if len(idx) >= 2 else 'corner-wall',
                                 'points': P[i][idx].tolist(), 'normal': wallN[w].tolist(), 'gap': float(max(g, 0))})
            else:
                near.append({'i': i, 'j': None, 'wall': w, 'gap': g})
    C = sq[:, :2]
    NRM = _edge_normals_all(P)
    pi, pj = near_pairs(C, math.sqrt(2) + NEAR_MAX)
    pg = _sat_batch(P, NRM, pi, pj)[0] if len(pi) else np.zeros(0)
    for m in np.where((pg >= CONTACT_TOL) & (pg < NEAR_MAX))[0]:
        near.append({'i': int(pi[m]), 'j': int(pj[m]), 'gap': float(pg[m])})
    for m in np.where(pg < CONTACT_TOL)[0]:
            i, j = int(pi[m]), int(pj[m]); g = float(pg[m])
            # classify: which axes are (near) maximal
            cands = []
            for ow, Pp, Nn in ((0, P[i], NRM[i]), (1, P[j], NRM[j])):
                Q = P[j] if ow == 0 else P[i]
                for k, nrm in enumerate(Nn):
                    gg = np.min(Q @ nrm) - np.max(Pp @ nrm)
                    if gg > g - CONTACT_TOL: cands.append((ow, k, nrm, gg))
            # vertices of the "other" body touching the face
            def touching(ow, nrm):
                Pp, Q = (P[i], P[j]) if ow == 0 else (P[j], P[i])
                proj = Q @ nrm; return np.where(proj < proj.min() + CONTACT_TOL)[0], Q
            # group candidate axes by direction (n and -n are the same axis)
            axes = []
            for ow, k, nrm, gg in cands:
                for ax in axes:
                    if abs(abs(np.dot(ax['n'], nrm)) - 1) < 1e-9: ax['members'].append((ow, k, nrm)); break
                else: axes.append({'n': nrm, 'members': [(ow, k, nrm)]})
            def features(ow, nrm):
                """touching vertices of the other body, overlap interval along the face, endpoints"""
                Pp, Q = (P[i], P[j]) if ow == 0 else (P[j], P[i])
                proj = Q @ nrm; idx = np.where(proj < proj.min() + CONTACT_TOL)[0]
                t = np.array([-nrm[1], nrm[0]])
                fidx = np.where(Pp @ nrm > (Pp @ nrm).max() - CONTACT_TOL)[0]
                lo = max((Pp[fidx] @ t).min(), (Q[idx] @ t).min()); hi = min((Pp[fidx] @ t).max(), (Q[idx] @ t).max())
                base = Q[idx[0]]; b0 = base @ t
                pts = [(base + (lo - b0) * t).tolist(), (base + (hi - b0) * t).tolist()]
                return idx, Q, lo, hi, pts
            best_axis, best_ov = None, -1
            for ax in axes:
                ow, k, nrm = ax['members'][0]
                idx, Q, lo, hi, pts = features(ow, nrm)
                ax['ov'] = hi - lo; ax['idx'] = idx; ax['pts'] = pts; ax['ow'] = ow; ax['nrm'] = nrm; ax['Q'] = Q
                if hi - lo > best_ov: best_ov, best_axis = hi - lo, ax
            if best_ov > CONTACT_TOL or len(axes) == 1:
                ax = best_axis; ow, nrm = ax['ow'], ax['nrm']
                nrm_i = -nrm if ow == 0 else nrm
                if ax['ov'] > CONTACT_TOL:
                    contacts.append({'i': i, 'j': int(j), 'type': 'edge-edge', 'points': ax['pts'], 'normal': nrm_i.tolist(),
                                     'face': int(ow), 'gap': float(max(g, 0)), 'overlap': float(ax['ov'])})
                else:
                    contacts.append({'i': i, 'j': int(j), 'type': 'corner-edge', 'points': ax['Q'][ax['idx']][:1].tolist(),
                                     'normal': nrm_i.tolist(), 'face': int(ow), 'gap': float(max(g, 0))})
            else:
                alts = [{'normal': (-ax['nrm'] if ax['ow'] == 0 else ax['nrm']).tolist(), 'points': ax['pts'][:1], 'face': int(ax['ow'])} for ax in axes]
                contacts.append({'i': i, 'j': int(j), 'type': 'corner-corner', 'points': alts[0]['points'],
                                 'normal': alts[0]['normal'], 'alts': alts, 'gap': float(max(g, 0))})
    # refine tiny gaps in mpmath: exact contact or tiny near miss?
    for c in contacts:
        if c['j'] is not None and c['gap'] > 1e-13:
            c['gap_mp'] = float(mp_gap(squares_str[c['i']], squares_str[c['j']]))
    for c in contacts:
        c['exact'] = (c['gap_mp'] < 1e-20) if 'gap_mp' in c else True   # no gap_mp => float gap <= 1e-13 => exact
    for m in near:
        if m['j'] is not None and m['gap'] < 1e-6:
            m['gap_mp'] = float(mp_gap(squares_str[m['i']], squares_str[m['j']]))
    near.sort(key=lambda m: m['gap'])

    # ---- first-order rigidity ----
    # Linearised non-penetration at every contact: n . (vel_i(p) - vel_j(p)) >= 0 over x = (vx_i, vy_i, w_i)_i.
    # A corner-corner contact is not one half-space but a union of them (the motion has to clear the
    # corner on one side or the other), which is not an LP constraint.  Dropping them outright -- the
    # obvious relaxation -- is far too loose: a grid of axis-aligned squares touches corner-to-corner
    # along every diagonal, and without those the grid shears freely, so hundreds of squares come out
    # "first-order mobile but nothing verifies".  Instead they are enforced by branch and cut: solve,
    # look for corner-corner contacts the answer drives into, pin each to the side it is closest to
    # satisfying, and re-solve.  Every solution is then a genuine first-order motion.
    import scipy.sparse as sp
    def row_entries(ii, jj, p, nrm):
        ent = []
        rx, ry = p[0] - sq[ii][0], p[1] - sq[ii][1]
        ent += [(3*ii, nrm[0]), (3*ii+1, nrm[1]), (3*ii+2, -nrm[0]*ry + nrm[1]*rx)]
        if jj is not None:
            rx, ry = p[0] - sq[jj][0], p[1] - sq[jj][1]
            ent += [(3*jj, -nrm[0]), (3*jj+1, -nrm[1]), (3*jj+2, -(-nrm[0]*ry + nrm[1]*rx))]
        return ent
    rows_i, cols_i, vals_i, nrows = [], [], [], 0
    row_i, row_j = [], []
    crows, ccols, cvals, ncr = [], [], [], 0
    crow_i, crow_j, cgroups = [], [], []
    for c in contacts:
        if c['type'] == 'corner-corner':
            grp = []
            for alt in c.get('alts', ()):
                for col, v in row_entries(c['i'], c['j'], alt['points'][0], np.array(alt['normal'])):
                    crows.append(ncr); ccols.append(col); cvals.append(v)
                crow_i.append(c['i']); crow_j.append(-1 if c['j'] is None else c['j'])
                grp.append(ncr); ncr += 1
            if len(grp) > 1: cgroups.append(grp)
            continue
        for p in c['points']:
            for col, v in row_entries(c['i'], c['j'], p, np.array(c['normal'])):
                rows_i.append(nrows); cols_i.append(col); vals_i.append(v)
            row_i.append(c['i']); row_j.append(-1 if c['j'] is None else c['j'])
            nrows += 1
    A = sp.csr_matrix((vals_i, (rows_i, cols_i)), shape=(nrows, 3*n)) if nrows else None
    ent_r = np.array(rows_i, int); ent_c = np.array(cols_i, int); ent_v = np.array(vals_i)
    row_i = np.array(row_i, int); row_j = np.array(row_j, int)
    AC = sp.csr_matrix((cvals, (crows, ccols)), shape=(ncr, 3*n)) if ncr else None
    cent_r = np.array(crows, int); cent_c = np.array(ccols, int); cent_v = np.array(cvals)
    crow_i = np.array(crow_i, int); crow_j = np.array(crow_j, int)
    def cc_cuts(Cl, groups, x):
        """Corner-corner contacts this motion drives into; pin each to its least violated side."""
        if Cl is None or not len(groups): return []
        v = Cl @ x; tol = 1e-9 * max(1.0, float(np.abs(x).max()))
        out = []
        for g in groups:
            k = max(g, key=lambda t: v[t])
            if v[k] < -tol: out.append(k)
        return out
    def branch_cut(run, Cl, groups):
        """Iterate `run(extra rows)` until no corner-corner contact is violated."""
        add = []
        for _ in range(12):
            x = run(add)
            if x is None: return None
            cuts = cc_cuts(Cl, groups, x)
            if not cuts: return x
            add = sorted(set(add) | set(cuts))
        return None
    Aub = -A if A is not None else None; bub = np.zeros(nrows) if nrows else None
    bounds = [(-1, 1)] * (3*n)
    verified = np.zeros(n, bool); firstorder = np.zeros(n, bool); motion_dirs = {}
    state0 = np.column_stack([sq[:, 0], sq[:, 1], np.radians(sq[:, 2])])
    allpairs = near_pairs(state0[:, :2], PAIR_R)
    def mark(x):
        for i in np.where(np.abs(x.reshape(n, 3)).max(axis=1) > 1e-7)[0]: firstorder[i] = True
    def try_step(x, target=None):
        """Cheap test: does stepping the whole configuration along x stay feasible?
        Only motions that are not tangent to any constraint pass this."""
        x = x.reshape(n, 3); amp = np.abs(x).max() if target is None else np.abs(x[target]).max()
        if amp <= 1e-7 or (target is not None and np.abs(x).max() > 500 * amp): return False
        mark(x)
        for eps in (1e-2, 1e-3):
            if config_feasible(sq, None, s, x * (eps / amp)):
                accept(x * (eps / amp)); return True
        return False
    def accept(dx):
        for i in np.where(np.abs(dx).max(axis=1) >= MIN_DISP)[0]:
            if not verified[i]: verified[i] = True; motion_dirs[int(i)] = dx[i].tolist()
    tried = set()
    def try_motion(x, target=None):
        """Full test: project onto the true (nonlinear) feasible set near q0 + eps*x.
        Neighbouring squares keep proposing the same mechanism, so directions are memoised."""
        x = x.reshape(n, 3); amp = np.abs(x).max()
        if amp <= 1e-7: return False
        mark(x)
        key = np.round(x / amp, 9).tobytes() + (b'' if target is None else b'#%d' % target)
        if key in tried: return False
        tried.add(key)
        dx = verify_motion(sq, s, x, allpairs, state0, target)
        if dx is None or not (np.abs(dx).max(axis=1) >= MIN_DISP).any(): return False
        accept(dx); return True
    def solve(cobj):
        def run(add):
            Au = -sp.vstack([A, AC[add]], format='csr') if add else Aub
            r = linprog(cobj, A_ub=Au, b_ub=np.zeros(Au.shape[0]) if Au is not None else None,
                        bounds=bounds, method='highs')
            return r.x if r.status == 0 and -r.fun > 1e-9 else None
        return branch_cut(run, AC, cgroups)
    # An infinitesimal motion of square i can only involve its connected component of the contact
    # graph, so the LP for square i is exact on that component.  For big packings the component is
    # the whole packing and a 3n-variable LP per square is hopeless, so it is truncated to the
    # CLUSTER_CAP squares nearest i in the contact graph (everything outside is held fixed) --
    # an approximation only for packings larger than the cap.
    adj = [[] for _ in range(n)]
    for c in contacts:
        if c['j'] is not None: adj[c['i']].append(c['j']); adj[int(c['j'])].append(c['i'])
    def cluster(i, cap):
        seen = {i}; frontier = [i]
        while frontier and len(seen) < cap:
            nxt = []
            for u in frontier:
                for v in adj[u]:
                    if v not in seen:
                        seen.add(v); nxt.append(v)
                        if len(seen) >= cap: break
                if len(seen) >= cap: break
            frontier = nxt
        return np.fromiter(sorted(seen), int, len(seen))
    def local_system(cl):
        """Rows that constrain the cluster, restricted to its columns: the hard ones, and the
        corner-corner alternatives with their grouping renumbered."""
        inc = np.zeros(n, bool); inc[cl] = True
        cmap = np.full(3*n, -1); cmap[(3*cl[:, None] + np.arange(3)).ravel()] = np.arange(3*len(cl))
        def cut(rmask, tot, er, ec, ev):
            nr = int(rmask.sum()); rmap = np.full(tot, -1); rmap[rmask] = np.arange(nr)
            k = (rmap[er] >= 0) & (cmap[ec] >= 0)
            return sp.csr_matrix((ev[k], (rmap[er[k]], cmap[ec[k]])), shape=(nr, 3*len(cl))), rmap
        Al = cut(inc[row_i] | ((row_j >= 0) & inc[row_j]), nrows, ent_r, ent_c, ent_v)[0]
        if AC is None: return Al, None, []
        cmask = inc[crow_i] | ((crow_j >= 0) & inc[crow_j])
        Cl, crmap = cut(cmask, ncr, cent_r, cent_c, cent_v)
        return Al, Cl, [[crmap[t] for t in g] for g in cgroups if crmap[g[0]] >= 0]
    def local_solve(Al, Cl, groups, cobj):
        m3 = Al.shape[1]
        def run(add):
            Au = sp.vstack([Al, Cl[add]], format='csr') if add else Al
            r = linprog(cobj, A_ub=-Au if Au.shape[0] else None,
                        b_ub=np.zeros(Au.shape[0]) if Au.shape[0] else None,
                        bounds=[(-1, 1)]*m3, method='highs')
            return r.x if r.status == 0 and -r.fun > 1e-9 else None
        return branch_cut(run, Cl, groups)
    def local_sparse(Al, Cl, groups, cobj):
        """Sparsest motion with unit progress on the target: min sum|x| s.t. A x >= 0, cobj.x = -1.
        A plain LP vertex drags along every square it legally can, and such a sprawling compound
        motion is nearly always blocked at second order somewhere; the sparse one is a single
        mechanism the nonlinear projection can follow."""
        m3 = Al.shape[1]; I = sp.eye(m3, format='csr')
        obj = np.concatenate([np.zeros(m3), np.ones(m3)])
        Aeq = sp.csr_matrix(np.concatenate([cobj, np.zeros(m3)]).reshape(1, -1))
        def run(add):
            Au = sp.vstack([Al, Cl[add]], format='csr') if add else Al
            Ab = sp.bmat([[-Au, None], [I, -I], [-I, -I]], format='csr')
            r = linprog(obj, A_ub=Ab, b_ub=np.zeros(Ab.shape[0]), A_eq=Aeq, b_eq=[-1.],
                        bounds=[(-1e3, 1e3)]*m3 + [(0, 1e3)]*m3, method='highs')
            return r.x[:m3] if r.status == 0 else None
        return branch_cut(run, Cl, groups)
    rng = np.random.default_rng(0)
    # Which squares have first-order freedom at all?  Rounds of one global LP: maximise a random
    # objective supported on the squares not yet known to move.  The optimal vertex of the motion
    # cone drags along essentially everything that can move, so a couple of rounds suffice; and
    # when both +c and -c give zero, every feasible x has c.x_R = 0, which for a generic c means
    # x_R = 0 -- the remaining squares are first-order rigid.  (n LPs of 3n variables, one per
    # square, is what made this unusable on the big packings.)
    movable = np.zeros(n, bool)
    for _ in range(12):
        R = ~movable
        if A is None or not R.any(): break
        c = rng.normal(size=3*n) * np.repeat(R, 3)
        got = False
        for sgn in (1, -1):
            x = solve(-sgn * c)
            if x is None: continue
            sup = np.abs(x.reshape(n, 3)).max(axis=1) > 1e-7
            got |= bool((sup & R).any())
            movable |= sup
            try_step(x)                     # collective motions that need no correction at all
        if not got: break
    firstorder |= movable
    # Every square the LP says can move but that no verified motion has moved yet needs its own
    # search.  Usually there are a handful; but in a packing that is second-order rigid over a wide
    # region there can be hundreds, none of which will ever verify, and a thorough search of each is
    # wasted.  So: a cheap pass over all of them first, and a thorough one (larger neighbourhood,
    # more directions) only over what is left, and only if the cheap pass was paying off.
    def search(i, cap, tries):
        cl = cluster(i, cap); Al, Cl, ccg = local_system(cl)
        loc = int(np.where(cl == i)[0][0])
        for _ in range(tries):
            c = rng.normal(size=3); dead = True
            for sgn in (1, -1):
                cobj = np.zeros(3*len(cl)); cobj[3*loc:3*loc+3] = -sgn * c
                xl = local_solve(Al, Cl, ccg, cobj)
                if xl is None: continue     # nothing in this half-space
                dead = False
                xs = local_sparse(Al, Cl, ccg, cobj)
                for xc in ((xs, xl) if xs is not None else (xl,)):
                    x = np.zeros((n, 3)); x[cl] = xc.reshape(len(cl), 3)
                    try_step(x, i) or try_motion(x)
                    if verified[i]: break   # a motion may verify its neighbours but not i itself
                if verified[i]: break
            if verified[i] or dead: break   # neither +c nor -c: the motion cone at i is trivial
    nmov = max(int((movable & ~verified).sum()), 1)
    attempts = wins = 0
    for cap, tries in ((80, 2), (min(n, CLUSTER_CAP, max(80, CLUSTER_BUDGET // nmov)), 3)):
        if attempts and wins * 5 < attempts: break
        for i in np.where(movable)[0]:
            if verified[i]: continue
            attempts += 1
            search(i, cap, tries)
            wins += bool(verified[i])

    # alone: only square i's variables, everything else fixed (2 LPs each; the sampled region verifies)
    free_alone = np.zeros(n, bool); alone_dirs = {}
    Acsc = A.tocsc() if A is not None else None
    ACcsc = AC.tocsc() if AC is not None else None
    for i in range(n):
        if A is None: free_alone[i] = True; continue
        Ai = Acsc[:, 3*i:3*i+3].toarray(); Ai = Ai[np.abs(Ai).sum(axis=1) > 0]
        Ci = np.zeros((0, 3)); gi = []
        if ACcsc is not None:
            Cf = ACcsc[:, 3*i:3*i+3].toarray()
            keep = np.abs(Cf).sum(axis=1) > 0
            ren = np.full(ncr, -1); ren[keep] = np.arange(int(keep.sum()))
            Ci = Cf[keep]; gi = [[ren[t] for t in g] for g in cgroups if ren[g[0]] >= 0]
        c = rng.normal(size=3)
        for sgn in (1, -1):
            def run(add, Ai=Ai, Ci=Ci, sgn=sgn, c=c):
                Au = np.vstack([Ai, Ci[add]]) if add else Ai
                r = linprog(-sgn*c, A_ub=-Au if len(Au) else None,
                            b_ub=np.zeros(len(Au)) if len(Au) else None, bounds=[(-1, 1)]*3, method='highs')
                return r.x if r.status == 0 and -r.fun > 1e-9 else None
            x = branch_cut(run, Ci, gi)
            if x is not None:
                free_alone[i] = True; alone_dirs.setdefault(i, []).append(x)
    # ---- free regions for squares free alone (sampled) ----
    regions = {}
    if want_regions:
        for i in np.where(free_alone)[0]:
            regions[int(i)] = free_region(i, sq, P, s, seeds=alone_dirs.get(int(i), []))
    truly_free = [int(i) for i in np.where(free_alone)[0] if regions.get(int(i), {}).get('free', True)]
    for i in truly_free: verified[i] = True
    # The sampler only walks translations at fixed theta plus a pure spin about the centre, so it
    # misses a square that can only roll (turn and slide together).  Re-test those nonlinearly,
    # with the square alone: it then counts as mobile, though not as 'free' with a drawn region.
    for i in np.where(free_alone)[0]:
        if verified[i]: continue
        for d3 in alone_dirs.get(int(i), []):
            x = np.zeros((n, 3)); x[i] = d3
            if try_motion(x): break
    wedged = [int(i) for i in range(n) if not verified[i] and (free_alone[i] or firstorder[i])]
    mobile = verified.copy()
    # ---- symmetry ----
    sym = symmetry(sq, s)
    # ---- category ----
    nz = [g for g in groups if g['theta'] > 1e-9]
    if not nz: cat = 'trivial'
    elif all(abs(g['theta'] - 45) < 1e-9 for g in nz): cat = 'diagonal'
    else: cat = 'tilted'
    return {
        'n': n, 's': s, 'waste': s*s - n,
        'theta': th.round(12).tolist(), 'tilt': tilt.round(12).tolist(),
        'angle_groups': groups, 'n_angles': len(groups), 'n_rotated': int(np.sum(th > 1e-9)),
        'category': cat, 'symmetry': sym,
        'contacts': contacts, 'near_misses': near[:40],
        'mobile': mobile.tolist(), 'free_alone': free_alone.tolist(), 'motion_dirs': motion_dirs,
        'rigid': not mobile.any(), 'regions': regions, 'wedged': wedged, 'free': truly_free, 'first_order_only': [int(i) for i in np.where(firstorder & ~verified)[0]],
        'n_corner_corner': sum(1 for c in contacts if c['type'] == 'corner-corner'),
    }

def config_feasible(sq, P, s, dx, tol=1e-12):
    """Is the packing still valid after moving every square i by (dx[i,0], dx[i,1]) and rotating by dx[i,2] rad?
    Penetration up to `tol` is tolerated (P is unused, kept for the old call signature)."""
    n = len(sq)
    st = np.column_stack([sq[:, 0] + dx[:, 0], sq[:, 1] + dx[:, 1], np.radians(sq[:, 2]) + dx[:, 2]])
    Q = _corners_all(st)
    if Q[:,:,0].min() < -tol or Q[:,:,1].min() < -tol or Q[:,:,0].max() > s + tol or Q[:,:,1].max() > s + tol: return False
    C = st[:, :2]
    ii, jj = near_pairs(C, PAIR_R)
    if not len(ii): return True
    N = _edge_normals_all(Q)
    return bool(_sat_batch(Q, N, ii, jj)[0].min() >= -tol)

def feasible_batch(Qs, P, s, neigh):
    """Qs: (M,4,2) candidate corner sets. Returns bool (M,) feasibility vs container and neighbours."""
    ok = (Qs[:,:,0].min(1) >= -1e-12) & (Qs[:,:,1].min(1) >= -1e-12) & (Qs[:,:,0].max(1) <= s+1e-12) & (Qs[:,:,1].max(1) <= s+1e-12)
    for j in neigh:
        Pj = P[j]; sep = np.full(len(Qs), -np.inf)
        for nrm in normals(Pj):                      # axes of the fixed neighbour
            sep = np.maximum(sep, (Qs @ nrm).min(1) - (Pj @ nrm).max())
        NQ = np.stack([Qs[:,1]-Qs[:,0], Qs[:,3]-Qs[:,0]], 1)   # two edge directions of the moving square
        for e in range(2):
            nrm = np.stack([NQ[:,e,1], -NQ[:,e,0]], 1); nrm /= np.linalg.norm(nrm, axis=1)[:,None]
            pq = np.einsum('mkd,md->mk', Qs, nrm); pp = np.einsum('kd,md->mk', Pj, nrm)
            sep = np.maximum(sep, np.maximum(pp.min(1) - pq.max(1), pq.min(1) - pp.max(1)))
        ok &= sep >= -1e-12
    return ok

def _slide_scan(i, P, s, neigh, dirs, max_t, K=64, rounds=3):
    """For each unit direction, how far square i can be pushed before it hits something.
    One batched feasibility call marches all directions out to max_t together; each further round
    subdivides the bracket K-fold, so `rounds` extra calls give a resolution of max_t / K**(rounds+1).
    Marching (rather than bisecting from max_t) keeps the answer inside the component holding t = 0."""
    D = np.asarray(dirs, float)                       # (nd, 2)
    nd = len(D)
    lo = np.zeros(nd); hi = np.full(nd, max_t)
    ts = np.concatenate([[SLIDE_TOL], np.linspace(0, max_t, K + 1)[1:]])
    for r in range(rounds + 1):
        if r: ts = np.linspace(0, 1, K + 1)[1:]        # fractions of the current bracket
        off = D[:, None, :] * (lo[:, None] + (hi - lo)[:, None] * ts[None, :])[:, :, None] if r \
              else D[:, None, :] * ts[None, :, None]   # (nd, K, 2)
        ok = feasible_batch((P[i][None, None] + off[:, :, None, :]).reshape(-1, 4, 2), P, s, neigh).reshape(nd, -1)
        first_bad = np.where(ok.all(1), ok.shape[1], ok.argmin(1))
        vals = lo[:, None] + (hi - lo)[:, None] * ts[None, :] if r else np.broadcast_to(ts, (nd, len(ts))).copy()
        newlo = np.where(first_bad > 0, vals[np.arange(nd), np.maximum(first_bad - 1, 0)], 0.0)
        newhi = np.where(first_bad < ok.shape[1], vals[np.arange(nd), np.minimum(first_bad, ok.shape[1] - 1)], hi)
        lo, hi = newlo, np.maximum(newhi, newlo)
        if r == 0: lo = np.where(ok[:, 0], lo, 0.0); hi = np.where(ok[:, 0], hi, 0.0)
    return lo

def _slide_dirs(i, P, neigh, seeds, sweep=72):
    """Unit directions worth a line search.  A square can be free to move while the feasible set of
    its translations is a *segment* -- it slides along a channel of exactly zero width -- and such a
    set almost never meets a node of the sampling grid.  A zero-width channel runs parallel to the
    faces that form it, so every edge direction of the square and of its neighbours is a candidate;
    the LP's first-order directions and a coarse sweep cover everything else."""
    V = []
    for Q in [P[i]] + [P[j] for j in neigh]:
        for k in range(4):
            e = Q[(k+1) % 4] - Q[k]
            V += [e, -e, np.array([e[1], -e[0]]), np.array([-e[1], e[0]])]
    V += [np.array([1., 0.]), np.array([-1., 0.]), np.array([0., 1.]), np.array([0., -1.])]
    for d in seeds: V.append(np.asarray(d, float)[:2])
    V += [np.array([math.cos(a), math.sin(a)]) for a in np.linspace(0, 2*math.pi, sweep, endpoint=False)]
    U = []
    for v in V:
        L = float(np.linalg.norm(v))
        if L < 1e-12: continue
        u = v / L
        if not any(float(np.dot(u, w)) > 1 - 1e-12 for w in U): U.append(u)
    return U

def free_region(i, sq, P, s, R=0.25, N=61, max_R=4.0, seeds=()):
    C = sq[:, :2]
    while True:
        neigh = [j for j in range(len(sq)) if j != i and np.sum((C[j]-C[i])**2) < (math.sqrt(2) + R*1.5)**2]
        xs = np.linspace(-R, R, N); h = xs[1] - xs[0]
        gx, gy = np.meshgrid(xs, xs)                    # gy rows, gx cols
        offs = np.stack([gx.ravel(), gy.ravel()], 1)
        Qs = P[i][None] + offs[:, None, :]
        mask = feasible_batch(Qs, P, s, neigh).reshape(N, N)
        c0 = N // 2; comp = np.zeros_like(mask); stack = [(c0, c0)]
        if not mask[c0, c0]:
            return {'note': 'origin infeasible (numerical)', 'R': R, 'free': False}
        while stack:
            b, a = stack.pop()
            if comp[b, a] or not mask[b, a]: continue
            comp[b, a] = True
            for db, da in ((1,0),(-1,0),(0,1),(0,-1)):
                bb, aa = b+db, a+da
                if 0 <= bb < N and 0 <= aa < N: stack.append((bb, aa))
        # line searches: exact reach along each axis, and slides the grid cannot see (see _slide_dirs)
        dirs = _slide_dirs(i, P, neigh, seeds)
        reach = _slide_scan(i, P, s, neigh, dirs, R)
        # Grow the window if the square reaches its edge -- whether the grid component runs off it or
        # a line search is still unobstructed at R.  Without the second test a long slide would be
        # reported as exactly R.
        touches = comp[0].any() or comp[-1].any() or comp[:,0].any() or comp[:,-1].any()
        if (touches or reach.max() >= R * (1 - 1e-9)) and R < max_R: R *= 2; continue
        ys, xs_ = np.where(comp)
        rot = rotation_range(i, sq, P, s, neigh)
        k = int(np.argmax(reach)); best_t, best_u = float(reach[k]), dirs[k]
        ext = {}
        for name, v in (('+x', (1, 0)), ('-x', (-1, 0)), ('+y', (0, 1)), ('-y', (0, -1))):
            for u, t in zip(dirs, reach):
                if abs(u[0]-v[0]) < 1e-12 and abs(u[1]-v[1]) < 1e-12: ext[name] = float(t)
        area_free = bool(comp.sum() > 1)
        out = {'R': R, 'N': N, 'cell': float(h), 'area': float(comp.sum() * h * h) if area_free else 0.0,
               'dx': [-ext.get('-x', 0.0), ext.get('+x', 0.0)], 'dy': [-ext.get('-y', 0.0), ext.get('+y', 0.0)],
               'rot': rot, 'mask': [''.join('1' if v else '0' for v in rowv) for rowv in comp]}
        out['free'] = bool(area_free or best_t > SLIDE_TOL or max(abs(r) for r in rot) > 1e-3)
        out['kind'] = 'region' if area_free else ('slide' if best_t > SLIDE_TOL else ('turn' if out['free'] else 'none'))
        if best_t > SLIDE_TOL: out['slide'] = {'dir': [float(best_u[0]), float(best_u[1])], 'len': best_t}
        return out

def rotation_range(i, sq, P, s, neigh):
    out = []
    for sgn in (1, -1):
        lo, hi = 0.0, 45.0
        def ok(d):
            Q = corners(sq[i][0], sq[i][1], sq[i][2] + sgn*d)
            if Q[:,0].min() < -1e-12 or Q[:,1].min() < -1e-12 or Q[:,0].max() > s+1e-12 or Q[:,1].max() > s+1e-12: return False
            return all(sat(Q, P[j])[0] >= -1e-12 for j in neigh)
        if not ok(1e-2) and not ok(1e-3): out.append(0.0); continue
        if ok(hi): out.append(sgn*45.0); continue
        for _ in range(40):
            m = (lo + hi) / 2
            if ok(m): lo = m
            else: hi = m
        out.append(sgn * lo)
    return out   # [ccw max, cw max] in degrees (cw negative)

def symmetry(sq, s):
    C = sq[:, :2] - s/2; th = sq[:, 2] % 90
    def key(pts, ths):
        return sorted((round(x, 7), round(y, 7), round(t % 90, 6) if (round(t % 90, 6) < 89.9999995) else 0.0) for (x, y), t in zip(pts, ths))
    base = key(C, th)
    ops = {'rot90': (lambda p: np.stack([-p[:,1], p[:,0]], 1), lambda t: t),
           'rot180': (lambda p: -p, lambda t: t),
           'rot270': (lambda p: np.stack([p[:,1], -p[:,0]], 1), lambda t: t),
           'mirror_v': (lambda p: np.stack([-p[:,0], p[:,1]], 1), lambda t: -t),   # left-right flip
           'mirror_h': (lambda p: np.stack([p[:,0], -p[:,1]], 1), lambda t: -t),
           'mirror_d': (lambda p: np.stack([p[:,1], p[:,0]], 1), lambda t: 90 - t),
           'mirror_a': (lambda p: np.stack([-p[:,1], -p[:,0]], 1), lambda t: 90 - t)}
    have = [name for name, (fp, ft) in ops.items() if key(fp(C), ft(th)) == base]
    mirrors = sum(1 for h in have if h.startswith('mirror'))
    rots = sum(1 for h in have if h.startswith('rot'))
    if mirrors == 4: cls = 'D4 (four mirror axes)'
    elif mirrors == 2: cls = 'two mirror axes'
    elif mirrors == 1: cls = 'one mirror axis'
    elif rots == 3: cls = 'rotational (90°) only'
    elif rots == 1: cls = 'rotational (180°) only'
    else: cls = 'none'
    return {'elements': have, 'class': cls}
