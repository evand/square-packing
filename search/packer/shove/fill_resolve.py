"""Fill C (local re-solve, 10-16; notes and every test run: shove/NOTES_resolve.md).

Release a band of jumble squares next to the opened gap and re-place band + m new squares there, the rest fixed.
The opened gap follows the block's staircase, so the region is measured from the block itself (Chebyshev distance to
the block squares <= r + 1) as well as from the strip guide (the wall strips).

Default variants (RESOLVE_VARIANTS='keep:1,kick:1,kick:1,kick:1,kick:1'); each is screen-quenched once (fq SCREEN,
loosen 1.0) and the best RESOLVE_SCREEN = 3 distinct screened sides are returned (RESOLVE_SCREEN=0: every variant,
unscreened, in order):
  keep:r     band = non-block squares within r + 0.5 of the block (or r of the guide), kept in place; add clean cells
             on the block's own lattice in the gap (nearest the block first), then the rest greedily (clean
             max-contact poses with angles 0, 45 and the band's own, else least penetration), then relax band + new
             only (L-BFGS on sum pen^2 + wall^2, the rest fixed).
  kick:r     keep, then a gaussian kick on band + new (sigma 0.03, angle 0.6 deg) and the local relax again.
Kept for the record (tested, not in the default; NOTES_resolve.md says why):
  refill:r   remove the band and re-place band + m from scratch (greedy): loses squares against the jumble.
  slide:r    slide band squares toward the gap, then insert.  tilt:r  keep without lattice cells.
  ext:r      keep with extra cells continuing the band's own local lattices.
  squeeze:r  keep at d = 1, then re-shove at d' < 1 with the same added squares and relax the band locally (bisection
             on d'): the band never absorbs any misfit.
"""
from __future__ import annotations
import math, os
import numpy as np
import harness as H

TAU = 1e-9


# ------------------------------------------------------------------ vectorised SAT for unit squares

def _hw(a, th):
    """half-width of a unit square at angle a (rad) projected on the axis at angle th (rad)"""
    t = a - th
    return 0.5 * (np.abs(np.cos(t)) + np.abs(np.sin(t)))


def pen_matrix(P, Q, grow=1.0):
    """(k, q) penetration depth of poses P (k, 3; deg; side `grow`) against unit squares Q (q, 3); 0 = disjoint"""
    P = np.asarray(P, float); Q = np.asarray(Q, float)
    if len(P) == 0 or len(Q) == 0:
        return np.zeros((len(P), len(Q)))
    ap = np.radians(P[:, 2])[:, None]; aq = np.radians(Q[:, 2])[None, :]
    dx = Q[None, :, 0] - P[:, None, 0]; dy = Q[None, :, 1] - P[:, None, 1]
    best = np.full(dx.shape, np.inf)
    for th in (ap, ap + math.pi / 2, aq, aq + math.pi / 2):
        th = np.broadcast_to(th, dx.shape)
        c, s = np.cos(th), np.sin(th)
        ov = grow * _hw(ap, th) + _hw(aq, th) - np.abs(dx * c + dy * s)
        best = np.minimum(best, ov)
    far = dx * dx + dy * dy > (0.5 * math.sqrt(2) * (1 + grow)) ** 2
    best[far] = 0.0
    return np.maximum(best, 0.0)


def wall_pen(P, S, grow=1.0):
    P = np.asarray(P, float)
    r = np.radians(P[:, 2] % 90.0)
    h = 0.5 * grow * (np.abs(np.cos(r)) + np.abs(np.sin(r)))
    return np.maximum.reduce([np.zeros(len(P)), h - P[:, 0], h - P[:, 1], P[:, 0] + h - S, P[:, 1] + h - S])


class Board:
    """fixed squares + placed squares with a cell grid for neighbour queries"""
    def __init__(self, S, sq):
        self.S = S; self.sq = [tuple(q) for q in np.asarray(sq, float)]
        self.cell = {}
        for i, q in enumerate(self.sq):
            self.cell.setdefault((int(q[0]), int(q[1])), []).append(i)

    def add(self, q):
        self.sq.append(tuple(q)); self.cell.setdefault((int(q[0]), int(q[1])), []).append(len(self.sq) - 1)

    def near(self, x0, y0, x1, y1):
        out = []
        for i in range(int(x0) - 2, int(x1) + 3):
            for j in range(int(y0) - 2, int(y1) + 3):
                out.extend(self.cell.get((i, j), ()))
        return out

    def score(self, P, grow=1.12):
        """(max pen incl. walls, contact) for poses P: contact = summed overlap of the pose grown by `grow` (more = snugger)"""
        P = np.asarray(P, float)
        idx = self.near(P[:, 0].min(), P[:, 1].min(), P[:, 0].max(), P[:, 1].max())
        Q = np.array([self.sq[i] for i in idx]) if idx else np.zeros((0, 3))
        pm = pen_matrix(P, Q)
        pen = np.maximum(pm.max(axis=1) if len(Q) else 0.0, wall_pen(P, self.S))
        cm = pen_matrix(P, Q, grow)
        contact = cm.sum(axis=1) + wall_pen(P, self.S, grow) * 2
        return pen, contact


# ------------------------------------------------------------------ region helpers

def dist_rects(x, y, rects):
    d = np.full(np.shape(x), np.inf)
    for x0, y0, x1, y1 in rects:
        dx = np.maximum(np.maximum(x0 - x, x - x1), 0); dy = np.maximum(np.maximum(y0 - y, y - y1), 0)
        d = np.minimum(d, np.hypot(dx, dy))
    return d


def angle_set(sq_band, extra=(0.0, 45.0)):
    a = sorted(set([round(v % 90.0, 1) for v in np.asarray(sq_band, float)[:, 2]] + list(extra))) if len(sq_band) else list(extra)
    out = []
    for v in a:                                            # merge angles within 1 degree
        if not out or v - out[-1] > 1.0:
            out.append(v)
    return out


def greedy_place(B, count, region_pts, angles, rng, samples=2500, rounds_clean=None, local=True, spread=1):
    """place up to `count` squares on board B: clean poses first (pen < TAU, max contact), then least penetration.
    region_pts (p, 2): candidate centres.  -> (list of poses, number placed cleanly)"""
    placed, nclean = [], 0
    angles = np.asarray(angles, float)
    pts = np.asarray(region_pts, float)
    if len(pts) == 0:
        return placed, 0
    clean_ok = True
    while len(placed) < count:
        k = min(samples, len(pts) * len(angles))
        ii = np.array([rng.randrange(len(pts)) for _ in range(k)])
        aa = angles[np.array([rng.randrange(len(angles)) for _ in range(k)])]
        jit = np.array([[rng.uniform(-0.03, 0.03), rng.uniform(-0.03, 0.03)] for _ in range(k)])
        P = np.column_stack([pts[ii] + jit, aa])
        pen, con = B.score(P)
        if clean_ok:
            ok = pen < TAU
            if ok.any():
                j = np.nonzero(ok)[0][np.argmax(con[ok])]
                pose = P[j]
                if local:                                   # nudge it snug: small moves that keep it clean, raise contact
                    for step in (0.05, 0.02, 0.01):
                        for _ in range(6):
                            cand = pose + np.array([[step, 0, 0], [-step, 0, 0], [0, step, 0], [0, -step, 0],
                                                    [0, 0, 20 * step], [0, 0, -20 * step]])
                            pn, cn = B.score(cand)
                            good = pn < TAU
                            if not good.any():
                                break
                            jj = np.nonzero(good)[0][np.argmax(cn[good])]
                            if cn[jj] <= B.score(pose[None])[1][0] + 1e-12:
                                break
                            pose = cand[jj]
                B.add(pose); placed.append(tuple(pose)); nclean += 1
                continue
            clean_ok = False
        key = pen - 1e-3 * con
        if spread > 1:                                     # random pick among the `spread` best (restarts differ)
            j = int(np.argsort(key)[rng.randrange(min(spread, len(key)))])
        else:
            j = int(np.argmin(key))
        B.add(P[j]); placed.append(tuple(P[j]))
    return placed, nclean


def free_points(S, sq, rects, rmax, h=0.05):
    xs, ys, free = H.free_mask(S, sq, h=h)
    X, Y = np.meshgrid(xs, ys)
    m = free & (dist_rects(X, Y, rects) <= rmax)
    return np.column_stack([X[m], Y[m]])


# ------------------------------------------------------------------ variants

def v_refill(S, sq, meta, m, rng, r):
    sq = np.asarray(sq, float)
    blk = np.zeros(len(sq), bool); blk[meta['block']] = True
    rects = meta['strip']
    dist = dist_rects(sq[:, 0], sq[:, 1], rects)
    band = (~blk) & (dist <= r)
    keep = sq[~band]
    nb = int(band.sum())
    angs = angle_set(sq[band])
    pts = free_points(S, keep, rects, r + 0.8)
    B = Board(S, keep)
    placed, nc = greedy_place(B, nb + m, pts, angs, rng)
    return np.vstack([keep, np.array(placed)]) if placed else keep, dict(nb=nb, clean=nc, target=nb + m)


def v_slide(S, sq, meta, m, rng, r, rounds=3):
    sq = np.asarray(sq, float).copy()
    blk = np.zeros(len(sq), bool); blk[meta['block']] = True
    rects = meta['strip']
    dist = dist_rects(sq[:, 0], sq[:, 1], rects)
    band = np.nonzero((~blk) & (dist <= r))[0]
    for _ in range(rounds):
        for i in band[np.argsort(dist[band])]:
            x, y, a = sq[i]
            # direction toward the nearest strip point
            best = None
            for x0, y0, x1, y1 in rects:
                px, py = min(max(x, x0), x1), min(max(y, y0), y1)
                dd = math.hypot(px - x, py - y)
                if best is None or dd < best[0]:
                    best = (dd, px - x, py - y)
            if best is None or best[0] < 1e-6:
                continue
            ux, uy = best[1] / best[0], best[2] / best[0]
            others = np.delete(sq, i, axis=0)
            Bo = Board(S, others)
            lo, hi = 0.0, min(best[0] + 0.5, 1.0)
            if Bo.score(np.array([[x + ux * hi, y + uy * hi, a]]))[0][0] < TAU:
                lo = hi
            else:
                for _ in range(20):
                    mid = 0.5 * (lo + hi)
                    if Bo.score(np.array([[x + ux * mid, y + uy * mid, a]]))[0][0] < TAU:
                        lo = mid
                    else:
                        hi = mid
            sq[i, 0] += ux * lo; sq[i, 1] += uy * lo
        dist = dist_rects(sq[:, 0], sq[:, 1], rects)
    angs = angle_set(sq[band])
    pts = free_points(S, sq, rects, r + 1.5)
    B = Board(S, sq)
    placed, nc = greedy_place(B, m, pts, angs, rng)
    return np.vstack([sq, np.array(placed)]) if placed else sq, dict(nb=len(band), clean=nc, target=m)



# ------------------------------------------------------------------ local relaxation (only `mov` squares move)

def pair_sep(P, Q):
    """for pose pairs P[i], Q[i] (k, 3): (pen, ux, uy) with u the min-overlap SAT axis oriented from P to Q"""
    ap = np.radians(P[:, 2]); aq = np.radians(Q[:, 2])
    dx = Q[:, 0] - P[:, 0]; dy = Q[:, 1] - P[:, 1]
    best = np.full(len(P), np.inf); bu = np.zeros((len(P), 2))
    for th in (ap, ap + math.pi / 2, aq, aq + math.pi / 2):
        c, s = np.cos(th), np.sin(th)
        proj = dx * c + dy * s
        ov = _hw(ap, th) + _hw(aq, th) - np.abs(proj)
        sg = np.where(proj >= 0, 1.0, -1.0)
        better = ov < best
        best = np.where(better, ov, best)
        bu[better] = np.column_stack([c * sg, s * sg])[better]
    return np.maximum(best, 0.0), bu[:, 0], bu[:, 1]


def local_relax(S, X, mov, iters=400, step=0.6, rot=True, rng=None, tol=1e-7):
    """Jacobi push relaxation of squares `mov` (indices into X) against everything else (fixed) and the walls, side S.
    Each overlapping pair pushes its moving member(s) apart along the min SAT axis; random small rotations accepted when
    they reduce a square's own overlap.  -> (X, max residual penetration among pairs touching mov)"""
    X = np.asarray(X, float).copy(); mov = np.asarray(sorted(set(int(i) for i in mov)), int)
    ismov = np.zeros(len(X), bool); ismov[mov] = True
    worst = np.inf
    for it in range(iters):
        M = X[mov]
        dx = X[None, :, 0] - M[:, None, 0]; dy = X[None, :, 1] - M[:, None, 1]
        close = (dx * dx + dy * dy < 2.0 + 1e-9)
        close[np.arange(len(mov)), mov] = False
        ii, jj = np.nonzero(close)
        disp = np.zeros((len(mov), 2))
        worst = 0.0
        if len(ii):
            pen, ux, uy = pair_sep(M[ii], X[jj])
            hit = pen > 0
            worst = float(pen.max()) if len(pen) else 0.0
            f = np.where(ismov[jj], 0.5, 1.0) * pen
            np.add.at(disp, ii[hit], -np.column_stack([ux, uy])[hit] * f[hit, None])
        r = np.radians(M[:, 2] % 90.0); h = 0.5 * (np.abs(np.cos(r)) + np.abs(np.sin(r)))
        wl = np.maximum(h - M[:, 0], 0); wr = np.maximum(M[:, 0] + h - S, 0)
        wb = np.maximum(h - M[:, 1], 0); wt = np.maximum(M[:, 1] + h - S, 0)
        disp[:, 0] += wl - wr; disp[:, 1] += wb - wt
        worst = max(worst, float(np.max(np.maximum.reduce([wl, wr, wb, wt]))) if len(M) else 0.0)
        if worst < tol:
            break
        X[mov, :2] += step * disp
        if rot and rng is not None and it % 5 == 4:
            for k in np.nonzero(np.abs(disp).sum(axis=1) > 1e-9)[0]:
                i = mov[k]
                cur = _own_pen(S, X, i)
                da = rng.gauss(0, 3.0)
                X[i, 2] += da
                if _own_pen(S, X, i) >= cur:
                    X[i, 2] -= da
    return X, worst


def _own_pen(S, X, i):
    d = np.hypot(X[:, 0] - X[i, 0], X[:, 1] - X[i, 1])
    nb = np.nonzero((d < math.sqrt(2) + 1e-9) & (np.arange(len(X)) != i))[0]
    p = 0.0
    if len(nb):
        pen, _, _ = pair_sep(np.repeat(X[i:i + 1], len(nb), axis=0), X[nb])
        p = float((pen ** 2).sum())
    return p + float(wall_pen(X[i:i + 1], S)[0] ** 2)


def _own(S, Xm, Xall, ii, jj, mm, w_mm):
    """per moving square: sum over its pairs (ii -> local mov index, jj -> global) of weight * pen^2, + wall^2"""
    e = np.zeros(len(Xm))
    if len(ii):
        pen, _, _ = pair_sep(Xm[ii], Xall[jj])
        np.add.at(e, ii, np.where(mm, w_mm, 1.0) * pen ** 2)
    return e + wall_pen(Xm, S) ** 2


def lbfgs_relax(S, X, mov, rounds=6, maxiter=300, h=1e-6, tol=1e-9):
    """overlap-energy minimisation (sum pen^2 + wall^2) over the poses of `mov`, the rest fixed; L-BFGS with a
    finite-difference gradient of each square's own energy (vectorised).  Neighbour lists rebuilt each round."""
    from scipy.optimize import minimize
    X = np.asarray(X, float).copy(); mov = np.asarray(sorted(set(int(i) for i in mov)), int)
    ismov = np.zeros(len(X), bool); ismov[mov] = True
    scale = np.array([1.0, 1.0, 180 / math.pi])        # optimise angle in radians
    worst = np.inf
    for _ in range(rounds):
        M = X[mov]
        dx = X[None, :, 0] - M[:, None, 0]; dy = X[None, :, 1] - M[:, None, 1]
        close = dx * dx + dy * dy < (math.sqrt(2) + 0.4) ** 2
        close[np.arange(len(mov)), mov] = False
        ii, jj = np.nonzero(close)
        mm = ismov[jj]
        pos = {g: k for k, g in enumerate(mov)}

        def f(v):
            Xm = (v.reshape(-1, 3)) * scale
            Xall = X.copy(); Xall[mov] = Xm
            E = _own(S, Xm, Xall, ii, jj, mm, 0.5).sum()
            g = np.zeros_like(Xm)
            for c in range(3):
                d = np.zeros(3); d[c] = h * scale[c]
                ep = _own(S, Xm + d, Xall, ii, jj, mm, 1.0); em = _own(S, Xm - d, Xall, ii, jj, mm, 1.0)
                g[:, c] = (ep - em) / (2 * h)
            return E, g.ravel()
        r = minimize(f, (M / scale).ravel(), jac=True, method='L-BFGS-B', options=dict(maxiter=maxiter, gtol=1e-14, ftol=1e-16))
        X[mov] = r.x.reshape(-1, 3) * scale
        Xm = X[mov]
        pen, _, _ = pair_sep(Xm[ii], X[jj]) if len(ii) else (np.zeros(0), 0, 0)
        worst = max(float(pen.max()) if len(pen) else 0.0, float(wall_pen(Xm, S).max()))
        if worst < tol:
            break
    return X, worst


def block_dist(x, y, Bk):
    """Chebyshev distance from points to the union of the (upright) block squares Bk"""
    x = np.asarray(x, float); y = np.asarray(y, float)
    d = np.full(x.shape, np.inf)
    for bx, by, _ in Bk:
        d = np.minimum(d, np.maximum(np.maximum(np.abs(x - bx), np.abs(y - by)) - 0.5, 0))
    return d


def region(S, sq, meta, r):
    """(band indices, region test) : band = non-block squares whose centre lies within r + 0.5 of the block (Chebyshev)
    or within r of the strip guide (the wall strips); the region = points within r + 1 of the block or r + 0.8 of the
    strip guide (the opened gap follows the block's staircase, not the guide's bounding-box rectangles)"""
    sq = np.asarray(sq, float)
    blk = np.zeros(len(sq), bool); blk[meta['block']] = True
    Bk = sq[blk]
    db = block_dist(sq[:, 0], sq[:, 1], Bk); ds = dist_rects(sq[:, 0], sq[:, 1], meta['strip'])
    band = np.nonzero((~blk) & ((db <= r + 0.5) | (ds <= r)))[0]
    test = lambda X, Y: (block_dist(X, Y, Bk) <= r + 1.0) | (dist_rects(X, Y, meta['strip']) <= r + 0.8)
    return band, test


def region_points(S, sq, test, h=0.05):
    xs, ys, free = H.free_mask(S, sq, h=h)
    X, Y = np.meshgrid(xs, ys)
    msk = free & test(X, Y)
    return np.column_stack([X[msk], Y[msk]])


def lattice_cells(S, sq, meta, test, cap):
    """clean upright poses on the block's lattice inside the region, nearest the block first (greedy)"""
    sq = np.asarray(sq, float)
    Bk = sq[meta['block']]
    ox, oy = Bk[0, 0] % 1.0, Bk[0, 1] % 1.0
    xs = np.arange(ox, S, 1.0); ys = np.arange(oy, S, 1.0)
    X, Y = np.meshgrid(xs, ys); X = X.ravel(); Y = Y.ravel()
    ok = test(X, Y) & (X >= 0.5 - 1e-9) & (X <= S - 0.5 + 1e-9) & (Y >= 0.5 - 1e-9) & (Y <= S - 0.5 + 1e-9)
    X, Y = X[ok], Y[ok]
    order = np.argsort(block_dist(X, Y, Bk), kind='stable')
    B = Board(S, sq); add = []
    for k in order:
        if len(add) >= cap:
            break
        P = np.array([[X[k], Y[k], 0.0]])
        if B.score(P)[0][0] < TAU:
            B.add(P[0]); add.append(tuple(P[0]))
    return add


def ext_cells(S, cur, seeds, test, cap):
    """clean poses continuing the jumble's own local lattices: each seed square translated by +-1 along its two
    frame axes (same angle), repeatedly (new squares become seeds); max contact first.  -> list of poses"""
    B = Board(S, cur); add = []
    seeds = [tuple(q) for q in np.asarray(seeds, float)]
    def cands(qs):
        out = []
        for x, y, a in qs:
            r = math.radians(a); c, s_ = math.cos(r), math.sin(r)
            for ex, ey in ((c, s_), (-c, -s_), (-s_, c), (s_, -c)):
                out.append((x + ex, y + ey, a))
        return out
    C = np.array(cands(seeds)) if seeds else np.zeros((0, 3))
    while len(add) < cap and len(C):
        C = C[test(C[:, 0], C[:, 1])]
        if not len(C):
            break
        pen, con = B.score(C)
        ok = pen < TAU
        if not ok.any():
            break
        j = np.nonzero(ok)[0][np.argmax(con[ok])]
        B.add(C[j]); add.append(tuple(C[j]))
        C = np.vstack([np.delete(C, j, axis=0), np.array(cands([tuple(C[j])]))])
    return add


def v_keep(S, sq, meta, m, rng, r, relax='lbfgs', lattice=True, tries=None):
    """keep the band, add clean block-lattice cells, the rest greedily (clean max-contact, then least penetration) in
    the region, then relax band + new only (rest fixed) by overlap-energy minimisation.  When the clean placement falls
    short of m, the greedy + relax is restarted `tries` times (RESOLVE_TRIES, default 6; restarts > 0 pick the
    overlapping squares at random among the 30 least-penetrating samples) and the lowest residual wins."""
    sq = np.asarray(sq, float)
    if tries is None:                                      # m = 2k + 2 never left the grid in dev runs: 2 restarts
        tries = int(os.environ.get('RESOLVE_TRIES', '6')) if m <= 2 * meta['k0'] + 1 else 2
    band, test = region(S, sq, meta, r)
    add = lattice_cells(S, sq, meta, test, m) if lattice in (True, 'block', 'both') else []
    if lattice in ('ext', 'both') and len(add) < m:
        base = np.vstack([sq] + ([np.array(add)] if add else []))
        tl = band[np.minimum(sq[band, 2] % 90, 90 - sq[band, 2] % 90) > 1.0] if lattice == 'ext' else band
        add = add + ext_cells(S, base, sq[tl], test, m - len(add))
    cur = np.vstack([sq] + ([np.array(add)] if add else []))
    angs = angle_set(sq[band])
    pts = region_points(S, cur, test)
    best = None
    for t in range(max(1, tries)):
        B = Board(S, cur)
        more, nc = greedy_place(B, m - len(add), pts, angs, rng, spread=1 if t == 0 else 30) if len(add) < m else ([], 0)
        X = np.vstack([cur] + ([np.array(more)] if more else []))
        mov = list(band) + list(range(len(sq), len(X)))
        if relax == 'none':
            w = None
        else:
            X, w = local_relax(S, X, mov, rng=rng) if relax == 'push' else lbfgs_relax(S, X, mov)
        if best is None or (w is not None and w < best[1]):
            best = (X, w, nc)
        if w is None or w < 1e-9 or nc >= m - len(add):
            break
    X, w, nc = best
    return X, dict(nb=len(band), clean=len(add) + nc, lat=len(add), target=m, resid=w)

def _offsets(corner, d):
    bx = d if corner[1] == 'R' else 0.0; by = d if corner[0] == 'U' else 0.0
    return bx, by, d - bx, d - by


def unshove(sq, meta):
    """the source packing (side s0) from the shoved one"""
    sq = np.asarray(sq, float).copy()
    blk = np.zeros(len(sq), bool); blk[meta['block']] = True
    bx, by, ox, oy = _offsets(meta['corner'], meta['d'])
    sq[blk, 0] -= bx; sq[blk, 1] -= by; sq[~blk, 0] -= ox; sq[~blk, 1] -= oy
    return sq


def v_squeeze(S, sq, meta, m, rng, r, dlo=0.9, steps=5, tol=1e-6, extra=0):
    """the keep fill (lattice cells first) built at d = 1, then re-shoved at smaller d' with the same added squares
    (carried in the block's frame) and band + new relaxed (rest fixed) in the box s0 + d'; bisection on d' for the
    smallest d' whose local relaxation reaches max penetration <= tol.  -> the feasible state at the best d'."""
    src = unshove(sq, meta)
    s0, corner = meta['s0'], meta['corner']
    S1, sq1, meta1 = H.shove(s0, src, corner, 1.0, block=meta['block'])
    X1, info = v_keep(S1, sq1, meta1, m + extra, rng, r)
    n0 = len(src)
    new1 = X1[n0:]
    bx1, by1, _, _ = _offsets(corner, 1.0)

    def build(dp):
        Sd, sqd, metad = H.shove(s0, src, corner, dp, block=meta['block'])
        bx, by, _, _ = _offsets(corner, dp)
        nw = new1.copy(); nw[:, 0] += bx - bx1; nw[:, 1] += by - by1
        band, _ = region(Sd, sqd, metad, r)
        X = np.vstack([sqd, nw])
        mov = list(band) + list(range(n0, len(X)))
        X, w = lbfgs_relax(Sd, X, mov)
        return Sd, X, w

    best = (S1, X1, info['resid']) if info['resid'] <= tol else None
    lo, hi = dlo, 1.0
    Sd, X, w = build(lo)
    if w <= tol:
        best = (Sd, X, w); hi = lo
    else:
        for _ in range(steps):
            mid = 0.5 * (lo + hi)
            Sd, X, w = build(mid)
            if w <= tol:
                best = (Sd, X, w); hi = mid
            else:
                lo = mid
    info = dict(info, dbest=None if best is None else best[0] - s0, resid=None if best is None else best[2])
    return (X1 if best is None else best[1]), info


def v_ext(S, sq, meta, m, rng, r):
    """keep, but block-lattice cells then cells continuing the band's own lattices (tilted rows) before the greedy"""
    return v_keep(S, sq, meta, m, rng, r, lattice='both')


def v_norelax(S, sq, meta, m, rng, r):
    """ablation: keep without the local relaxation"""
    return v_keep(S, sq, meta, m, rng, r, relax='none')


def v_tilt(S, sq, meta, m, rng, r):
    """keep without the lattice cells: every added square by the greedy (clean max-contact first), band relaxed"""
    return v_keep(S, sq, meta, m, rng, r, lattice=False)


def v_kick(S, sq, meta, m, rng, r, sig=0.03):
    """keep, then release the band: gaussian kick on band + new (x, y: sig; angle: 20 sig deg) and local relax again"""
    X, info = v_keep(S, sq, meta, m, rng, r)
    band, _ = region(S, np.asarray(sq, float), meta, r)
    mov = np.array(list(band) + list(range(len(sq), len(X))), int)
    X[mov, 0] += [rng.gauss(0, sig) for _ in mov]; X[mov, 1] += [rng.gauss(0, sig) for _ in mov]
    X[mov, 2] += [rng.gauss(0, 20 * sig) for _ in mov]
    X, w = lbfgs_relax(S, X, mov)
    return X, dict(info, resid=w)


VARIANTS = os.environ.get('RESOLVE_VARIANTS', 'keep:1,kick:1,kick:1,kick:1,kick:1')
SCREEN_KEEP = int(os.environ.get('RESOLVE_SCREEN', '3'))     # > 0: screen-quench every variant, return the best this many
LOG = os.environ.get('RESOLVE_LOG')


VFUN = dict(refill=v_refill, slide=v_slide, keep=v_keep, squeeze=v_squeeze, tilt=v_tilt, kick=v_kick, ext=v_ext, norelax=v_norelax)


def _screen(S, X):
    """one fq SCREEN quench at loosen 1.0 -> side (inf if infeasible)"""
    import hop, tempfile
    from pk.engine import SCREEN
    with tempfile.TemporaryDirectory(prefix='resolve_') as tmp:
        q = hop.quench(S, [tuple(map(float, v)) for v in X], tmp, opt=SCREEN.but(loosen=1.0))
    return q[0] if q else math.inf


def fill(S, sq, meta, m, rng):
    out = []
    for spec in VARIANTS.split(','):
        name, r = spec.split(':'); r = float(r)
        X, info = VFUN[name](S, sq, meta, m, rng, r)
        sc = _screen(S, X) if SCREEN_KEEP else None
        if LOG:
            with open(LOG, 'a') as f:
                f.write(f"n0 {meta['n0']} {meta['corner']} m {m} {spec}: band {info['nb']} target {info['target']} "
                        f"clean {info['clean']} resid {info.get('resid')} cells {len(meta['cells'])} screen {sc}\n")
        out.append((sc, len(out), X))
    if SCREEN_KEEP:
        # distinct screened sides first (a kick that lands in keep's basin adds nothing), best first
        out.sort(key=lambda t: (t[0], t[1]))
        pick, seen = [], []
        for sc, i, X in out:
            if all(abs(sc - v) > 1e-7 for v in seen):
                pick.append(X); seen.append(sc)
        for sc, i, X in out:
            if len(pick) >= SCREEN_KEEP:
                break
            if not any(X is P for P in pick):
                pick.append(X)
        return pick[:SCREEN_KEEP]
    return [X for _, _, X in out]
