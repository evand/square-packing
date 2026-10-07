#!/usr/bin/env python3
"""floatsep.py -- a CONSERVATIVE float separation oracle for the cutting-plane loop (task K).

Why.  `search/branch.py` and `search/tighten.py` call the exact i128 verifier every round only to
be told *which poses the current weights under-cover*; the answer becomes LP rows.  Measured on the
live runs that call was 80-99 % of a 3-5 h round (`search/LOOPSPEED.md`).  Exactness is needed
once, for the file that is finally claimed -- not to steer the LP.  This module answers the same
question with a dense float scan in numpy, in seconds.

**Nothing is claimed on the strength of this module.**  It only chooses rows.  Every row it emits
is a valid constraint by construction (see "Rows are valid" below), so even a badly wrong scan can
only make the LP slower or weaker, never unsound; and the loop still ends with the unchanged
`verify` + `xcheck.py` step on the final file.

What it computes.  For each angle bin `k` of the verifier's own `N`-net that we choose to look at:

  * `theta_k = 2 arctan(k/N)`, `h = sigma_k/2` and the admissible centre box `[L, s-L]^2`,
    `L = w_min(k)/2`, are taken from the verifier's own integer formulas (`bin_geometry`), so a
    row emitted here has exactly the shape of a verifier witness and lands in the same
    `--row-grid` bucket;
  * the captured atom weight is evaluated on a regular grid of centres in the bin's rotated frame
    (pitch `pitch`), by the sliding-window/prefix-sum scan of `lp_search.scan_angle`;
  * every anchor-clique image is credited its weight on the poses it *provably* holds, with the
    same conservative predicates the verifier uses on a cell -- `contains` by the sigma_k-square
    (the verifier uses the larger exact bin core) and `meets` by the hexagon `A + [-h,h]^2`
    including the segment's own normal -- each clique counted once;
  * for a branch certificate the region threshold `1 + lambda_j` is subtracted where the centre
    lies in corner box `j`.

Conservative direction.  Clique membership here is a SUBSET of what the verifier credits (sigma_k
square vs exact bin core; `TOL` shaved off every inequality in the safe direction), so the float
oracle never claims coverage the exact sweep would deny: it can only report a violation the exact
sweep would not, never miss one for the reason of over-crediting.  It *can* miss a violation
because the grid is finite -- which is why `branch.py --sep float` still runs the exact sweep every
`--exact-every` rounds and at convergence, and treats "the float oracle is clean" as a hint, never
as a verdict.

Rows are valid.  A row `(c, theta_k, sigma_k/2)` asserts "the closed square of side sigma_k at
centre `c`, angle `theta_k`, captures >= 1".  For any centre in the bin's admissible box there is
an angle `theta` in `[theta_k, theta_{k+1}]` for which the unit square at `(c, theta)` lies in the
container, and that unit square CONTAINS the concentric sigma_k-square at angle `theta_k` (this is
exactly the verifier's bin argument).  So the row is implied by the true covering constraint at an
admissible pose, whatever this module computed.

Usage (as a library; `branch.py --sep float` drives it)
    import floatsep
    mv, rows, info = floatsep.separate(P, w, s, N, bins, topk=6, pitch=0.004,
                                       cliques=[(weight, images), ...], lam=[...], r=1.0, pool=pool)
"""
import math, os
import numpy as np

TOL = 1e-9          # a pose within TOL of a piece boundary is NOT credited (safe direction)

SCALE = 10 ** 6     # the verifier's rounding scale for sigma_k and w_min


# ------------------------------------------------------------------ the verifier's bin geometry
def bin_geometry(N, k):
    """(cos, sin, sigma, w_min) of bin k of the N-net, with the verifier's own integer rounding:
    sigma_k rounded DOWN to 1e-6 (a smaller square) and w_min rounded DOWN (a bigger centre box).
    Mirrors `bin_geometry` / the `main` loop of verify/src/main.rs.

    `int()` is not decoration: g0*g1*SCALE is ~6e19 at N = 2000, k near N, which overflows int64.
    Called with a numpy integer (as any caller iterating over an array of bins does) the result
    silently wraps and both sigma_k and the admissible box come out garbage."""
    N = int(N); k = int(k)
    c0, s0, g0 = N * N - k * k, 2 * k * N, N * N + k * k
    k2 = k + 1
    c1, s1, g1 = N * N - k2 * k2, 2 * k2 * N, N * N + k2 * k2
    cd = c0 * c1 + s0 * s1
    sd = c0 * s1 - s0 * c1
    sg = (g0 * g1 * SCALE) // (cd + sd) / SCALE
    wm = min(((c0 + s0) * SCALE) // g0, ((c1 + s1) * SCALE) // g1) / SCALE
    return c0 / g0, s0 / g0, sg, wm


def bin_list(N, sym, dtheta):
    """the bins to scan: every bin whose angle is a multiple of `dtheta` radians, over [0,45deg]
    when the certificate has the verifier's D4 reduction and [0,90deg) otherwise.  Bin 0 and the
    last bin are always included (the extreme angles are where covering is tightest)."""
    # the verifier sweeps k = 0 .. kk-1 (bin k covers [theta_k, theta_{k+1}]), so the last bin is
    # kk-1 in both cases
    if sym:
        kk = 0
        while (kk + N) * (kk + N) < 2 * N * N: kk += 1
        kk -= 1
    else:
        kk = N - 1
    ks = []
    last = -1e9
    for k in range(kk + 1):
        th = 2.0 * math.atan(k / N)
        if k == 0 or k == kk or th - last >= dtheta - 1e-12:
            ks.append(k); last = th
    return ks


# ------------------------------------------------------------------ clique geometry, float
def prep_cliques(cliques):
    """[(weight, images)] -> a light, fork-friendly description:
    [(weight, [ [ (Apts, [Fpts, ...]), ... ] ]) ...] with Apts an (n,2) float array (n = 1 for a
    point anchor, 2 for a segment).  One entry per IMAGE (the certificate writes one file clique
    per image and the verifier credits each file clique once, so the LP coefficient is the number
    of images containing the pose -- exactly `anchorclique.coeff`)."""
    out = []
    for (w, imgs) in cliques:
        if w <= 0: continue
        for cl in imgs:
            anchors, pieces = cl
            ps = []
            for (a, filt) in pieces:
                A = np.array(_pts(anchors[a]), dtype=float)
                F = [np.array(_pts(anchors[f]), dtype=float) for f in filt]
                ps.append((A, F))
            out.append((float(w), ps))
    return out


def _pts(a):
    return [(float(a[1]), float(a[2]))] if a[0] == 'P' else [(float(a[1]), float(a[2])), (float(a[3]), float(a[4]))]


def _rect_contains(A, ct, st, h):
    """u-space rectangle of { centres whose sigma_k-square CONTAINS the anchor A } (conservative:
    the verifier's exact bin core is larger, so this is a subset)"""
    q0 = A[:, 0] * ct + A[:, 1] * st
    q1 = -A[:, 0] * st + A[:, 1] * ct
    return (q0.max() - h + TOL, q0.min() + h - TOL, q1.max() - h + TOL, q1.min() + h - TOL)


def _rect_meets(F, ct, st, h):
    """u-space bounding rectangle of { centres whose sigma_k-square MEETS the anchor F }, plus the
    segment's own normal half-planes as (n0, n1, base0, base1, bound) or None for a point anchor"""
    q0 = F[:, 0] * ct + F[:, 1] * st
    q1 = -F[:, 0] * st + F[:, 1] * ct
    rect = (q0.min() - h + TOL, q0.max() + h - TOL, q1.min() - h + TOL, q1.max() + h - TOL)
    if len(F) == 1:
        return (q0.max() - h + TOL, q0.min() + h - TOL, q1.max() - h + TOL, q1.min() + h - TOL), None
    nx, ny = -(F[1, 1] - F[0, 1]), (F[1, 0] - F[0, 0])       # normal in container coordinates
    n0 = nx * ct + ny * st
    n1 = -nx * st + ny * ct
    return rect, (n0, n1, q0[0], q1[0], h * (abs(n0) + abs(n1)) - TOL)


def clique_credit(CQ, ct, st, h, g0, g1, out):
    """add every clique's weight to `out` (shape (len(g0), len(g1))) on the grid CELLS it holds.
    `out` is indexed [i, j] = the grid cell around centre (g0[i], g1[j]) in the bin's rotated
    frame -- with the SAME half-side the atoms use, so that a pose is credited exactly what the
    verifier credits the cell that pose sits in when the cell is small (see `scan_bin` on why the
    half-side must not be shrunk)."""
    if not CQ: return
    n0, n1 = len(g0), len(g1)
    d0 = g0[1] - g0[0] if n0 > 1 else 1.0
    d1 = g1[1] - g1[0] if n1 > 1 else 1.0
    for (w, pieces) in CQ:
        acc = None; i0g = j0g = 0
        for (A, F) in pieces:
            a0, b0, a1, b1 = _rect_contains(A, ct, st, h)
            slabs = []
            for Fp in F:
                (r0, r1, s0_, s1_), nrm = _rect_meets(Fp, ct, st, h)
                a0 = max(a0, r0); b0 = min(b0, r1); a1 = max(a1, s0_); b1 = min(b1, s1_)
                if nrm is not None: slabs.append(nrm)
            if a0 > b0 or a1 > b1: continue
            i0 = max(0, int(math.ceil((a0 - g0[0]) / d0 - 1e-12)))
            i1 = min(n0, int(math.floor((b0 - g0[0]) / d0 + 1e-12)) + 1)
            j0 = max(0, int(math.ceil((a1 - g1[0]) / d1 - 1e-12)))
            j1 = min(n1, int(math.floor((b1 - g1[0]) / d1 + 1e-12)) + 1)
            if i1 <= i0 or j1 <= j0: continue
            U0 = g0[i0:i1][:, None]; U1 = g1[j0:j1][None, :]
            m = np.ones((i1 - i0, j1 - j0), dtype=bool)
            for (nn0, nn1, p0, p1, bnd) in slabs:
                m &= np.abs((U0 - p0) * nn0 + (U1 - p1) * nn1) <= bnd
            if not m.any(): continue
            if acc is None:
                acc = np.zeros((n0, n1), dtype=bool)
            acc[i0:i1, j0:j1] |= m
        if acc is not None:
            out += w * acc


# ------------------------------------------------------------------ the region threshold
def region_required(cx, cy, s, r, lam4, band):
    """The weight a pose must capture beyond 1, conservatively w.r.t. the verifier's CELL rule.

    `verify/src/main.rs`, `required(may, inside)`: a cell entirely inside corner box `j` must
    capture `1 + lambda_j`; a cell that merely MAY MEET one or more boxes must capture
    `1 + max_j max(lambda_j, 0)`; a cell meeting none must capture `1`.  So the `1 + lambda`
    threshold applies in a band OUTSIDE each box as wide as the cell, and a negative `lambda`
    (a box switched off, `lambda = -2`) is a discount only for a cell wholly inside that box --
    never for one straddling its boundary.

    Deciding the threshold by "is the centre in a box?", as this module first did, is optimistic
    twice over: it charges nothing in the outside band, and it hands the -2 discount to poses
    within a cell-width of the free box's boundary.  `band` must cover the cell, so pass at least
    the grid pitch.  Returns (required, flag), flag = box index + 1 for a pose inside a box."""
    boxes = ((0.0, 0.0), (s - r, 0.0), (0.0, s - r), (s - r, s - r))
    req = np.zeros(cx.shape); flag = np.zeros(cx.shape, dtype=np.int8)
    for j, (bx, by) in enumerate(boxes):
        lj = lam4[j] if j < len(lam4) else lam4[0]
        dx = np.maximum(np.maximum(bx - cx, cx - (bx + r)), 0.0)
        dy = np.maximum(np.maximum(by - cy, cy - (by + r)), 0.0)
        d = np.maximum(dx, dy)                                  # inf-norm distance, 0 inside
        req = np.maximum(req, np.where(d <= band, max(lj, 0.0), 0.0))
        inside = d <= 0.0
        flag = np.where(inside, j + 1, flag)
        if lj < 0.0:            # the discount, only where a whole cell is certainly inside
            deep = inside & (np.minimum(np.minimum(cx - bx, bx + r - cx),
                                        np.minimum(cy - by, by + r - cy)) > band)
            req = np.where(deep, lj, req)
    return req, flag


# ------------------------------------------------------------------ the scan
def bin_frame(s, N, k):
    """(cos, sin, h, theta, lo, hi) of bin k, or None if its admissible centre box is empty"""
    ct, st, sg, wm = bin_geometry(N, k)
    lo, hi = wm / 2.0, s - wm / 2.0
    if hi <= lo: return None
    return ct, st, sg / 2.0, 2.0 * math.atan(int(k) / int(N)), lo, hi


def bin_atoms(P, w, ct, st):
    """the bin's atoms in u0 and u1 order, the only per-bin O(m log m) work"""
    q0 = P[:, 0] * ct + P[:, 1] * st
    q1 = -P[:, 0] * st + P[:, 1] * ct
    o = np.argsort(q0, kind='stable'); q0s = q0[o]; q1o = q1[o]; wo = w[o]
    ordy = np.argsort(q1o, kind='stable')
    return q0s, q1o[ordy], wo[ordy], np.arange(len(o))[ordy]


def grid_values(P, w, s, ct, st, h, hc, g0, g1, CQ, lam4, r, band, lo, hi, pre=None):
    """value of every cell of the grid g0 x g1 in the bin's rotated frame: the weight all of its
    poses capture, plus the cliques that hold all of them, minus the region requirement.
    `hc` is the shrunk half-side (see `scan_bin`); `h` itself is only the admissible box."""
    q0s, q1s, wy, idx_in_q0 = pre if pre is not None else bin_atoms(P, w, ct, st)
    a = np.searchsorted(q0s, g0 - hc + TOL, 'left')
    b = np.searchsorted(q0s, g0 + hc - TOL, 'right')
    lo1 = np.searchsorted(q1s, g1 - hc + TOL, 'left')
    hi1 = np.searchsorted(q1s, g1 + hc - TOL, 'right')
    vals = np.zeros((len(g0), len(g1)))
    for i in range(len(g0)):
        if b[i] <= a[i]: continue
        mask = (idx_in_q0 >= a[i]) & (idx_in_q0 < b[i])
        cw = np.concatenate([[0.0], np.cumsum(np.where(mask, wy, 0.0))])
        vals[i, :] = cw[hi1] - cw[lo1]
    clique_credit(CQ, ct, st, hc, g0, g1, vals)
    U0 = g0[:, None]; U1 = g1[None, :]
    cx = U0 * ct - U1 * st
    cy = U0 * st + U1 * ct
    ok = (cx >= lo - 1e-12) & (cx <= hi + 1e-12) & (cy >= lo - 1e-12) & (cy <= hi + 1e-12)
    flag = np.zeros(vals.shape, dtype=np.int8)
    if lam4 is not None:
        req, flag = region_required(cx, cy, s, r, lam4, band)
        vals = vals - req
    return np.where(ok, vals, 9e9), flag


def scan_bin(P, w, s, N, k, pitch, CQ, lam, r, topk, block, thr, refine=0, rfac=4):
    """one bin: returns (min value seen, [(value, cx, cy, theta_k, h, flag), ...] worst first).

    The half-side is `h`, NOT `h` shrunk to the grid cell.  Shrinking looks like the conservative
    thing to do -- the verifier minimises over what EVERY pose of an arrangement cell captures --
    but the cutting-plane loop puts its atoms exactly ON the tight rows' square boundaries, so
    shrinking destroys the signal instead of tightening it: measured on `branch_L1110f_probe.txt`,
    the pose (3.48, 0.5) at theta = 0 captures 42 atoms of total weight 2.499999 at half-side
    `h = 0.499501` and 2 atoms of weight 0.013433 at `h - 0.002`.  With the shrink every converged
    certificate reports `float_min ~ -lambda` and the oracle is useless.  The pose-versus-cell gap
    is real but small (on that checkpoint, 375 of the exact sweep's 6000 witnesses capture >= 1 at
    the pose while their cell does not) and it is not what made the oracle optimistic in
    production -- the angle net was (`search/LOOPSPEED.md`).

    `refine`: after the grid pass, re-scan the worst `refine * topk` cells on a local grid at
    `pitch / rfac`.  OFF by default: measured on `branch_L1110f_probe.txt` it cost 2.5x the scan
    (24 s -> 60 s coarse, 162 s -> 403 s over the whole net) and moved the minimum in neither case.
    `descend`, seeded from last round's rows, is the local search that actually pays."""
    fr = bin_frame(s, N, k)
    if fr is None: return None, []
    ct, st, h, th, lo, hi = fr
    lam4 = None if lam is None else (list(lam) if len(lam) == 4 else [lam[0]] * 4)
    corners = np.array([[lo, lo], [hi, lo], [hi, hi], [lo, hi]])
    cu0 = corners[:, 0] * ct + corners[:, 1] * st
    cu1 = -corners[:, 0] * st + corners[:, 1] * ct
    g0 = np.arange(cu0.min(), cu0.max() + pitch * 0.5, pitch)
    g1 = np.arange(cu1.min(), cu1.max() + pitch * 0.5, pitch)
    if len(g0) < 2 or len(g1) < 2: return None, []
    hc = h
    pre = bin_atoms(P, w, ct, st)
    v, flag = grid_values(P, w, s, ct, st, h, hc, g0, g1, CQ, lam4, r, pitch, lo, hi, pre)
    mn = float(v.min())
    n0, n1 = v.shape
    p0 = (-n0) % block; p1 = (-n1) % block
    VP = np.pad(v, ((0, p0), (0, p1)), constant_values=9e9)
    bl = VP.reshape(VP.shape[0] // block, block, VP.shape[1] // block, block).transpose(0, 2, 1, 3).reshape(-1, block * block)
    bm = bl.min(axis=1); am = bl.argmin(axis=1)
    nb1 = VP.shape[1] // block
    want = max(topk, refine * topk)
    order = np.argsort(bm)[:want]
    sb = order[bm[order] < 9e8]
    bi, bj = np.divmod(sb, nb1); ii, jj = np.divmod(am[sb], block)
    r0 = bi * block + ii; r1 = bj * block + jj
    keep = (r0 < n0) & (r1 < n1); r0 = r0[keep]; r1 = r1[keep]
    cand = [(float(v[i, j]), float(g0[i]), float(g1[j]), int(flag[i, j])) for i, j in zip(r0, r1)]
    # local descent around the worst cells
    if refine and rfac > 1 and cand:
        fp = pitch / rfac; hc2 = h
        for (_, u0c, u1c, _) in list(cand):
            lg0 = np.arange(u0c - pitch, u0c + pitch + fp * 0.5, fp)
            lg1 = np.arange(u1c - pitch, u1c + pitch + fp * 0.5, fp)
            if len(lg0) < 2 or len(lg1) < 2: continue
            lv, lf = grid_values(P, w, s, ct, st, h, hc2, lg0, lg1, CQ, lam4, r, fp, lo, hi, pre)
            m2 = float(lv.min())
            if m2 < mn: mn = m2
            i2, j2 = np.unravel_index(int(np.argmin(lv)), lv.shape)
            if lv[i2, j2] < 9e8:
                cand.append((float(lv[i2, j2]), float(lg0[i2]), float(lg1[j2]), int(lf[i2, j2])))
    cand.sort()
    out = []
    for (val, u0c, u1c, fl) in cand:
        if val >= thr or len(out) >= topk: break
        out.append((val, float(u0c * ct - u1c * st), float(u0c * st + u1c * ct), th, h, fl))
    return mn, out


def descend(P, w, s, N, seeds, pitch, CQ, lam, r, thr, levels=3, rfac=4, half=2.0):
    """Local descent from given poses -- last round's rows, i.e. the poses that were worst then.

    A grid at `--sep-pitch` over a tenth of the angle net will not, on a converged leaf, land on
    the pose the exact sweep finds: on `runs/branch_L1110f_probe.txt` (exact minimum 0.636114) the
    grid alone reports 0.950 at 0.4 deg and still 0.688 with every bin scanned, while evaluating
    the exact sweep's own worst pose gives 0.636114 to six decimals.  The worst pose barely moves
    between rounds, so re-descending from it costs one small window per seed -- the atoms of a bin
    are sorted once and each window is a few hundred cells -- and it recovers the resolution that
    a global grid could only buy at O(pitch^-2).

    `seeds` are (cx, cy, theta, h, flag) rows; the bin is recovered from theta."""
    if len(seeds) == 0: return 9e9, []
    CQ = CQ or []
    byk = {}
    for row in seeds:                       # (cx, cy, theta, h, flag[, credit]) -- see branch.read_witnesses
        cx, cy, th = row[0], row[1], row[2]
        k = int(round(N * math.tan(th / 2.0)))
        byk.setdefault(k, []).append((cx, cy))
    mn = 9e9; out = []
    for k, ps in byk.items():
        fr = bin_frame(s, N, k)
        if fr is None: continue
        ct, st, h, th, lo, hi = fr
        lam4 = None if lam is None else (list(lam) if len(lam) == 4 else [lam[0]] * 4)
        pre = bin_atoms(P, w, ct, st)
        for (cx, cy) in ps:
            u0 = cx * ct + cy * st; u1 = -cx * st + cy * ct
            wid = pitch * half; fp = pitch; best = None
            for _ in range(levels):
                fp = fp / rfac
                g0 = np.arange(u0 - wid, u0 + wid + fp * 0.5, fp)
                g1 = np.arange(u1 - wid, u1 + wid + fp * 0.5, fp)
                if len(g0) < 2 or len(g1) < 2: break
                lv, lf = grid_values(P, w, s, ct, st, h, h, g0, g1, CQ, lam4, r, fp, lo, hi, pre)
                i2, j2 = np.unravel_index(int(np.argmin(lv)), lv.shape)
                if lv[i2, j2] >= 9e8: break
                u0 = float(g0[i2]); u1 = float(g1[j2]); wid = fp
                val = float(lv[i2, j2])
                if val < mn: mn = val
                best = (val, u0, u1, int(lf[i2, j2]))
            # one row per seed -- the descent's own end point.  Emitting every level as well
            # multiplied the row count by three for no new information.
            if best is not None and best[0] < thr:
                val, u0, u1, fl = best
                out.append((val, float(u0 * ct - u1 * st), float(u0 * st + u1 * ct), th, h, fl))
    out.sort()
    return mn, out


def _chunk(args):
    P, w, s, N, ks, pitch, CQ, lam, r, topk, block, thr, refine = args
    mn = 9e9; rows = []
    for k in ks:
        m, o = scan_bin(P, w, s, N, k, pitch, CQ, lam, r, topk, block, thr, refine=refine)
        if m is not None and m < mn: mn = m
        rows += o
    return mn, rows


def separate(P, w, s, N, ks, topk=6, pitch=0.004, cliques=None, lam=None, r=1.0,
             block=8, thr=1.0, pool=None, nproc=1, refine=0, keep_value=False, seeds=None):
    """dense float separation over the bins `ks`.  Returns (min value, rows, info) with
    rows = [(cx, cy, theta_k, h, flag), ...] worst first."""
    CQ = prep_cliques(cliques or [])
    P = np.ascontiguousarray(P, dtype=float); w = np.ascontiguousarray(w, dtype=float)
    if pool is not None and len(ks) > 8:
        nk = getattr(pool, '_processes', nproc)
        chunks = [list(ks[i::nk]) for i in range(nk)]
        outs = pool.map(_chunk, [(P, w, s, N, ch, pitch, CQ, lam, r, topk, block, thr, refine) for ch in chunks if ch])
    else:
        outs = [_chunk((P, w, s, N, list(ks), pitch, CQ, lam, r, topk, block, thr, refine))]
    mn = min(o[0] for o in outs)
    rows = sum((o[1] for o in outs), [])
    nseed = 0
    if seeds is not None and len(seeds):
        smn, srows = descend(P, w, s, N, seeds, pitch, CQ, lam, r, thr)
        nseed = len(srows)
        if smn < mn: mn = smn
        rows = rows + srows
    rows.sort()
    info = dict(bins=len(ks), pitch=pitch, viol=len(rows), cliques=len(CQ), seeded=nseed)
    if keep_value: return mn, rows, info
    return mn, [(cx, cy, th, h, fl) for (v, cx, cy, th, h, fl) in rows], info


# ------------------------------------------------------------------ certificate reader (self-test)
def read_cert(path):
    """a whole certificate -- atoms, anchor cliques, region -- in the form `separate` wants.
    Returns (s, P, w, cliques, lam, r, sym) with one `cliques` entry per FILE clique (which is how
    the verifier credits them) and `lam` in weight units (or None for a plain certificate)."""
    from fractions import Fraction as Fr
    t = open(path).read().split()
    i = 0
    def nx():
        nonlocal i
        v = t[i]; i += 1; return v
    sn = int(nx()); sd = int(nx()); D = int(nx()); WD = int(nx()); m = int(nx())
    s = sn / sd
    A = np.empty((m, 2)); W = np.empty(m)
    for j in range(m):
        A[j, 0] = int(nx()) / D; A[j, 1] = int(nx()) / D; W[j] = int(nx()) / WD
    cliques = []
    if i < len(t) and t[i] == 'cliques':
        raise SystemExit("floatsep: box-clique blocks are not supported (anchor cliques only)")
    if i < len(t) and t[i] == 'anchors':
        nx()
        na = int(nx()); nc = int(nx())
        anc = []
        for _ in range(na):
            kind = nx()
            if kind == 'anchorP':
                x = int(nx()); y = int(nx()); d = int(nx()); anc.append(('P', Fr(x, d), Fr(y, d)))
            else:
                x0 = int(nx()); y0 = int(nx()); x1 = int(nx()); y1 = int(nx()); d = int(nx())
                anc.append(('S', Fr(x0, d), Fr(y0, d), Fr(x1, d), Fr(y1, d)))
        for _ in range(nc):
            wc = int(nx()) / WD; npc = int(nx()); pieces = []
            for _ in range(npc):
                assert nx() == 'piece'
                a = int(nx()); rr = int(nx()); pieces.append((a, tuple(int(nx()) for _ in range(rr))))
            cliques.append((wc, [(tuple(anc), tuple(pieces))]))
    lam = None; r = 1.0
    if i < len(t) and t[i] == 'region':
        nx(); assert nx() == 'corner'
        rn = int(nx()); rd = int(nx()); r = rn / rd
        assert nx() == 'lambda'
        ls = []
        while t[i] != 'k': ls.append(int(nx()) / WD)
        nx()
        ks = [int(v) for v in t[i:]]
        lam = ls * 4 if len(ls) == 1 else ls
    # D4 symmetry of the atom multiset (the verifier's own test), and no cliques
    sd_ = np.round(A * D).astype(np.int64); K = int(round(s * D))
    key = lambda X: sorted(map(tuple, np.c_[X, W]))
    sym = (key(sd_) == key(np.c_[K - sd_[:, 0], sd_[:, 1]]) and key(sd_) == key(np.c_[sd_[:, 0], K - sd_[:, 1]])
           and key(sd_) == key(np.c_[sd_[:, 1], sd_[:, 0]]))
    sym = sym and not cliques and (lam is None or len(set(lam)) == 1)
    return s, A, W, cliques, lam, r, sym


def main():
    import argparse, time, multiprocessing as mp
    ap = argparse.ArgumentParser(description="dense float separation on a certificate (self-test / benchmark)")
    ap.add_argument('cert'); ap.add_argument('--N', type=int, default=2000)
    ap.add_argument('--pitch', type=float, default=0.004)
    ap.add_argument('--dtheta', type=float, default=0.4, help='angle pitch in degrees')
    ap.add_argument('--topk', type=int, default=6); ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--refine', type=int, default=0, help='re-scan the worst refine*topk grid cells at pitch/4 (default 0 = off: it cost 2.5x and moved nothing; --seeds is the local search that pays)')
    ap.add_argument('--seeds', default=None, help="a witness file (v theta cx cy flag) to descend from as well -- in the loop these are last round's rows")
    ap.add_argument('--out', default=None, help='write the rows as a witness file (v theta cx cy flag)')
    a = ap.parse_args()
    t0 = time.time()
    s, P, W, cliques, lam, r, sym = read_cert(a.cert)
    ks = bin_list(a.N, sym, math.radians(a.dtheta))
    t1 = time.time()
    pool = mp.get_context('fork').Pool(a.threads) if a.threads > 1 else None
    seeds = []
    if a.seeds:
        for line in open(a.seeds):
            q = line.split()
            if len(q) >= 4: seeds.append((float(q[2]), float(q[3]), float(q[1]), 0.0, int(q[4]) if len(q) > 4 else 0))
    mv, rows, info = separate(P, W, s, a.N, ks, topk=a.topk, pitch=a.pitch, cliques=cliques,
                              lam=lam, r=r, pool=pool, nproc=a.threads, refine=a.refine,
                              keep_value=bool(a.out), seeds=seeds)
    if pool: pool.close(); pool.join()
    print(f"{a.cert}: atoms={len(P)} cliques={len(cliques)} sym={sym} bins={len(ks)}/{a.N} pitch={a.pitch} dtheta={a.dtheta}deg "
          f"-> float_min={mv:.7f} rows={len(rows)} seeds={len(seeds)} seeded_rows={info.get('seeded', 0)} read={t1-t0:.1f}s scan={time.time()-t1:.1f}s")
    if a.out:
        with open(a.out, 'w') as f:
            for (v, cx, cy, th, h, fl) in rows: f.write(f"{v} {th} {cx} {cy} {fl}\n")


if __name__ == '__main__':
    main()
