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
    Mirrors `bin_geometry` / the `main` loop of verify/src/main.rs."""
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
    """add every clique's weight to `out` (shape (len(g0), len(g1))) on the grid poses it holds.
    `out` is indexed [i, j] = centre (g0[i], g1[j]) in the bin's rotated frame."""
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


# ------------------------------------------------------------------ the scan
def scan_bin(P, w, s, N, k, pitch, CQ, lam, r, topk, block, thr):
    """one bin: returns (min value seen, [(value, cx, cy, theta_k, h, flag), ...] worst first).
    `value` is captured weight minus the region requirement, so a violation is `value < thr`."""
    ct, st, sg, wm = bin_geometry(N, k)
    h = sg / 2.0
    lo, hi = wm / 2.0, s - wm / 2.0
    if hi <= lo: return None, []
    th = 2.0 * math.atan(k / N)
    # the admissible centre box, rotated into the bin's frame
    corners = np.array([[lo, lo], [hi, lo], [hi, hi], [lo, hi]])
    cu0 = corners[:, 0] * ct + corners[:, 1] * st
    cu1 = -corners[:, 0] * st + corners[:, 1] * ct
    g0 = np.arange(cu0.min(), cu0.max() + pitch * 0.5, pitch)
    g1 = np.arange(cu1.min(), cu1.max() + pitch * 0.5, pitch)
    if len(g0) < 2 or len(g1) < 2: return None, []
    # captured atom weight at every grid centre (prefix sums over a sliding window in u0)
    q0 = P[:, 0] * ct + P[:, 1] * st
    q1 = -P[:, 0] * st + P[:, 1] * ct
    o = np.argsort(q0, kind='stable'); q0s = q0[o]; q1o = q1[o]; wo = w[o]
    ordy = np.argsort(q1o, kind='stable'); q1s = q1o[ordy]; wy = wo[ordy]
    idx_in_q0 = np.arange(len(o))[ordy]
    a = np.searchsorted(q0s, g0 - h + TOL, 'left')
    b = np.searchsorted(q0s, g0 + h - TOL, 'right')
    lo1 = np.searchsorted(q1s, g1 - h + TOL, 'left')
    hi1 = np.searchsorted(q1s, g1 + h - TOL, 'right')
    vals = np.zeros((len(g0), len(g1)))
    for i in range(len(g0)):
        if b[i] <= a[i]: continue
        mask = (idx_in_q0 >= a[i]) & (idx_in_q0 < b[i])
        cw = np.concatenate([[0.0], np.cumsum(np.where(mask, wy, 0.0))])
        vals[i, :] = cw[hi1] - cw[lo1]
    # anchor cliques
    clique_credit(CQ, ct, st, h, g0, g1, vals)
    # container coordinates of every grid centre, admissibility, and the region requirement
    U0 = g0[:, None]; U1 = g1[None, :]
    cx = U0 * ct - U1 * st
    cy = U0 * st + U1 * ct
    ok = (cx >= lo - 1e-12) & (cx <= hi + 1e-12) & (cy >= lo - 1e-12) & (cy <= hi + 1e-12)
    flag = np.zeros(vals.shape, dtype=np.int8)
    if lam is not None:
        loX = cx <= r + 1e-9; hiX = cx >= s - r - 1e-9
        loY = cy <= r + 1e-9; hiY = cy >= s - r - 1e-9
        flag = np.where(loX & loY, 1, np.where(hiX & loY, 2, np.where(loX & hiY, 3, np.where(hiX & hiY, 4, 0)))).astype(np.int8)
        req = np.zeros(vals.shape)
        for j in range(4): req = np.where(flag == j + 1, lam[j] if j < len(lam) else lam[0], req)
        vals = vals - req
    v = np.where(ok, vals, 9e9)
    mn = float(v.min())
    bad = v < thr
    if not bad.any(): return mn, []
    # spread the witnesses out: the worst cell of each block x block tile, best `topk` tiles
    n0, n1 = v.shape
    p0 = (-n0) % block; p1 = (-n1) % block
    VP = np.pad(v, ((0, p0), (0, p1)), constant_values=9e9)
    bl = VP.reshape(VP.shape[0] // block, block, VP.shape[1] // block, block).transpose(0, 2, 1, 3).reshape(-1, block * block)
    bm = bl.min(axis=1); am = bl.argmin(axis=1)
    sb = np.nonzero(bm < thr)[0]
    if len(sb) > topk: sb = sb[np.argsort(bm[sb])[:topk]]
    nb1 = VP.shape[1] // block
    bi, bj = np.divmod(sb, nb1); ii, jj = np.divmod(am[sb], block)
    r0 = bi * block + ii; r1 = bj * block + jj
    keep = (r0 < n0) & (r1 < n1); r0 = r0[keep]; r1 = r1[keep]
    out = []
    for i, j in zip(r0, r1):
        out.append((float(v[i, j]), float(g0[i] * ct - g1[j] * st), float(g0[i] * st + g1[j] * ct), th, h, int(flag[i, j])))
    out.sort()
    return mn, out


def _chunk(args):
    P, w, s, N, ks, pitch, CQ, lam, r, topk, block, thr = args
    mn = 9e9; rows = []
    for k in ks:
        m, o = scan_bin(P, w, s, N, k, pitch, CQ, lam, r, topk, block, thr)
        if m is not None and m < mn: mn = m
        rows += o
    return mn, rows


def separate(P, w, s, N, ks, topk=6, pitch=0.004, cliques=None, lam=None, r=1.0,
             block=8, thr=1.0, pool=None, nproc=1):
    """dense float separation over the bins `ks`.  Returns (min value, rows, info) with
    rows = [(cx, cy, theta_k, h, flag), ...] worst first."""
    CQ = prep_cliques(cliques or [])
    P = np.ascontiguousarray(P, dtype=float); w = np.ascontiguousarray(w, dtype=float)
    if pool is not None and len(ks) > 8:
        nk = getattr(pool, '_processes', nproc)
        chunks = [list(ks[i::nk]) for i in range(nk)]
        outs = pool.map(_chunk, [(P, w, s, N, ch, pitch, CQ, lam, r, topk, block, thr) for ch in chunks if ch])
    else:
        outs = [_chunk((P, w, s, N, list(ks), pitch, CQ, lam, r, topk, block, thr))]
    mn = min(o[0] for o in outs)
    rows = sum((o[1] for o in outs), [])
    rows.sort()
    info = dict(bins=len(ks), pitch=pitch, viol=len(rows), cliques=len(CQ))
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
    ap.add_argument('--out', default=None, help='write the rows as a witness file (v theta cx cy flag)')
    a = ap.parse_args()
    t0 = time.time()
    s, P, W, cliques, lam, r, sym = read_cert(a.cert)
    ks = bin_list(a.N, sym, math.radians(a.dtheta))
    t1 = time.time()
    pool = mp.get_context('fork').Pool(a.threads) if a.threads > 1 else None
    mv, rows, info = separate(P, W, s, a.N, ks, topk=a.topk, pitch=a.pitch, cliques=cliques,
                              lam=lam, r=r, pool=pool, nproc=a.threads)
    if pool: pool.close(); pool.join()
    print(f"{a.cert}: atoms={len(P)} cliques={len(cliques)} sym={sym} bins={len(ks)}/{a.N} pitch={a.pitch} dtheta={a.dtheta}deg "
          f"-> float_min={mv:.7f} rows={len(rows)} read={t1-t0:.1f}s scan={time.time()-t1:.1f}s")
    if a.out:
        with open(a.out, 'w') as f:
            for (cx, cy, th, h, fl) in rows: f.write(f"0 {th} {cx} {cy} {fl}\n")


if __name__ == '__main__':
    main()
