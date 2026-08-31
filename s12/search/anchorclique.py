#!/usr/bin/env python3
"""Anchor-clique columns for the cover LP (task I).

Object: the family of `notes/clique-family.md`,

    K(p, A) = { S : p in S and S meets A }  union  { S : A subseteq S }        (Lemma 0)

written in a certificate as an `anchors` block (`certificates/FORMAT.md`, "Anchor cliques") and
checked by `verify/` and `xcheck.py`.  This module supplies the LP side:

  * `kpa` builds `K(p, A)` for a point `p` within distance 1 of a wall, with `A` the perpendicular
    segment of Lemma 2 (`rho >= rho* = eps p_x / sqrt(1 - p_x^2)`, so that `K(p, A)` contains the
    whole point clique of `p` and the clique column **dominates** the point column of `p`);
  * `images` gives the D4 orbit of a clique (one LP column of cost = number of distinct images,
    exactly as a point orbit);
  * `coeff` decides membership of a *row* of the LP (a pose `(cx, cy, theta_k, h = sigma_k/2)`,
    which is what the verifier's witnesses are) with the SAME conservative predicates the verifier
    uses on a cell -- `contains` by the sigma_k-square (the verifier uses the exact bin core, which
    is larger, so the LP never over-credits itself) and `meets` by the hexagon `A + [-h,h]^2`
    including the segment's own normal (the axis whose omission is the near miss of
    `notes/clique-family.md` 7);
  * `block` writes the certificate block.

Everything that ends up in a certificate is an exact `Fraction` over the certificate's coordinate
denominator; floats are used only to choose candidates and to fill the LP matrix.
"""
import math
from fractions import Fraction as F
import numpy as np

TOL = 1e-9          # a pose within TOL of the boundary of a piece is NOT counted (safe direction)


# ---------------------------------------------------------------------------- the objects
# An anchor is ('P', x, y) or ('S', x0, y0, x1, y1) with Fraction coordinates.
# A clique is (anchors, pieces) with pieces = ((a, (f, ...)), ...).

def _map(g, a):
    if a[0] == 'P':
        return ('P',) + g(a[1], a[2])
    return ('S',) + g(a[1], a[2]) + g(a[3], a[4])


def images(s, cl):
    """the distinct images of a clique under the dihedral group of the container [0,s]^2"""
    gs = [lambda x, y: (x, y), lambda x, y: (s - x, y), lambda x, y: (x, s - y), lambda x, y: (s - x, s - y),
          lambda x, y: (y, x), lambda x, y: (s - y, x), lambda x, y: (y, s - x), lambda x, y: (s - y, s - x)]
    out = []
    for g in gs:
        im = (tuple(_map(g, a) for a in cl[0]), cl[1])
        if im not in out: out.append(im)
    return out


def rho_star(d, eps):
    """Lemma 2: every admissible unit square containing p meets the perpendicular segment at
    offset eps and half-length rho iff rho >= rho* = eps*d/sqrt(1-d^2), d = the wall distance."""
    if not (0 < d < 1): return None
    return eps * d / math.sqrt(1.0 - d * d)


def kpa(s, D, X, Y, wall, eps_n, rho_n):
    """K(p, A) with p = (X/D, Y/D) and A the segment perpendicular to `wall`
    (0 = x low, 1 = x high, 2 = y low, 3 = y high) at offset eps_n/D, half-length rho_n/D.
    Returns (anchors, pieces) or None if the anchor would leave the container."""
    px, py = F(X, D), F(Y, D)
    e, r = F(eps_n, D), F(rho_n, D)
    if wall == 0:   ax = px + e; A = ('S', ax, py - r, ax, py + r)
    elif wall == 1: ax = px - e; A = ('S', ax, py - r, ax, py + r)
    elif wall == 2: ay = py + e; A = ('S', px - r, ay, px + r, ay)
    else:           ay = py - e; A = ('S', px - r, ay, px + r, ay)
    for v in A[1:]:
        if v < 0 or v > s: return None
    # piece 0 = {S : p in S, S meets A}, piece 1 = {S : A subseteq S}; piece 0 filters anchor 1,
    # which is what Lemma 0 needs for the pair (they do not intersect)
    return ((('P', px, py), A), ((0, (1,)), (1, ())))


def cand_params(s, D, X, Y, frac=0.95, slack=0.10, pad=0.002, wallpad=0.0005):
    """the (wall, eps_n, rho_n) of the widest Lemma-2 anchor at p = (X/D, Y/D): eps a fraction
    `frac` of its cap 1 - d, rho = (1+slack) rho* + pad.

    Both pads exist because the verifier's `meets` test is a shade stricter than Lemma 2, which is
    about the *unit* square: it uses the concentric sigma_k-square, and it tests centres down to
    `w_lo/2` (the bin's smallest admissible bound), i.e. a sliver of slightly inadmissible poses.
    `pad` covers that in the rho direction.  `wallpad` covers it in the eps direction, and it is the
    one that matters: with `eps` at its cap `1 - d` the anchor sits exactly on the far edge of the
    unit square of a pose jammed against the wall, so the sigma_k-square of that pose misses it by
    ~1e-4 and the clique then does NOT contain the whole point clique of `p` -- which is the one
    property (Lemma 2) that makes the column dominate the point column.  Measured on the k = 4
    leaf's dual: with `eps` at the cap, one wall-jammed row of dual mass 0.06 fell out of every
    clique -- more, there, than the extra mass the clique gained, and the LP never used a clique
    column.  `wallpad` is deliberately small (it only bites when `frac` is essentially 1): the
    honest way to choose `eps` is to price several of them, since `ybar(K)` already accounts for
    both what a longer reach gains and what it loses, and that is what `anchorsep.separate` does."""
    px, py = X / D, Y / D; sf = float(s)
    out = []
    for wall, d in ((0, px), (1, sf - px), (2, py), (3, sf - py)):
        if not (0.05 < d < 1.0): continue
        eps = min(frac * (1.0 - d), (1.0 - d) - wallpad)
        if eps <= 0: continue
        rho = (1.0 + slack) * rho_star(d, eps) + pad
        if rho >= 0.5: continue                       # a unit square cannot contain the segment
        en = int(math.floor(eps * D)); rn = int(math.ceil(rho * D))
        if en <= 0 or rn <= 0: continue
        out.append((wall, en, rn, d))
    return out


# ---------------------------------------------------------------------------- LP membership
def _rot(rows):
    ct = np.cos(rows[:, 2]); st = np.sin(rows[:, 2])
    return ct, st, rows[:, 0] * ct + rows[:, 1] * st, -rows[:, 0] * st + rows[:, 1] * ct


def _pts(a):
    return [(float(a[1]), float(a[2]))] if a[0] == 'P' else [(float(a[1]), float(a[2])), (float(a[3]), float(a[4]))]


def contains_mask(a, ct, st, u0, u1, h):
    """rows whose sigma_k-square (centre (u0,u1) in the bin frame, half-side h) contains the anchor
    -- conservative: the verifier uses the larger exact bin core, so this never over-credits"""
    ok = np.ones(len(h), dtype=bool)
    for (zx, zy) in _pts(a):
        q0 = zx * ct + zy * st; q1 = -zx * st + zy * ct
        ok &= np.maximum(np.abs(q0 - u0), np.abs(q1 - u1)) <= h - TOL
    return ok


def meets_mask(a, ct, st, u0, u1, h):
    """rows whose sigma_k-square meets the anchor: the cell (here a point) must lie in
    A + [-h,h]^2, i.e. inside the anchor's bounding box grown by h AND between the two
    half-planes of the segment's own normal"""
    P = _pts(a)
    if len(P) == 1:
        return contains_mask(a, ct, st, u0, u1, h)
    (x0, y0), (x1, y1) = P
    q00 = x0 * ct + y0 * st; q01 = -x0 * st + y0 * ct
    q10 = x1 * ct + y1 * st; q11 = -x1 * st + y1 * ct
    ok = (u0 <= np.maximum(q00, q10) + h - TOL) & (u0 >= np.minimum(q00, q10) - h + TOL)
    ok &= (u1 <= np.maximum(q01, q11) + h - TOL) & (u1 >= np.minimum(q01, q11) - h + TOL)
    nx, ny = -(y1 - y0), (x1 - x0)                      # normal in container coordinates
    n0 = nx * ct + ny * st; n1 = -nx * st + ny * ct     # ... and in the bin frame
    ok &= np.abs((u0 - q00) * n0 + (u1 - q01) * n1) <= h * (np.abs(n0) + np.abs(n1)) - TOL
    return ok


def member(cl, rows, rot=None):
    """boolean mask: which rows (poses (cx, cy, theta, h)) are in the clique `cl`"""
    if len(rows) == 0: return np.zeros(0, dtype=bool)
    ct, st, u0, u1 = rot if rot is not None else _rot(rows)
    h = rows[:, 3]
    anchors, pieces = cl
    out = np.zeros(len(h), dtype=bool)
    for (a, filt) in pieces:
        m = contains_mask(anchors[a], ct, st, u0, u1, h)
        for f in filt:
            if not m.any(): break
            m &= meets_mask(anchors[f], ct, st, u0, u1, h)
        out |= m
    return out


def coeff(imgs, rows, rot=None):
    """LP column of a clique orbit: the number of images of the clique containing each row"""
    if len(rows) == 0: return np.zeros(0)
    if rot is None: rot = _rot(rows)
    tot = np.zeros(len(rows))
    for cl in imgs: tot += member(cl, rows, rot)
    return tot


# ---------------------------------------------------------------------------- certificate block
def block(cliques, weights, D):
    """the `anchors` block for a list of cliques (each already an image, i.e. one clique of the
    file) with integer weight numerators over the file's W.  Coordinates are exact over `D`."""
    anchors = []; lines = []
    for cl, w in zip(cliques, weights):
        base = len(anchors)
        for a in cl[0]:
            if a[0] == 'P':
                x, y = F(a[1]) * D, F(a[2]) * D
                assert x.denominator == 1 and y.denominator == 1, "anchor is not a multiple of 1/D"
                anchors.append(f"anchorP {int(x)} {int(y)} {D}")
            else:
                v = [F(t) * D for t in a[1:]]
                assert all(t.denominator == 1 for t in v), "anchor is not a multiple of 1/D"
                anchors.append("anchorS " + " ".join(str(int(t)) for t in v) + f" {D}")
        lines.append((w, [(base + a, tuple(base + f for f in filt)) for (a, filt) in cl[1]]))
    out = [f"anchors {len(anchors)} {len(lines)}"]
    out += anchors
    for w, pieces in lines:
        out.append(f"{w} {len(pieces)}")
        for a, filt in pieces:
            out.append(f"piece {a} {len(filt)}" + ("".join(f" {f}" for f in filt)))
    return "\n".join(out) + "\n"


def key(cl):
    """hashable exact key of a clique"""
    return (tuple((a[0],) + tuple((F(v).numerator, F(v).denominator) for v in a[1:]) for a in cl[0]), cl[1])
