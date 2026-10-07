#!/usr/bin/env python3
"""s60_prep.py -- inputs for the m = 8 line-cover LP (search/S60_COVER.md).  HEURISTIC helpers: nothing here is a proof.

  transplant IN M OUT
      grow a D4-symmetric mixed (v1) cover of [0,m0]^2 (m0 odd) to [0,M]^2 (any M >= m0) by replicating its centre
      column/row of unit cells (h0 = (m0-1)/2, centre strip [h0, h0+1]).  A coordinate X (units 1/D) of a point, or the
      midpoint of a segment piece (along or across its line), maps to
          X                 if X < h0*D   (or X == h0*D)
          X + k*D, k=0..ext if h0*D < X < (h0+1)*D          (ext = M - m0)
          X + ext*D         if X > (h0+1)*D
      and a segment's line coordinate A == (h0+1)*D (the right edge of the centre strip) goes to every A + k*D.
      So corner and wall-band cells are carried over unchanged and the centre cell tiles the new interior.  The map
      commutes with the reflection X -> m0*D - X / M*D - X, so D4 invariance is kept (checked).  Masses unchanged:
      the result is a warm start / column set, NOT a cover (the seams are not checked).
  ring FILE [FILE ...]
      exact (Fraction) savings profile of mixed covers of [0,k]^2: mu of the closed / half-open boxes
      [w, k-w]^2 for w = 0, 1/2, 1, ..., and per-ring deficits (area - mu) with an explicit grid-line convention
      (see S60_COVER.md sec 4).
"""
import sys, os, math
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mixed_cover as MC


def transplant(inp, M, out):
    cv = MC.load(inp); MC.validate(cv)
    D = cv['D']; assert cv['s_den'] == 1, "integer container only"
    m0 = cv['s_num']; assert m0 % 2 == 1, "odd m0 only"
    h0 = (m0 - 1) // 2; ext = M - m0; assert ext >= 0
    lo, hi = h0 * D, (h0 + 1) * D

    def mapc(X2):                    # X2 = 2*coordinate (so segment midpoints stay integral)
        if X2 <= 2 * lo: return [0]
        if X2 >= 2 * hi: return [ext]
        return list(range(ext + 1))

    pts = {}
    for X, Y, W in cv['points']:
        for a in mapc(2 * X):
            for b in mapc(2 * Y):
                k = (X + a * D, Y + b * D); pts[k] = pts.get(k, 0) + W
    segs = {}
    for X0, Y0, X1, Y1, W in cv['segments']:
        if X0 == X1:                 # vertical: line x = X0, along y
            A, T0, T1, o = X0, min(Y0, Y1), max(Y0, Y1), 0
        else:
            A, T0, T1, o = Y0, min(X0, X1), max(X0, X1), 1
        if A <= lo: As = [A]
        elif A >= hi and A != hi: As = [A + ext * D]
        else: As = [A + k * D for k in range(ext + 1)]      # A == hi, or strictly inside the centre strip
        for A2 in As:
            for b in mapc(T0 + T1):
                t0, t1 = T0 + b * D, T1 + b * D
                key = (A2, t0, A2, t1) if o == 0 else (t0, A2, t1, A2)
                segs[key] = segs.get(key, 0) + W
    K = M * D
    new = dict(s_num=M, s_den=1, D=D, W=cv['W'], points=[(x, y, w) for (x, y), w in sorted(pts.items())],
               segments=[(a, b, c, d, w) for (a, b, c, d), w in sorted(segs.items())], polygons=[])
    MC.validate(new)
    # D4 check (rotation by 90 deg: (x,y) -> (K-y, x))
    P = {(x, y): w for x, y, w in new['points']}
    assert all(P.get((K - y, x)) == w for (x, y), w in P.items()), "points not C4-invariant"
    S = {}
    for a, b, c, d, w in new['segments']:
        S[(min(a, c), min(b, d), max(a, c), max(b, d))] = w
    rot = lambda a, b, c, d: (min(K - b, K - d), min(a, c), max(K - b, K - d), max(a, c))
    assert all(S.get(rot(*k)) == w for k, w in S.items()), "segments not C4-invariant"
    refl = lambda a, b, c, d: (K - c, b, K - a, d)
    assert all(S.get(refl(*k)) == w for k, w in S.items()), "segments not reflection-invariant"
    MC.write(out, new, f"s60_prep.py transplant of {os.path.basename(inp)} (s={m0}) to s={M} (centre cell replicated)")
    print(f"transplant {inp} (s={m0}, {len(cv['points'])} pts, {len(cv['segments'])} segs, total {float(MC.total(cv)):.6f}) -> "
          f"{out} (s={M}, {len(new['points'])} pts, {len(new['segments'])} segs, total {float(MC.total(new)):.6f})")


# ------------------------------------------------------------------------------------------------ savings profile
def _fac(v, A, B):
    """weight of coordinate v (integer, units 1/D) for the box side interval [A, B]: 1 strictly inside, 1/2 on an
    end, 0 outside (the 'half-boundary' convention: grid-line mass on a box boundary is split equally between the
    two sides, so ring deficits add up exactly to the total saving)."""
    if A < v < B: return Fraction(1)
    if v == A or v == B: return Fraction(1, 2)
    return Fraction(0)


def box_mass(cv, a, b, mode='half'):
    """exact mu of the axis box [a,b]^2 (Fractions).  mode 'half' (boundary mass x 1/2 per side, corners 1/4),
    'closed' (boundary counts fully), 'open' (boundary excluded).  A segment piece lying across the box counts by
    the length fraction of its part inside (endpoints have no mass)."""
    D, W = cv['D'], cv['W']; A, B = a * D, b * D
    f = _fac
    if mode == 'closed': f = lambda v, A, B: Fraction(1) if A <= v <= B else Fraction(0)
    if mode == 'open': f = lambda v, A, B: Fraction(1) if A < v < B else Fraction(0)
    tot = Fraction(0)
    for X, Y, w in cv['points']:
        g = f(X, A, B)
        if g: tot += g * f(Y, A, B) * w
    for X0, Y0, X1, Y1, w in cv['segments']:
        if X0 == X1: L, t0, t1 = X0, min(Y0, Y1), max(Y0, Y1)
        else: L, t0, t1 = Y0, min(X0, X1), max(X0, X1)
        g = f(L, A, B)
        if not g: continue
        ov = min(B, t1) - max(A, t0)
        if ov > 0: tot += g * Fraction(w * ov, t1 - t0)
    return tot / W


def profile(path):
    """returns dict: k, total, saving, rows [(w, area, mu_half, mu_closed)] for w = 0, 1/2, 1, ..."""
    cv = MC.load(path); k = Fraction(cv['s_num'], cv['s_den']); tot = MC.total(cv)
    rows = []
    for i in range(0, int(k) + 1):
        w = Fraction(i, 2)
        if 2 * w >= k: break
        a, b = w, k - w
        rows.append((w, (b - a) ** 2, box_mass(cv, a, b, 'half'), box_mass(cv, a, b, 'closed')))
    return dict(k=k, total=tot, saving=k * k - tot, rows=rows)


def ring(files):
    for f in files:
        P = profile(f); k = P['k']; R = P['rows']
        print(f"== {f}: k = {k}, total {float(P['total']):.6f}, saving k^2 - total = {float(P['saving']):.6f}")
        print("   box [w,k-w]^2: w, area, deficit (half-boundary), deficit (closed box)")
        for w, ar, mh, mc in R:
            print(f"     w={str(w):4s} area {float(ar):6.2f}  def_half {float(ar - mh):+9.5f}  def_closed {float(ar - mc):+9.5f}")
        D = {w: ar - mh for w, ar, mh, mc in R}
        print("   unit rings [w, w+1] (deficit of [w,k-w]^2 minus that of [w+1,k-w-1]^2; half-boundary convention):")
        for w in sorted(D):
            if w + 1 in D: print(f"     ring [{w},{w + 1}]: {float(D[w] - D[w + 1]):+9.5f}")
            else: print(f"     core [{w},{k - w}]: {float(D[w]):+9.5f}")


def ring_table(specs):
    """specs FILE[:F] -- one line per cover scaled by F (default 1): saving k^2 - F*total, and the half-boundary
    deficits (area - F*mu) of the wall ring [0,1], the interior unit rings [w,w+1] for integer w >= 1, and the
    centre (core [w,k-w] with k-2w in {1,2}).  They add up exactly to the saving."""
    print("k  cover                         F         F*total    saving | wall ring [0,1]  per unit wall | interior rings [1,2],[2,3],... | core | interior sum")
    for sp_ in specs:
        f, F = (sp_.rsplit(':', 1) + ['1'])[:2]; F = Fraction(F)
        P = profile(f); k = P['k']
        dd = {w: ar - F * mh for w, ar, mh, mc in P['rows'] if w.denominator == 1}
        ws = sorted(dd); rings = [dd[w] - dd[w + 1] for w in ws if w + 1 in dd]; core = dd[ws[-1]]
        sav = k * k - F * P['total']
        inner = sum(rings[1:]) + core
        print(f"{k}  {os.path.basename(f):28s} {float(F):.7f} {float(F * P['total']):9.4f}  {float(sav):7.4f} | {float(rings[0]):+8.4f}  "
              f"{float(rings[0]) / (4 * (float(k) - 1)):+.4f} | " + ' '.join(f"{float(r):+.4f}" for r in rings[1:]) +
              f" | {float(core):+.4f} | {float(inner):+.4f}")


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'transplant': transplant(sys.argv[2], int(sys.argv[3]), sys.argv[4])
    elif cmd == 'ring': ring(sys.argv[2:])
    elif cmd == 'ring-table': ring_table(sys.argv[2:])
    else: print(__doc__); sys.exit(1)
