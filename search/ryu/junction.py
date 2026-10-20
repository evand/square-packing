"""Filled-end L-shape (r1prime.py) with the top band flipped so that both thick ends meet at the junction corner,
and one merged axis fill there.  Exact rationals; same frames as r1prime.py.

Layout (S = m + w, m = k - b): grid [0, m]^2; right band in [m, S] x [0, S] as in r1prime; the top band is the
band of length m mapped by the rotation P(x', y') = (y', S - x') (r1prime used the transposition (y', m + x')), so
its top end, thickest at band-frame x' = 0, now sits under the ceiling at x = m, next to the right band's top end
(thickest at x = m too).  Fills: the two bottom ends as in r1prime (the top band's lands at the left wall), and a
merged corner fill replacing the two separate top fills:
  * rows r = 0..nr-1, [X, X+1] x [S-1-r, S-r], X = start + q, crossing x = m; the part right of x = m stays above
    the right band's last row (V1 >= R1 - sig^2 at the corner (X+1, Y)), the part left of it above the top band's
    (V2 >= R2 - sig^2 at the corner (X, Y), V2 = -(S - y) sig + x gam);
  * columns X = m-1-c below the rows, [X, X+1] x [Y, Y+1], m <= Y, Y + 1 <= S - nr, V2 >= R2 - sig^2 at (X, Y).
nr and the row start (merged, or x = m, i.e. the separate fills) are chosen to maximize the count.

  python3 search/ryu/junction.py table KMIN KMAX        best b per k, old c (r1prime) vs new c
  python3 search/ryu/junction.py check K [B]             build, exact pairwise check of the whole packing
  python3 search/ryu/junction.py localcheck K [B]        exact check of every pair involving a corner-fill square
"""
import sys, math
from fractions import Fraction as F
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from r1prime import trig, band_counts, count as count_old, corners, separated

floor = math.floor


def params(k, b, t=None):
    if t is None:
        t = F(1, b) + F(1, 10 ** 12)
    gam, sig = trig(t)
    w = b * gam + sig
    m = k - b
    return t, gam, sig, w, m, m + w


def corner_fill(k, b, t=None):
    """best merged corner fill: list of lower-left corners (X, Y) of axis squares."""
    t, gam, sig, w, m, S = params(k, b, t)
    R1 = band_counts(b, S, t)[0]
    R2 = band_counts(b, F(m), t)[0]
    c1 = R1 - sig * sig
    c2 = R2 - sig * sig
    best = None
    for merged in (True, False):
        for nr in range(0, 6):
            sq = []
            for r in range(nr):
                Y = S - 1 - r
                if Y < m:
                    break
                rb = min(S, m + (Y * gam - c1) / sig)          # X + 1 <= rb
                xl = (c2 + (S - Y) * sig) / gam                 # X >= xl if the square reaches x < m
                start = max(min(xl, F(m)), F(0)) if merged else F(m)
                n = floor(rb - start)
                sq += [(start + q, Y) for q in range(max(n, 0))]
            top = S - nr
            c = 0
            while m - 1 - c >= 0:
                X = F(m - 1 - c)
                yl = max(S - (X * gam - c2) / sig, F(m))
                n = floor(top - yl)
                if n <= 0:
                    break
                sq += [(X, top - 1 - i) for i in range(n)]
                c += 1
            if best is None or len(sq) > len(best):
                best = sq
    return best


def separate_tops(k, b, t=None):
    """r1prime's top-fill count for both bands (what the corner fill replaces)."""
    t, gam, sig, w, m, S = params(k, b, t)
    return band_counts(b, S, t)[2] + band_counts(b, F(m), t)[2]


def count_new(k, b, t=None):
    r = count_old(k, b, t)
    if r is None:
        return None
    return r[0] - separate_tops(k, b, t) + len(corner_fill(k, b, t))


def best_b(k, fn):
    res = None
    for b in range(2, k - 3):
        N = fn(k, b)
        if N and (res is None or N > res[1]):
            res = (b, N)
    return res


def squares(k, b, t=None):
    """(cx, cy, c, s) for the whole packing, container [0, S]^2."""
    t, gam, sig, w, m, S = params(k, b, t)
    out = [(F(2 * i + 1, 2), F(2 * j + 1, 2), F(1), F(0)) for i in range(m) for j in range(m)]

    def rows(L):
        R, nb, _, _ = band_counts(b, L, t)
        sq = []
        for j in range(R):
            for i in range(b):
                a = i + F(1, 2) + j * (sig / gam)
                c2 = j + F(1, 2)
                sq.append((sig + a * gam - c2 * sig, a * sig + c2 * gam, gam, sig))
        for q in range(nb):                                   # bottom fill, as r1prime
            sq.append((w - 1 - q + F(1, 2), F(1, 2), F(1), F(0)))
        return sq

    for (x, y, c, s) in rows(S):                              # right band
        out.append((m + x, y, c, s))
    for (x, y, c, s) in rows(F(m)):                           # top band, rotated: (x', y') -> (y', S - x')
        out.append((y, S - x, s, -c))
    fill = [(X + F(1, 2), Y + F(1, 2), F(1), F(0)) for (X, Y) in corner_fill(k, b, t)]
    return out, fill, S


def check_pairs(sq, S, focus=None):
    for q in sq + (focus or []):
        for (X, Y) in corners(q):
            assert 0 <= X <= S and 0 <= Y <= S, ('wall', q)
    allsq = sq + (focus or [])
    if focus is None:
        targets = range(len(allsq))
    else:
        targets = range(len(sq), len(allsq))
    # spatial hash on centres (cell 1): overlapping unit squares have centres within sqrt(2)
    cells = {}
    for idx, q in enumerate(allsq):
        cells.setdefault((floor(q[0]), floor(q[1])), []).append(idx)
    for i in targets:
        p = allsq[i]
        cx, cy = floor(p[0]), floor(p[1])
        for dx in (-2, -1, 0, 1, 2):
            for dy in (-2, -1, 0, 1, 2):
                for j in cells.get((cx + dx, cy + dy), ()):
                    if j == i or (focus is None and j < i):
                        continue
                    q = allsq[j]
                    if (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 >= 2:
                        continue
                    assert separated(p, q), ('pair', p, q)
    return len(allsq)


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'table':
        lo, hi = int(sys.argv[2]), int(sys.argv[3])
        for k in range(lo, hi + 1):
            bo, No = best_b(k, lambda k, b: (count_old(k, b) or (None,))[0])
            bn, Nn = best_b(k, count_new)
            print(k, bo, k * k - No, bn, k * k - Nn, (k * k - No) - (k * k - Nn), flush=True)
    else:
        k = int(sys.argv[2])
        b = int(sys.argv[3]) if len(sys.argv) > 3 else best_b(k, count_new)[0]
        base, fill, S = squares(k, b)
        assert S < k
        if cmd == 'check':
            n = check_pairs(base + fill, S)
        else:
            n = check_pairs(base, S, focus=fill)
        assert n == count_new(k, b), (n, count_new(k, b))
        print(f'k={k} b={b}: {n} squares, side {float(S):.12f} < {k}, c = {k*k - n}: VALID ({cmd}, exact)')
