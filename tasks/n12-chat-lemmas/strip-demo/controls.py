"""Negative controls for check.py (all exact unless marked 'float').

    python3 controls.py certs/h113.json

 N1  exact counterexample at h = 12081/10000 (= h* + 0.00099): four pairwise-disjoint closed unit
     squares in [0, 39995/10000]^2 with all centre heights <= h.  For these four poses the chain
     functional, maximised over EVERY choice of separation heights in a fine menu, is < 4, so every
     certificate at this h must fail (K5): leaf bounds are <= true values.
 N2  the given certificate with all separation heights set to a (the single-line argument): REJECT.
 N3  the given certificate with one box deleted: REJECT (coverage).
 N4  the given certificate with h raised by 1/1000: REJECT (coverage).
 N5  the given certificate with a = 0.915 > sqrt2 - 1/2: REJECT (parameters).
 N6  a deliberately wrong leaf: one claim raised above its proven bound: REJECT (leaf); and (float)
     a sampled pose in the box where the true node value is below the wrong claim.
"""
import copy, json, math, random, sys
from fractions import Fraction as Fr
import check
from leaf import CS
import sanity


def rej(cert_dict, label):
    try:
        check.check(cert_dict, quiet=True)
        print(f"  {label}: ACCEPTED  <-- control FAILED")
        return False
    except check.CertError as e:
        print(f"  {label}: REJECT ({e})")
        return True


# ---------- N1: exact counterexample just above h* --------------------------------------------
def poly_square(cx, cy, t):
    C, S = CS(t)
    h = Fr(1, 2)
    return [(cx + sx * h * C - sy * h * S, cy + sx * h * S + sy * h * C) for sx, sy in ((1, 1), (-1, 1), (-1, -1), (1, -1))]


def separated(P, Q):
    """Exact separating-axis test for convex polygons: True iff some edge normal separates strictly."""
    for poly in (P, Q):
        n = len(poly)
        for i in range(n):
            (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % n]
            nx, ny = y2 - y1, x1 - x2
            a = [nx * x + ny * y for x, y in P]
            b = [nx * x + ny * y for x, y in Q]
            if max(a) < min(b) or max(b) < min(a):
                return True
    return False


def half_chords(t, c, z):
    """Exact (lL, lR) for a pose with rational t in (0,1), centre height c, line y = z."""
    C, S = CS(t)
    e = z - c
    if S == 0:
        return Fr(1, 2), Fr(1, 2)
    lR = min((Fr(1, 2) - e * S) / C, (e * C + Fr(1, 2)) / S)
    lL = min((Fr(1, 2) + e * S) / C, (Fr(1, 2) - e * C) / S)
    return lL, lR


def N1():
    h = Fr(12081, 10000)
    hstar_plus = (1 + math.sqrt(2)) / 2 + 0.001
    L = Fr(39995, 10000)
    g, d = Fr(999, 1000), Fr(1, 10000)
    tD = Fr(41421, 100000)                      # theta ~ 44.9997 deg
    xB, xC = 1 + d, 2 + d + g                   # left edges of the 2nd and 3rd floor squares
    xD = (2 + d + xC) / 2                       # diamond above the gap
    # floor squares: axis-parallel t = 0, centres at height 1/2
    sqs = [('A', Fr(1, 2), Fr(1, 2), Fr(0)), ('B', xB + Fr(1, 2), Fr(1, 2), Fr(0)), ('C', xC + Fr(1, 2), Fr(1, 2), Fr(0))]
    # lowest diamond height that clears both corners, found by bisection on exact tests, then checked
    lo, hi = Fr(1), Fr(13, 10)
    for _ in range(60):
        mid = (lo + hi) / 2
        P = poly_square(xD, mid, tD)
        ok = all(separated(P, poly_square(x, y, t)) for _, x, y, t in sqs)
        lo, hi = (lo, mid) if ok else (mid, hi)
    cD = Fr(math.ceil(hi * 10 ** 6), 10 ** 6)
    sqs.append(('D', xD, cD, tD))
    polys = {n: poly_square(x, y, t) for n, x, y, t in sqs}
    names = [s[0] for s in sqs]
    disjoint = all(separated(polys[u], polys[v]) for i, u in enumerate(names) for v in names[i + 1:])
    inside = all(0 <= x <= L and 0 <= y <= L for P in polys.values() for x, y in P)
    heights_ok = all(y <= h for _, _, y, _ in sqs)
    print(f"N1  counterexample at h = {h} = {float(h)} (h* + 0.001 = {hstar_plus:.7f}), container side L = {L}:")
    print(f"    diamond centre height {cD} = {float(cD):.6f}; pairwise disjoint (exact SAT): {disjoint}; "
          f"inside [0,L]^2: {inside}; all centre heights <= h: {heights_ok}")
    assert disjoint and inside and heights_ok
    # chain functional for the left-to-right order A, B, D, C, over every height assignment
    order = [sqs[0], sqs[1], sqs[3], sqs[2]]
    menu = [Fr(k, 1000) for k in range(1, 1000) if h - Fr(1, 2) < Fr(k, 1000) < 1]
    pvals = [(CS(t)[0] + CS(t)[1]) / 2 for _, _, _, t in order]
    sep = []
    for k in range(3):
        _, _, c1, t1 = order[k]
        _, _, c2, t2 = order[k + 1]
        best = max(half_chords(t1, c1, z)[1] + half_chords(t2, c2, z)[0] for z in menu)
        sep.append(best)
    T = pvals[0] + sum(sep) + pvals[3]
    print(f"    max over all separation heights z in (h-1/2, 1) (step 1/1000) of the chain value = {float(T):.6f} < 4: {T < 4}")
    assert T < 4
    return True


def main(path):
    base = json.load(open(path))
    print("certificate under test:", path)
    check.check(base)
    N1()
    ok = True
    c = copy.deepcopy(base)
    ia = [str(z) for z in c['menu']].index(c['a'])
    c['rule'] = [[ia] * len(r) for r in c['rule']]
    # the one-line rule needs (a,a) leaves; supply the exact bounds for them
    from leaf import Box
    for i, key in enumerate(c['boxes']):
        eb = Box(*[Fr(x) for x in key])
        A = Fr(c['a'])
        extra = []
        for l, r in ((ia, ia), (-1, ia), (ia, -1)):
            v = eb.V(None if l == -1 else A, None if r == -1 else A)
            if v is not None:
                extra.append([l, r, str(v)])
        c['leaves'][i] = extra
    print("N2  single-line rule (all separations on y = a), exact leaves:")
    ok &= rej(c, "N2")
    c = copy.deepcopy(base)
    k = max(range(len(c['boxes'])), key=lambda i: Fr(c['boxes'][i][3]) - Fr(c['boxes'][i][2]))
    del c['boxes'][k]; del c['leaves'][k]
    del c['rule'][k]
    for r in c['rule']:
        del r[k]
    print("N3  delete one box:")
    ok &= rej(c, "N3")
    c = copy.deepcopy(base)
    c['h'] = str(Fr(c['h']) + Fr(1, 1000))
    print("N4  raise h by 1/1000 with the same boxes:")
    ok &= rej(c, "N4")
    c = copy.deepcopy(base)
    c['a'] = '183/200'
    c['menu'] = [('183/200' if z == base['a'] else z) for z in c['menu']]
    print("N5  a = 183/200 = 0.915 > sqrt2 - 1/2:")
    ok &= rej(c, "N5")
    # N6: wrong leaf.  Pick the H box's (z, z) leaf with the largest index; set the claim to
    # (smallest sampled true value) + 1/100, which is false at the sampled pose.
    c = copy.deepcopy(base)
    from leaf import Box
    flags = [Box(*[Fr(x) for x in k]).may_be_H(Fr(c['a'])) for k in c['boxes']]
    bi = flags.index(True)
    claims = c['leaves'][bi]
    li = max(range(len(claims)), key=lambda j: (claims[j][0] != -1 and claims[j][1] != -1, j))
    l, r, q = claims[li]
    M = [float(Fr(z)) for z in c['menu']]
    Z = lambda i: None if i == -1 else M[i]
    rng = random.Random(7)
    samples = sanity.sample_box([Fr(x) for x in c['boxes'][bi]], 2000, rng)
    vals = [(sanity.node_value(th, cc, Z(l), Z(r)), th, cc) for th, cc in samples]
    vmin, th, cc = min(vals)
    wrong = Fr(math.floor(vmin * 10 ** 6), 10 ** 6) + Fr(1, 100)
    claims[li] = [l, r, str(wrong)]
    print(f"N6  wrong leaf: box {bi} (may hold H), V({Z(l)},{Z(r)}) claimed {float(wrong):.6f} "
          f"(honest claim {float(Fr(q)):.6f}).  Float witness: pose theta = {math.degrees(th):.4f} deg, "
          f"c = {cc:.6f} has true value {vmin:.6f} < claim.")
    ok &= rej(c, "N6")
    print("all controls behaved as expected" if ok else "SOME CONTROL FAILED")
    return ok


if __name__ == '__main__':
    sys.exit(0 if main(sys.argv[1]) else 1)
