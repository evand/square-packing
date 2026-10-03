"""Independent checks of the k2m4 record (leaves reviewer, 2026-10-03).  No code imported from public/.
1. own parse of the box file: total, D4 invariance (as a measure: per-line density profiles), duplicates.
2. root set == the 16200 D4 root boxes tiling [0,9/2]^2 x [0,1/2]; per-root stats == leaf counts; totals == .out.
3. per root: the dumped leaves + top-gap slabs form a guillotine partition of the root (own reconstruction,
   split planes taken from leaf faces, not from the checker's rule); every gap slab and every EMPTY leaf has no
   admissible pose (own exact test with h(u) = (cos+sin)/2, rational in u); kind-specific sanity
   (AXIS: u1 == 0; EXACT0: u0 == 0; SYM: u0 >= sqrt2-1; EXACT45: straddles; nothing else certifies past u0 > sqrt2-1).
"""
import json, sys, hashlib
from fractions import Fraction as F
from collections import Counter, defaultdict

PUB = '/home/evand/math/square-packing/public/s12'
BOX = PUB + '/runs/qx2_k4x_k008/sol_exact_box9.txt'
REC = PUB + '/runs/qx2_k4x_k008/qxzm_full.jsonl'
M = F(9)


def parse(path):
    toks = []
    for ln in open(path):
        ln = ln.split('#')[0]
        toks += ln.split()
    assert toks[0] == 'mixed' and toks[1] == '1'
    it = iter(toks[2:])
    nx = lambda: int(next(it))
    sn, sd = nx(), nx(); D = nx(); W = nx(); npt = nx()
    pts = [(F(nx(), D), F(nx(), D), F(nx(), W)) for _ in range(npt)]
    ns = nx()
    segs = [(F(nx(), D), F(nx(), D), F(nx(), D), F(nx(), D), F(nx(), W)) for _ in range(ns)]
    npg = nx()
    polys = []
    for _ in range(npg):
        k = nx(); w = F(nx(), W)
        V = [(F(nx(), D), F(nx(), D)) for _ in range(k)]
        polys.append((w, V))
    rest = list(it)
    assert not rest, rest[:5]
    return F(sn, sd), pts, segs, polys


def line_profiles(segs):
    """key ('H', y) / ('V', x) -> list of (t0, t1, density)"""
    L = defaultdict(list)
    for x0, y0, x1, y1, w in segs:
        if y0 == y1:
            t0, t1 = sorted((x0, x1)); key = ('H', y0)
        elif x0 == x1:
            t0, t1 = sorted((y0, y1)); key = ('V', x0)
        else:
            raise ValueError('non-axis segment')
        assert t1 > t0
        L[key].append((t0, t1, w / (t1 - t0)))
    return L


def canon(prof):
    """piecewise-constant density as a canonical tuple of (b_i, b_{i+1}, rho) with merged equal neighbours"""
    bs = sorted(set([a for a, _, _ in prof] + [b for _, b, _ in prof]))
    out = []
    for lo, hi in zip(bs, bs[1:]):
        r = sum((d for a, b, d in prof if a <= lo and b >= hi), F(0))
        if r == 0: continue
        if out and out[-1][1] == lo and out[-1][2] == r: out[-1] = (out[-1][0], hi, r)
        else: out.append((lo, hi, r))
    return tuple(out)


def part1():
    m, pts, segs, polys = parse(BOX)
    print('container', m, 'points', len(pts), 'segments', len(segs), 'polygons', len(polys))
    assert m == M and not pts and len(polys) == 1
    tot = sum(s[4] for s in segs) + polys[0][0]
    D = F(214770225571, 200000000000)
    print('total', tot, float(tot), ' == 81 - 4D:', tot == 81 - 4 * D, ' < 77:', tot < 77)
    w, V = polys[0]
    xs = sorted(set(v[0] for v in V)); ys = sorted(set(v[1] for v in V))
    print('polygon', V, 'mass', w, 'area', (xs[1] - xs[0]) * (ys[1] - ys[0]))
    assert w == (xs[1] - xs[0]) ** 2 and xs == ys == [F(14, 5), F(31, 5)]
    # orientation ccw
    area2 = sum(V[i][0] * V[(i + 1) % 4][1] - V[(i + 1) % 4][0] * V[i][1] for i in range(4))
    print('polygon signed 2*area', area2)
    assert all(s[4] > 0 for s in segs), 'nonpositive mass'
    for x0, y0, x1, y1, _ in segs:
        assert all(0 <= v <= m for v in (x0, y0, x1, y1))
    c = Counter((s[0], s[1], s[2], s[3]) for s in segs)
    cn = Counter(tuple(sorted([(s[0], s[1]), (s[2], s[3])])) for s in segs)
    dup = {k: v for k, v in cn.items() if v > 1}
    print('segments listed more than once (as point sets):', len(dup))
    byline = Counter()
    for (p, q), v in dup.items():
        key = ('H', p[1]) if p[1] == q[1] else ('V', p[0])
        byline[key] += 1
    print('  by line:', dict(byline))
    # unit segments?
    lens = Counter(max(abs(s[2] - s[0]), abs(s[3] - s[1])) for s in segs)
    print('segment lengths:', dict(lens))
    L = line_profiles(segs)
    P = {k: canon(v) for k, v in L.items()}
    print('lines', len(P), sorted(set(k[1] for k in P)))
    # D4: the 8 maps of [0,9]^2; check each maps the line measure onto itself and U onto itself
    def img(key, prof, g):
        o, p = key
        # g: (swap, fx, fy) applied as: (x, y) -> swap? (y, x) ; then x -> 9-x if fx, y -> 9-y if fy
        out = []
        for a, b, r in prof:
            if o == 'H':
                pts_ = [(a, p), (b, p)]
            else:
                pts_ = [(p, a), (p, b)]
            q = []
            for x, y in pts_:
                if g[0]: x, y = y, x
                if g[1]: x = M - x
                if g[2]: y = M - y
                q.append((x, y))
            if q[0][1] == q[1][1]:
                k2 = ('H', q[0][1]); t = sorted((q[0][0], q[1][0]))
            else:
                k2 = ('V', q[0][0]); t = sorted((q[0][1], q[1][1]))
            out.append((k2, (t[0], t[1], r)))
        return out
    ok = True
    for g in [(s, fx, fy) for s in (0, 1) for fx in (0, 1) for fy in (0, 1)]:
        Q = defaultdict(list)
        for key, prof in P.items():
            for k2, piece in img(key, prof, g): Q[k2].append(piece)
        Qc = {k: canon(v) for k, v in Q.items()}
        if Qc != P: ok = False; print('D4 FAIL for', g)
    print('D4 invariance of the line measure (8 maps):', ok, '; U = [14/5,31/5]^2 symmetric:', F(14, 5) + F(31, 5) == M)
    # where the doubled segments sit and the summed density there
    for key in sorted(byline):
        print('  line', key, 'profile', [(str(a), str(b), float(r)) for a, b, r in P[key]][:12])
    return P


# ------------------------------------------------------------------ part 2/3: the record
def h(u):
    """half of cos + sin at theta = 2 atan u"""
    return (1 - u * u + 2 * u) / (2 * (1 + u * u))


def U45(u):  # sign of u^2 + 2u - 1 (u vs tan 22.5)
    v = u * u + 2 * u - 1
    return (v > 0) - (v < 0)


def gap_inadmissible(x0, x1, y0, y1, ua, ub, open_lo=True):
    """no admissible pose with centre in [x0,x1]x[y0,y1] (all <= 9/2) and u in (ua, ub] (or [ua, ub])."""
    Mx = min(x1, y1)
    assert max(x1, y1) <= M / 2
    # admissible at u iff h(u) <= min(cx, cy) for some centre, i.e. h(u) <= Mx (cx <= 9 - h automatic for cx <= 9/2)
    # h is increasing on [0, sqrt2-1], decreasing after; inf over the interval at its ends
    if not (h(ub) > Mx): return False
    if open_lo:
        if h(ua) > Mx: return True
        if h(ua) == Mx and U45(ua) < 0: return True
        return False
    return h(ua) > Mx


class Fail(Exception): pass


def check_root(root, leaves, gaps):
    """leaves: list of (box, kind).  Guillotine reconstruction."""
    def rec(node, ls):
        x0, x1, y0, y1, u0, u1 = node
        if not ls:
            # no leaf: the closed node must hold no admissible pose (degenerate EMPTY leaves are dropped below)
            if not gap_inadmissible(x0, x1, y0, y1, u0, u1, open_lo=False): raise Fail(('empty node admissible', node))
            gaps.append((node, None)); return
        for b, k in ls:
            if not (x0 <= b[0] < b[1] <= x1 and y0 <= b[2] < b[3] <= y1 and u0 <= b[4] <= b[5] <= u1):
                raise Fail(('leaf outside node', node, b))
        top = max(b[5] for b, k in ls)
        if top < u1:
            if not gap_inadmissible(x0, x1, y0, y1, top, u1): raise Fail(('gap admissible', node, top))
            gaps.append((node, top))
            u1 = top; node = (x0, x1, y0, y1, u0, u1)
        if len(ls) == 1:
            b, k = ls[0]
            if tuple(b) != node: raise Fail(('single leaf != node', node, b))
            return
        # find a guillotine plane
        for ax in (0, 1, 2):
            lo, hi = node[2 * ax], node[2 * ax + 1]
            cands = sorted(set(b[2 * ax] for b, k in ls if lo < b[2 * ax] < hi))
            mid = (lo + hi) / 2
            if mid in cands: cands.remove(mid); cands.insert(0, mid)
            for v in cands:
                A = [(b, k) for b, k in ls if b[2 * ax + 1] <= v]
                B = [(b, k) for b, k in ls if b[2 * ax] >= v]
                if len(A) + len(B) == len(ls):
                    na = list(node); na[2 * ax + 1] = v
                    nb = list(node); nb[2 * ax] = v
                    rec(tuple(na), A); rec(tuple(nb), B)
                    return
        raise Fail(('no guillotine plane', node, len(ls)))
    sys.setrecursionlimit(10000)
    # degenerate EMPTY leaves (u0 == u1) certify nothing (their slice is inadmissible, checked separately); drop them
    # and require every leaf-free region to be inadmissible instead
    rec(root, [(b, k) for b, k in leaves if not (k == 'EMPTY' and b[4] == b[5])])


def part23():
    fh = open(REC)
    hdr = json.loads(fh.readline())
    print('header shas', hdr['sha256'])
    print('argv', ' '.join(hdr['argv']))
    roots = []; tot = Counter(); kinds = Counter(); bad = []; ngaps = 0; vol = F(0)
    cpu = {}; maxd = 0
    sanity = Counter()
    for ln in fh:
        d = json.loads(ln)
        root = tuple(F(v) for v in d['root'])
        roots.append(root)
        st = d['st']
        for k, v in st.items():
            if k == 'maxdepth': maxd = max(maxd, v)
            else: tot[k] += v
        cpu[root] = st['cpu']
        assert d['unc'] == []
        leaves = [(tuple(F(v) for v in b), k) for b, k in d['leaves']]
        kc = Counter(k for b, k in leaves)
        for k in kc:
            if kc[k] != st.get(k, -1): bad.append(('count', root, k, kc[k], st.get(k)))
        nleaf = len(leaves)
        if st['boxes'] != 2 * nleaf - 1: bad.append(('boxes != 2 leaves - 1', root, st['boxes'], nleaf))
        for b, k in leaves:
            kinds[k] += 1
            x0, x1, y0, y1, u0, u1 = b
            if k == 'AXIS' and u1 != 0: bad.append(('AXIS u1', root, b))
            if u1 == 0 and k not in ('AXIS', 'EMPTY'): bad.append(('u1=0 non-AXIS', root, b, k))
            if k == 'EXACT0' and u0 != 0: bad.append(('EXACT0 u0', root, b))
            if k == 'SYM' and not U45(u0) >= 0: bad.append(('SYM below 45', root, b))
            if k == 'EXACT45' and not (U45(u1) > 0 and U45(u0) < 0): bad.append(('EXACT45', root, b))
            if k == 'EXACT' and not (u0 > 0 and U45(u1) <= 0): bad.append(('EXACT', root, b))
            if k == 'EMPTY':
                if not gap_inadmissible(x0, x1, y0, y1, u0, u1, open_lo=False): bad.append(('EMPTY admissible', root, b))
            if u1 == u0: sanity['degenerate u ' + k] += 1
            if U45(u0) >= 0 and k not in ('SYM', 'EMPTY'): sanity['past45 ' + k] += 1
            vol += (x1 - x0) * (y1 - y0) * (u1 - u0)
        gaps = []
        try:
            check_root(root, leaves, gaps)
        except Fail as e:
            bad.append(('cover', root, e.args[0]))
        ngaps += len(gaps)
    print('roots', len(roots), 'distinct', len(set(roots)))
    exp = [(F(i, 10), F(i + 1, 10), F(j, 10), F(j + 1, 10), F(k, 16), F(k + 1, 16))
           for i in range(45) for j in range(45) for k in range(8)]
    print('root set == 45x45x8 grid on [0,9/2]^2 x [0,1/2]:', set(roots) == set(exp), len(exp))
    print('totals', {k: tot[k] for k in ('LEB', 'CAP', 'EXACT', 'EXACT0', 'EXACT45', 'AXIS', 'SYM', 'PIECE', 'ADM',
                                          'P1', 'MIX', 'SPLIT', 'EMPTY', 'UNCERT', 'boxes')}, 'maxdepth', maxd,
          'cpu', round(sum(cpu.values())))
    print('leaf kinds from dumps', dict(kinds))
    print('gap slabs verified inadmissible:', ngaps)
    print('sanity', dict(sanity))
    print('BAD', len(bad))
    for b in bad[:30]: print('  ', b)
    json.dump(sorted(([str(v) for v in r], c) for r, c in cpu.items()), open('root_cpu.json', 'w'))


if __name__ == '__main__':
    part1()
    part23()
