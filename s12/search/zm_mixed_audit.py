#!/usr/bin/env python3
"""Adversarial audit tests for search/zm_mixed.py (search/ZM_MIXED_AUDIT.md), 2026-09-27.

Written by the auditor.  Independent of zm_mixed.py where it matters: the file checks and the exact mass oracle
below read the raw integer file and do their own arithmetic (they do not use zm_mixed.Cover / exact_mass).
The component tests call zm_mixed's internals (piece_bound, cert_split with diag) and compare every certified
COMPONENT bound with the exact mass of that component at rational poses -- a much sharper test than the
leaf stress test (which only looks at the total mass, and cannot see a checker bug on a cover with margin).

  python3 search/zm_mixed_audit.py file FILE                  total, container, D4 invariance (independent)
  python3 search/zm_mixed_audit.py roots JSONL [--m 5 --pitch 1/20 --ubins 16]   root coverage of a run
  python3 search/zm_mixed_audit.py oracle FILE [--n 300]      zm_mixed.exact_mass vs the independent oracle
  python3 search/zm_mixed_audit.py leaves FILE DUMP [--per 40] component test of a leaf dump
  python3 search/zm_mixed_audit.py boxes FILE [--n 200]       piece bound vs exact piece mass, adversarial boxes
"""
import sys, os, math, random, argparse, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from collections import Counter

HALF = F(1, 2)


# ============================================================================ independent raw-file reading / checks
def read_raw(path):
    tok = []
    for line in open(path):
        tok += line.split('#', 1)[0].split()
    assert tok[0] == 'mixed' and tok[1] == '1', 'not a mixed v1 file'
    it = iter(tok[2:])
    nx = lambda: int(next(it))
    sn, sd, D, W = nx(), nx(), nx(), nx()
    pts = [(nx(), nx(), nx()) for _ in range(nx())]
    segs = [(nx(), nx(), nx(), nx(), nx()) for _ in range(nx())]
    polys = []
    for _ in range(nx()):
        k = nx(); w = nx(); polys.append((w, [(nx(), nx()) for _ in range(k)]))
    rest = list(it)
    assert not rest, 'trailing tokens'
    return dict(sn=sn, sd=sd, D=D, W=W, pts=pts, segs=segs, polys=polys)


def file_check(path):
    r = read_raw(path)
    D, W = r['D'], r['W']
    assert (r['sn'] * D) % r['sd'] == 0
    S = r['sn'] * D // r['sd']
    tot = sum(w for *_, w in r['pts']) + sum(s[4] for s in r['segs']) + sum(p[0] for p in r['polys'])
    total = F(tot, W)
    ok = all(0 <= X <= S and 0 <= Y <= S and w >= 0 for X, Y, w in r['pts'])
    ok &= all(0 <= a <= S for s in r['segs'] for a in s[:4]) and all(s[4] >= 0 for s in r['segs'])
    ok &= all((s[0], s[1]) != (s[2], s[3]) for s in r['segs'])
    print(f"container [0,{F(r['sn'], r['sd'])}]^2, D={D}, W={W}: {len(r['pts'])} points, {len(r['segs'])} segments, "
          f"{len(r['polys'])} polygons; all in container, w >= 0, no degenerate segment: {ok}")
    print(f"total = {total} = {float(total):.12f}  (< 21: {total < 21})")
    # D4 invariance of the MEASURE: point measure = coordinate -> summed weight; segment measure: a uniform
    # segment is determined by its unordered endpoints and weight, but the same line measure can be cut into
    # segments differently; we compare the multiset of (unordered endpoints, weight) -- sufficient, not necessary.
    def canon(g):
        P = Counter()
        for X, Y, w in r['pts']:
            if w: P[g(X, Y)] += w
        Sg = Counter()
        for X0, Y0, X1, Y1, w in r['segs']:
            if w: Sg[(min(g(X0, Y0), g(X1, Y1)), max(g(X0, Y0), g(X1, Y1)), w)] += 1
        return P, Sg
    c0 = canon(lambda X, Y: (X, Y))
    gens = {'x->S-x': lambda X, Y: (S - X, Y), 'x<->y': lambda X, Y: (Y, X), 'y->S-y': lambda X, Y: (X, S - Y),
            'rot90': lambda X, Y: (S - Y, X)}
    for name, g in gens.items():
        print(f"  invariant under {name}: {canon(g) == c0}")
    # points exactly on segment-carrying lines (Corollary T' material)
    Vx = {s[0] for s in r['segs'] if s[0] == s[2]}; Hy = {s[1] for s in r['segs'] if s[1] == s[3]}
    onl = [(X, Y) for X, Y, w in r['pts'] if w and (X in Vx or Y in Hy)]
    print(f"  points lying on a segment line: {len(onl)}; lines V {sorted(F(x, D) for x in Vx)}, "
          f"H {sorted(F(y, D) for y in Hy)}")
    return r


# ============================================================================ independent exact mass oracle
def trig(u):
    d = 1 + u * u
    return (1 - u * u) / d, 2 * u / d


def oracle_mass(r, cx, cy, u, parts=False):
    """exact mu(Q) of the raw cover r at the pose (cx, cy, theta = 2 atan u): closed square."""
    D, W = r['D'], r['W']
    c, s = trig(u)
    pm = F(0)
    fx, fy = float(cx), float(cy)
    for X, Y, w in r['pts']:
        if not w: continue
        if abs(X / D - fx) > 0.71 or abs(Y / D - fy) > 0.71: continue
        a = F(X, D) - cx; b = F(Y, D) - cy
        x = a * c + b * s; y = -a * s + b * c
        if -HALF <= x <= HALF and -HALF <= y <= HALF: pm += F(w, W)
    sm = F(0)
    for X0, Y0, X1, Y1, w in r['segs']:
        if not w: continue
        if min(X0, X1) / D > fx + 0.71 or max(X0, X1) / D < fx - 0.71: continue
        if min(Y0, Y1) / D > fy + 0.71 or max(Y0, Y1) / D < fy - 0.71: continue
        ax, ay = F(X0, D) - cx, F(Y0, D) - cy
        dx, dy = F(X1 - X0, D), F(Y1 - Y0, D)
        lo, hi = F(0), F(1)                      # parameter s in [0, 1]
        for (p, q) in ((ax * c + ay * s, dx * c + dy * s), (-ax * s + ay * c, -dx * s + dy * c)):
            # -1/2 <= p + q t <= 1/2
            if q == 0:
                if not (-HALF <= p <= HALF): lo, hi = F(1), F(0)
            else:
                t1, t2 = (-HALF - p) / q, (HALF - p) / q
                if t1 > t2: t1, t2 = t2, t1
                lo, hi = max(lo, t1), min(hi, t2)
        if hi > lo: sm += F(w, W) * (hi - lo)
    assert not r['polys'], 'oracle: polygons not implemented (none in the m5 cover)'
    return (pm, sm) if parts else pm + sm


def admissible(m, cx, cy, u):
    c, s = trig(u); w = c + s
    return w / 2 <= cx <= m - w / 2 and w / 2 <= cy <= m - w / 2


# ============================================================================ run coverage
def roots_check(jsonl, m, pitch, ubins):
    m = F(m); pitch = F(pitch)
    n = int(m / 2 / pitch)
    want = set()
    for i in range(n):
        for j in range(n):
            for k in range(ubins):
                want.add((i * pitch, (i + 1) * pitch, j * pitch, (j + 1) * pitch, F(k, 2 * ubins), F(k + 1, 2 * ubins)))
    got = Counter(); boxes = 0; unc = 0; labels = Counter(); stsum = Counter()
    for ln in open(jsonl):
        if not ln.strip(): continue
        d = json.loads(ln)
        got[tuple(F(v) for v in d['root'])] += 1
        labels[d['label']] += 1
        boxes += d['st']['boxes']; unc += d['st']['UNCERT'] + len(d['unc'])
        for k, v in d['st'].items():
            if k not in ('maxdepth', 'cpu'): stsum[k] += v
    dup = [k for k, v in got.items() if v > 1]
    miss = want - set(got); extra = set(got) - want
    # geometric tiling check (independent of the enumeration): volumes add up and boxes are grid cells
    vol = sum((b[1] - b[0]) * (b[3] - b[2]) * (b[5] - b[4]) for b in got)
    print(f"{sum(got.values())} lines, {len(got)} distinct roots, labels {dict(labels)}; wanted {len(want)}; "
          f"missing {len(miss)}, extra {len(extra)}, duplicated {len(dup)}")
    print(f"volume of the distinct roots {vol} vs region (m/2)^2 * 1/2 = {(m / 2) ** 2 / 2}")
    print(f"boxes {boxes}, uncertified {unc}; leaf census {dict(stsum)}")
    return not miss and not dup and unc == 0


# ============================================================================ helpers on zm_mixed
def load_Z(path):
    import mixed_cover as MC
    import zm_mixed as Z
    cv = MC.load(path)
    return Z, Z.Cover(cv)


def poses_in_box(box, m, n, lines_v=(), lines_h=(), rng=random):
    """admissible rational poses of the box, adversarial: the 8 corners (clamped into admissibility),
    face points, tiny u (1e-9 .. 1e-3 above u0), exact germ/pivot poses c = line + 1/(2 cos) (+- 1e-12),
    square edge exactly on a line at theta = 0, the wall w/2 exactly, and random."""
    cx0, cx1, cy0, cy1, u0, u1 = box
    out = []
    def add(cx, cy, u):
        if not (u0 <= u <= u1): return
        c, s = trig(u); w = c + s
        # clamp into admissibility (stays in the box if the box has admissible poses at this u)
        cx = min(max(cx, w / 2), m - w / 2); cy = min(max(cy, w / 2), m - w / 2)
        if cx0 <= cx <= cx1 and cy0 <= cy <= cy1: out.append((cx, cy, u))
    us = [u0, u1, (u0 + u1) / 2]
    for e in (F(1, 10 ** 9), F(1, 10 ** 12), F(1, 10 ** 6), F(1, 1000)):
        if u0 + e <= u1: us.append(u0 + e)
    for u in us:
        for cx in (cx0, cx1, (cx0 + cx1) / 2):
            for cy in (cy0, cy1, (cy0 + cy1) / 2):
                add(cx, cy, u)
    for _ in range(n):
        q = rng.random()
        u = u0 + (u1 - u0) * F(rng.randint(0, 10 ** 6), 10 ** 6)
        if q < 0.15: u = rng.choice(us)
        c, s = trig(u)
        cx = cx0 + (cx1 - cx0) * F(rng.randint(0, 10 ** 6), 10 ** 6)
        cy = cy0 + (cy1 - cy0) * F(rng.randint(0, 10 ** 6), 10 ** 6)
        t = rng.random()
        if t < 0.25 and lines_v:            # pivot of a vertical germ pair: left edge touches x = xi at y = ...
            xi = rng.choice(lines_v)
            cx = xi + 1 / (2 * c) + rng.choice([0, 0, F(1, 10 ** 12), -F(1, 10 ** 12), F(1, 10 ** 7)])
            cx = xi - 1 / (2 * c) if rng.random() < 0.3 else cx
        elif t < 0.5 and lines_h:
            eta = rng.choice(lines_h)
            cy = eta + 1 / (2 * c) + rng.choice([0, 0, F(1, 10 ** 12), -F(1, 10 ** 12), F(1, 10 ** 7)])
            cy = eta - 1 / (2 * c) if rng.random() < 0.3 else cy
        elif t < 0.6 and lines_v:           # a vertex of Q exactly on a vertical line: cx +- (c + s)/2 = xi
            xi = rng.choice(lines_v); cx = xi + rng.choice([-1, 1]) * (c + s) / 2
        elif t < 0.7 and lines_h:
            eta = rng.choice(lines_h); cy = eta + rng.choice([-1, 1]) * (c + s) / 2
        elif t < 0.75:                      # wall exactly
            cx = (c + s) / 2 if rng.random() < 0.5 else cx
            cy = (c + s) / 2 if rng.random() < 0.5 else cy
        add(cx, cy, u)
    return out


def float_min_piece(r, box, m, starts=6, rng=random):
    """a crude float local search for the minimum PIECE mass over the admissible poses of the box; returns
    candidate rational poses (to be evaluated exactly)."""
    import numpy as np
    cx0, cx1, cy0, cy1, u0, u1 = [float(v) for v in box]
    D, W = r['D'], r['W']
    segs = [(X0 / D, Y0 / D, X1 / D, Y1 / D, w / W) for X0, Y0, X1, Y1, w in r['segs']
            if not (max(X0, X1) / D < cx0 - 0.8 or min(X0, X1) / D > cx1 + 0.8 or
                    max(Y0, Y1) / D < cy0 - 0.8 or min(Y0, Y1) / D > cy1 + 0.8)]
    if not segs: return []
    A = np.array(segs)
    def f(cx, cy, u):
        c, s = (1 - u * u) / (1 + u * u), 2 * u / (1 + u * u)
        ax, ay = A[:, 0] - cx, A[:, 1] - cy; dx, dy = A[:, 2] - A[:, 0], A[:, 3] - A[:, 1]
        lo = np.zeros(len(A)); hi = np.ones(len(A))
        for p, q in ((ax * c + ay * s, dx * c + dy * s), (-ax * s + ay * c, -dx * s + dy * c)):
            with np.errstate(divide='ignore', invalid='ignore'):
                t1 = (-0.5 - p) / q; t2 = (0.5 - p) / q
            z = q == 0
            t1 = np.where(z, np.where(np.abs(p) <= 0.5, -1e9, 1e9), t1)
            t2 = np.where(z, np.where(np.abs(p) <= 0.5, 1e9, -1e9), t2)
            lo = np.maximum(lo, np.minimum(t1, t2)); hi = np.minimum(hi, np.maximum(t1, t2))
        return float((A[:, 4] * np.clip(hi - lo, 0, None)).sum())
    def clampad(cx, cy, u):
        u = min(max(u, u0), u1)
        w = ((1 - u * u) + 2 * u) / (1 + u * u)
        cx = min(max(cx, cx0, w / 2), cx1, m - w / 2); cy = min(max(cy, cy0, w / 2), cy1, m - w / 2)
        return cx, cy, u
    res = []
    for _ in range(starts):
        x = clampad(rng.uniform(cx0, cx1), rng.uniform(cy0, cy1), rng.uniform(u0, u1))
        fx = f(*x); step = [(cx1 - cx0) / 4, (cy1 - cy0) / 4, (u1 - u0) / 4]
        for _ in range(60):
            imp = False
            for i in range(3):
                for sg in (1, -1):
                    y = list(x); y[i] += sg * step[i]; y = clampad(*y); fy = f(*y)
                    if fy < fx: x, fx, imp = y, fy, True
            if not imp: step = [v / 2 for v in step]
        res.append(x)
    out = []
    for cx, cy, u in res:
        out.append((F(cx).limit_denominator(10 ** 12), F(cy).limit_denominator(10 ** 12), F(u).limit_denominator(10 ** 12)))
    return out


# ============================================================================ component tests
def prep(Z, chk, box):
    import zeromargin as zm
    cx0, cx1, cy0, cy1, u0, u1 = box
    if u1 > u0:
        cu1 = zm.clip_bin(box, chk.m)
        if cu1 < u1: box = (cx0, cx1, cy0, cy1, u0, cu1)
    return box, zm.bin_data(box[4], box[5])


def g_exact(kind, px, py, cx, cy, u):
    """zeromargin's violation polynomial G (<= 0 iff inequality `kind` holds), evaluated exactly -- written
    here from the geometry (X = (aC + bS)/N, Y = (-aS + bC)/N), not copied."""
    a, b = px - cx, py - cy
    C, S, N = 1 - u * u, 2 * u, 1 + u * u
    X2 = 2 * (a * C + b * S); Y2 = 2 * (-a * S + b * C)
    return (X2 - N, -X2 - N, Y2 - N, -Y2 - N)[kind]


class Tally:
    def __init__(self): self.n = 0; self.viol = []; self.minslack = {}
    def check(self, tag, lhs, rhs, info):
        """record a check lhs >= rhs"""
        self.n += 1
        d = lhs - rhs
        if tag not in self.minslack or d < self.minslack[tag][0]: self.minslack[tag] = (d, info)
        if d < 0: self.viol.append((tag, float(d), info))


def split_diag(chk, box, B, L, lparts):
    diag = []
    chk.cert_split(box, B, None, L, lparts, diag=diag)
    return diag


def check_box_components(Z, chk, r, box, poses, tally, Lclaim=None, kind=None, mode='full'):
    """At each pose: exact piece mass >= L (this box's piece bound, and >= Lclaim if given: the dumped phantom
    weight, which includes the parent's inherited bound); for leaves certified by point primitives also
    exact point mass + Lclaim >= 1; for SPLIT data every region's claimed points are in Q, its piece bound is
    <= the exact piece mass, a region proved EMPTY contains no admissible pose, and points + region bound >= 1."""
    m = chk.m
    bx, B = prep(Z, chk, box)
    L, _ = chk.piece_bound(bx, B)
    lparts = chk._lparts
    diag = None
    if mode == 'full' and lparts is not None and bx[4] > 0:
        diag = split_diag(chk, bx, B, max(L, Lclaim or 0), lparts)
    P = chk.zc.P
    for (cx, cy, u) in poses:
        if not (bx[0] <= cx <= bx[1] and bx[2] <= cy <= bx[3] and bx[4] <= u <= bx[5]): continue
        if not admissible(m, cx, cy, u): continue
        pm, sm = oracle_mass(r, cx, cy, u, parts=True)
        info = (str(box), str(cx), str(cy), str(u))
        tally.check('piece>=L', sm, L, info)
        if Lclaim is not None:
            tally.check('total>=1 in a certified leaf', pm + sm, 1, info)
            tally.check('piece>=Lclaim', sm, Lclaim, info)
            if kind in ('ADM', 'MIX', 'P1', 'CHAIN', 'SPLIT', 'PIECE'):
                tally.check('pts+Lclaim>=1 (' + kind + ')', pm + Lclaim, 1, info) if kind != 'SPLIT' else None
        if diag:
            for dk in diag:
                kd, chain = dk['kind'], dk['chain']
                rr = 0
                for j, q in enumerate(chain):
                    if g_exact(kd, P[q][0], P[q][1], cx, cy, u) <= 0: rr = j + 1
                reg = [x for x in dk['regions'] if x['r'] == rr]
                if not reg: continue
                reg = reg[0]
                c, s = trig(u)
                for k in reg['pts']:
                    a, b = P[k][0] - cx, P[k][1] - cy
                    x = a * c + b * s; y = -a * s + b * c
                    inq = -HALF <= x <= HALF and -HALF <= y <= HALF
                    tally.check('SPLIT region point in Q', 1 if inq else 0, 1, info + (k, rr))
                if reg['empty']:
                    tally.check('SPLIT EMPTY region has no pose', 0, 1, info + (rr,))
                else:
                    tally.check('SPLIT region pieces <= exact', sm, reg['pieces'], info + (rr,))
                    if kind == 'SPLIT' and dk.get('ok'):
                        tally.check('pts + region pieces >= 1 (SPLIT ok chain)', pm + reg['pieces'], 1, info + (rr,))


_ME = {}
def float_min_total(path, box, m, n=6, top=3):
    """lowest-total poses of a 6x6x6 float grid over the box (then evaluated exactly by the caller)."""
    import numpy as np
    import mixed_cover as MC
    if path not in _ME: _ME[path] = MC.MixedEval.from_cover(MC.load(path), tol=0.0)
    E = _ME[path]
    cx0, cx1, cy0, cy1, u0, u1 = [float(v) for v in box]
    X, Y, U = np.meshgrid(np.linspace(cx0, cx1, n), np.linspace(cy0, cy1, n), np.linspace(u0, u1, n), indexing='ij')
    X, Y, U = X.ravel(), Y.ravel(), U.ravel()
    T = 2 * np.arctan(U); w = np.cos(T) + np.sin(T)
    ok = (X >= w / 2) & (X <= m - w / 2) & (Y >= w / 2) & (Y <= m - w / 2)
    if not ok.any(): return []
    X, Y, U, T = X[ok], Y[ok], U[ok], T[ok]
    v = E.mass(X, Y, T)
    out = []
    for i in np.argsort(v)[:top]:
        out.append((F(X[i]).limit_denominator(10 ** 12), F(Y[i]).limit_denominator(10 ** 12), F(U[i]).limit_denominator(10 ** 12)))
    return out


def parse_dump(path):
    out = []
    for ln in open(path):
        if ln.startswith('#') or not ln.strip(): continue
        lab, bs, kind, wit = [x.strip() for x in ln.split(' ; ', 3)]
        box = tuple(F(v) for v in bs.split())
        Lc = None
        if kind == 'PIECE': Lc = F(wit)
        elif kind not in ('EMPTY', 'UNCERT'):
            Lc = F(wit.split("'")[1])
        out.append((box, kind, Lc))
    return out


def leaves_test(path, dump, per, seed, maxleaves, kinds=None):
    random.seed(seed)
    Z, cov = load_Z(path)
    r = read_raw(path)
    chk = Z.MixedChecker(cov, max_depth=24, use_chain=True, chain_from=0)
    lv = sorted(set(cov.V)); lh = sorted(set(cov.H))
    L = [x for x in parse_dump(dump) if x[1] not in ('EMPTY', 'UNCERT') and (kinds is None or x[1] in kinds)]
    random.shuffle(L)
    L = L[:maxleaves]
    tally = Tally(); t0 = time.time()
    for i, (box, kind, Lc) in enumerate(L):
        ps = poses_in_box(box, chk.m, per, lv, lh)
        ps += [p for p in float_min_piece(r, box, float(chk.m), starts=3)]
        ps += float_min_total(path, box, float(chk.m))
        check_box_components(Z, chk, r, box, ps, tally, Lclaim=Lc, kind=kind,
                             mode='full' if kind == 'SPLIT' else 'lite')
        if (i + 1) % 200 == 0:
            print(f"  {i+1}/{len(L)} leaves, {tally.n} checks, {len(tally.viol)} violations, {time.time()-t0:.0f}s", flush=True)
    report(tally, f"leaves {dump}: {len(L)} leaves")
    return tally


def report(tally, head):
    print(f"{head}: {tally.n} checks, {len(tally.viol)} violations")
    for tag, (d, info) in sorted(tally.minslack.items()):
        print(f"   min slack {tag:45s} {float(d): .3e}   at {info}")
    for v in tally.viol[:20]:
        print("   VIOLATION", v)


def adversarial_boxes(n, rng, m=5):
    """small boxes placed on the geometry the lemmas are delicate about (m = 5, lines at 1..4)."""
    out = []
    for _ in range(n):
        t = rng.random()
        s = F(1, rng.choice([20, 40, 80, 160, 320, 640, 1280, 2560]))
        su = F(1, rng.choice([32, 64, 128, 256, 1024, 4096, 2 ** 16]))
        def al(v, h): return (F(math.floor(v / h))) * h
        if t < 0.3:      # tile / wall germ: centre at (i + 1/2, j + 1/2), u from 0
            cx = F(rng.randint(0, 4)) + HALF + rng.choice([0, -s, -s / 2, s / 3]); cy = F(rng.randint(0, 4)) + HALF + rng.choice([0, -s, -s / 2])
            u0 = F(0) if rng.random() < 0.7 else su * rng.randint(0, 3)
        elif t < 0.45:   # wall band: cx near w/2 (the corner_choices split), all u
            u0 = su * rng.randint(0, int(F(1, 2) / su) - 1)
            c, sn = trig(u0 + su / 2)
            cx = F(float((c + sn) / 2)).limit_denominator(10 ** 6) - s / 2
            cy = F(rng.randint(500, 2500), 1000)
        elif t < 0.65:   # a vertex of Q crossing a grid line (Lemma L' / V)
            u0 = su * rng.randint(1, int(F(1, 2) / su) - 1)
            c, sn = trig(u0)
            cx = F(rng.randint(1, 4)) + rng.choice([-1, 1]) * F(float((c + sn) / 2)).limit_denominator(10 ** 6)
            cy = F(rng.randint(500, 2500), 1000)
        else:            # anywhere, any size
            cx = F(rng.randint(500, 2500), 1000); cy = F(rng.randint(500, 2500), 1000)
            u0 = su * rng.randint(0, int(F(1, 2) / su) - 1)
        box = (cx, cx + s, cy, cy + s, u0, u0 + su)
        if rng.random() < 0.5: box = (box[2], box[3], box[0], box[1], box[4], box[5])
        out.append(box)
    return out


def boxes_test(path, n, per, seed):
    rng = random.Random(seed); random.seed(seed)
    Z, cov = load_Z(path)
    r = read_raw(path)
    chk = Z.MixedChecker(cov, max_depth=24, use_chain=True, chain_from=0)
    lv = sorted(set(cov.V)); lh = sorted(set(cov.H))
    tally = Tally(); t0 = time.time(); used = 0
    for i, box in enumerate(adversarial_boxes(n, rng, int(chk.m))):
        ps = poses_in_box(box, chk.m, per, lv, lh)
        if not ps: continue
        used += 1
        ps += float_min_piece(r, box, float(chk.m), starts=4)
        check_box_components(Z, chk, r, box, ps, tally, mode='full')
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{n} boxes ({used} with poses), {tally.n} checks, {len(tally.viol)} violations, {time.time()-t0:.0f}s", flush=True)
    report(tally, f"adversarial boxes on {path}: {used} boxes")
    return tally


def oracle_test(path, n, seed):
    random.seed(seed)
    Z, cov = load_Z(path)
    r = read_raw(path)
    m = cov.m
    worst = 0; t0 = time.time()
    for i in range(n):
        u = F(random.randint(0, 10 ** 6), 2 * 10 ** 6)
        if i % 5 == 0: u = F(0)
        if i % 7 == 0: u = F(1, 10 ** random.randint(6, 12))
        c, s = trig(u); w = c + s
        cx = w / 2 + (m - w) * F(random.randint(0, 10 ** 6), 10 ** 6)
        cy = w / 2 + (m - w) * F(random.randint(0, 10 ** 6), 10 ** 6)
        if i % 3 == 0: cx = F(random.randint(1, 4)) + HALF + F(random.randint(-3, 3), 10 ** 9)
        if i % 4 == 0: cy = F(random.randint(1, 4)) + 1 / (2 * c)
        if i % 11 == 0: cx = w / 2
        if not admissible(m, cx, cy, u): continue
        a = oracle_mass(r, cx, cy, u); b = Z.exact_mass(cov, cx, cy, u)
        if a != b:
            print("MISMATCH", cx, cy, u, float(a), float(b)); worst += 1
    print(f"oracle vs zm_mixed.exact_mass: {n} poses, {worst} mismatches, {time.time()-t0:.0f}s")



def tprime_test(n, per, seed):
    """Corollary T' (points on germ lines, withheld from the point primitives): a grid cover with points AT the
    grid crossings (on two group lines at once) and along the lines.  Check: the with-points piece bound L2 <=
    exact segment mass + exact mass of the moved points, at adversarial poses; and no point is moved twice."""
    import mixed_cover as MC
    import zm_mixed as Z
    import zm_mixed_test as ZT
    rng = random.Random(seed); random.seed(seed)
    D, W = 1000, 10 ** 6
    pts = []
    for i in range(1, 4):
        for j in range(1, 4):
            pts.append((i * D, j * D, 60000))                        # crossing: on V line i and H line j
        for k in range(0, 4 * 7 + 1):
            t = k * D // 7
            pts.append((i * D, t, 7000)); pts.append((t, i * D, 7000))
    cv = ZT.grid_cover(4, F(1, 2), 10, D, W, pts)
    MC.validate(cv)
    cov = Z.Cover(cv)
    r = dict(sn=4, sd=1, D=D, W=W, pts=cv['points'], segs=cv['segments'], polys=[])
    chk = Z.MixedChecker(cov, max_depth=20, use_chain=True)
    tally = Tally(); nmoved = 0
    for _ in range(n):
        s = F(1, rng.choice([20, 80, 320, 1280]))
        su = F(1, rng.choice([64, 256, 4096, 2 ** 16]))
        cx = F(rng.randint(1, 3)) + HALF + rng.choice([0, -s, -s / 2]); cy = F(rng.randint(0, 3)) + HALF + rng.choice([0, -s, -s / 2])
        u0 = F(0) if rng.random() < 0.7 else su * rng.randint(0, 5)
        box = (cx, cx + s, cy, cy + s, u0, u0 + su)
        bx, B = prep(Z, chk, box)
        moved = []
        L2, _ = chk.piece_bound(bx, B, with_pts=True, moved=moved)
        if len(set(moved)) != len(moved): tally.check('moved twice', 0, 1, str(box))
        nmoved += len(moved)
        ps = poses_in_box(bx, chk.m, per, [F(1), F(2), F(3)], [F(1), F(2), F(3)])
        for (x, y, u) in ps:
            if not admissible(chk.m, x, y, u): continue
            pm, sm = oracle_mass(r, x, y, u, parts=True)
            c, sn = trig(u); mm = F(0)
            for k in moved:
                a, b = cov.points[k][0] - x, cov.points[k][1] - y
                X = a * c + b * sn; Y = -a * sn + b * c
                if -HALF <= X <= HALF and -HALF <= Y <= HALF: mm += cov.points[k][2]
            tally.check("L2 <= segs + moved pts (Cor. T')", sm + mm, L2, (str(box), str(x), str(y), str(u)))
    report(tally, f"Corollary T' test: {n} boxes, {nmoved} point moves")
    return tally


def find_hole(path, cxr, cyr, thr, n=40, top=8, seed=1):
    """float grid search + local refinement of the minimum of mu over a region (cx range, cy range, theta range in
    degrees); the best poses are snapped to rationals and evaluated EXACTLY with the independent oracle."""
    import numpy as np
    import mixed_cover as MC
    rng = random.Random(seed)
    cv = MC.load(path); r = read_raw(path); m = float(F(r['sn'], r['sd']))
    E = MC.MixedEval.from_cover(cv, tol=0.0)
    X, Y, T = np.meshgrid(np.linspace(*cxr, n), np.linspace(*cyr, n), np.radians(np.linspace(*thr, n)), indexing='ij')
    X, Y, T = X.ravel(), Y.ravel(), T.ravel()
    w = np.cos(T) + np.sin(T)
    ok = (X >= w / 2) & (X <= m - w / 2) & (Y >= w / 2) & (Y <= m - w / 2)
    X, Y, T = X[ok], Y[ok], T[ok]
    v = E.mass(X, Y, T)
    idx = np.argsort(v)[:top]
    res = []
    for i in idx:
        x = [X[i], Y[i], T[i]]; fx = v[i]
        st = [(cxr[1] - cxr[0]) / n, (cyr[1] - cyr[0]) / n, math.radians(thr[1] - thr[0]) / n]
        for _ in range(80):
            imp = False
            for k in range(3):
                for sg in (1, -1):
                    y = list(x); y[k] += sg * st[k]
                    ww = math.cos(y[2]) + math.sin(y[2])
                    if not (ww / 2 <= y[0] <= m - ww / 2 and ww / 2 <= y[1] <= m - ww / 2): continue
                    fy = E.mass([y[0]], [y[1]], [y[2]])[0]
                    if fy < fx: x, fx, imp = y, fy, True
            if not imp: st = [s_ / 2 for s_ in st]
        u = F(math.tan(x[2] / 2)).limit_denominator(10 ** 9)
        cx = F(x[0]).limit_denominator(10 ** 9); cy = F(x[1]).limit_denominator(10 ** 9)
        mm = F(r['sn'], r['sd'])
        if admissible(mm, cx, cy, u):
            ex = oracle_mass(r, cx, cy, u)
            res.append((ex, cx, cy, u, fx))
    res.sort()
    for ex, cx, cy, u, fx in res[:top]:
        print(f"  exact {float(ex):.9f} (float {fx:.9f}) at cx={cx} cy={cy} u={u}  (theta {math.degrees(2*math.atan(float(u))):.4f} deg)")
    return res


def locate_in_dump(dump, cx, cy, u):
    hits = []
    for box, kind, Lc in parse_dump(dump):
        if box[0] <= cx <= box[1] and box[2] <= cy <= box[3] and box[4] <= u <= box[5]: hits.append((kind, box))
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('what'); ap.add_argument('path'); ap.add_argument('extra', nargs='?')
    ap.add_argument('--m', default='5'); ap.add_argument('--pitch', default='1/20'); ap.add_argument('--ubins', type=int, default=16)
    ap.add_argument('--n', type=int, default=200); ap.add_argument('--per', type=int, default=30)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--max', type=int, default=10 ** 9)
    ap.add_argument('--kinds', default=None)
    a = ap.parse_args()
    if a.what == 'file': file_check(a.path)
    elif a.what == 'roots': sys.exit(0 if roots_check(a.path, a.m, a.pitch, a.ubins) else 1)
    elif a.what == 'oracle': oracle_test(a.path, a.n, a.seed)
    elif a.what == 'leaves':
        t = leaves_test(a.path, a.extra, a.per, a.seed, a.max, a.kinds.split(',') if a.kinds else None)
        sys.exit(1 if t.viol else 0)
    elif a.what == 'hole':
        cx0, cx1, cy0, cy1, t0, t1 = [float(v) for v in a.extra.split(',')]
        find_hole(a.path, (cx0, cx1), (cy0, cy1), (t0, t1), n=a.n)
    elif a.what == 'locate':
        cx, cy, u = [F(v) for v in a.extra.split(',')]
        for kind, box in locate_in_dump(a.path, cx, cy, u): print('  ', kind, [str(v) for v in box])
    elif a.what == 'tprime':
        t = tprime_test(a.n, a.per, a.seed); sys.exit(1 if t.viol else 0)
    elif a.what == 'boxes':
        t = boxes_test(a.path, a.n, a.per, a.seed); sys.exit(1 if t.viol else 0)


if __name__ == '__main__':
    main()
