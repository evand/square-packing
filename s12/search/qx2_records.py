#!/usr/bin/env python3
"""Re-check the shipped records of the s(k^2 - 3) = k bundle (certificates/k2m3/README.md) from scratch.
Independent of qx2_zm.py's control flow: it reads only the files, enumerates the root region itself and re-derives
the admissibility condition; nothing is imported from the checkers.  Exact (Fraction) arithmetic throughout.

  python3 search/qx2_records.py cover BOX N
      own parser of the mixed v1 box cover: well-formed, no points, every segment axis-parallel and non-degenerate,
      every piece in [0,m]^2, masses >= 0; exactly one polygon, the square [a, m-a]^2 with mass = area (Lebesgue);
      exact total < N; the measure is invariant under x -> m-x and x <-> y (so under D4), as multisets of pieces.

  python3 search/qx2_records.py record REC.jsonl.gz REC.out CHECKERDIR BOX
      REC = the qx2_zm.py --dump-leaves record (run V3).  Checks:
      * header sha256 of qx2_zm.py, zm_mixed.py, zeromargin.py, mixed_cover.py = the files in CHECKERDIR, and of the
        input = BOX; the argv has the certificate settings (--depth 18, --exact-umax 1/2, --exact-from 3,
        --dump-leaves, default pitch and u-bins, no region restriction, EXACT on); the .out repeats these shas;
      * the roots are exactly the D4 grid [0,m/2]^2 x u in [0,1/2], pitch 1/10, 8 u-bins, each once;
      * per root: UNCERT 0, no uncertified box; every leaf kind is a certifying kind; the census counts equal the
        leaf list; boxes = 2 leaves - 1 (a binary tree); max depth <= 18; every leaf lies in its root; the labels are
        consistent (AXIS: u1 = 0; SYM: u0^2 + 2u0 - 1 >= 0, i.e. theta0 >= 45 deg; EXACT0: u0 = 0; EXACT45:
        u1^2 + 2u1 - 1 > 0; EXACT: u0 > 0 and theta1 <= 45 deg; EMPTY: the closed leaf holds no admissible pose);
      * COVERAGE, per root: the positive-volume leaves have pairwise disjoint interiors, and every part of the root
        they leave uncovered (the slabs removed by clip_bin) contains no admissible pose with u > 0, exactly; so
        leaves + slabs tile the root.  (u = 0 is Lemma Z's; a pose on the bottom face u = u0 > 0 of an uncovered part
        is on the closed top face of the part or leaf below it, which is checked there.)
        Admissible: the closed unit square at centre c, angle theta = 2 atan u lies in [0,m]^2, i.e.
        w(u)/2 <= cx <= m - w(u)/2 and the same for cy, with w(u) = cos theta + sin theta = (1 + 2u - u^2)/(1 + u^2);
        w increases on [0, sqrt2 - 1] and decreases after, so on a box it is decided at the end points of the bin.
      * the census totals (boxes, max depth, every leaf kind) equal the .out's, and the .out says VERIFIED-D4.
      NOT checked: the mass inequality itself at any PIECE / EXACT / EXACT0 / EXACT45 / LEB / CAP leaf.  Those labels
      are accepted on their intervals and the coverage above; re-proving them needs a re-run of qx2_zm.py
      (certificates/k2m3/verify.sh --full).  So a clean record is a complete, well-formed proof skeleton, not a
      fresh geometric check and not a second checker.
Exit 0 iff clean.
"""
import sys, os, json, gzip, hashlib, re
from fractions import Fraction as F
from collections import Counter


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def fail(msg):
    print("NOT CLEAN: " + msg)
    sys.exit(1)


# ------------------------------------------------------------------------------------------------ the cover
def read_cover(path):
    tok = []
    for line in open(path):
        tok += line.split('#', 1)[0].split()
    if tok[:2] != ['mixed', '1']: fail("not a mixed v1 file")
    v = [int(t) for t in tok[2:]]
    i = 0

    def nx(k=1):
        nonlocal i
        if i + k > len(v): fail("file too short")
        i += k
        return v[i - k:i] if k > 1 else v[i - 1]
    sn, sd, D, W = nx(), nx(), nx(), nx()
    pts = [tuple(nx(3)) for _ in range(nx())]
    segs = [tuple(nx(5)) for _ in range(nx())]
    polys = []
    for _ in range(nx()):
        k, w = nx(), nx()
        polys.append((w, [tuple(nx(2)) for _ in range(k)]))
    if i != len(v): fail("trailing tokens")
    return F(sn, sd), D, W, pts, segs, polys


def cover_check(path, n):
    m, D, W, pts, segs, polys = read_cover(path)
    if m.denominator != 1 or m < 1: fail(f"container side {m} is not a positive integer")
    m = int(m); S = m * D
    if pts: fail("points present (the qx2 box covers have none)")
    for X0, Y0, X1, Y1, w in segs:
        if not all(0 <= c <= S for c in (X0, Y0, X1, Y1)) or w < 0: fail("segment outside the container or w < 0")
        if (X0, Y0) == (X1, Y1): fail("degenerate segment")
        if X0 != X1 and Y0 != Y1: fail("segment not axis-parallel")
    if len(polys) != 1: fail(f"{len(polys)} polygons (expected exactly one, the Lebesgue square)")
    pw, vs = polys[0]
    xs = sorted(set(x for x, _ in vs)); ys = sorted(set(y for _, y in vs))
    if len(vs) != 4 or len(xs) != 2 or xs != ys or set(vs) != {(x, y) for x in xs for y in xs}:
        fail("the polygon is not an axis-parallel square")
    A, B = xs
    if A + B != S or not 0 < A < B: fail("the polygon square is not [a, m-a]^2")
    # counter-clockwise (FORMAT.md): positive shoelace area
    area2 = sum(vs[j][0] * vs[(j + 1) % 4][1] - vs[(j + 1) % 4][0] * vs[j][1] for j in range(4))
    if area2 <= 0: fail("polygon not counter-clockwise")
    area = F(area2, 2 * D * D)
    if F(pw, W) != area: fail(f"polygon mass {F(pw, W)} != its area {area} (not Lebesgue)")
    total = F(sum(s_[4] for s_ in segs) + pw, W)

    def canon(g):
        Sg = Counter()
        for X0, Y0, X1, Y1, w in segs:
            if w: Sg[(frozenset((g(X0, Y0), g(X1, Y1))), w)] += 1
        Pg = Counter((frozenset(g(x, y) for x, y in vs_), w) for w, vs_ in polys)
        return Sg, Pg
    ident = canon(lambda X, Y: (X, Y))
    d4 = canon(lambda X, Y: (S - X, Y)) == ident and canon(lambda X, Y: (Y, X)) == ident
    print(f"cover: s = {m}, D = {D}, W = {W}: 0 points, {len(segs)} axis-parallel segments on "
          f"{len(set(('H', s_[1]) if s_[1] == s_[3] else ('V', s_[0]) for s_ in segs))} lines, 1 polygon = "
          f"[{F(A, D)}, {F(B, D)}]^2 with mass = area = {area}; total = {total} = {float(total):.12f}; "
          f"< {n}: {total < n}; D4-invariant: {d4}")
    if not total < n: fail(f"total {total} >= {n}")
    if not d4: fail("measure not D4-invariant")
    print("cover: CLEAN")


# ------------------------------------------------------------------------------------------------ the record
GOOD = {'LEB', 'CAP', 'EXACT', 'EXACT0', 'EXACT45', 'AXIS', 'SYM', 'PIECE', 'ADM', 'P1', 'MIX', 'SPLIT', 'EMPTY'}
ALLKEYS = GOOD | {'UNCERT', 'CORE', 'CHAIN', 'TRI', 'THR', 'LIN', 'TPTS', 'boxes', 'maxdepth', 'cpu'}
SETTINGS = {'--depth': '18', '--exact-umax': '1/2', '--exact-from': '3'}
FREE = {'--nproc', '--progress', '--resume', '--unc-out'}          # do not change what is certified
FLAGS = {'--dump-leaves'}


def w(u):
    return (1 + 2 * u - u * u) / (1 + u * u)


def admissible_M(b, m):
    """2 min(cx1, m - cx0, cy1, m - cy0): a pose of the box's centre range with angle u is admissible iff w(u) <= M
    (and the centre range meets [w/2, m - w/2] in both coordinates, which is what M expresses)."""
    x0, x1, y0, y1 = b[:4]
    return 2 * min(x1, m - x0, y1, m - y0)


def has_adm_halfopen(b, m):
    """does [x0,x1] x [y0,y1] x (u0, u1] contain an admissible pose?  w is strictly increasing then strictly
    decreasing, so inf over (u0, u1] is min(w(u0), w(u1)), attained at u1 or approached from above at u0."""
    M = admissible_M(b, m)
    return w(b[5]) <= M or w(b[4]) < M


def has_adm_closed(b, m):
    M = admissible_M(b, m)
    return min(w(b[4]), w(b[5])) <= M


def meets(a, b):          # positive-volume intersection of two boxes
    return all(max(a[2 * i], b[2 * i]) < min(a[2 * i + 1], b[2 * i + 1]) for i in range(3))


def contains(a, b):       # a contains b
    return all(a[2 * i] <= b[2 * i] and b[2 * i + 1] <= a[2 * i + 1] for i in range(3))


def vol(b):
    return (b[1] - b[0]) * (b[3] - b[2]) * (b[5] - b[4])


def coverage(root, leaves, m):
    """k-d partition of the root along leaf faces.  Returns (covered volume, uncovered volume, #uncovered parts,
    problems).  A terminal part is either inside exactly one leaf or meets none (then: no admissible pose)."""
    cov_v = F(0); unc_v = F(0); nunc = 0; bad = []
    stack = [(root, [L for L in leaves if meets(L, root)])]
    while stack:
        reg, Ls = stack.pop()
        if not Ls:
            unc_v += vol(reg); nunc += 1
            if has_adm_halfopen(reg, m): bad.append(('admissible pose in an uncovered part', reg))
            continue
        cut = None
        for L in Ls:
            for i in range(3):
                for c in (L[2 * i], L[2 * i + 1]):
                    if reg[2 * i] < c < reg[2 * i + 1]: cut = (i, c); break
                if cut: break
            if cut: break
        if cut is None:             # every leaf here contains the part
            if len(Ls) > 1: bad.append(('overlapping leaves', reg))
            cov_v += vol(reg)
            continue
        i, c = cut
        lo = list(reg); lo[2 * i + 1] = c; lo = tuple(lo)
        hi = list(reg); hi[2 * i] = c; hi = tuple(hi)
        stack.append((lo, [L for L in Ls if meets(L, lo)]))
        stack.append((hi, [L for L in Ls if meets(L, hi)]))
    return cov_v, unc_v, nunc, bad


def record_check(rec, out, chkdir, box):
    m, *_ = read_cover(box)
    m = int(m)
    L = [json.loads(l) for l in gzip.open(rec, 'rt') if l.strip()]
    hdr, recs = L[0], L[1:]
    if hdr.get('kind') != 'header': fail("first line is not a header")
    want = {n: sha(os.path.join(chkdir, n)) for n in ('qx2_zm.py', 'zm_mixed.py', 'zeromargin.py', 'mixed_cover.py')}
    want['input'] = sha(box)
    if hdr['sha256'] != want:
        fail("header sha256 != the checker files / box cover: " + str({k: (hdr['sha256'].get(k), v) for k, v in want.items()}))
    print("header: sha256 of qx2_zm.py " + want['qx2_zm.py'][:8] + "…, zm_mixed.py " + want['zm_mixed.py'][:8] +
          "…, zeromargin.py " + want['zeromargin.py'][:8] + "…, mixed_cover.py " + want['mixed_cover.py'][:8] +
          "…, input " + want['input'][:8] + "… = the shipped files")
    argv = hdr['argv'][2:]
    if not hdr['argv'][0].endswith('qx2_zm.py') or hdr['argv'][1].startswith('-') or hdr['argv'][1] == 'axis':
        fail("argv is not a qx2_zm.py certificate run: " + str(hdr['argv']))
    seen = {}; j = 0
    while j < len(argv):
        a = argv[j]
        if a in FLAGS: seen[a] = True; j += 1
        elif a in SETTINGS or a in FREE: seen[a] = argv[j + 1]; j += 2
        else: fail(f"argv option {a} not allowed in a certificate run (settings: {SETTINGS}, flags {FLAGS})")
    for k, v in SETTINGS.items():
        if seen.get(k) != v: fail(f"argv {k} = {seen.get(k)}, want {v}")
    if '--dump-leaves' not in seen: fail("no --dump-leaves")
    print(f"argv: {' '.join(hdr['argv'][:2])} " + ' '.join(f"{k} {v}" for k, v in SETTINGS.items()) +
          " --dump-leaves (pitch, u-bins default 1/10, 8; whole region; EXACT on)")
    # the roots
    pitch, ub = F(1, 10), 8
    n = int(F(m, 2) / pitch)
    grid = {(i * pitch, (i + 1) * pitch, j * pitch, (j + 1) * pitch, F(k, 2 * ub), F(k + 1, 2 * ub))
            for i in range(n) for j in range(n) for k in range(ub)}
    roots = [tuple(F(v) for v in r['root']) for r in recs]
    if len(roots) != len(set(roots)): fail("a root appears twice")
    if set(roots) != grid: fail("the roots are not the D4 grid")
    print(f"roots: exactly the {len(grid)} boxes of [0,{F(m, 2)}]^2 x u in [0,1/2] (pitch 1/10, 8 u-bins), each once")
    tot = Counter(); maxdepth = 0; nleaves = 0; ndeg = 0; kinds = Counter()
    cv = unv = F(0); nunc = 0; nbad = 0
    for r, root in zip(recs, roots):
        st = r['st']
        if set(st) - ALLKEYS: fail(f"unknown census key {set(st) - ALLKEYS} at root {r['root']}")
        if st.get('UNCERT', 0) != 0 or r['unc']: fail(f"uncertified box at root {r['root']}")
        lv = [(tuple(F(v) for v in b), k) for b, k in r['leaves']]
        c = Counter(k for _, k in lv)
        if set(c) - GOOD: fail(f"non-certifying leaf kind {set(c) - GOOD} at root {r['root']}")
        for k in GOOD | {'UNCERT'}:
            if st.get(k, 0) != c.get(k, 0): fail(f"census {k} = {st.get(k)} != {c.get(k, 0)} leaves at root {r['root']}")
        if st['boxes'] != 2 * len(lv) - 1: fail(f"boxes {st['boxes']} != 2 * {len(lv)} - 1 at root {r['root']}")
        if st['maxdepth'] > 18: fail(f"depth {st['maxdepth']} > 18 at root {r['root']}")
        for b, k in lv:
            if not contains(root, b) or b[0] >= b[1] or b[2] >= b[3] or b[4] > b[5]:
                fail(f"leaf {b} not a box inside its root {r['root']}")
            u0, u1 = b[4], b[5]
            ok = {'AXIS': u1 == 0, 'SYM': u0 * u0 + 2 * u0 - 1 >= 0, 'EXACT0': u0 == 0 and u1 > 0,
                  'EXACT45': u1 * u1 + 2 * u1 - 1 > 0 and u0 < u1, 'EMPTY': not has_adm_closed(b, m),
                  'EXACT': u0 > 0 and u1 * u1 + 2 * u1 - 1 <= 0 and u0 < u1}.get(k, u0 < u1)
            if not ok: fail(f"leaf {[str(x) for x in b]} labelled {k} fails its label condition")
        pos = [b for b, _ in lv if vol(b) > 0]
        ndeg += len(lv) - len(pos)
        a, bnc, nu, bad = coverage(root, pos, m)
        if bad or a + bnc != vol(root) or a != sum(vol(b) for b in pos):
            nbad += 1
            print(f"  root {r['root']}: {bad[:3]}")
            continue
        cv += a; unv += bnc; nunc += nu
        for k in st:
            if k == 'maxdepth': maxdepth = max(maxdepth, st[k])
            elif k != 'cpu': tot[k] += st[k]
        tot['cpu'] += st['cpu']
        nleaves += len(lv); kinds += c
    if nbad: fail(f"{nbad} roots not covered")
    print(f"leaves: {nleaves} ({ndeg} of zero volume: AXIS/EMPTY at a single angle), all certifying kinds, labels "
          f"consistent, census = leaf lists, boxes = 2 leaves - 1 per root, UNCERT 0")
    print(f"coverage: in every root the positive-volume leaves are interior-disjoint and the {nunc} uncovered parts "
          f"(the clip_bin slabs) hold no admissible pose with u > 0; volume: leaves {float(cv):.9f} + slabs "
          f"{float(unv):.9f} = {cv + unv} (exact) = the region's (m/2)^2 / 2: {cv + unv == F(m * m, 4) / 2}")
    if cv + unv != F(m * m, 4) / 2: fail("volumes do not add up to the region")
    # the .out
    txt = open(out).read()
    for k, v in want.items():
        name = 'input' if k == 'input' else k
        if not re.search(rf"^sha256 {v}  {re.escape(name)}$", txt, re.M): fail(f".out does not state sha256 {k} = {v}")
    mo = re.search(r"done in \d+s: boxes (\d+), max depth (\d+), CPU (\d+) s", txt)
    ml = re.search(r"^  leaves: (.*)$", txt, re.M)
    if not mo or not ml: fail(".out has no census")
    outc = {k: int(v) for k, v in re.findall(r"(\w+) (\d+)", ml.group(1))}
    mine = {k: tot.get(k, 0) for k in outc}
    if int(mo.group(1)) != tot['boxes'] or int(mo.group(2)) != maxdepth or mine != outc:
        fail(f".out census {mo.groups()} {outc} != records {tot['boxes']}, {maxdepth}, {mine}")
    if abs(int(mo.group(3)) - tot['cpu']) > 1: fail(".out CPU != sum of the records' cpu")
    if not re.search(r"^VERIFIED-D4 ", txt, re.M) or 'NOT VERIFIED' in txt: fail(".out does not say VERIFIED-D4")
    print(f"census = .out: boxes {tot['boxes']}, max depth {maxdepth}, CPU {tot['cpu']:.0f} s; " +
          ", ".join(f"{k} {tot.get(k, 0)}" for k in outc if tot.get(k, 0)) + "; VERIFIED-D4")
    print("record: CLEAN")


if __name__ == '__main__':
    a = sys.argv[1:]
    if len(a) == 3 and a[0] == 'cover': cover_check(a[1], int(a[2]))
    elif len(a) == 5 and a[0] == 'record': record_check(*a[1:])
    else:
        print(__doc__); sys.exit(2)
