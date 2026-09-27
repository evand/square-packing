#!/usr/bin/env python3
"""Re-summarise the shipped s(21) run records from scratch (certificates/s21/README.md).

  python3 search/s21_records.py cover COVER            own parser: well-formed, all pieces in [0,5]^2, w >= 0, no
                                                      degenerate segment, no polygon; exact total < 21; the measure
                                                      is invariant under x -> 5-x and x <-> y (so under D4).
  python3 search/s21_records.py zm_mixed DIR COVER     DIR = certificates/s21/zm_mixed_d4 (roots.jsonl, manifest.json,
                                                      checker/): header shas = the shipped checker files and the
                                                      cover; settings = D4, cert mode, whole region; the roots are
                                                      exactly the 40,000 boxes of [0,5/2]^2 x u in [0,1/2] (pitch
                                                      1/20, 16 u-bins), each once, none with an uncertified box.
  python3 search/s21_records.py zmx2 LOG COVER MODE    LOG = certificates/s21/zmx2_{d4,full}/roots.log, MODE = d4|full:
                                                      the header's mode/region; every root of the region (2,500
                                                      resp. 2 passes x 10,000) once, 0 uncertified, 0 capped.

Independent of the checkers: it enumerates the region itself and reads only the records.  Exit 0 iff clean.
"""
import sys, os, json, hashlib, re
from fractions import Fraction as F


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def fail(msg):
    print("NOT CLEAN: " + msg)
    sys.exit(1)


def cover_check(path):
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
    npg = nx()
    if i != len(v): fail("trailing tokens")
    if npg != 0: fail("polygons present (the s(21) certificate has none)")
    if F(sn, sd) != 5: fail(f"container side {sn}/{sd} != 5")
    S = 5 * D
    ok = all(0 <= X <= S and 0 <= Y <= S and w >= 0 for X, Y, w in pts)
    ok &= all(0 <= c <= S for s_ in segs for c in s_[:4]) and all(s_[4] >= 0 for s_ in segs)
    if not ok: fail("a piece outside [0,5]^2 or a negative mass")
    if any((X0, Y0) == (X1, Y1) for X0, Y0, X1, Y1, _ in segs): fail("degenerate segment")
    total = F(sum(w for *_, w in pts) + sum(s_[4] for s_ in segs), W)
    from collections import Counter
    def canon(g):
        P = Counter(); Sg = Counter()
        for X, Y, w in pts:
            if w: P[g(X, Y)] += w
        for X0, Y0, X1, Y1, w in segs:
            if w: Sg[(frozenset((g(X0, Y0), g(X1, Y1))), w)] += 1
        return P, Sg
    ident = canon(lambda X, Y: (X, Y))
    inv_r = canon(lambda X, Y: (S - X, Y)) == ident
    inv_s = canon(lambda X, Y: (Y, X)) == ident
    axis = sum(1 for X0, Y0, X1, Y1, _ in segs if X0 == X1 or Y0 == Y1)
    print(f"cover: s = 5, D = {D}, W = {W}: {len(pts)} points, {len(segs)} segments ({axis} axis-parallel), 0 polygons; "
          f"all in [0,5]^2, masses >= 0")
    print(f"total = {total} = {float(total):.12f}  (< 21: {total < 21})")
    print(f"measure invariant under x -> 5-x: {inv_r}, under x <-> y: {inv_s} (these generate D4)")
    if not (total < 21 and inv_r and inv_s): fail("cover")
    print("COVER CLEAN")


def zm_mixed(d, cover):
    recs = os.path.join(d, 'roots.jsonl')
    with open(recs) as fh:
        lines = [json.loads(l) for l in fh if l.strip()]
    h = lines[0]
    if h.get('kind') != 'zm_mixed cert header': fail("first line of roots.jsonl is not a header")
    want = {'zm_mixed.py': sha(os.path.join(d, 'checker', 'zm_mixed.py')),
            'mixed_cover.py': sha(os.path.join(d, 'checker', 'mixed_cover.py')),
            'zeromargin.py': sha(os.path.join(d, 'checker', 'zeromargin.py')),
            'input': sha(cover)}
    for k, v in want.items():
        if h['sha256'].get(k) != v: fail(f"header sha256 of {k} {h['sha256'].get(k)} != {v}")
    s = h['settings']
    if s['mode'] != 'D4' or not s['cert_mode'] or s['tprime'] or any(v is not None for v in s['region'].values()):
        fail(f"settings are not a whole-region D4 cert-mode run: {s}")
    if F(s['pitch']) != F(1, 20) or s['ubins'] != 16: fail(f"pitch/ubins {s['pitch']}/{s['ubins']}")
    p, nb = F(1, 20), 16
    region = {(str(i * p), str((i + 1) * p), str(j * p), str((j + 1) * p), str(F(k, 2 * nb)), str(F(k + 1, 2 * nb)))
              for i in range(50) for j in range(50) for k in range(nb)}
    seen = set(); tot = {}; bad = 0
    for r in lines[1:]:
        if r['label'] != '': fail(f"record with label {r['label']!r} in a D4 run")
        key = tuple(str(F(v)) for v in r['root'])
        if key in seen: fail(f"root {key} recorded twice")
        seen.add(key)
        if r['st']['UNCERT'] != 0 or r['unc']: bad += 1
        for k, v in r['st'].items():
            tot[k] = max(tot.get(k, 0), v) if k == 'maxdepth' else tot.get(k, 0) + v
    missing, extra = region - seen, seen - region
    print(f"header: zm_mixed.py {h['sha256']['zm_mixed.py'][:16]}…, mixed_cover.py {h['sha256']['mixed_cover.py'][:16]}…, "
          f"zeromargin.py {h['sha256']['zeromargin.py'][:16]}… (= checker/), input {h['sha256']['input'][:16]}… (= cover); "
          f"settings D4, cert mode, depth {s['depth']}, CHAIN {s['chain']} from {s['chain_from']}, theta-bias {s['theta_bias']}, "
          f"SPLIT {s['split']}")
    print(f"region: {len(region)} roots = [0,5/2]^2 x u in [0,1/2]; present {len(seen & region)}, missing {len(missing)}, "
          f"extra {len(extra)}, with uncertified boxes {bad}")
    print(f"totals: boxes {tot['boxes']}, max depth {tot['maxdepth']}, leaves ADM {tot['ADM']} CHAIN {tot['CHAIN']} "
          f"SPLIT {tot['SPLIT']} PIECE {tot['PIECE']} P1 {tot['P1']} MIX {tot['MIX']} EMPTY {tot['EMPTY']} "
          f"UNCERTIFIED {tot['UNCERT']}; T' closures {tot['TPTS']}; CPU {tot['cpu'] / 3600:.2f} h")
    man = json.load(open(os.path.join(d, 'manifest.json')))
    if man['records']['sha256'] != sha(recs): fail("manifest's records sha256 != roots.jsonl")
    if man['header'] != h: fail("manifest header != roots.jsonl header")
    if man['result']['verdict'] != 'VERIFIED-D4': fail(f"manifest verdict {man['result']['verdict']}")
    if missing or extra or bad or tot['UNCERT']: fail("zm_mixed records")
    print("zm_mixed D4 RECORDS CLEAN: every root of the D4 region certified")


def zmx2(log, cover, mode):
    hdr = None; seen = set(); tot = {'boxes': 0, 'cert': 0, 'empty': 0, 'uncert': 0, 'capped': 0, 'ms': 0}
    maxd = 0; nunc = 0
    pat = re.compile(r'^ROOT (\d+) pass (\d) root (\S+) boxes (\d+) cert (\d+) empty (\d+) uncert (\d+) '
                     r'maxdepth (\d+) capped (\d+) ms (\d+)$')
    for ln in open(log):
        ln = ln.rstrip('\n')
        if ln.startswith('# zmx2 cert'):
            if hdr is None: hdr = ln
            elif ln != hdr: fail("two different headers in the log")
            continue
        if ln.startswith('UNCERT'): nunc += 1; continue
        m = pat.match(ln)
        if not m: fail(f"unparsed log line {ln[:80]!r}")
        key = (int(m[2]), tuple(str(F(v)) for v in m[3].split(',')))
        if key in seen: fail(f"root {key} recorded twice")
        seen.add(key)
        for k, g in (('boxes', 4), ('cert', 5), ('empty', 6), ('uncert', 7), ('capped', 9), ('ms', 10)):
            tot[k] += int(m[g])
        maxd = max(maxd, int(m[8]))
    if hdr is None: fail("no header")
    n, passes = (25, (0,)) if mode == 'd4' else (50, (0, 1))
    if f"mode={mode}" not in hdr or f"region=x0-{n-1},y0-{n-1},bins0-3" not in hdr: fail(f"header {hdr}")
    p = F(1, 10)
    region = {(ps, (str(i * p), str((i + 1) * p), str(j * p), str((j + 1) * p), str(F(k, 8)), str(F(k + 1, 8))))
              for ps in passes for i in range(n) for j in range(n) for k in range(4)}
    missing, extra = region - seen, seen - region
    what = "[0,5/2]^2 x u in [0,1/2]" if mode == 'd4' else "[0,5]^2 x u in [0,1/2], cover and its y -> 5-y mirror"
    print(f"header: {hdr[2:]}")
    print(f"region: {len(region)} roots = {what}; present {len(seen & region)}, missing {len(missing)}, extra {len(extra)}")
    print(f"totals: boxes {tot['boxes']}, certified {tot['cert']}, empty {tot['empty']}, uncertified {tot['uncert']} "
          f"(UNCERT lines {nunc}), capped {tot['capped']}, max depth {maxd}, CPU {tot['ms'] / 1000:.0f} s")
    if missing or extra or tot['uncert'] or nunc or tot['capped']: fail(f"zmx2 {mode} records")
    print(f"zmx2 {mode} RECORDS CLEAN: every root certified")


if __name__ == '__main__':
    if len(sys.argv) >= 3 and sys.argv[1] == 'cover':
        cover_check(sys.argv[2])
    elif len(sys.argv) >= 4 and sys.argv[1] == 'zm_mixed':
        zm_mixed(sys.argv[2], sys.argv[3])
    elif len(sys.argv) >= 5 and sys.argv[1] == 'zmx2':
        zmx2(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print(__doc__); sys.exit(2)
