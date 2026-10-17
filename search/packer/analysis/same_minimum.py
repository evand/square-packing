#!/usr/bin/env python3
"""Are two postings the same local minimum?  (10-17, Evan: "refinement" = the same minimum computed more exactly or
reached by descent; a different minimum is a new packing, "based on" the one it came from.)

For every lineage pair in site/data/lineage.yaml (posting -> {n: the posting it started from}), and every two postings
at one n (n <= 324) whose posted sides are within --close: the best stored packing
of each posting at n (sighting source pending:#N), exactsolve on each (Newton to the KKT point, 80 digits), then
  same      a certified continuous path (analysis/pathcert.py) joins the higher posting's core to the lower's with the
            side never above the higher's.  The core is every square that is a rattler in neither exact solution;
            rattlers carry no stress and enter no equation, so two placements of them are one minimum (10-17, Evan:
            "a packing gets credit for its refinements, even if the first posted improvement is a slightly different
            configuration").  Both ends are full packings; core_only marks pairs certified only with rattlers left out.
            Candidates: their exact sides agree to 1e-15 and at least one certifies (a KKT point with a certified minimum's side
            is that minimum; distinct minima essentially never share a side).  A posting that is not itself at a KKT
            point (f64 floats short of the minimum) is descended first (exact/descend.py) and the descended packing
            solved: a posting that descends to its successor's minimum makes the successor a refinement
  distinct  both certify and the exact sides differ; or the posted sides differ by more than GAP (no float posting of a
            minimum is that far above it), checked without solving
  same-side exact sides agree but no path was certified (treated as unknown)
  unknown   a packing missing from the store, or a solve that does not certify
Writes site/data/lineage_minima.json: {"<posting>|<n>": {verdict, root, S_ev, S_root, gap, why}}; read by
site/tools/provenance.py (a root is credited as the finder only on 'same').

  same_minimum.py [--procs 14] [--gap 1e-5] [--close 1e-7]
"""
import argparse, json, os, subprocess, sys, tempfile
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
PK = os.path.dirname(HERE)
ROOT = os.path.abspath(os.path.join(PK, '..', '..'))
sys.path.insert(0, PK)
EXACT = os.path.join(ROOT, 'search', 'exact', 'exactsolve.py')
DESCEND = os.path.join(ROOT, 'search', 'exact', 'descend.py')
OUT = os.path.join(PK, 'runs', 'same_minimum')
LINEAGE = os.path.join(ROOT, 'site', 'data', 'lineage.yaml')
RESULT = os.path.join(ROOT, 'site', 'data', 'lineage_minima.json')


def posting_packing(st, ev, n):
    """(store id, side) of the best feasible packing a posting carries at n, or None"""
    if not ev.startswith('#'):
        return None
    r = st.db.execute("SELECT p.id, p.s FROM packing p JOIN sighting g ON g.packing = p.id WHERE p.n=? AND g.source=? "
                      "AND p.status != 'infeasible' ORDER BY p.s LIMIT 1", (n, 'pending:' + ev)).fetchone()
    return (r['id'], r['s']) if r else None


def _report(d, stem):
    p = os.path.join(d, stem + '.json')
    if not os.path.exists(p):
        return None
    j = json.load(open(p))
    return dict(S=j.get('S_exact'), ok=bool(j.get('cert_valid') and j.get('newton_converged', True)))


def _floats(src, dst):
    """the store writes exact fractions; descend.py reads floats"""
    from fractions import Fraction
    rows = [l.split() for l in open(src) if l.strip() and not l.startswith('#')]
    out = [f'{rows[0][0]} {float(Fraction(rows[0][1]))!r}'] + [' '.join(repr(float(Fraction(x))) for x in t) for t in rows[1:]]
    open(dst, 'w').write('\n'.join(out) + '\n')


def solve(args):
    """exact side of a posting's packing: exactsolve on it as posted; when that finds no KKT point (a float posting
    short of its minimum: 'not jammed', a failed solve) or one that doesn't certify, exact/descend.py first (MILP kicks + slp2) and exactsolve on
    the descended packing (10-17: #432's 261 and 108 postings)"""
    pid, n = args
    d = os.path.join(OUT, f'n{n}_p{pid}')
    os.makedirs(d, exist_ok=True)
    inp = os.path.join(d, 'in.txt')
    if not os.path.exists(os.path.join(d, 'in.json')):
        subprocess.run([sys.executable, os.path.join(PK, 'pk.py'), 'get', str(pid), '--fmt', 'ours', '-o', inp],
                       cwd=PK, capture_output=True, timeout=300)
        subprocess.run([sys.executable, EXACT, inp, '--out', d, '-q'], capture_output=True, timeout=3600)
    r0 = _report(d, 'in')
    if r0 and r0['S'] and r0['ok']:
        return pid, dict(r0, via='posted')
    dd = os.path.join(d, 'descend')
    os.makedirs(dd, exist_ok=True)
    if not os.path.exists(os.path.join(dd, 'desc.json')):
        _floats(inp, os.path.join(dd, 'in_f.txt'))
        subprocess.run([sys.executable, DESCEND, os.path.join(dd, 'in_f.txt'), os.path.join(dd, 'desc.txt'), '--rounds', '40'],
                       capture_output=True, timeout=3600)
        if os.path.exists(os.path.join(dd, 'desc.txt')):
            subprocess.run([sys.executable, EXACT, os.path.join(dd, 'desc.txt'), '--out', dd, '-q'], capture_output=True, timeout=3600)
    r = _report(dd, 'desc')
    if r and r['S'] and (r['ok'] or not (r0 and r0['S'])):
        return pid, dict(r, via='descended')
    return pid, (dict(r0, via='posted') if r0 and r0['S'] else None)


def _exact_file(pid, n, via):
    d = os.path.join(OUT, f'n{n}_p{pid}')
    return os.path.join(d, 'in.exact.txt') if via == 'posted' else os.path.join(d, 'descend', 'desc.exact.txt')


def path_job(args):
    """analysis/pathcert.py on two postings whose exact sides agree: a certified path from the higher to the lower
    with the side never above the higher's (10-17, Evan: the real test of one basin)"""
    key, n, (pa, va), (pb, vb) = args
    out = os.path.join(OUT, 'paths', key.replace('|', '_').replace('#', '') + '.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if not os.path.exists(out):
        cmd = [sys.executable, '-W', 'ignore', os.path.join(HERE, 'pathcert.py'),
               os.path.join(OUT, f'n{n}_p{pa}', 'in.txt'), '--exact', _exact_file(pa, n, va),
               '--to', os.path.join(OUT, f'n{n}_p{pb}', 'in.txt'), '--exact-to', _exact_file(pb, n, vb), '--out', out]
        subprocess.run(cmd, capture_output=True, timeout=7200)
    if not os.path.exists(out):
        return key, None
    r = json.load(open(out))
    return key, dict(ok=r['ok'], ok_mod_rattlers=r.get('ok_mod_rattlers', r['ok']), L=r.get('L'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--procs', type=int, default=14); ap.add_argument('--gap', type=float, default=1e-5)
    ap.add_argument('--close', type=float, default=1e-7)
    a = ap.parse_args()
    import yaml
    from pk.store import Store
    st = Store()
    pairs = []                                       # (key, n, posting, root, A, B)
    for ev, m in (yaml.safe_load(open(LINEAGE)).get('lineage') or {}).items():
        for n, root in m.items():
            pairs.append((f'{ev}|{n}', int(n), str(ev), str(root), posting_packing(st, str(ev), int(n)),
                          posting_packing(st, str(root), int(n))))
    nl = len(pairs)
    # 10-17: every two postings at one n whose sides are within --close (provenance's SAME_TOL): side distance alone
    # can't tell them apart (146: #505 and #606 are distinct minima 1.4e-8 apart; #668's float postings sit 1.7e-8 above
    # their own minima).  Identical sides (< 1e-12) are copies and need no check.
    import itertools, collections
    by = collections.defaultdict(list)
    for n, src, _ in st.db.execute("SELECT n, source, s FROM frontier WHERE source LIKE 'pending:%' AND n <= 324"):
        P = posting_packing(st, src.split(':', 1)[1], n)
        if P:
            by[n].append((src.split(':', 1)[1], P))
    for n, L in sorted(by.items()):
        for (x, A), (y, B) in itertools.combinations(sorted(L), 2):
            if 1e-12 <= abs(A[1] - B[1]) <= a.close:
                pairs.append((f'{x}|{y}|{n}', n, x, y, A, B))
    todo = sorted({(p[0], n) for k, n, ev, root, A, B in pairs if A and B and abs(A[1] - B[1]) <= a.gap for p in (A, B)})
    print(f'{nl} lineage pairs + {len(pairs) - nl} close pairs; {len(todo)} packings to solve exactly', flush=True)
    with Pool(a.procs) as pool:
        S = dict(pool.imap_unordered(solve, todo))
    out = {}
    from fractions import Fraction
    for k, n, ev, root, A, B in pairs:
        if not A or not B:
            out[k] = dict(verdict='unknown', root=root, why='a posting has no stored packing at n'); continue
        gap = A[1] - B[1]
        if abs(gap) > a.gap:
            out[k] = dict(verdict='distinct', root=root, gap=gap, why=f'posted sides differ by {abs(gap):.2e} > {a.gap:g}')
            continue
        sa, sb = S.get(A[0]), S.get(B[0])
        if not sa or not sb or not (sa['ok'] or sb['ok']):
            out[k] = dict(verdict='unknown', root=root, gap=gap, why='no certified exact side for either posting',
                          S_ev=sa and sa['S'], S_root=sb and sb['S']); continue
        dS = float(Fraction(sa['S'][:60]) - Fraction(sb['S'][:60]))
        same = abs(dS) < 1e-15
        if not same and not (sa['ok'] and sb['ok']):
            out[k] = dict(verdict='unknown', root=root, gap=gap, why=f'exact sides differ by {dS:+.3e} but one does not certify',
                          S_ev=sa['S'][:40], S_root=sb['S'][:40]); continue
        via = ', '.join(f'{w} {x["via"]}' for w, x in ((ev, sa), (root, sb)) if x['via'] != 'posted')
        out[k] = dict(verdict='same' if same else 'distinct', root=root, gap=gap, S_ev=sa['S'][:40], S_root=sb['S'][:40],
                      dS_exact=dS, why='exact sides ' + ('agree' if same else f'differ by {dS:+.3e}') + (f' ({via})' if via else ''))
    # 'same' needs a certified path (pathcert.py) from the higher posting to the lower; exact sides agreeing without
    # one is 'same-side' (provenance treats it as unknown)
    PJ = []
    for k, n, ev, root, A, B in pairs:
        v = out.get(k)
        if v and v['verdict'] == 'same':
            hi, lo = (A, B) if A[1] >= B[1] else (B, A)
            PJ.append((k, n, (hi[0], S[hi[0]]['via']), (lo[0], S[lo[0]]['via'])))
    print(f'{len(PJ)} same verdicts: certifying paths', flush=True)
    with Pool(a.procs) as pool:
        for k, r in pool.imap_unordered(path_job, PJ):
            out[k]['path'] = r
            if r and r['ok']:
                out[k]['why'] += '; path certified at the higher posting\'s side'
            elif r and r.get('ok_mod_rattlers'):
                out[k]['core_only'] = True
                out[k]['why'] += ("; core path certified (rattlers excluded: they may sit in other cages, or no route "
                                  "was found for them)")
            else:
                out[k]['verdict'] = 'same-side'
                out[k]['why'] += '; no certified path'
    json.dump(dict(note='generated by search/packer/analysis/same_minimum.py; keys "posting|n" are lineage links '
                        '(root = the posting it started from), "a|b|n" two postings at n within --close',
                   gap=a.gap, close=a.close, pairs=out), open(RESULT, 'w'), indent=1, sort_keys=True)
    C = collections.Counter(v['verdict'] for v in out.values())
    for k, v in sorted(out.items(), key=lambda kv: (kv[1]['verdict'], kv[0])):
        if v['verdict'] != 'same' or k.count('|') == 1:
            print(f"{k:>16} ~ {v['root']:<10} {v['verdict']:<9} {v['why']}")
    print(dict(C))


if __name__ == '__main__':
    main()
