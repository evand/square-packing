#!/usr/bin/env python3
"""Refresh the live jlevy register (data only, CC BY 4.0, Joshua Levy, the squares project) from an extracted tarball.

  regsync.py PACKING_DIR [--date YYYY-MM-DD]

PACKING_DIR = <extracted squares-main>/packing.  Copies frontier/n-*.md and witnesses/known-best/n-*.yaml into regnow/
(the old regnow is kept as regnow_<old date>), rewrites inputs_live/n-<n>.txt (our text format: `n s` then `x y deg`)
and writes ../../packer/runs/register_live_<date>.json; prints every n whose upper bound changed vs the previous table.
"""
import argparse, datetime, glob, json, math, os, re, shutil, sys
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, '../../packer/runs')


def table_row(path):
    fm = yaml.safe_load(open(path).read().split('---')[1])['packing']
    rep = fm.get('reported_upper_bound') or {}
    lb = (fm.get('verified_lower_bound') or fm.get('reported_lower_bound') or {}).get('value')
    return dict(n=fm['n'], status=fm.get('status'), ub=str(rep.get('value')), by=rep.get('found_by') or [],
                year=rep.get('found_year'), opt=rep.get('analytically_optimized'), improved=rep.get('improved_by') or [],
                method=rep.get('construction_method'), reviewed=fm.get('source_reviewed'), lb=None if lb is None else str(lb))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('packing_dir'); ap.add_argument('--date', default=datetime.date.today().isoformat())
    a = ap.parse_args()
    src_f, src_w = f'{a.packing_dir}/frontier', f'{a.packing_dir}/witnesses/known-best'
    assert os.path.isdir(src_f) and os.path.isdir(src_w), 'not a packing dir'
    reg = f'{HERE}/regnow'
    old = None
    if os.path.exists(f'{reg}/FETCHED.txt'):
        m = re.search(r'\d{4}-\d{2}-\d{2}', open(f'{reg}/FETCHED.txt').read())
        old = m.group(0) if m else 'old'
        if old != a.date and not os.path.exists(f'{reg}_{old}'):
            shutil.copytree(reg, f'{reg}_{old}')
    for sub, src, pat in (('frontier', src_f, 'n-*.md'), ('known-best', src_w, 'n-*.yaml')):
        os.makedirs(f'{reg}/{sub}', exist_ok=True)
        for p in glob.glob(f'{src}/{pat}'):
            shutil.copy(p, f'{reg}/{sub}/')
    open(f'{reg}/FETCHED.txt', 'w').write(
        f'fetched {a.date} codeload.github.com/jlevy/squares main tarball (frontier n-*.md, witnesses/known-best n-*.yaml)\n')

    sys.path.insert(0, HERE)
    import register                                                       # reads regnow first
    T = {}
    for p in sorted(glob.glob(f'{reg}/frontier/n-*.md')):
        r = table_row(p); T[str(r['n'])] = r
    os.makedirs(f'{HERE}/inputs_live', exist_ok=True)
    bad = []
    for n in sorted(map(int, T)):
        try:
            s, sq, s_str = register.witness(n)
        except Exception as e:                                            # trivial grids may have no witness file
            bad.append((n, repr(e)[:60])); continue
        with open(f'{HERE}/inputs_live/n-{n}.txt', 'w') as f:
            if '/' in s_str:                                              # rational exact side -> 40-digit decimal (ceiling)
                from fractions import Fraction as F
                q = F(s_str); s_str = f'{-(-q.numerator * 10**40 // q.denominator) // 1}'; s_str = s_str[:-40] + '.' + s_str[-40:]
            f.write(f'{n} {s_str}\n')
            for x, y, th in sq:
                f.write(f'{x!r} {y!r} {math.degrees(th)!r}\n')
    out = f'{RUNS}/register_live_{a.date}.json'
    json.dump(T, open(out, 'w'), indent=0)
    print(f'{len(T)} cases -> {out}; witnesses failed: {bad[:10]}{" ..." if len(bad) > 10 else ""} ({len(bad)})')
    prev = sorted(glob.glob(f'{RUNS}/register_live_*.json'))
    prev = [p for p in prev if p != out]
    if prev:
        P = json.load(open(prev[-1]))
        ch = [(int(n), P[n]['ub'], r['ub'], r['by']) for n, r in T.items() if n in P and P[n]['ub'] != r['ub']]
        print(f'changed vs {os.path.basename(prev[-1])}: {len(ch)}')
        for n, u0, u1, by in sorted(ch):
            print(f'  {n:4d} {u0} -> {u1} {by}')


if __name__ == '__main__':
    main()
