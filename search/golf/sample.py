#!/usr/bin/env python3
"""Stratified sample of the 625 root cells (pitch 1/10, D4 region [0, 5/2]^2) for the Lean-tree census.

Strata are by the shipped zm_mixed census of the current cover (boxes per 1/10 cell, search/golf/zm_census.py):
the heavy cells are taken whole, the light ones sampled uniformly at random (seed fixed).  Writes
search/golf/data/sample_cells.txt (one "i,j" per line) and search/golf/data/strata.json.
Usage (from the repo root): python3 search/golf/sample.py
"""
import json
import random

STRATA = [  # (name, lo, hi, n sampled or None = all)
    ('S6', 11000, 10 ** 9, None),
    ('S5', 3000, 11000, None),
    ('S4', 1000, 3000, 10),
    ('S3', 300, 1000, 10),
    ('S2', 100, 300, 8),
    ('S1', 0, 100, 20),
]


def main():
    d = json.load(open('search/golf/data/zm_cells_current.json'))['cells']
    rng = random.Random(20260930)
    out, strata = [], {}
    for name, lo, hi, n in STRATA:
        cells = sorted(k for k, v in d.items() if lo <= v['boxes'] < hi)
        pick = cells if n is None else sorted(rng.sample(cells, n))
        strata[name] = dict(lo=lo, hi=hi, N=len(cells), cells=cells, sample=pick,
                            zm_boxes=sum(d[c]['boxes'] for c in cells))
        out += pick
    with open('search/golf/data/sample_cells.txt', 'w') as f:
        f.write('\n'.join(out) + '\n')
    json.dump(strata, open('search/golf/data/strata.json', 'w'), indent=1)
    for k, v in strata.items():
        print(k, v['N'], len(v['sample']), v['zm_boxes'])
    # S6 (the 8 heavy cells, 512 zm_mixed roots of pitch 1/20 x 1/32 in u) is sampled at root level instead
    from fractions import Fraction as F
    heavy = {tuple(map(int, c.split(','))) for c in strata['S6']['cells']}
    roots = []
    with open('certificates/s21/zm_mixed_d4/roots.jsonl') as f:
        next(f)
        for l in f:
            r = json.loads(l)
            x, y = F(r['root'][0]), F(r['root'][2])
            if (int(x * 10), int(y * 10)) in heavy:
                roots.append((r['st']['boxes'], r['root']))
    RST = [('R1', 0, 100, 6), ('R2', 100, 1000, 5), ('R3', 1000, 3000, 4), ('R4', 3000, 6000, 2),
           ('R5', 6000, 10 ** 9, 1)]
    rstrata, lines = {}, []
    for name, lo, hi, n in RST:
        rr = sorted((tuple(r), b) for b, r in roots if lo <= b < hi)
        pick = sorted(rng.sample(rr, n))
        rstrata[name] = dict(lo=lo, hi=hi, N=len(rr), zm_boxes=sum(b for _, b in rr),
                             sample=[dict(label=f"{name}_{i}", root=list(r), zm_boxes=b)
                                     for i, (r, b) in enumerate(pick)])
        for i, (r, b) in enumerate(pick):
            lines.append(f"{name}_{i} " + ' '.join(r))
        print(name, len(rr), n, rstrata[name]['zm_boxes'], [b for _, b in pick])
    json.dump(rstrata, open('search/golf/data/strata_roots.json', 'w'), indent=1)
    with open('search/golf/data/sample_roots.txt', 'w') as f:
        f.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
