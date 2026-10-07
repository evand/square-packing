#!/usr/bin/env python3
"""Aggregate a zm_mixed roots.jsonl census onto the 1/10 root cells that gen_zmmtree.py uses.

Usage (from the repo root):
  python3 search/golf/zm_census.py certificates/s21/zm_mixed_d4/roots.jsonl > search/golf/data/zm_cells_current.json
"""
import json
import sys
from collections import defaultdict
from fractions import Fraction as F

KEYS = ('boxes', 'ADM', 'CHAIN', 'SPLIT', 'PIECE', 'EMPTY', 'cpu', 'THR', 'LIN', 'UNCERT')


def main():
    agg = defaultdict(lambda: defaultdict(float))
    tot = defaultdict(float)
    with open(sys.argv[1]) as f:
        next(f)
        for line in f:
            r = json.loads(line)
            x0, y0 = F(r['root'][0]), F(r['root'][2])
            c = (int(x0 * 10), int(y0 * 10))
            st = r['st']
            for k in KEYS:
                agg[c][k] += st.get(k, 0)
                tot[k] += st.get(k, 0)
    out = {f"{i},{j}": dict(v) for (i, j), v in sorted(agg.items())}
    json.dump(dict(total=dict(tot), cells=out), sys.stdout, indent=0)
    print(json.dumps(dict(tot)), file=sys.stderr)


if __name__ == '__main__':
    main()
