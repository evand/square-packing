#!/usr/bin/env python3
"""tabulate runs/sos_probe/log.jsonl: lam(eps) per (system, family, orthant, order, options)."""
import json
import sys
from collections import defaultdict

import glob
path = sys.argv[1] if len(sys.argv) > 1 else 'runs/sos_probe/res'
filt = sys.argv[2:]
tab = defaultdict(dict)
meta = {}
for line in (l for f in sorted(glob.glob(path + '/*.json')) for l in open(f)):
    d = json.loads(line)
    key = (d['cfg'], ','.join(map(str, d['sub'])), d.get('classes', ''), d['orthant'], d.get('diag'), d.get('owner'),
           d.get('override', ''), d['order'], d['prod'], d['pdeg'], d.get('csp0', 0), d.get('bound', False), d.get('mult'))
    if any(f not in ' '.join(map(str, key)) for f in filt):
        continue
    e = 'bound' if d.get('bound') else round(3 - d['T'], 6)
    tab[key][e] = (d['lam'], d['status'], d.get('t_solve'))
    meta[key] = (d['n_mono'], d['gram'][:2], d['n_gram'], d['n_scal'])
for key in sorted(tab, key=str):
    row = tab[key]
    cells = []
    for e in sorted(row, key=lambda x: (isinstance(x, str), -x if not isinstance(x, str) else 0)):
        lam, st, ts = row[e]
        flag = '' if st == 'optimal' else '~' if st == 'optimal_inaccurate' else '!'
        cells.append(f"{e}:{lam:.2e}{flag}({ts:.0f}s)" if lam is not None else f"{e}:None!")
    print(' '.join(map(str, key)), '|', meta[key], '|', '  '.join(cells))
