#!/usr/bin/env python3
"""One-off: recompute contact 'exact' flags in exported JSONs (rule: exact unless gap_mp >= 1e-20)."""
import json, glob
for f in glob.glob('www/data/p/*.json'):
    d = json.load(open(f)); ch = False
    for c in d['analysis']['contacts']:
        e = (c['gap_mp'] < 1e-20) if 'gap_mp' in c else True
        if e != c.get('exact'): c['exact'] = e; ch = True
    if ch: json.dump(d, open(f, 'w'), separators=(',', ':'))
