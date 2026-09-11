#!/usr/bin/env python3
"""Parse + analyse a few named SVGs straight to www/data (no validation) for development."""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_svg import parse, fmt
from analysis import analyze
import export
HERE = os.path.dirname(os.path.abspath(__file__)); D = HERE + '/../data/'; W = HERE + '/../www/data/'
os.makedirs(W + 'p', exist_ok=True)
P = json.load(open(D + 'packings.json')) if os.path.exists(D + 'packings.json') else {}
for name in sys.argv[1:]:
    t = time.time(); s, sq = parse(D + 'ellsworth/svg/' + name)
    squares = [[fmt(x) for x in q] for q in sq]
    a = analyze(float(s), squares)
    P[name] = {'s': fmt(s), 'n': len(sq), 'squares': squares, 'errors': []}
    json.dump({'name': name, 's': fmt(s), 'n': len(sq), 'squares': [[float(x) for x in q] for q in squares], 'analysis': export.compact(a)},
              open(W + 'p/' + name.replace('.svg', '.json'), 'w'), separators=(',', ':'))
    print(name, 'n=%d free=%s wedged=%s mobile=%d rigid=%s sym=%s %.1fs' % (len(sq), a['free'], a['wedged'], sum(a['mobile']), a['rigid'], a['symmetry']['class'], time.time() - t))
export.build_index(P)
