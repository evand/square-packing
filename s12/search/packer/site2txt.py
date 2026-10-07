#!/usr/bin/env python3
"""Convert a site packing JSON (squares = [[cx, cy, deg], ...], 's') to packer's text format."""
import json, sys
d = json.load(open(sys.argv[1]))
sq = d['squares']
s = float(sys.argv[2]) if len(sys.argv) > 2 else float(d['s'])
print(len(sq), repr(s))
for x, y, a in sq:
    print(repr(float(x)), repr(float(y)), repr(float(a)))
