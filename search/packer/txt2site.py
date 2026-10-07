#!/usr/bin/env python3
"""Convert packer text output to site-style JSON on stdout."""
import json, sys
L = [l.split() for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
n, s = int(L[0][0]), float(L[0][1])
print(json.dumps({'n': n, 's': repr(s), 'squares': [[float(a), float(b), float(c)] for a, b, c, *_ in L[1:n + 1]]}))
