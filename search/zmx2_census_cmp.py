#!/usr/bin/env python3
"""zmx2_census_cmp.py A.log B.log -- compare two zmx2 cert logs root for root (ZMX2.md sec 12).

The census of a root is (root box, boxes, cert, empty, uncert, maxdepth, capped), keyed by (id, pass);
the 'ms' field is ignored.  Exit 0 iff the two logs have the same roots with identical censuses.
"""
import sys


def census(p):  # ROOT id pass P root R boxes b cert c empty e uncert u maxdepth m capped k ms t
    out = {}
    for l in open(p):
        if l.startswith('ROOT '):
            f = l.split()
            out[(f[1], f[3])] = tuple(f[i] for i in (5, 7, 9, 11, 13, 15, 17))
    return out


a, b = census(sys.argv[1]), census(sys.argv[2])
diff = [k for k in set(a) | set(b) if a.get(k) != b.get(k)]
print('%s vs %s: %d / %d roots, %d differ' % (sys.argv[1], sys.argv[2], len(a), len(b), len(diff)))
sys.exit(1 if diff or not a else 0)
