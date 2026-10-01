#!/usr/bin/env python3
"""List zmcheck --d4 roots (side 8: 40x40 cells x 4 u-bins, order x-cell, y-cell, bin) not yet in a
ZM_ROOTLOG stderr file, below the highest finished index: the roots in flight (or skipped).
    python3 search/s61_wand125/zmcheck_inflight.py search/s61_wand125/zmcheck_d4/roots.err"""
import re, sys
done = set()
for l in open(sys.argv[1]):
    m = re.match(r'ROOT x\[(\d+)/1000,\d+/1000\] y\[(\d+)/1000,\d+/1000\] u\[(\d+)/8,\d+/8\]', l)
    if m:
        done.add((int(m[1]) // 100, int(m[2]) // 100, int(m[3])))
order = [(i, j, k) for i in range(40) for j in range(40) for k in range(4)]
hi = max(order.index(d) for d in done)
print(f'{len(done)} roots logged; highest index {hi}')
for idx, r in enumerate(order[:hi + 1]):
    if r not in done:
        i, j, k = r
        print(f'  open: index {idx}  x[{i/10:.1f},{(i+1)/10:.1f}] y[{j/10:.1f},{(j+1)/10:.1f}] u[{k}/8,{k+1}/8]')
