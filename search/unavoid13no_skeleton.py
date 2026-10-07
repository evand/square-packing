#!/usr/bin/env python3
"""unavoid13-no: shape statistics of the k-sets a loop produced (runs/unavoid13no_<tag>/P_round*.npy).
Clusters the x- and the y-coordinates with gap > 0.3 and reports the column / row occupancy pattern,
the number of points within 0.12 of an interior lattice point, and the points within 0.12 of the
lines x = 1, 2, 3 or y = 1, 2, 3."""
import sys, glob, os
import numpy as np


def clusters(v, gap=0.3):
    v = np.sort(v); out = [[v[0]]]
    for a in v[1:]:
        if a - out[-1][-1] > gap: out.append([a])
        else: out[-1].append(a)
    return out


def main():
    for tag in sys.argv[1:]:
        files = sorted(glob.glob(f"runs/unavoid13no_{tag}/P_round*.npy"))
        print(f"== {tag}: {len(files)} sets")
        for f in files:
            P = np.load(f)
            cx = clusters(P[:, 0]); cy = clusters(P[:, 1])
            colpat = "+".join(str(len(c)) for c in cx); rowpat = "+".join(str(len(c)) for c in cy)
            colpos = ",".join(f"{np.mean(c):.2f}" for c in cx); rowpos = ",".join(f"{np.mean(c):.2f}" for c in cy)
            lat = sum(1 for x, y in P if min(abs(x - i) for i in (1, 2, 3)) < 0.12 and min(abs(y - j) for j in (1, 2, 3)) < 0.12)
            onx = sum(1 for x, y in P if min(abs(x - i) for i in (1, 2, 3)) < 0.12)
            ony = sum(1 for x, y in P if min(abs(y - j) for j in (1, 2, 3)) < 0.12)
            print(f"  {os.path.basename(f)[:-4]}: cols {colpat} at x={colpos}; rows {rowpat} at y={rowpos}; near-lattice {lat}, near x∈{{1,2,3}} {onx}, near y∈{{1,2,3}} {ony}")


if __name__ == '__main__':
    main()
