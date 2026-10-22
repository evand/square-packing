#!/usr/bin/env python3
"""Grid DP of the chain relaxation vs the conjectured frontier, many H, n <= NMAX."""
import sys
from concurrent.futures import ProcessPoolExecutor
from dp import dp
from targets import target
NMAX = 30
def one(H):
    out, _ = dp(H, NMAX, nt=91, ny=41)
    return H, out
if __name__ == '__main__':
    Hs = [1.5, 1.7, 1.82, 1.83, 1.835, 1.84, 1.845, 1.85, 1.86, 1.87, 1.88, 1.89, 1.9, 1.92, 1.94, 1.96, 1.98, 1.995]
    with ProcessPoolExecutor(6) as ex:
        for H, out in ex.map(one, Hs):
            gaps = {n: out[n] - target(n, H) for n in out}
            print(f"H={H:.3f} min gap {min(gaps.values()):+.5f} (n={min(gaps, key=gaps.get)}) max gap {max(gaps.values()):+.5f} | "
                  + " ".join(f"{n}:{gaps[n]:+.4f}" for n in (3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30)), flush=True)
