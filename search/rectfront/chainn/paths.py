#!/usr/bin/env python3
"""Optimal DP chains: tilt (deg) and floor/ceiling/mid of each middle square."""
import sys, math
import numpy as np
from dp import dp
from targets import target
def path(H, n, nt=91, ny=41):
    out, (T, Y, Sg, St, En, back) = dp(H, n, nt, ny)
    V = back[n - 3]                     # values with n-2 middle squares
    j = int(np.argmin(V + En)); seq = [j]
    for r in range(n - 4, -1, -1):
        j = int(np.argmin(back[r] + Sg[:, j])); seq.append(j)
    seq.reverse()
    def lab(i):
        k = (math.cos(T[i]) + math.sin(T[i])) / 2
        return 'F' if abs(Y[i] - k) < 1e-9 else ('C' if abs(Y[i] - (H - k)) < 1e-9 else f'm{(Y[i]-k)/(H-2*k):.2f}')
    return out[n], [f"{math.degrees(T[i]):.0f}{lab(i)}" for i in seq]
if __name__ == "__main__":
  for H, n in [(1.85, 9), (1.85, 10), (1.9, 7), (1.9, 8), (1.96, 11), (1.96, 12), (1.835, 30), (1.84, 25), (1.99, 20)]:
      v, p = path(H, n)
      print(f"H={H} n={n} DP={v:.5f} target={target(n,H):.5f}  ", ' '.join(p))
