#!/usr/bin/env python3
"""mu_k for k = 6..14 from the family file (own builder): total == k^2 - 4D, D4 invariance, and mu_k == mu_7 pattern:
also check mu_k agrees with the quadrant measure mu_Q on the closed window [0, k-R]^2 (segments clipped to it)."""
from fractions import Fraction as F
from collections import Counter
B = '/home/evand/math/square-packing/public/s12/search/qx2_data/'
R = 2; w = 2; a = F(9, 5)
prof = []; mod = []
for line in open(B + 'L4_k02_family.txt'):
    if line.startswith('#') or not line.strip(): continue
    t = line.split(); (prof if t[0] in 'hv' else mod).append((t[0],) + tuple(F(x) for x in t[1:]))
D = F(423621306389, 500000000000)
def norm(A, Bp, m): return (min(A, Bp), max(A, Bp), m)
def box(k):
    segs = []
    refl = [lambda x, y: (x, y), lambda x, y: (k - x, y), lambda x, y: (x, k - y), lambda x, y: (k - x, k - y)]
    for (o, c, t0, t1, m) in mod:
        A, Bp = ((t0, c), (t1, c)) if o == 'H' else ((c, t0), (c, t1))
        for f in refl: segs.append(norm(f(*A), f(*Bp), m))
    rot = [lambda x, y: (x, y), lambda x, y: (k - y, x), lambda x, y: (k - x, k - y), lambda x, y: (y, k - x)]
    for j in range(-1, k + 1):
        for p in prof:
            A, Bp = ((j + p[2], p[1]), (j + p[3], p[1])) if p[0] == 'h' else ((j + p[1], p[2]), (j + p[1], p[3]))
            if min(A[0], Bp[0]) >= R and max(A[0], Bp[0]) <= k - R:
                for f in rot: segs.append(norm(f(*A), f(*Bp), p[-1]))
    return segs
def quad(N):
    """mu_Q singular part: nu + band pi on [R, N] x [0,w] + diagonal image; (Lebesgue on U handled separately)"""
    segs = []
    for (o, c, t0, t1, m) in mod:
        segs.append(norm((t0, c), (t1, c), m) if o == 'H' else norm((c, t0), (c, t1), m))
    for j in range(R, N + 1):
        for p in prof:
            A, Bp = ((j + p[2], p[1]), (j + p[3], p[1])) if p[0] == 'h' else ((j + p[1], p[2]), (j + p[1], p[3]))
            if min(A[0], Bp[0]) >= R:
                segs.append(norm(A, Bp, p[-1])); segs.append(norm(A[::-1], Bp[::-1], p[-1]))
    return segs
def clip(segs, X):
    """restrict to closed window [0,X]^2, positive-length parts only (masses rescaled)"""
    out = Counter()
    for (P, Q, m) in segs:
        if P[1] == Q[1]:
            if not (0 <= P[1] <= X): continue
            lo, hi = max(P[0], 0), min(Q[0], X)
            if hi > lo: out[((lo, P[1]), (hi, P[1]))] += m * (hi - lo) / (Q[0] - P[0])
        else:
            if not (0 <= P[0] <= X): continue
            lo, hi = max(P[1], 0), min(Q[1], X)
            if hi > lo: out[((P[0], lo), (P[0], hi))] += m * (hi - lo) / (Q[1] - P[1])
    return out
for k in range(6, 15):
    S = box(k)
    tot = sum(s[2] for s in S) + (k - 2 * a) ** 2
    maps = [lambda x, y: (k - y, x), lambda x, y: (k - x, y), lambda x, y: (y, x)]
    C0 = Counter(S); d4 = all(Counter(norm(f(*P), f(*Q), m) for P, Q, m in S) == C0 for f in maps)
    agree = clip(S, k - R) == clip(quad(k + 2), k - R)
    print(k, 'segments', len(S), 'total == k^2-4D:', tot == k * k - 4 * D, 'D4:', d4, 'mu_k == mu_Q on [0,k-R]^2 (singular part):', agree)
