#!/usr/bin/env python3
"""t3_chain_hform: the closed forms of the H-lemma, and the checks that they are exact.

Three objects, all for the disjunctive LP of S6_SKELETON.md sec 3.1 with the walls carrying
delta (BANDCUT_K.md sec 1.1):

  H(T, L, t)   the value of the explicit H-shaped Farkas dual: one wall-to-wall chain of T
               squares at common tilt t, with L transverse links (a below the low end, c above
               the high end, L = a + c) and the six wall rows that close it.  Weights:
                   chain links (T-1)      1
                   lo-ax(s_1), hi-ax(s_T) 1/cos t
                   leg links (L)          tan t
                   lo-tr(b_1), hi-tr(d_c) sin t
                   hi-ax(b_1), lo-ax(d_c) sin t tan t
               Exact cancellation of every centre coordinate is checked in `check`.

  omega(t)     THE UNIFORM-TILT IDENTITY.  If every square has tilt t (so every P_i = u/2,
               u = cos t + sin t, and every m_ij = 1), then EVERY normalised Farkas dual obeys
                   delta  <=  omega * (T + 2 - u)  -  1 ,
               omega the total weight on the `lo` wall rows (= that on the `hi` rows).  So the
               whole question at uniform tilt is how small omega can be made, and the H dual has
                   omega  =  M / ( (T-1) + L tan t + 2 M ) ,   M = (1 + sin^2 t + sin t cos t)/cos t.

  Lstar(T,t)   the least L for which H(T,L,t) <= 0:
                   L*  =  ( M (T - u) - (T-1) ) / tan t ,      L* -> T - 2 as t -> 0 .
"""
import math
import sys


def M_of(t):
    C, S = math.cos(t), math.sin(t)
    return (1.0 + S * S + S * C) / C


def H(T, L, t):
    """the H dual's bound, in closed form."""
    C, S = math.cos(t), math.sin(t)
    u = C + S
    Mv = M_of(t)
    W = (T - 1) + L * S / C + 2 * Mv
    N = -(T - 1) - L * S / C + (T - u) * Mv
    return N / W


def omega_H(T, L, t):
    C, S = math.cos(t), math.sin(t)
    Mv = M_of(t)
    return Mv / ((T - 1) + L * S / C + 2 * Mv)


def Lstar(T, t):
    if t <= 0:
        return float(T - 2)
    C, S = math.cos(t), math.sin(t)
    u = C + S
    return (M_of(t) * (T - u) - (T - 1)) * C / S


def check():
    """the H weights cancel every centre coordinate exactly, and H = omega(T+2-u) - 1."""
    import numpy as np
    bad = 0
    for T in (3, 4, 5, 11):
        for a in (0, 1, 2):
            for c in (0, 1, 2):
                for td in (0.5, 5, 20, 45):
                    t = math.radians(td)
                    C, S = math.cos(t), math.sin(t)
                    u = C + S
                    # squares: chain s_0..s_{T-1}; bottom leg b_0..b_{a-1} (b_0 on the wall,
                    # b_{a-1} -> s_0); top leg d_0..d_{c-1} (s_{T-1} -> d_0, d_{c-1} on the wall)
                    nsq = T + a + c
                    ch = list(range(a, a + T))
                    bl = list(range(a))                       # b_0 .. b_{a-1}
                    tl = list(range(a + T, a + T + c))
                    A = []          # (weight, coeff dict square -> (dx, dy), b)
                    e = (C, S)                                # main-chain normal
                    f = (-S, C)                               # transverse normal
                    P = u / 2.0
                    for s in range(T - 1):
                        A.append((1.0, {ch[s + 1]: e, ch[s]: (-e[0], -e[1])}, 1.0))
                    A.append((1.0 / C, {ch[0]: (1.0, 0.0)}, P))
                    A.append((1.0 / C, {ch[-1]: (-1.0, 0.0)}, P - T))
                    low = bl + [ch[0]]
                    for s in range(a):
                        A.append((S / C, {low[s + 1]: f, low[s]: (-f[0], -f[1])}, 1.0))
                    up = [ch[-1]] + tl
                    for s in range(c):
                        A.append((S / C, {up[s + 1]: f, up[s]: (-f[0], -f[1])}, 1.0))
                    if a:
                        A.append((S, {bl[0]: (0.0, 1.0)}, P))
                        A.append((S * S / C, {bl[0]: (-1.0, 0.0)}, P - T))
                    else:
                        A.append((S, {ch[0]: (0.0, 1.0)}, P))
                        A.append((S * S / C, {ch[0]: (-1.0, 0.0)}, P - T))
                    if c:
                        A.append((S, {tl[-1]: (0.0, -1.0)}, P - T))
                        A.append((S * S / C, {tl[-1]: (1.0, 0.0)}, P))
                    else:
                        A.append((S, {ch[-1]: (0.0, -1.0)}, P - T))
                        A.append((S * S / C, {ch[-1]: (1.0, 0.0)}, P))
                    co = np.zeros((nsq, 2))
                    tw = 0.0
                    wb = 0.0
                    for (w, d, bb) in A:
                        for q, v in d.items():
                            co[q] += w * np.array(v)
                        tw += w
                        wb += w * bb
                    resid = float(np.max(np.abs(co)))
                    val = -wb / tw
                    closed = H(T, a + c, t)
                    om = omega_H(T, a + c, t)
                    om_id = om * (T + 2 - u) - 1.0
                    if resid > 1e-12 or abs(val - closed) > 1e-12 or abs(val - om_id) > 1e-12:
                        bad += 1
                        print(f'  FAIL T={T} a={a} c={c} t={td} resid={resid:.2e} '
                              f'val={val:.12f} closed={closed:.12f} omega-form={om_id:.12f}')
    print('check: %s (dual feasibility, closed form, omega identity over T in {3,4,5,11}, '
          'a,c in {0,1,2}, t in {0.5,5,20,45} deg)' % ('OK' if bad == 0 else '%d FAILURES' % bad))


def tables():
    print('\n# H(T, L, t): the H dual bound.  <= 0 is a kill.   T = 3, n = 6, L <= n - T = 3')
    print('%7s' % 't deg', ''.join('%12s' % ('L=%d' % L) for L in range(5)), '%10s' % 'L*(3,t)')
    for td in (0, 0.5, 1, 2, 5, 10, 15, 20, 25, 28, 29, 30, 35, 40, 45):
        t = math.radians(td)
        print('%7.1f' % td, ''.join('%12.6f' % H(3, L, t) for L in range(5)),
              '%10.4f' % Lstar(3, t))
    print('\n# L*(T, t) = least number of transverse links that makes the H dual kill.')
    print('# L_max = n - T = T(T-2) is all the squares there are.')
    print('%5s %8s' % ('T', 'L_max'), ''.join('%9s' % ('%g deg' % d)
                                              for d in (0, 1, 5, 10, 20, 30, 45)))
    for T in (2, 3, 4, 5, 6, 11, 12, 17):
        print('%5d %8d' % (T, T * (T - 2)),
              ''.join('%9.3f' % Lstar(T, math.radians(d))
                      for d in (0, 1, 5, 10, 20, 30, 45)))
    print('\n# and the same as a FRACTION L*/L_max -- what share of the squares the H must use')
    print('%5s' % 'T', ''.join('%9s' % ('%g deg' % d) for d in (0, 1, 5, 10, 20, 30, 45)))
    for T in (3, 4, 5, 6, 11, 12, 17):
        print('%5d' % T, ''.join('%9.3f' % (Lstar(T, math.radians(d)) / (T * (T - 2)))
                                 for d in (0, 1, 5, 10, 20, 30, 45)))
    print('\n# the uniform-tilt threshold: delta <= 0 iff some dual has omega <= 1/(T+2-u)')
    print('%7s %10s %12s %12s' % ('t deg', 'u', '1/(T+2-u)', 'omega_H(3,3,t)'))
    for td in (0, 5, 15, 30, 45):
        t = math.radians(td)
        u = math.cos(t) + math.sin(t)
        print('%7.1f %10.6f %12.6f %12.6f' % (td, u, 1.0 / (3 + 2 - u), omega_H(3, 3, t)))


if __name__ == '__main__':
    check()
    if 'notables' not in sys.argv:
        tables()
