#!/usr/bin/env python3
"""Numeric checks for notes/bandcut-cost.md.

Five independent checks, all cheap:

  A  W(D) = |cos D| + |sin D|:  the mixed-tilt link penalty is first order in D.
  B  The mixed-tilt / mixed-normal chain bound MT, and that it reduces to the
     common-tilt bound B_T of search/ARCH_TLEDGER.md sec 2.1.
  C  MT evaluated on the measured T = 4 band optimum
     (runs/bandcut_scan_pT_T4.json): the axis chain and the tilted stack.
  D  The stack slack s(p,t) = p - (p cos t + sin t): s(p,t_p) = 0 and
     ds/dt = 1 at t = t_p for every p  (exact, Fractions).
  E  The waste table of the published squeezable rectangles against the fit
     waste = min(a,b) + 2, the budget identity waste(A)+waste(B) = T, and the
     12x12 decomposition into two (4,8) blocks.

Usage:  python3 search/bandcut_cost.py [A|B|C|D|E|all]
"""
import json
import math
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def W(d):
    return abs(math.cos(d)) + abs(math.sin(d))


# ---------------------------------------------------------------- A

def check_A():
    print("A. link half-width m(D) = 1/2 + W(D)/2, W(D) = |cos D| + |sin D|")
    print(f"{'D(deg)':>8} {'W(D)':>12} {'W-1':>12} {'(W-1)/2':>12} {'D - D^2/2':>12}")
    for ddeg in (0, 0.5, 1, 2, 5, 10, 20, 28.0725, 36.8699, 45):
        d = math.radians(ddeg)
        print(f"{ddeg:8.4f} {W(d):12.8f} {W(d)-1:12.8f} {(W(d)-1)/2:12.8f} "
              f"{d - d*d/2:12.8f}")
    # W >= 1 with equality only at multiples of 90 deg
    worst = min(W(math.radians(k * 90.0 / 2000)) for k in range(2001))
    print(f"   min over 0..90 deg of W = {worst:.12f}   (theory 1)")
    print("   => every tilt-mismatched link costs (W-1)/2 ~ |D|/2 in the numerator,")
    print("      i.e. FIRST order in the tilt difference.")


# ---------------------------------------------------------------- B

TAU_MAX = (1 + math.sqrt(2)) / 2      # max transverse offset of a max-surplus link


def MT(T, k, thetas, phi, R, u1, uk, beta=0.0):
    """Mixed-tilt wall-to-wall chain bound (MT).

    thetas: the k tilts (rad).  phi: the COMMON separation-normal angle (rad),
    cos phi >= 0.  R: total rise (c_k - c_1)_y.  u1, uk: cos+sin of the two
    wall squares.  beta: half-spread of the link normals about phi (0 = common
    normal, the exact case).  Returns the upper bound on delta.

        delta * (cos(beta)(k-1) + 2 cos(phi))
            <=  cos(phi)(T - ubar) + sin(phi) R
                - cos(beta) * M  +  sin(beta) * (k-1) * TAU_MAX
        M = sum_i [ 1/2 + W(theta_{i+1}-theta_i)/2 ].
    """
    M = sum(0.5 + W(thetas[i + 1] - thetas[i]) / 2 for i in range(k - 1))
    ubar = (u1 + uk) / 2
    num = (math.cos(phi) * (T - ubar) + math.sin(phi) * R
           - math.cos(beta) * M + math.sin(beta) * (k - 1) * TAU_MAX)
    den = math.cos(beta) * (k - 1) + 2 * math.cos(phi)
    return num / den


def B_T(T, R, t):
    u = math.cos(t) + math.sin(t)
    return ((T - u) * math.cos(t) + R * math.sin(t)) / (T - 1) - 1


def check_B():
    print("B. MT vs the common-tilt bound B_T of ARCH_TLEDGER sec 2.1")
    print("   (common tilt t, common normal phi = t, k = T, rise R)")
    print(f"{'T':>3} {'t(deg)':>8} {'R':>8} {'B_T':>14} {'MT':>14} "
          f"{'num equal?':>11} {'sign same?':>11}")
    ok = True
    for T in (3, 4, 5, 12):
        for tdeg in (0.8, 3.2, 12.0, 28.0725):
            t = math.radians(tdeg)
            for R in (0.0, 1.0, 1.0 - t, (T - 1) * math.sin(t), 2.0):
                u = math.cos(t) + math.sin(t)
                b = B_T(T, R, t)
                m = MT(T, T, [t] * T, t, R, u, u)
                # same numerator, denominators (T-1) vs (T-1+2cos t)
                numb = b * (T - 1)
                numm = m * ((T - 1) + 2 * math.cos(t))
                eq = abs(numb - numm) < 1e-12
                sg = (b > 0) == (m > 0)
                ok &= eq and sg
            print(f"{T:3d} {tdeg:8.4f} {R:8.4f} {b:14.9f} {m:14.9f} "
                  f"{str(eq):>11} {str(sg):>11}")
    print(f"   identical numerators everywhere: {ok}")
    print("   MT denominator (k-1+2 cos phi) > (T-1): MT is the sharper form")
    print("   (it charges delta at the two walls as well as at the k-1 links).")

    # the straight-stack specialisation: MT > 0  <=>  k cos t + sin t < T
    print("\n   straight stack of k at tilt t, R = (k-1) sin t, walls at tilt t:")
    print(f"{'T':>3} {'k':>3} {'t(deg)':>9} {'k cos+sin':>11} {'MT':>14} {'sign ok':>8}")
    for (T, k) in ((4, 4), (4, 3), (12, 4), (12, 12)):
        for tdeg in (0.0, 10.0, 28.0725, 30.0, 36.8699, 45.0):
            t = math.radians(tdeg)
            u = math.cos(t) + math.sin(t)
            R = (k - 1) * math.sin(t)
            m = MT(T, k, [t] * k, t, R, u, u)
            span = k * math.cos(t) + math.sin(t)
            good = (m > 1e-12) == (span < T - 1e-12)
            print(f"{T:3d} {k:3d} {tdeg:9.4f} {span:11.6f} {m:14.9f} {str(good):>8}")
    print("   MT > 0 exactly when the straight stack's span k cos t + sin t < T.")

    # how fast a spread of link NORMALS destroys the bound
    print("\n   effect of a spread beta of the link normals (T = k, tilt 0, R = 0):")
    print(f"{'T':>3} {'beta(deg)':>10} {'MT':>14}")
    for T in (4, 12):
        for bdeg in (0.0, 0.1, 0.5, 1.0, 5.0):
            m = MT(T, T, [0.0] * T, 0.0, 0.0, 1.0, 1.0, beta=math.radians(bdeg))
            print(f"{T:3d} {bdeg:10.4f} {m:14.9f}")


# ---------------------------------------------------------------- C

def check_C():
    path = os.path.join(ROOT, "runs", "bandcut_scan_pT_T4.json")
    d = json.load(open(path))
    T = d["T"]
    sq = d["sq"]
    print(f"C. MT on the measured T = {T} band optimum, delta = {d['delta']}")
    print(f"   p = {d['p']}, tilt = {d['tdeg']:.6f} deg, {len(sq)} squares")

    def u(th):
        return abs(math.cos(th)) + abs(math.sin(th))

    # chain 1: the four axis-parallel squares of the bottom row
    bot = sorted([s for s in sq if s[2] == 0.0 and abs(s[1] - 0.5) < 1e-9])
    th = [s[2] for s in bot]
    R = bot[-1][1] - bot[0][1]
    m = MT(T, len(bot), th, 0.0, R, u(th[0]), u(th[-1]))
    print(f"   axis chain of {len(bot)}: rise R = {R:+.6f}, MT bound on delta = {m:+.9e}")

    # chain 2: the four tilted squares of the (4,1) stack
    band = sorted([s for s in sq if s[2] != 0.0])
    th = [s[2] for s in band]
    R = band[-1][1] - band[0][1]
    t = th[0]
    m = MT(T, len(band), th, t, R, u(th[0]), u(th[-1]))
    span = len(band) * math.cos(t) + math.sin(t)
    print(f"   band chain of {len(band)}: rise R = {R:+.6f}, MT bound on delta = {m:+.9e}")
    print(f"      (k cos t + sin t = {span:.9f} vs T = {T})")
    print("   => the band chain's own bound is slack; the AXIS chain is what")
    print("      certifies delta <= 0 at T = 4.")

    # what a mixed-tilt chain would cost: replace the band by alternating tilts
    print("\n   cost of mixing tilts on the same geometry (same R, same normals):")
    print(f"{'spread(deg)':>12} {'sum(W-1)/2':>12} {'MT bound':>14}")
    for spread in (0.0, 1.0, 5.0, 10.0, 28.0725):
        th = [ (0.0 if i % 2 == 0 else math.radians(spread)) for i in range(4) ]
        t = math.radians(spread / 2)
        mm = MT(T, 4, th, t, 1.0, u(th[0]), u(th[-1]))
        pen = sum((W(th[i+1] - th[i]) - 1) / 2 for i in range(3))
        print(f"{spread:12.4f} {pen:12.8f} {mm:14.9f}")


# ---------------------------------------------------------------- D

def check_D():
    print("D. stack slack s(p,t) = p - (p cos t + sin t), exact at t = t_p = 2 atan(1/p)")
    print(f"{'p':>3} {'cos t_p':>12} {'sin t_p':>12} {'s(p,t_p)':>10} "
          f"{'ds/dt at t_p':>14} {'height C+pS':>12}")
    ok = True
    for p in range(2, 13):
        # tan(t/2) = 1/p  =>  cos = (p^2-1)/(p^2+1), sin = 2p/(p^2+1)   exact
        C = Fraction(p * p - 1, p * p + 1)
        S = Fraction(2 * p, p * p + 1)
        assert C * C + S * S == 1
        s = p - (p * C + S)
        dsdt = p * S - C          # d/dt [p - p cos t - sin t]
        ok &= (s == 0) and (dsdt == 1)
        print(f"{p:3d} {str(C):>12} {str(S):>12} {str(s):>10} {str(dsdt):>14} "
              f"{float(C + p * S):12.6f}")
    print(f"   s(p,t_p) = 0 and ds/dt = 1 for every p tested: {ok}")
    print("   => a band chain of ANY length gains slack eps (rad) per unit of excess")
    print("      tilt, to first order.  The gain does not depend on p or on T.")


# ---------------------------------------------------------------- E

PUB = [   # (a, b, squares, waste)  -- [pub] Table 2 / Lemma 1 of the paper
    (4, 8, 26, 6),
    (5, 9, 38, 7),
    (5, 10, 43, 7),
    (6, 11, 58, 8),
    (6, 12, 64, 8),
]


def check_E():
    print("E. waste of the published squeezable rectangles vs the fit min(a,b) + 2")
    print(f"{'a':>3} {'b':>3} {'squares':>8} {'area':>5} {'waste':>6} "
          f"{'min+2':>6} {'fit?':>5}")
    ok = True
    for a, b, m, w in PUB:
        assert a * b - m == w
        f = min(a, b) + 2
        ok &= (w == f)
        print(f"{a:3d} {b:3d} {m:8d} {a*b:5d} {w:6d} {f:6d} {str(w == f):>5}")
    for k in range(3, 8):
        a, b, m = 2 * k, 2 * k + 4, 4 * k * k + 6 * k - 2
        w = a * b - m
        f = min(a, b) + 2
        ok &= (w == f == 2 * k + 2)
        print(f"{a:3d} {b:3d} {m:8d} {a*b:5d} {w:6d} {f:6d} {str(w == f):>5}"
              "   (Lemma 1)")
    print(f"   waste = min(a,b) + 2 on every published ingredient: {ok}")

    print("\n   the paper's recipes, read through the fit:")
    for parity, (A, B) in (
            ("even n>=14", (lambda n: (12, 6), lambda n: (n - 10, n - 6))),
            ("odd  n>=13", (lambda n: (10, 5), lambda n: (n - 9, n - 5)))):
        for n in (13, 14, 20, 30):
            a1, b1 = A(n)
            a2, b2 = B(n)
            w1, w2 = min(a1, b1) + 2, min(a2, b2) + 2
            print(f"   {parity}  n={n:3d}  A={(a1,b1)} w={w1:3d}  "
                  f"B={(a2,b2)} w={w2:3d}  sum={w1+w2:3d}  T={n:3d}  "
                  f"{'OK' if w1 + w2 == n else 'MISMATCH'}")

    print("\n   budget identity waste(A)+waste(B) = T with waste = min+2 and min >= 4:")
    print("      T = (mA + 2) + (mB + 2) = mA + mB + 4,  mA, mB >= 4  =>  T >= 12")
    print("   the 12x12 witness (two (4,8) blocks + two integer blocks):")
    tot = 26 + 26 + 4 * 4 + 8 * 8
    print(f"      A=(4,8):26  B=(8,4):26  C=(4,4):16  D=(8,8):64  "
          f"total={tot}  n=T^2-T={12*12-12}  {'OK' if tot == 132 else 'BAD'}")
    print("   at T = 11: waste(A)+waste(B) = 11 with both >= 6 is infeasible.")
    print("   BUT T = 11 is a verified positive (Cantrell) -- see sec 6 of the note.")

    print("\n   filter s(5) = 3:  n = T^2 - 4 = 5 at T = 3.")
    print("      waste of the 3x3 container at n=5 is 4;  an absolute waste bound")
    print("      W0 >= 5 would 'prove' s(5) = 3, which is false: s(5) = 2+1/sqrt2")
    s5 = 2 + 1 / math.sqrt(2)
    print(f"      = {s5:.7f}, whose own waste is s^2 - 5 = {s5*s5-5:.7f}.")
    print(f"      => no absolute-constant waste bound with W0 > {s5*s5-5:.4f} is true.")


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    fns = {"A": check_A, "B": check_B, "C": check_C, "D": check_D, "E": check_E}
    keys = list(fns) if which == "all" else [which.upper()]
    for k in keys:
        fns[k]()
        print()


if __name__ == "__main__":
    main()
