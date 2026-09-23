#!/usr/bin/env python3
"""t3_exist_pigeon: the centre-pigeonhole mechanism -- the far field of (E3) at T = 3 contains no
packing at all, and why the same argument is empty at T = 4.

Task tasks/t3-existence/README.md.  Nothing in search/ is modified.

MECHANISM.  In a packing of unit squares in [0,T]^2 with margin delta >= 0:
  * the centre of square i lies in [u_i/2 + delta, T - u_i/2 - delta]^2, a square of side
    T - u_i - 2 delta, with u_i = |cos th_i| + |sin th_i| in [1, sqrt2]; every one of those boxes
    sits inside the largest of them, of side T - u_min - 2 delta;
  * every separated pair has a unit normal n with n.(c_j - c_i) >= m_ij + delta >= 1 + delta,
    so |c_j - c_i|_2 >= 1 + delta.
So the n centres are n points at pairwise Euclidean distance >= 1 + delta inside a square of side
S = T - min_i u_i - 2 delta.  If d_n is the largest possible minimum pairwise distance of n points
in the UNIT square, this forces  S d_n >= 1 + delta.  d_6 = sqrt(13)/6 (Graham 1963, proved), so
at T = 3, n = 6 and delta = 0:

    3 - u_min  >=  6/sqrt(13) = 1.6641006   <=>   u_min <= 1.3359   <=>  SOME tilt <= 25.8431 deg.

Hence: if EVERY square has tilt in [25.8431 deg, 45 deg] there is no packing of six unit squares
in [0,3]^2 -- the entire regime in which Lemma H fails (tilt >= 28.7959 deg, t3-chain.md sec 2.3)
is empty of packings.  At T = 4, n = 12 the same computation needs 4 - sqrt2 = 2.5858 >= 1/d_12 =
2.5725: TRUE for every tilt, so the argument never fires.

usage:  python3 search/t3_exist_pigeon.py
"""
import math

# largest minimum pairwise distance of n points in the unit square (Erich Friedman, "Points in
# Squares"; n <= 9 are proved optimal, the rest are the best known PACKINGS, i.e. lower bounds on
# d_n -- which is the safe direction for the "mechanism is empty at T = 4" claim).
D = {2: math.sqrt(2), 3: math.sqrt(6) - math.sqrt(2), 4: 1.0, 5: math.sqrt(2) / 2,
     6: math.sqrt(13) / 6, 7: 2 * (2 - math.sqrt(3)), 8: (math.sqrt(6) - math.sqrt(2)) / 2,
     9: 0.5, 10: 0.421279543983903, 11: 0.398207310236844, 12: 0.388730126323020,
     20: 0.286611652351681, 30: 0.226882699666606, 110: 0.102734, 132: 0.093755}
PROVED = {2, 3, 4, 5, 6, 7, 8, 9, 14, 16, 25, 36}


def tilt_threshold(T, n):
    """the largest u for which n points at distance 1 still fit: u <= T - 1/d_n."""
    if n not in D:
        return None
    umax = T - 1.0 / D[n]
    if umax >= math.sqrt(2):
        return math.inf            # never fires
    if umax <= 1.0:
        return 0.0                 # fires everywhere
    # u = sqrt2 sin(t + 45 deg)  ->  t = asin(u/sqrt2) - 45 deg
    return math.degrees(math.asin(umax / math.sqrt(2))) - 45.0


def main():
    print(__doc__.split('usage:')[0])
    print(f"{'T':>4} {'n':>5} {'d_n':>12} {'1/d_n':>10} {'T-sqrt2':>10} {"u_min":>9} "
          f"{'tilt >=':>10}  verdict")
    for (T, n) in ((2, 2), (3, 5), (3, 6), (4, 12), (4, 15), (11, 110), (12, 132)):
        if n not in D:
            continue
        th = tilt_threshold(T, n)
        umax = T - 1.0 / D[n]
        ok = '(proved d_n)' if n in PROVED else '(d_n is a lower bound: safe for "never fires")'
        if th == math.inf:
            v = f"never fires  {ok}"
        elif n not in PROVED:
            v = f"would fire above {th:.4f} deg IF d_n were optimal -- NOT proved"
        elif th == 0.0:
            v = f"fires at every angle  {ok}"
        else:
            v = f"fires for all tilts >= {th:.4f} deg  {ok}"
        print(f"{T:>4} {n:>5} {D[n]:12.9f} {1.0/D[n]:10.6f} {T-math.sqrt(2):10.6f} "
              f"{umax:9.5f} {'' if th in (0.0, math.inf) else f'{th:9.4f}':>10}  {v}")

    print("\nthe T = 3, n = 6 threshold in full:")
    print(f"  6/sqrt(13)           = {6/math.sqrt(13):.10f}")
    print(f"  u_min = 3 - 6/sqrt13 = {3 - 6/math.sqrt(13):.10f}")
    t = math.degrees(math.asin((3 - 6 / math.sqrt(13)) / math.sqrt(2))) - 45.0
    print(f"  tilt threshold       = {t:.10f} deg")
    print(f"  Lemma H fails above  = 28.7959070 deg   (t3-chain.md sec 2.3)")
    print(f"  overlap              = [{t:.4f}, 28.7959] deg  -> the two mechanisms cover the")
    print("                          whole uniform-tilt line with room to spare.")
    print("\nmargin of the pigeonhole at uniform tilt t (positive = no packing):")
    print(f"{'t deg':>8} {'u':>9} {'S = 3-u':>9} {'S d_6':>9} {'1 - S d_6':>11}")
    for td in (0, 10, 20, 25, 25.8431, 28.7959, 30, 35, 40, 45):
        t = math.radians(td)
        u = math.cos(t) + math.sin(t)
        S = 3 - u
        print(f"{td:8.3f} {u:9.6f} {S:9.6f} {S*D[6]:9.6f} {1-S*D[6]:+11.6f}")
    print("\nthe same quantity at T = 4, n = 12 (never positive):")
    print(f"{'t deg':>8} {'u':>9} {'S = 4-u':>9} {'S d_12':>9} {'1 - S d_12':>11}")
    for td in (0, 20, 30, 45):
        t = math.radians(td)
        u = math.cos(t) + math.sin(t)
        S = 4 - u
        print(f"{td:8.3f} {u:9.6f} {S:9.6f} {S*D[12]:9.6f} {1-S*D[12]:+11.6f}")
    print("\nfilters:")
    print("  s(5) = 2.7071 < 3 : five squares, S d_5 = (3-u)/sqrt2 >= (3-sqrt2)/sqrt2 = 1.1213")
    print("                      >= 1 at every tilt, so the mechanism never fires at n = 5 -- it")
    print("                      cannot prove the false statement s(5) = 3.  PASSED.")
    print("  T = 2, n = 2      : S d_2 = (2-u) sqrt2 >= 1 iff u <= 2 - 1/sqrt2 = 1.2929, i.e.")
    print("                      tilt <= 21.09 deg: the mechanism proves that two unit squares at")
    print("                      a common tilt >= 21.09 deg do not fit in [0,2]^2 (true), and is")
    print("                      silent at tilt 0 (also true: they do fit).  PASSED.")


if __name__ == '__main__':
    main()
