#!/usr/bin/env python3
"""arch_chain.py -- the wall-to-wall chain identity behind c_T  (2026-09-20).

THE IDENTITY (exact, elementary).  Let Q_1, ..., Q_k be unit squares at a COMMON tilt t whose
consecutive pairs are separated along the common edge normal e = (cos t, sin t) with gap >= 1 + d
(this is the pair row of `S6_SKELETON.md` §3.1: for equal tilts the width m = 1).  Summing the k-1
row inequalities telescopes:

        (c_k - c_1) . e  >=  (k-1) (1 + d).

If the chain is *wall-to-wall in x* -- Q_1 touches x = 0 and Q_k touches x = T, so
(c_1)_x = P and (c_k)_x = T - P with P = u/2, u = cos t + sin t -- then (c_k - c_1)_x = T - u and

        d  <=  B_T(R, t) := [ (T - u) cos t + R sin t ] / (k - 1)  -  1 ,     R := (c_k - c_1)_y

with k = T for a chain of T squares.  R is the chain's total RISE and is the only free quantity.

        B_T(R, t) = [ (R - 1) t + (1 - T/2) t^2 + O(t^3) ] / (T - 1)

so writing R = 1 + r1 t + O(t^2),

        B_T = - c(r1, T) t^2 + O(t^3),        c(r1, T) = ( T/2 - 1 - r1 ) / (T - 1).

Commands:
    analyse   read a configuration (x y tilt per line, or `s6local --show` output) and report every
              tight wall-to-wall chain, its rise R, and the bound B_T it certifies;
    table     the c(r1, T) table and the conjectured c_T;
    fit       given measured delta*/t^2 at several t, extrapolate to t = 0 and back out r1.

    python3 search/arch_chain.py table
    python3 search/arch_chain.py analyse --T 4 --tdeg 0.8 --file runs/arch_tledger_T4_show.txt
    python3 search/arch_chain.py fit --T 5 --tdeg 0.1,0.2,0.4,0.8 --vals ...
"""
import argparse
import math
import re
import sys

import numpy as np
import sympy as sp


# --------------------------------------------------------------------------- symbolic
def series_bound(T, R_expr, order=5):
    t = sp.symbols('t', positive=True)
    u = sp.cos(t) + sp.sin(t)
    B = ((sp.nsimplify(T) - u) * sp.cos(t) + R_expr * sp.sin(t)) / (sp.nsimplify(T) - 1) - 1
    return t, sp.simplify(B), sp.series(B, t, 0, order).removeO().expand()


def cmd_table(a):
    t = sp.symbols('t', positive=True)
    Tv, Rv, r1 = sp.symbols('T R r1')
    u = sp.cos(t) + sp.sin(t)
    B = ((Tv - u) * sp.cos(t) + Rv * sp.sin(t)) / (Tv - 1) - 1
    ser = sp.series(B.subs(Rv, 1 + r1 * t), t, 0, 4).removeO().expand()
    print("# B_T(R,t) = [ (T-u) cos t + R sin t ] / (T-1) - 1,   u = cos t + sin t")
    print("# with R = 1 + r1 t:")
    print("#   B =", sp.simplify(sp.collect(ser, t)))
    c = -sp.simplify(ser.coeff(t, 2))
    print("#   => c(r1,T) = -[t^2] B =", sp.simplify(c), "=", sp.factor(sp.simplify(c)))
    print()
    print("  T    c(r1=0,T) = (T-2)/(2(T-1))   c(r1=-1,T) = T/(2(T-1))    measured c_T")
    meas = {2: '1  (exact, -(u-1)^2)', 3: '3/4  (exact on leaf)', 4: '1/3  (exact on leaf)',
            5: '3/8  (exact on leaf, ARCH_TLEDGER.md sec 3)',
            6: '2/5  (search extrapolates to 0.3999999)'}
    for T in (2, 3, 4, 5, 6, 7, 10, 17, 100):
        c0 = sp.Rational(T - 2, 2 * (T - 1))
        c1 = sp.Rational(T, 2 * (T - 1))
        print(f" {T:3d}      {str(c0):>8s} = {float(c0):.6f}        {str(c1):>8s} = {float(c1):.6f}"
              f"      {meas.get(T, '')}")
    print()
    print("# limits:  c(0,T) -> 1/2  and  c(-1,T) -> 1/2  as T -> infinity;  both are > 0 for")
    print("# every T >= 3, and c(0,T) is strictly increasing in T from 1/3 at T = 4.")
    print("# c(0,T) = 0 only at T = 2;  c(0,T) < 0 never for T >= 2.  NO SIGN CHANGE.")
    print()
    print("# exact closed forms of B_T(R,t) on the two observed branches:")
    tt = sp.symbols('t', positive=True)
    for T in (2, 3, 4, 5, 6, 17):
        for r1, nm in ((0, 'r1= 0'), (-1, 'r1=-1')):
            _t, B, ser = series_bound(T, 1 + r1 * tt)
            print(f"  T={T:3d} {nm}:  B = {sp.nsimplify(sp.simplify(ser))}")


# --------------------------------------------------------------------------- configuration
def read_config(path):
    """accept `sqNN  x=..  y=..  tilt=..` (s6local --show / arch_ct --show) or plain `x y tilt`."""
    pts = []
    for line in open(path):
        m = re.search(r'x=\s*([-\d.eE+]+)\s+y=\s*([-\d.eE+]+)\s+tilt=\s*([-\d.eE+]+)', line)
        if m:
            pts.append((float(m.group(1)), float(m.group(2)), math.radians(float(m.group(3)))))
            continue
        f = line.split()
        if len(f) == 3 and not line.lstrip().startswith('#'):
            try:
                pts.append((float(f[0]), float(f[1]), math.radians(float(f[2]))))
            except ValueError:
                pass
    return pts


def analyse(pts, T, tol=1e-7, verbose=True):
    n = len(pts)
    X = np.array([p[0] for p in pts]); Y = np.array([p[1] for p in pts])
    th = np.array([p[2] for p in pts])
    C, S = np.cos(th), np.sin(th)
    P = 0.5 * (np.abs(C) + np.abs(S))
    # the geometric margin: d = min over pairs of (best normal gap - m) and over walls
    best = 1e18
    gaps = {}
    for i in range(n):
        for j in range(i + 1, n):
            dx, dy = X[j] - X[i], Y[j] - Y[i]
            m = 0.5 + 0.5 * (abs(math.cos(th[j] - th[i])) + abs(math.sin(th[j] - th[i])))
            g = -1e18
            axis = None
            for o in (i, j):
                for k, (d0, d1) in enumerate(((C[o], S[o]), (-S[o], C[o]))):
                    pr = abs(d0 * dx + d1 * dy)
                    if pr > g:
                        g, axis = pr, (o, k)
            gaps[(i, j)] = (g - m, axis)
            best = min(best, g - m)
    wall = min(list(X - P) + list(T - P - X) + list(Y - P) + list(T - P - Y))
    delta = min(best, 1e18)
    if verbose:
        print(f"# n={n} T={T}  geometric pair margin delta = {delta:+.12e}   wall slack = {wall:+.3e}")
    # tight links, split by which normal family they use (k = 0: the "x-ish" normal e)
    links = {0: [], 1: []}
    for (i, j), (g, (o, k)) in gaps.items():
        if g <= delta + tol:
            links[k].append((i, j, o))
    if verbose:
        print(f"# tight links: {len(links[0])} on normal e=(cos,sin) (x-chains), "
              f"{len(links[1])} on e'=(-sin,cos) (y-chains)")
    out = {}
    for k, axname, coord in ((0, 'x', X), (1, 'y', Y)):
        adj = {i: set() for i in range(n)}
        for (i, j, _o) in links[k]:
            adj[i].add(j); adj[j].add(i)
        # wall-to-wall chains: start at a square touching the low wall, end at the high wall
        lo = [i for i in range(n) if coord[i] - P[i] <= wall + tol]
        hi = set(i for i in range(n) if T - P[i] - coord[i] <= wall + tol)
        chains = []

        def walk(path):
            i = path[-1]
            if i in hi and len(path) >= 2:
                chains.append(list(path))
            if len(path) > T + 2:
                return
            for j in adj[i]:
                if j in path:
                    continue
                if coord[j] > coord[i] + 1e-9:
                    walk(path + [j])
        for i in lo:
            walk([i])
        # the strongest certificate of each length: the one with the SMALLEST signed rise R
        # (the rise enters the bound with a + sign along e and a - sign along e').
        best_len = {}
        for ch in chains:
            kk = len(ch)
            dperp = (Y[ch[-1]] - Y[ch[0]]) if k == 0 else (X[ch[-1]] - X[ch[0]])
            dpar = (X[ch[-1]] - X[ch[0]]) if k == 0 else (Y[ch[-1]] - Y[ch[0]])
            R = dperp if k == 0 else -dperp
            t0 = th[ch[0]]
            B = (dpar * math.cos(t0) + R * math.sin(t0)) / (kk - 1) - 1
            if kk not in best_len or B < best_len[kk][0]:
                best_len[kk] = (B, R, dpar, ch)
        out[axname] = best_len
        if verbose:
            for kk in sorted(best_len):
                B, R, dpar, ch = best_len[kk]
                t0 = th[ch[0]]
                u = math.cos(t0) + math.sin(t0)
                Bwall = ((T - u) * math.cos(t0) + R * math.sin(t0)) / (kk - 1) - 1
                print(f"  {axname}-chain of {kk}: span={dpar:+.9f} (T-u={T - u:+.9f})  "
                      f"rise R={R:+.9f}  B={B:+.6e}  B(wall-to-wall)={Bwall:+.6e}")
                if kk == T:
                    r1 = (R - 1) / t0
                    print(f"      -> R = 1 + r1*t with r1 = {r1:+.6f};  "
                          f"c(r1,T) = {(T / 2 - 1 - r1) / (T - 1):.6f}   "
                          f"[delta/t^2 = {delta / t0**2:+.6f}]")
                    print(f"      -> chain squares (x, y): "
                          + ' '.join(f"({X[i]:.6f},{Y[i]:.6f})" for i in ch))
    return delta, out


def cmd_analyse(a):
    pts = read_config(a.file)
    if a.tdeg is not None:
        pts = [(x, y, math.radians(a.tdeg)) for (x, y, _t) in pts]
    analyse(pts, a.T, tol=a.tol)


def cmd_leaf(a):
    """the exact uniform-tilt leaf bounds of S6_LOCAL.md sec 3 / ARCH_TLEDGER.md sec 3, checked for
    sign over the whole tilt range.  A2c asks for an explicit radius; these have none -- they are
    <= 0 for every t."""
    t = sp.symbols('t', positive=True)
    u = sp.cos(t) + sp.sin(t)
    leaves = {
        2: -(u - 1) ** 2,
        3: -3 * (u - 1) ** 2 / (u ** 2 + 3),
        4: ((32 * sp.sqrt(2) * sp.sin(t) ** 4 * sp.cos(t + sp.pi / 4) + 128 * sp.sin(t) ** 4
             + 52 * sp.sin(t) + 16 * sp.sin(2 * t) - 21 * sp.sin(3 * t) + 8 * sp.sin(4 * t)
             - sp.sin(5 * t) - 48 * sp.sqrt(2) * sp.sin(t + sp.pi / 4) - 48 * sp.cos(t) ** 5
             + 14 * sp.cos(t) + 112 * sp.cos(2 * t) - 14 * sp.cos(3 * t) - 16)
            / (4 * (3 * sp.sin(t) + sp.sin(3 * t) + 11 * sp.cos(t) + 7 * sp.cos(3 * t)))),
        5: ((-2 * (1 - sp.cos(2 * t)) ** 4 - 4 * (1 - sp.cos(2 * t)) ** 3 * sp.sin(2 * t)
             - 28 * (1 - sp.cos(2 * t)) ** 3 - 10 * (1 - sp.cos(2 * t)) ** 2 * sp.sin(t)
             + 36 * (1 - sp.cos(2 * t)) ** 2 * sp.sin(2 * t)
             - 10 * sp.sqrt(2) * (1 - sp.cos(2 * t)) ** 2 * sp.sin(3 * t + sp.pi / 4)
             - 30 * (1 - sp.cos(2 * t)) ** 2 * sp.cos(t) - 120 * (1 - sp.cos(2 * t)) ** 2
             + 120 * sp.sin(t) - 207 * sp.sin(2 * t) + 200 * sp.sin(3 * t) - 160 * sp.sin(4 * t)
             + 80 * sp.sin(5 * t) - 11 * sp.sin(6 * t) + 440 * sp.cos(t) - 634 * sp.cos(2 * t)
             + 180 * sp.cos(3 * t) + 20 * sp.cos(5 * t) + 10 * sp.cos(6 * t) - 16)
            / (2 * (-(1 - sp.cos(2 * t)) ** 2 + 6 * sp.sin(2 * t) + 10 * sp.cos(2 * t) + 6) ** 2)),
    }
    print("# exact fixed-assignment bounds on the uniform-tilt leaf (S6_LOCAL sec 3; T=5 new here).")
    print("# delta*(t,...,t) <= leaf_T(t) holds only where that leaf's dual stays FEASIBLE; the")
    print("# closed form continues past that point and is then meaningless.  The interval of")
    print("# interest is (0, 45] deg, since delta*(t) = delta*(90-t) (reflect in y = x).")
    for T, ex in sorted(leaves.items()):
        f = sp.lambdify(t, ex, 'math')
        ser = sp.series(ex, t, 0, 5).removeO().expand()
        vals, worst45, cross = [], -1e18, None
        for i in range(1, 9000):
            d = i / 100.0
            try:
                v = f(math.radians(d))
            except Exception:
                continue
            if d <= 45.0 and v > worst45:
                worst45 = v
            if cross is None and v > 0:
                cross = d
            if abs(d - round(d)) < 1e-9 and round(d) in (5, 15, 28, 45):
                vals.append((round(d), v))
        print(f"\n  T = {T}:  series  {sp.nsimplify(ser)}")
        print("    " + "  ".join(f"{d}deg: {v:+.6f}" for d, v in vals))
        print(f"    max over (0, 45] deg: {worst45:+.6f}   "
              f"first sign change of the formula: {cross if cross else 'none below 90'} deg")
    print()
    print("# On (0, 45] every leaf bound is strictly negative, so A2c needs NO radius on these")
    print("# leaves -- this is the 'shape B' form of the local theorem (proof-architecture sec 2').")
    print("# The T = 4, 5 formulas do turn positive above ~50 deg; that is the fixed assignment")
    print("# ceasing to be the optimal one, not delta* turning positive.")


def cmd_threshold(a):
    """B_T(R,t) <= 0  <=>  R <= R*(T,t) = T tan(t/2) + cos t - sin t."""
    t, T, R = sp.symbols('t T R', positive=True)
    u = sp.cos(t) + sp.sin(t)
    B = ((T - u) * sp.cos(t) + R * sp.sin(t)) / (T - 1) - 1
    Rstar = sp.solve(sp.Eq(B, 0), R)[0]
    Rstar = sp.simplify(Rstar)
    print("# B_T(R,t) = [ (T-u) cos t + R sin t ] / (T-1) - 1 ,   u = cos t + sin t")
    print("# B_T(R,t) <= 0  <=>  R <= R*(T,t), with")
    print("#   R*(T,t) =", Rstar)
    cand = T * sp.tan(t / 2) + sp.cos(t) - sp.sin(t)
    print("#   claimed  R* = T tan(t/2) + cos t - sin t ;  difference simplifies to",
          sp.simplify(sp.expand_trig(sp.simplify(Rstar - cand))))
    print("#   series   R* =", sp.series(cand, t, 0, 3))
    print()
    print("# THREE SPECIALISATIONS (all exact):")
    print("#  (a) R = 1  (the tiling / pinwheel: the chain climbs exactly one cell).")
    print("#      B_T(1,t) <= 0  <=>  T tan(t/2) + cos t - sin t >= 1  <=>  T >= 1 + cos t + sin t")
    e = sp.simplify(sp.expand_trig(cand - 1))
    print("#      check: R* - 1 =", e, " -> T >= 1 + u, and max u = sqrt(2), so")
    print("#      B_T(1,t) <= 0 for EVERY t as soon as T >= 1 + sqrt2 = 2.41421..., i.e. T >= 3,")
    print("#      with equality only at t = 0 mod 90 deg.   At T = 2 it FAILS for every t in (0,90).")
    print("#  (b) R = (T-1) sin t  (a STRAIGHT row of T at tilt t: consecutive centres differ by")
    print("#      the edge vector (cos t, sin t)).  Then")
    Bstraight = sp.simplify(B.subs(R, (T - 1) * sp.sin(t)))
    print("#      B_T = ", Bstraight, " = (T-u) cos t/(T-1) - cos^2 t,")
    print("#      and B_T > 0  <=>  T cos t + sin t < T  <=>  tan(t/2) > 1/T   (§3.2's window).")
    print("#  (c) R = 0  (a chain with no rise): B_T = (T-u) cos t/(T-1) - 1 < 0 for all t in (0,90).")
    print()
    print("# numeric check of (a) and (b) on a grid:")
    bad_a = bad_b = 0
    for Tv in (2, 3, 4, 5, 6, 17):
        for i in range(1, 9000):
            tv = math.radians(i / 100.0)
            uu = math.cos(tv) + math.sin(tv)
            Bv = ((Tv - uu) * math.cos(tv) + 1 * math.sin(tv)) / (Tv - 1) - 1
            if (Bv <= 1e-15) != (Tv >= 1 + uu - 1e-15):
                bad_a += 1
            Rs = (Tv - 1) * math.sin(tv)
            Bs = ((Tv - uu) * math.cos(tv) + Rs * math.sin(tv)) / (Tv - 1) - 1
            if (Bs > 1e-15) != (math.tan(tv / 2) > 1.0 / Tv + 1e-15):
                bad_b += 1
    print(f"  (a) mismatches: {bad_a} / {6 * 8999};   (b) mismatches: {bad_b} / {6 * 8999}")
    print()
    print("# how much room the rise has, r1 <= T/2 - 1 (R = 1 + r1 t):")
    print("  T    R*(T,t) - 1 at t = t_T = 2 atan(1/T)   R* slope T/2 - 1   B_T(1, t_T)  "
          " depth -c_T t_T^2")
    for Tv in (3, 4, 5, 6, 10, 17, 100):
        tT = 2 * math.atan(1.0 / Tv)
        uu = math.cos(tT) + math.sin(tT)
        Rs = Tv * math.tan(tT / 2) + math.cos(tT) - math.sin(tT)
        Bv = ((Tv - uu) * math.cos(tT) + math.sin(tT)) / (Tv - 1) - 1
        cT = (Tv - 2) / (2.0 * (Tv - 1))
        print(f" {Tv:3d}       {Rs - 1:+.6f}                    {Tv / 2 - 1:8.2f}     "
              f"{Bv:+.6e}   {-cT * tT**2:+.6e}")


def cmd_fit(a):
    """delta*/t^2 measured at several t -> linear extrapolation to t = 0 -> c_T and r1."""
    ts = [math.radians(float(s)) for s in a.tdeg.split(',')]
    vs = [float(s) for s in a.vals.split(',')]
    assert len(ts) == len(vs)
    deg = min(2, len(ts) - 1)
    co = np.polyfit(ts, vs, deg)
    c0 = co[-1]
    print(f"# T={a.T}:  fit of delta*/t^2 in t (degree {deg}) -> value at t=0: {c0:+.8f}")
    print(f"#   => c_{a.T} <= {-c0:.8f}")
    for q, nm in ((sp.Rational(a.T - 2, 2 * (a.T - 1)), '(T-2)/(2(T-1))'),
                  (sp.Rational(a.T, 2 * (a.T - 1)), 'T/(2(T-1))'),
                  (sp.Rational(3, 4), '3/4'), (sp.Rational(1), '1')):
        print(f"#      {nm:>16s} = {str(q):>6s} = {float(q):.8f}   "
              f"(diff {abs(float(q) + c0):.2e})")
    r1 = a.T / 2 - 1 + c0 * (a.T - 1)
    print(f"#   implied chain rise slope r1 = T/2 - 1 + c0*(T-1) = {r1:+.6f}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('table'); p.set_defaults(fn=cmd_table)
    p = sub.add_parser('analyse')
    p.add_argument('--T', type=float, required=True)
    p.add_argument('--file', required=True)
    p.add_argument('--tdeg', type=float, default=None)
    p.add_argument('--tol', type=float, default=1e-7)
    p.set_defaults(fn=cmd_analyse)
    p = sub.add_parser('leaf'); p.set_defaults(fn=cmd_leaf)
    p = sub.add_parser('threshold'); p.set_defaults(fn=cmd_threshold)
    p = sub.add_parser('fit')
    p.add_argument('--T', type=int, required=True)
    p.add_argument('--tdeg', required=True)
    p.add_argument('--vals', required=True)
    p.set_defaults(fn=cmd_fit)
    a = ap.parse_args()
    a.fn(a)


if __name__ == '__main__':
    main()
