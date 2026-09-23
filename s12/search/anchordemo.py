#!/usr/bin/env python3
"""Build the anchor-clique demonstration certificate (task I).

Takes a plain certificate, DELETES one point and puts the anchor clique `K(p, A)` of Lemma 2 at
that point with the same weight.  `K(p, A)` contains the whole point clique of `p`
(`notes/clique-family.md` Lemma 2), so the covering is at least as strong and the total weight is
unchanged -- and the clique is *load-bearing*: with its weight set to 0 the file no longer covers.
That is what the box-clique demonstration could not be (a box clique can never contain a point
clique), and it is the point of the family.

    python3 search/anchordemo.py CERT OUT --i 0 [--frac 0.95] [--slack 0.1] [--pad 0.002] [--D 1000000]
"""
import sys, os, math, argparse
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import anchorclique as AC


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cert'); ap.add_argument('out')
    ap.add_argument('--i', type=int, default=0, help='index of the point to replace')
    ap.add_argument('--frac', type=float, default=0.95); ap.add_argument('--slack', type=float, default=0.10)
    ap.add_argument('--pad', type=float, default=0.002); ap.add_argument('--D', type=int, default=1000000)
    a = ap.parse_args()
    t = open(a.cert).read().split()
    sn, sd, D, WD, m = (int(v) for v in t[:5])
    pts = [(int(t[5+3*i]), int(t[6+3*i]), int(t[7+3*i])) for i in range(m)]
    assert len(t) == 5 + 3*m, "this tool only handles plain certificates"
    s = F(sn, sd)
    X, Y, w = pts[a.i]
    px, py = X / D, Y / D
    cands = AC.cand_params(F(int(round(float(s) * a.D)), a.D), a.D, int(round(px*a.D)), int(round(py*a.D)),
                           frac=a.frac, slack=a.slack, pad=a.pad)
    assert cands, f"point {a.i} = ({px:.4f}, {py:.4f}) is not within distance 1 of a wall"
    wall, en, rn, d = cands[0]
    Xa, Ya = int(round(px*a.D)), int(round(py*a.D))
    cl = AC.kpa(F(int(round(float(s)*a.D)), a.D), a.D, Xa, Ya, wall, en, rn)
    assert cl is not None, "the anchor would leave the container"
    rest = [p for j, p in enumerate(pts) if j != a.i]
    # Keep the anchor's own points as ZERO-WEIGHT atoms.  They carry no weight, but the sweep's
    # cell boundaries are the atoms' breakpoints, and a cell is credited only if it lies WHOLLY in
    # a piece: without them the cells straddle the boundary of {S : p in S} and of
    # {S : A subseteq S} and a whole band around each loses the credit.
    for (vx, vy) in [(Xa, Ya)] + [(int(F(v) * a.D), int(F(u) * a.D)) for (v, u) in
                                  [(cl[0][1][1], cl[0][1][2]), (cl[0][1][3], cl[0][1][4])]]:
        gx = F(vx, a.D) * D; gy = F(vy, a.D) * D
        if gx.denominator == 1 and gy.denominator == 1:
            rest.append((int(gx), int(gy), 0))
    body = "\n".join(f"{x} {y} {ww}" for x, y, ww in rest)
    blk = AC.block([cl], [w], a.D)
    with open(a.out, 'w') as f:
        f.write(f"{sn} {sd}\n{D}\n{WD}\n{len(rest)}\n{body}\n{blk}")
    tot = sum(p[2] for p in rest) + w
    print(f"{a.out}: {len(rest)} points + 1 anchor clique of weight {w}/{WD} = {w/WD:.7f}")
    print(f"  p = ({px:.6f}, {py:.6f}) (wall {wall}, distance {d:.4f}), eps = {en/a.D:.6f}, rho = {rn/a.D:.6f}, "
          f"rho* = {AC.rho_star(d, en/a.D):.6f}")
    print(f"  total weight {tot}/{WD} = {tot/WD:.7f} (unchanged)")


if __name__ == '__main__':
    main()
