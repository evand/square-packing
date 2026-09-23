#!/usr/bin/env python3
"""Check that `verify/` and `xcheck.py` accept an INTERIOR general-segment anchor clique
(`anchorclique.kseg`), i.e. that the certificate format never needed the wall restriction.

Builds two files from a shipped certificate: the same points plus one interior K(p, A) at weight
0 (must VERIFY with the same minimum -- no weight leak, and the Lemma-0 pairwise check must pass)
and at a small positive weight (must VERIFY with the total increased by it).
"""
import sys, os, subprocess, math
from fractions import Fraction as F
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, 'search'))
import anchorclique as AC

src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'certificates/s12_lower_3.931795_sparse.txt')
N = sys.argv[2] if len(sys.argv) > 2 else '2000'
L = open(src).read().rstrip('\n').split('\n')
K, Dn = (int(v) for v in L[0].split())
D = int(L[1]); WD = int(L[2]); m = int(L[3])
pts = L[4:4 + m]
tail = L[4 + m:]
s = F(K, Dn)
print(f"# {src}: s = {s} = {float(s):.6f}, D = {D}, WD = {WD}, {m} points, tail {tail}")

# an interior anchor point and a segment anchor of length ~0.37 at offset ~0.22, in the middle of
# the container -- the shape RECONCILE.md 2 found violated
px = F(int(0.52 * float(s) * D), D)
py = F(int(0.55 * float(s) * D), D)
ex, ey = 0.16, 0.16
rho = 0.185
mx, my = float(px) + ex, float(py) + ey
dx, dy = ex / math.hypot(ex, ey), ey / math.hypot(ex, ey)
v0 = (mx + rho * dy, my - rho * dx)
v1 = (mx - rho * dy, my + rho * dx)
cl = AC.kseg(s, D, int(px * D), int(py * D),
             int(round(v0[0] * D)), int(round(v0[1] * D)),
             int(round(v1[0] * D)), int(round(v1[1] * D)))
assert cl is not None
print(f"# clique: p = ({float(px):.5f}, {float(py):.5f}), wall distance "
      f"{min(float(px), float(py), float(s)-float(px), float(s)-float(py)):.4f}, "
      f"A = ({v0[0]:.5f},{v0[1]:.5f})-({v1[0]:.5f},{v1[1]:.5f}), |A| = {2*rho:.4f}")
imgs = AC.images(s, cl)
print(f"# {len(imgs)} D4 images")

VER = os.path.join(HERE, 'verify/target/release/verify')
for w in (0, 10):
    extra = set()
    for im in imgs:
        for anc in im[0]:
            vs = [(anc[1], anc[2])] if anc[0] == 'P' else [(anc[1], anc[2]), (anc[3], anc[4])]
            for (vx, vy) in vs:
                gx, gy = F(vx) * D, F(vy) * D
                assert gx.denominator == 1 and gy.denominator == 1
                extra.add((int(gx), int(gy)))
    lines = list(pts) + [f"{X} {Y} 0" for (X, Y) in sorted(extra)]
    out = os.path.join(HERE, 'runs', f'kseg_w{w}.txt')
    with open(out, 'w') as f:
        f.write(f"{K} {Dn}\n{D}\n{WD}\n{len(lines)}\n" + "\n".join(lines) + "\n")
        f.write(AC.block(imgs, [w] * len(imgs), D))
        if tail:
            f.write("\n".join(tail) + "\n")
    r = subprocess.run([VER, out, '12', N, '4', '0'], capture_output=True, text=True)
    tot = [l for l in r.stdout.split('\n') if 'total' in l.lower() or 'weight' in l.lower()]
    mn = [l for l in r.stdout.split('\n') if l.startswith('min covered')]
    vd = [l for l in r.stdout.split('\n') if 'VERIFIED' in l or 'REJECT' in l or 'ERROR' in l]
    print(f"\n== clique weight {w}/{WD} ==")
    for l in tot + mn + vd:
        print("   ", l)
