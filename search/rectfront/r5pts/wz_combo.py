"""Witness-combination analysis for pair lemmas P1 or P2 on cell Z (1<=x<=a/2, 1<=y<=b/2, y >= 2-k optional).
Witnesses of 'not In' for offset (dx,dy): u>1/2, u<-1/2, v>1/2, v<-1/2 (u=dx c+dy s, v=-dx s+dy c).
For each pair of witnesses: max over samples of min(w1, w2) - 1/2 (>0 means the combo occurs: no contradiction),
and for multipliers (l1,l2) in {(1,1),(c,s),(s,c)}: max of l1 w1 + l2 w2 - (l1+l2)/2 (<= 0 => certificate)."""
import numpy as np, sys
R2 = np.sqrt(2); C0 = 0.91
rng = np.random.default_rng(0); n = 6_000_000
a = rng.uniform(1+R2, 4-R2, n); b = 4+R2-a
th = rng.uniform(0, np.pi/2, n); c, s = np.cos(th), np.sin(th); k = (c+s)/2
x = rng.uniform(1, 1.3, n); y = rng.uniform(1, 1.5, n)
ok = (x <= a/2) & (y <= b/2)
if 'Z2' in sys.argv: ok &= (y >= 2-k)
P = dict(r=(C0, 1+0*a), s=(1+0*a, 1+0*a), m=(1+0*a, b/2), mp=(a-1, b/2), sh=(1+0*a, b-1))
def wit(p):
    dx, dy = P[p][0]-x, P[p][1]-y
    u = dx*c + dy*s; v = -dx*s + dy*c
    return {'u+': u, 'u-': -u, 'v+': v, 'v-': -v}
for p1, p2 in [(sys.argv[1], sys.argv[2])]:
    W1, W2 = wit(p1), wit(p2)
    for k1, w1 in W1.items():
        for k2, w2 in W2.items():
            occ = np.minimum(w1, w2)[ok].max() - 0.5
            if occ > 0:
                print(f"{p1}:{k1} {p2}:{k2}  OCCURS ({occ:+.4f})"); continue
            m1 = w1[ok].max()-0.5; m2 = w2[ok].max()-0.5
            if m1 <= 0 or m2 <= 0:
                print(f"{p1}:{k1} {p2}:{k2}  single impossible ({p1}:{k1} max {m1:+.4f}, {p2}:{k2} max {m2:+.4f})"); continue
            certs = []
            for nm, l1, l2 in [('1,1', 1, 1), ('c,s', c, s), ('s,c', s, c)]:
                v = (l1*w1 + l2*w2 - (l1+l2)/2)[ok].max()
                certs.append(f"{nm}:{v:+.4f}")
            print(f"{p1}:{k1} {p2}:{k2}  combo excluded (min max {occ:+.4f}); certs {' '.join(certs)}")
