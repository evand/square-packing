import numpy as np
from wz_sub import *
rng = np.random.default_rng(5); n = 10_000_000
a = rng.uniform(1+R2, 4-R2, n); b = 4+R2-a
th = rng.uniform(0, np.pi/2, n); c, s = np.cos(th), np.sin(th); k = (c+s)/2
x = rng.uniform(1, 1.3, n); y = rng.uniform(1, 1.5, n)
Z = (x <= a/2) & (y <= b/2); Z2 = Z & (y >= 2-k)
g = pts(a, x, y, c, s)
for xs in [1.06, 1.08, 1.09, 1.1, 1.11, 1.117, 1.13]:
    L = Z2 & (x <= xs); Rr = Z2 & (x >= xs)
    print(f"x*={xs}: left r|rh {np.minimum(g['r'],g['rh'])[L].max():+.4f}  right m&mp {M(g['m'],g['mp'])[Rr].max():+.4f}")
print("Z2 m:", g['m'][Z2].max(), " Z s|sh:", np.minimum(g['s'], g['sh'])[Z].max())
