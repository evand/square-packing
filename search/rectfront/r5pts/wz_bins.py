import numpy as np
from wz_sub import *
rng = np.random.default_rng(3); n = 8_000_000
a = rng.uniform(1+R2, 4-R2, n); b = 4+R2-a
th = rng.uniform(0, np.pi/2, n); c, s = np.cos(th), np.sin(th); k = (c+s)/2
t = rng.uniform(0, 1, n); x = 1 + t*(a/2-1); y = rng.uniform(1, 1.5, n)
ok = (y <= b/2) & (y >= 2-k)
g = pts(a, x, y, c, s)
rs = M(g['r'], g['s']); mm = M(g['m'], g['mp'])
for lo in np.arange(0, 1, 0.1):
    sel = ok & (t >= lo) & (t < lo+0.1)
    print(f"t in [{lo:.1f},{lo+.1:.1f}): max(r&s) {rs[sel].max():+.4f}  max(m&mp) {mm[sel].max():+.4f}  max(min) {np.minimum(rs,mm)[sel].max():+.4f}")
