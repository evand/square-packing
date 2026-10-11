import numpy as np, sys
from feas import attempt
from rectpack import gap
def feasible(n, W, H, starts, rng):
    for s in range(starts):
        v = attempt(n, W, H, rng)
        if gap(v, n, W, H) > -1e-10:
            return True, v
    return False, None
n, H, lo, hi, starts = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
rng = np.random.default_rng(int(sys.argv[6]) if len(sys.argv) > 6 else 0)
ok, v = feasible(n, hi, H, starts, rng)
if not ok:
    print(f"n={n} H={H}: not found even at W={hi}"); sys.exit()
while hi - lo > 2e-3:
    mid = (lo+hi)/2
    ok, vv = feasible(n, mid, H, starts, rng)
    if ok: hi, v = mid, vv
    else: lo = mid
print(f"n={n} H={H}: feasible at W={hi:.4f}, none found at W={lo:.4f}")
np.save(f"cfg_n{n}_H{H}.npy", v)
