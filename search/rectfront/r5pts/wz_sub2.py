from wz_sub import *
def band(rng, a, b, k, n):
    x = rng.uniform(0.5, 1, n); y = rng.uniform(1, 2, n)
    return x, y, (y <= b-1) & (x >= k)
def band2(rng, a, b, k, n):
    x = rng.uniform(0.5, 1.2, n); y = rng.uniform(1, 2, n)
    return x, y, (y <= b-1) & (x >= k)
run('band x in [k,1], y in [1,b-1]: r or rh', band, lambda g: m_(g['r'], g['rh']))
run('band x in [k,1.2]: r or rh', band2, lambda g: m_(g['r'], g['rh']))
