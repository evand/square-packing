from wz_sub import *
def Z1(rng, a, b, k, n):
    x = rng.uniform(1, 1.3, n); y = rng.uniform(1, 1.5, n)
    return x, y, (x <= a/2) & (y <= b/2) & (y <= 2-k)
def Z2(rng, a, b, k, n):
    x = rng.uniform(1, 1.3, n); y = rng.uniform(1, 1.5, n)
    return x, y, (x <= a/2) & (y <= b/2) & (y >= 2-k)
def Z(rng, a, b, k, n):
    x = rng.uniform(1, 1.3, n); y = rng.uniform(1, 1.5, n)
    return x, y, (x <= a/2) & (y <= b/2)
run('Z1: r', Z1, lambda g: g['r'])
run('Z1: s', Z1, lambda g: g['s'])
run('Z2: m', Z2, lambda g: g['m'])
run('Z2: mp', Z2, lambda g: g['mp'])
run('Z2: s or sh', Z2, lambda g: m_(g['s'], g['sh']))
run('Z: m and mp', Z, lambda g: M(g['m'], g['mp']))
run('Z: s or sh', Z, lambda g: m_(g['s'], g['sh']))
run('Z: (r&s) or (m&mp)', Z, lambda g: m_(M(g['r'], g['s']), M(g['m'], g['mp'])))
