"""Test sub-claims: worst value of a claim (<= 0 means holds) over a sampled cell."""
import numpy as np
R2 = np.sqrt(2); C0 = 0.91
def G(px, py, x, y, c, s):
    dx, dy = px - x, py - y
    return np.maximum(abs(dx*c + dy*s), abs(-dx*s + dy*c)) - 0.5
def pts(a, x, y, c, s):
    b = 4 + R2 - a
    P = dict(r=(C0, 1), s=(1, 1), m=(1, b/2), rh=(C0, b-1), sh=(1, b-1), mp=(a-1, b/2), rp=(a-C0, 1), sp=(a-1, 1))
    return {k: G(v[0], v[1], x, y, c, s) for k, v in P.items()}
def run(name, cell, claim, n=4_000_000, seed=0):
    rng = np.random.default_rng(seed)
    a = rng.uniform(1+R2, 4-R2, n); b = 4+R2-a
    th = rng.uniform(0, np.pi/2, n); c, s = np.cos(th), np.sin(th); k = (c+s)/2
    x, y, ok = cell(rng, a, b, k, n)
    g = pts(a, x, y, c, s)
    v = claim(g)[ok]
    i = np.argmax(v)
    print(f"{name:40s} worst {v[i]:+.5f} at a={a[ok][i]:.4f} th={np.degrees(th[ok][i]):6.2f} x={x[ok][i]:.4f} y={y[ok][i]:.4f}  (n={ok.sum()})")
M = np.maximum; m_ = np.minimum
def W2b(rng, a, b, k, n):
    x = rng.uniform(C0, 1, n); y = rng.uniform(1, 1.5, n)
    return x, y, (y >= 2-k) & (y <= b/2) & (x >= k)
def W2a(rng, a, b, k, n):
    x = rng.uniform(0.5, C0, n); y = rng.uniform(1, 1.5, n)
    return x, y, (y >= 2-k) & (y <= b/2) & (x >= k)
if __name__ == '__main__':
    run('W2b: r', W2b, lambda g: g['r'])
    run('W2b: r or rh', W2b, lambda g: m_(g['r'], g['rh']))
    run('W2b: m', W2b, lambda g: g['m'])
    run('W2b: s', W2b, lambda g: g['s'])
    run('W2a: r or rh', W2a, lambda g: m_(g['r'], g['rh']))
    run('W2a: m', W2a, lambda g: g['m'])
