"""Census of local minima of the width of 4 unit squares in a W x H box (H fixed).
Signature: squares sorted bottom to top; each = flat/tilted + wall contacts (L,R,B,T); canonical under the box's 4 symmetries."""
import numpy as np, sys, json
from scipy.optimize import minimize
from rectpack import penalty, gap, half
H = float(sys.argv[1]); N = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3]))
def minW(z, W):
    zz = np.r_[z, W]
    for mu in (1e2, 1e4, 1e6, 1e8, 1e10):
        zz = minimize(lambda q: q[-1] + mu*penalty(q[:-1], 4, q[-1], H), zz, method='L-BFGS-B', options={'ftol':1e-25,'gtol':1e-16,'maxiter':20000}).x
    z, W = zz[:-1], zz[-1]
    for k in range(400):
        z = minimize(lambda q: penalty(q, 4, W, H), z, method='L-BFGS-B', options={'ftol':1e-25,'gtol':1e-16}).x
        if gap(z, 4, W, H) > -1e-10: return z, W
        W += 1e-4
    return None, None
def sig(z, W):
    x, y, t = z[:4], z[4:8], z[8:]
    def one(x, y, t):
        hx = half(t); a = np.degrees(t) % 90; a = min(a, 90-a)
        s = ('F' if a < 1.0 else 'T')
        s += ''.join(c for c, ok in (('L', x-hx < 2e-3), ('R', W-x-hx < 2e-3), ('B', y-hx < 2e-3), ('T', H-y-hx < 2e-3)) if ok)
        return s
    def word(x, y, t):
        o = np.argsort(y); return '|'.join(one(x[k], y[k], t[k]) for k in o)
    sw = lambda s: s.translate(str.maketrans('LRBT', 'RLTB'))
    cands = []
    for fx in (0, 1):
        for fy in (0, 1):
            xx = W - x if fx else x; yy = H - y if fy else y
            cands.append(word(xx, yy, t))
    return min(cands)
out = []
for s in range(N):
    W0 = rng.uniform(1.9, 2.3)
    z = np.r_[rng.uniform(0.5, W0-0.5, 4), rng.uniform(0.5, H-0.5, 4), rng.uniform(0, np.pi/2, 4)]
    if rng.random() < 0.5: z[8:] = np.where(rng.random(4) < 0.5, 0.0, z[8:])
    z, W = minW(z, W0)
    if z is None: continue
    a = np.degrees(z[8:]) % 90; a = np.minimum(a, 90-a)
    out.append(dict(W=float(W), sig=sig(z, W), tilts=sorted(np.round(a, 1).tolist()), z=z.tolist()))
json.dump(out, open(f'census4_H{H}_s{sys.argv[3]}.json', 'w'))
print(len(out), 'minima')
