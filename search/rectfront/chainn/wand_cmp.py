# 2026-10-25: our reimplementation (from the paper text only) of wand125 minimum-widths v2 G_N, vs our Morandi family.
# Our own implementation of wand125's G_N (formulas from their paper text) vs our Morandi family F.
import math, sys
sys.path.insert(0, '.')
from targets import F
from scipy.optimize import minimize_scalar
def G(N, th, a):
    c, s = math.cos(th), math.sin(th); h = c + s; b = 1.5*a - 1
    P = (2*h - a)/(h*h - 1); Q = 1 + (2*h - b)/(h*h - 1); A = (2 - a)/4*(math.tan(th) - 1/math.tan(th))
    return (N-4)*P + 2*Q - 2*abs(A)
lo, hi = 2*math.atan(4339/12800), 2*math.atan(4553/8800)
worst = 0
for N in [12, 14, 16, 20, 30, 50]:
    for a in [1.851, 1.86, 1.88, 1.9, 1.92506, 1.95, 1.99]:
        r = minimize_scalar(lambda t: G(N, t, a), bounds=(math.pi/4 + 1e-9, hi), method='bounded', options={'xatol': 1e-13})
        g = min(r.fun, G(N, math.pi/4, a)); f = F(N, a); odd = 4 + 2*(N-2)*math.sqrt(2) - (N-1)*a
        worst = max(worst, abs(g - f))
        print(f"N={N:3d} a={a:.5f}  G_min={g:.10f}  ours={f:.10f}  diff={g-f:+.2e}  45deg={odd:.6f}  th={math.degrees(r.x):.3f}")
print("max |diff|", worst)
# odd check
for N in [11, 13, 21]:
    for a in [1.851, 1.9]:
        print(N, a, 4 + 2*(N-2)*math.sqrt(2) - (N-1)*a, F(N, a))
