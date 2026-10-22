"""(S): is there theta(t) in [0,1] with kap >= theta G + (1-theta) Zb for all slopes?  (max over theta of the
min over vertices; concave in theta -> golden search).  Prints the feasible theta interval per t."""
import math, numpy as np
from vertex import kap, vertices
import vertexS
def worst(t, H, th, V):
    b = 1 - H / 2; Z = 2 + 2 * math.sqrt(2) - 2 * H
    return min(kap(t, H, m, mp) - th * (1 + b * abs(abs(m) - abs(mp))) - (1 - th) * (Z - b * (abs(m) + abs(mp)))
               for m, mp in V)
for H in (1.8285, 1.85, 1.9):
    print(f"H={H}")
    for td in (0, 0.5, 1, 2, 3, 5, 8, 12, 16, 20, 25, 29, 35, 45):
        t = math.radians(td)
        V = [v for v in vertices(vertexS.lines_S(t, H)) if abs(v[0]) < 50 and abs(v[1]) < 50]
        ths = np.linspace(0, 1, 2001)
        w = np.array([worst(t, H, th, V) for th in ths])
        ok = ths[w >= -1e-12]
        print(f"  t={td:5.1f}: best margin {w.max():+.2e} at theta={ths[w.argmax()]:.3f}; feasible theta "
              + (f"[{ok.min():.3f}, {ok.max():.3f}]" if len(ok) else "none"))
