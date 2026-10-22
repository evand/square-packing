"""(S): for each tilt t, does one branch hold for all slopes?  min over vertices of kap - G and kap - Zb."""
import math, numpy as np
from vertex import kap, vertices
from vertexS import lines_S
for H in (1.8285, 1.85, 1.87, 1.9, 1.95, 1.99):
    b = 1 - H / 2; Z = 2 + 2 * math.sqrt(2) - 2 * H
    bad = []; grid_ok = []; zig_ok = []
    for td in np.linspace(0, 45, 451):
        t = math.radians(td)
        if math.cos(t) + math.sin(t) > H: continue
        V = [v for v in vertices(lines_S(t, H)) if abs(v[0]) < 50 and abs(v[1]) < 50]
        g = min(kap(t, H, m, mp) - 1 - b * abs(abs(m) - abs(mp)) for m, mp in V)
        z = min(kap(t, H, m, mp) - Z + b * (abs(m) + abs(mp)) for m, mp in V)
        if g >= -1e-12: grid_ok.append(td)
        if z >= -1e-12: zig_ok.append(td)
        if g < -1e-12 and z < -1e-12: bad.append(td)
    rng = lambda L: f"[{min(L):.1f}, {max(L):.1f}]" if L else "none"
    print(f"H={H}: grid branch for all slopes at t in {rng(grid_ok)} ({len(grid_ok)} pts); zig branch at {rng(zig_ok)} ({len(zig_ok)}); neither: {rng(bad)}")
