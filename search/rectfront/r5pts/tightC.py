"""Tight families of cover C (c0 = 0.9): probes with capture exactly 1, grouped by the centre cell (folded to the
lower-left quadrant) and the caught subset.  Labels: r = (c0,1) w 2/3, s = (1,1) w 1/3, m = (1,b/2) w 1/3, with
' for the mirror x -> a - x and ^ for y -> b - y."""
import sys, math, collections, numpy as np, ptcover as pc
from c10 import cover, line_b
LAB = ["r", "r'", "s", "s'"]
def labels(a, b, c0):
    P, W = cover(a, b, c0)
    names = []
    for X in ("r", "r'", "s", "s'"):
        for y in ("", "^"): names.append(X + y)
    names += ["m", "m'"]
    return P, W, names
def probes(P, a, b, th):
    L, k = pc.lines_at(P, th, a, b)
    N = np.array([l[0] for l in L]); D = np.array([l[1] for l in L])
    i, j = np.triu_indices(len(L), 1)
    det = N[i, 0]*N[j, 1] - N[i, 1]*N[j, 0]; ok = np.abs(det) > 1e-12
    i, j, det = i[ok], j[ok], det[ok]
    V = np.stack([(D[i]*N[j, 1] - D[j]*N[i, 1])/det, (N[i, 0]*D[j] - N[j, 0]*D[i])/det], 1)
    ang = np.arange(16)*np.pi/8 + 0.1
    M = (V[:, None, :] + 1e-7*np.stack([np.cos(ang), np.sin(ang)], 1)[None]).reshape(-1, 2)
    return M[(M[:, 0] >= k) & (M[:, 0] <= a - k) & (M[:, 1] >= k) & (M[:, 1] <= b - k)], k
def cell(x, y, a, b, c0):
    x = min(x, a - x); y = min(y, b - y)
    cx = "x<c0" if x < c0 else ("c0<x<1" if x < 1 else "1<x")
    cy = "y<1" if y < 1 else "1<y"
    return cx + "," + cy
if __name__ == '__main__':
    c0 = 0.9
    agg = collections.defaultdict(lambda: [999, -1, 0])
    for a in [1 + math.sqrt(2), 2.47, 2.53, 4 - math.sqrt(2)]:
        b = line_b(a); P, W, names = labels(a, b, c0)
        for th in np.radians(np.arange(0, 90, 0.25)):
            M, k = probes(P, a, b, th)
            c, s = np.cos(th), np.sin(th)
            dx = P[None, :, 0] - M[:, None, 0]; dy = P[None, :, 1] - M[:, None, 1]
            ins = (np.abs(dx*c + dy*s) <= .5 + 1e-12) & (np.abs(-dx*s + dy*c) <= .5 + 1e-12)
            cap = ins @ W
            for r in np.where(cap < 1 + 1e-9)[0]:
                x, y = M[r]
                # fold names to the quadrant
                fx = x > a/2; fy = y > b/2
                nm = []
                for t in np.where(ins[r])[0]:
                    n = names[t]
                    if fx: n = n.replace("'", "!") if "'" in n else (n[0] + "'" + n[1:] if n[0] in "rsm" else n)
                    if fx: n = n.replace("!", "")
                    if fy: n = n.replace("^", "#") if "^" in n else n + "^"
                    if fy: n = n.replace("#", "")
                    nm.append(n)
                key = (cell(x, y, a, b, c0), " + ".join(sorted(nm)))
                d = math.degrees(th) % 90
                g = agg[key]; g[0] = min(g[0], d); g[1] = max(g[1], d); g[2] += 1
    for (cl, sub), (lo, hi, n) in sorted(agg.items()):
        print(f"{cl:14s} {sub:22s} angles {lo:5.2f}..{hi:5.2f}  probes {n}")
