#!/usr/bin/env python3
"""Where is the empty space in a record?  (10-20; Evan: "junctions between the two-bulb regions in the 270-273 region
look unsophisticated ... more slack than is apparent").

Rasterize the box (pixel h), mark covered pixels, give every empty pixel to its nearest covered pixel's square (EDT).
Squares are classed by tilt a mod 90: U upright (< 0.5 deg from axis), P (0.5, 45), M (45, 89.5) (the two bulb
orientations of the two-bulb records).  Each empty pixel is classed by the classes present among squares within R:
  junction   both P and M near (between the two bulbs)
  boundary   U and one of P / M (bulb against the upright field)
  wall       one class, within R of the box edge
  interior   one class, away from the walls (inside a bulb / the upright field)
Output per n: total waste s^2 - n (exact) and the raster split (in unit squares), interface waste per interface
square, and the raster's own error (raster empty - exact).

  python3 analysis/junction_waste.py 268 269 270 271 272 273 241 211 308 183 [--h 0.01]
"""
import argparse, math, os, sys
import numpy as np
from scipy import ndimage
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from pk.store import Store

R = 0.75


def tilt(a):
    t = a % 90.0
    return min(t, 90.0 - t)


def groups(T):
    order = np.argsort(T); g = np.zeros(len(T), int); cur = 0
    for a, b in zip(order[:-1], order[1:]):
        if T[b] - T[a] > 1.5:
            cur += 1
        g[b] = cur
    return g


def analyse(n, P, h):
    s = float(P.s); sq = np.asarray(P.sq, float)
    m = int(math.ceil(s / h)); xs = (np.arange(m) + 0.5) * h
    X, Y = np.meshgrid(xs, xs, indexing='ij')
    lab = np.full((m, m), -1, int)
    for i, (x, y, a) in enumerate(sq):
        c, d = math.cos(math.radians(a)), math.sin(math.radians(a))
        r = 0.75
        i0, i1 = max(0, int((x - r) / h)), min(m, int((x + r) / h) + 2)
        j0, j1 = max(0, int((y - r) / h)), min(m, int((y + r) / h) + 2)
        dx, dy = X[i0:i1, j0:j1] - x, Y[i0:i1, j0:j1] - y
        u, v = dx * c + dy * d, -dx * d + dy * c
        inside = (np.abs(u) <= 0.5) & (np.abs(v) <= 0.5)
        lab[i0:i1, j0:j1][inside] = i
    empty = lab < 0
    T = np.array([tilt(a) for a in sq[:, 2]])
    am = sq[:, 2] % 90.0
    G = np.where(T < 0.5, 0, np.where(am < 45, 1, 2))
    # groups present within R: dilate each group's mask
    k = int(round(R / h)); st = np.zeros((2 * k + 1, 2 * k + 1), bool)
    yy, xx = np.mgrid[-k:k + 1, -k:k + 1]; st[(xx ** 2 + yy ** 2) <= k * k] = True
    near = {}
    for g in (0, 1, 2):
        mask = np.isin(lab, np.where(G == g)[0])
        near[g] = ndimage.binary_dilation(mask, structure=st) if mask.any() else np.zeros_like(mask)
    wall = (np.minimum(np.minimum(X, s - X), np.minimum(Y, s - Y)) < R)
    A = h * h
    junc = empty & near[1] & near[2]
    bound = empty & ~junc & near[0] & (near[1] | near[2])
    inter = junc | bound
    wallc = empty & ~inter & wall
    inte = empty & ~inter & ~wall
    # interface squares: squares touching the interface zone
    _, ind = ndimage.distance_transform_edt(empty, return_indices=True)
    owner = lab[ind[0], ind[1]]
    isq = np.unique(owner[junc])
    up = (T < 0.5)
    return dict(n=n, s=s, waste=s * s - n, raster=empty.sum() * A, interface=junc.sum() * A, boundary=bound.sum() * A, wall=wallc.sum() * A,
                interior=inte.sum() * A, interior_tilted=(inte & ~up[owner]).sum() * A,
                interior_upright=(inte & up[owner]).sum() * A, n_iface_sq=len(isq), ngroups=int((G == 1).sum()), nM=int((G == 2).sum()),
                tilted=int((~up).sum()))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('ns', type=int, nargs='+'); ap.add_argument('--h', type=float, default=0.01)
    a = ap.parse_args()
    st = Store(); F = st.frontier()
    print(f"{'n':>4} {'s':>10} {'waste':>6} {'rastErr':>7} {'junct':>6} {'wall':>6} {'int.tilt':>8} {'int.up':>6} "
          f"{'#jsq':>5} {'j/sq':>6} {'bound':>6} {'P':>4} {'M':>4}")
    for n in a.ns:
        r = F.get(n)
        if not r or r.get('packing') is None:
            print(n, 'no frontier packing'); continue
        d = analyse(n, st.get(r['packing']), a.h)
        print(f"{n:4d} {d['s']:10.6f} {d['waste']:6.2f} {d['raster'] - d['waste']:+7.3f} {d['interface']:6.2f} "
              f"{d['wall']:6.2f} {d['interior_tilted']:8.2f} {d['interior_upright']:6.2f} {d['n_iface_sq']:5d} "
              f"{d['interface'] / max(1, d['n_iface_sq']):6.3f} {d['boundary']:6.2f} {d['ngroups']:4d} {d['nM']:4d}", flush=True)


if __name__ == '__main__':
    main()
