#!/usr/bin/env python3
"""t11_chains: every (near-)tight wall-to-wall chain of a verified positive packing, with its chain bound.

For a path Q_1..Q_k whose link s is separated along an x-type unit normal n_s = (cos tau_s, sin tau_s),
cos tau_s > 0, |tau_s| <= 45 deg (an edge normal of either square of the link):

    n_s.(c_{s+1} - c_s) >= m_s + delta,   m_s = 1/2 + W(theta_{s+1} - theta_s)/2,
    x_1 >= p_1 + delta,   x_k <= T - p_k - delta,   p = (|cos theta| + |sin theta|)/2.

Dividing each link row by cos tau_s and summing the x-components:

    delta * (2 + sum_s 1/cos tau_s)  <=  L + C,
    L = T - p_1 - p_k - sum_s m_s / cos tau_s        (rise-free length budget; = 0 for T axis-parallel squares)
    C = sum_s tan(tau_s) * Dy_s                       (rise credit: the only term that can be positive)

This is exact and valid for ANY path whose links are separated at level delta along the chosen normals.
At the actual configuration, (L + C)/D - delta = [wall slacks + sum_s (gap_s - delta)/cos tau_s]/D >= 0, so
"slack" below is that non-negative quantity and a slack-0 path is a tight chain.  (y-chains: swap x and y.)

Usage:  python3 search/t11_chains.py runs/bandcut_scan_exact110.json [--top 12] [--tol 1e-6]
"""
import argparse, heapq, json, math
from fractions import Fraction


def load(path):
    d = json.load(open(path))
    T = float(d['T'])
    delta = float(Fraction(int(d['delta_num']), int(d['delta_den'])))
    sq = [(float(Fraction(a)), float(Fraction(b)), math.atan2(float(Fraction(s)), float(Fraction(c))))
          for a, b, c, s in d['sq']]
    return T, delta, sq


def W(D):
    return abs(math.cos(D)) + abs(math.sin(D))


def tilt_deg(th):
    t = math.degrees(th) % 90.0
    return t if t <= 45.0 else t - 90.0


def build(T, delta, sq, axis, tol, same_row=True):
    """Directed links i -> j along `axis`-type normals, separated at level delta (gap >= delta - tol).
    Returns edges[i] = list of (j, w, tau, m, gap, dperp) with w = (gap - delta)/cos tau >= 0."""
    n = len(sq)
    edges = [[] for _ in range(n)]
    for i in range(n):
        xi, yi, ti = sq[i]
        for j in range(n):
            if i == j:
                continue
            xj, yj, tj = sq[j]
            dx, dy = xj - xi, yj - yi
            if axis == 'y':
                dx, dy = dy, dx
            if dx <= 0 or dx * dx + dy * dy > 9.0:          # only forward, only neighbours
                continue
            m = 0.5 + 0.5 * W(tj - ti)
            if same_row and abs(tilt_deg(ti)) < 1 and abs(tilt_deg(tj)) < 1 and abs(dy) > 0.5:
                continue        # axis-parallel pairs: keep only same-row links (tan tau = 0: L, C unchanged)
            best = None
            for th in (ti, tj):
                for k in range(4):
                    a = th + k * math.pi / 2
                    nx, ny = math.cos(a), math.sin(a)
                    if axis == 'y':
                        nx, ny = ny, nx                      # coordinates already swapped
                    if nx <= abs(ny) - 1e-12:                # not along-type, or pointing backwards
                        continue
                    gap = nx * dx + ny * dy - m
                    if gap < delta - tol:
                        continue
                    tau = math.atan2(ny, nx)
                    w = max(gap - delta, 0.0) / nx
                    if best is None or w < best[1]:
                        best = (j, w, tau, m, gap, dy)
            if best:
                edges[i].append(best)
    return edges


def walls(T, delta, sq, axis):
    lo, hi = [], []
    for (x, y, th) in sq:
        c = x if axis == 'x' else y
        p = 0.5 * W(th)
        lo.append(c - p - delta)
        hi.append(T - c - p - delta)
    return lo, hi


def k_best_paths(n, edges, lo, hi, top, cap):
    """Best-first enumeration of wall-to-wall paths by total slack (all weights >= 0), slack <= cap."""
    # exact remaining-cost heuristic by reverse Dijkstra
    h = list(hi)
    rev = [[] for _ in range(n)]
    for i in range(n):
        for e in edges[i]:
            rev[e[0]].append((i, e[1]))
    pq = [(h[i], i) for i in range(n)]
    heapq.heapify(pq)
    while pq:
        d, v = heapq.heappop(pq)
        if d > h[v] + 1e-15:
            continue
        for (u, w) in rev[v]:
            if d + w < h[u] - 1e-15:
                h[u] = d + w
                heapq.heappush(pq, (h[u], u))
    out, cnt = [], 0
    pq = [(lo[i] + h[i], lo[i], (i,), ()) for i in range(n) if lo[i] + h[i] <= cap]
    heapq.heapify(pq)
    while pq and len(out) < top:
        f, g, path, links = heapq.heappop(pq)
        v = path[-1]
        if links and links[-1] is None:                     # terminated
            out.append((g, path, links[:-1]))
            continue
        if g + hi[v] <= cap:
            heapq.heappush(pq, (g + hi[v], g + hi[v], path, links + (None,)))
        for e in edges[v]:
            j, w = e[0], e[1]
            if j in path:
                continue
            if g + w + h[j] <= cap:
                heapq.heappush(pq, (g + w + h[j], g + w, path + (j,), links + (e,)))
        cnt += 1
        if cnt > 2_000_000:
            break
    return out


def describe(T, delta, sq, axis, slack, path, links):
    p1, pk = 0.5 * W(sq[path[0]][2]), 0.5 * W(sq[path[-1]][2])
    L = T - p1 - pk - sum(e[3] / math.cos(e[2]) for e in links)
    C = sum(math.tan(e[2]) * e[5] for e in links)
    D = 2 + sum(1 / math.cos(e[2]) for e in links)
    perp = 1 if axis == 'x' else 0
    rise = sq[path[-1]][perp] - sq[path[0]][perp]
    taus = [math.degrees(e[2]) for e in links]
    tl = [tilt_deg(sq[i][2]) for i in path]
    mixed = sum(e[3] - 1 for e in links)                     # sum (W-1)/2: mixed-tilt cost
    dperp = [e[5] for e in links]
    net = [1 - e[3] / math.cos(e[2]) + math.tan(e[2]) * e[5] for e in links]
    return dict(dperp=dperp, net=net, k=len(path), slack=slack, L=L, C=C, D=D, bound=(L + C) / D, rise=rise,
                taus=taus, tilts=tl, mixed=mixed, path=path,
                n_tilted=sum(1 for t in tl if abs(t) > 1.0),
                spread=(max(taus) - min(taus)) if taus else 0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('inp')
    ap.add_argument('--top', type=int, default=12)
    ap.add_argument('--tol', type=float, default=1e-9)
    ap.add_argument('--cap', type=float, default=0.05, help='list paths with total slack <= cap')
    ap.add_argument('--json', default=None)
    ap.add_argument('--detail', type=int, default=6)
    a = ap.parse_args()
    T, delta, sq = load(a.inp)
    n = len(sq)
    tl = [abs(tilt_deg(s[2])) for s in sq]
    print('T = %g  n = %d  delta = %+.9e   near-axis (<1 deg): %d   tilted: %d'
          % (T, n, delta, sum(t < 1 for t in tl), sum(t >= 1 for t in tl)))
    res = {}
    for axis in 'xy':
        edges = build(T, delta, sq, axis, a.tol)
        lo, hi = walls(T, delta, sq, axis)
        paths = k_best_paths(n, edges, lo, hi, a.top, a.cap)
        print('\n=== %s-chains, wall to wall, by slack (slack 0 = every row tight); %d links in graph ==='
              % (axis, sum(len(e) for e in edges)))
        print('  #   k  tilted   slack        L (rise-free)   C (rise credit)   (L+C)/D      rise    mixed-tilt   normal spread')
        rows = []
        for r, (g, path, links) in enumerate(paths):
            d = describe(T, delta, sq, axis, g, path, links)
            rows.append(d)
            print(' %2d  %2d   %2d    %.3e   %+.6f       %+.6f        %+.3e   %+.4f   %.4f      %.2f deg'
                  % (r, d['k'], d['n_tilted'], d['slack'], d['L'], d['C'], d['bound'], d['rise'],
                     d['mixed'], d['spread']))
        seen = {}
        for r, d in enumerate(rows):
            key = tuple(i for i, t in zip(d['path'], d['tilts']) if abs(t) > 1.0)
            seen.setdefault(key, []).append(r)
        print('  distinct tilted sub-chains among the listed paths: %d' % len(seen))
        for key, rs in list(seen.items())[:a.detail]:
            d = rows[rs[0]]
            print('  rows %s: tilts   %s' % (rs[:6], ' '.join('%+.1f' % t for t in d['tilts'])))
            print('           normals   %s' % ' '.join('%+.1f' % t for t in d['taus']))
            print('           Dperp     %s' % ' '.join('%+.3f' % e for e in d['dperp']))
            print('           link net  %s   (= 1 - m/cos tau + tan tau * Dperp; 0 for an axis-parallel link)'
                  % ' '.join('%+.3f' % e for e in d['net']))
        res[axis] = rows
    if a.json:
        json.dump(res, open(a.json, 'w'))


if __name__ == '__main__':
    main()
