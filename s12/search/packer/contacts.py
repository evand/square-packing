#!/usr/bin/env python3
"""Feature-level contact graphs of packings.

Contact = pair (or square-wall) with separating gap < tol.  For a pair, the separating axes within tol of the best one are
found; each such axis belongs to a side of its owner square; the other square's corners within tol of the extreme on that
axis are its touching features.  Contact type: corner-side (one corner on a side) or side-side (two corners each way).
Features are labelled in each square's local frame: side k = outward normal at angle theta + 90k deg; corner k = between
sides k and k+1.  Angles are taken modulo 90 deg relative to a reference configuration (so a near-axis square flipping
between 89.99 and 0.01 keeps its labels).

  contacts.py ref.txt other1.txt ...    -> diff of each other's contact graph vs ref; and a count of distinct graphs
"""
import math, sys, collections

H = 0.5


def load(path):
    L = [l.split() for l in open(path) if l.strip()]
    n, s = int(L[0][0]), float(L[0][1])
    return s, [(float(a), float(b), math.radians(float(c))) for a, b, c, *_ in L[1:n + 1]]


def align(sq, ref):
    """Shift each angle by a multiple of 90 deg to be closest to the reference square's angle."""
    out = []
    for (x, y, t), (_, _, tr) in zip(sq, ref):
        k = round((tr - t) / (math.pi / 2))
        out.append((x, y, t + k * math.pi / 2))
    return out


def corners(q):
    x, y, t = q
    c, s = math.cos(t), math.sin(t)
    # corner k between side k (normal angle t + 90k) and side k+1
    return [(x + H * (c * a - s * b), y + H * (s * a + c * b)) for a, b in ((1, 1), (-1, 1), (-1, -1), (1, -1))]


def normals(q):
    t = q[2]
    return [(math.cos(t + k * math.pi / 2), math.sin(t + k * math.pi / 2)) for k in range(4)]


def pair_contact(i, j, qi, qj, tol):
    """Return frozenset of feature incidences if gap < tol, else None."""
    best, cand = -1e9, []
    for own, oth, oi, oj in ((qi, qj, i, j), (qj, qi, j, i)):
        co = corners(oth)
        for k, (nx, ny) in enumerate(normals(own)):
            face = nx * own[0] + ny * own[1] + H
            proj = [nx * p[0] + ny * p[1] for p in co]
            g = min(proj) - face
            cand.append((g, oi, k, oj, proj))
            best = max(best, g)
    if best > tol:
        return None
    feats = set()
    for g, oi, k, oj, proj in cand:
        if g >= best - tol and g > -tol:
            m = min(proj)
            touching = tuple(c for c in range(4) if proj[c] - m < tol)
            if len(touching) == 1:
                feats.add(('CS', oj, touching[0], oi, k))           # corner of oj on side k of oi
            else:
                a, b = sorted([(oi, k), (oj, (touching[0] + 1) % 4 if touching == (touching[0], (touching[0] + 1) % 4) else 0)])
                feats.add(('SS', min(oi, oj), max(oi, oj), (oi, k)))
    return frozenset(feats) if feats else None


def graph(s, sq, tol=1e-9):
    n = len(sq)
    E = set()
    for i in range(n):
        for k, (px, py) in enumerate(corners(sq[i])):
            for w, d in (('L', px), ('R', s - px), ('B', py), ('T', s - py)):
                if d < tol:
                    E.add(('CW', i, k, w))
    for i in range(n):
        for j in range(i + 1, n):
            if (sq[i][0] - sq[j][0]) ** 2 + (sq[i][1] - sq[j][1]) ** 2 < 2.0 + 1e-6:
                f = pair_contact(i, j, sq[i], sq[j], tol)
                if f:
                    E |= f
    return frozenset(E)


def describe(e):
    if e[0] == 'CW':
        return f'corner {e[2]} of #{e[1]} on wall {e[3]}'
    if e[0] == 'CS':
        return f'corner {e[2]} of #{e[1]} on side {e[4]} of #{e[3]}'
    return f'side-side #{e[1]}-#{e[2]} (side {e[3][1]} of #{e[3][0]})'


if __name__ == '__main__':
    s0, R = load(sys.argv[1])
    G0 = graph(s0, R)
    print(f'ref {sys.argv[1]}: s={s0:.10f}, {len(G0)} incidences '
          f'({sum(e[0] == "CW" for e in G0)} wall, {sum(e[0] == "CS" for e in G0)} corner-side, {sum(e[0] == "SS" for e in G0)} side-side)')
    graphs = collections.Counter()
    for p in sys.argv[2:]:
        s, sq = load(p)
        G = graph(s, align(sq, R))
        graphs[G] += 1
        if len(sys.argv) <= 6:
            print(f'\n{p}: s={s:.10f}, {len(G)} incidences; lost {len(G0 - G)}, gained {len(G - G0)}')
            for e in sorted(G0 - G, key=str):
                print('   - ' + describe(e))
            for e in sorted(G - G0, key=str):
                print('   + ' + describe(e))
    print(f'\n{len(sys.argv) - 2} configurations, {len(graphs)} distinct contact graphs')
