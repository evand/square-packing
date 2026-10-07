"""Independent floating-point sanity check of the leaf formulas (NOT part of the proof).

For random poses in every box, compute the node values directly from the square's polygon (no
half-chord formulas) and confirm that every leaf claim of a certificate is <= the sampled value.
Also used by controls.py to exhibit a pose that falsifies a deliberately wrong leaf.

    python3 sanity.py certs/h113.json [samples_per_box]
"""
import json, math, random, sys
from fractions import Fraction as Fr


def square(theta, c):
    C, S = math.cos(theta), math.sin(theta)
    return [(0.5 * (sx * C - sy * S), c + 0.5 * (sx * S + sy * C)) for sx, sy in ((1, 1), (-1, 1), (-1, -1), (1, -1))]


def chord(poly, z):
    xs = []
    n = len(poly)
    for i in range(n):
        (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % n]
        if (y1 - z) * (y2 - z) <= 0 and y1 != y2:
            xs.append(x1 + (z - y1) * (x2 - x1) / (y2 - y1))
        elif y1 == y2 == z:
            xs += [x1, x2]
    return (min(xs), max(xs)) if xs else None


def node_value(theta, c, z1, z2):
    """V(z1, z2) = lL(z1) + lR(z2) with END (None) meaning the bounding-box half-width."""
    P = square(theta, c)
    p = max(x for x, _ in P)
    if z1 is None:
        lL = p
    else:
        ch = chord(P, z1)
        lL = -ch[0]
    if z2 is None:
        lR = p
    else:
        ch = chord(P, z2)
        lR = ch[1]
    return lL + lR


def sample_box(key, n, rng):
    t0, t1, c0, c1 = [float(x) for x in key]
    out = []
    for i in range(n):
        if i < 4:
            t = (t0, t1)[i % 2]
        else:
            t = rng.uniform(t0, t1)
        th = 2 * math.atan(t)
        p = (math.cos(th) + math.sin(th)) / 2
        lo = max(c0, p)
        if lo > c1:
            continue
        c = (lo, c1)[(i // 2) % 2] if i < 4 else rng.uniform(lo, c1)
        out.append((th, c))
    return out


def run(path, n=400, seed=1):
    cert = json.load(open(path))
    M = [float(Fr(z)) for z in cert['menu']]
    Z = lambda i: None if i == -1 else M[i]
    rng = random.Random(seed)
    worst = math.inf
    for key, claims in zip(cert['boxes'], cert['leaves']):
        key = [Fr(x) for x in key]
        for th, c in sample_box(key, n, rng):
            for l, r, q in claims:
                v = node_value(th, c, Z(l), Z(r))
                worst = min(worst, v - float(Fr(q)))
    return worst


if __name__ == '__main__':
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    w = run(sys.argv[1], n)
    print(f"min over samples of (true node value - claimed bound) = {w:.3e}  ->  {'OK' if w >= -1e-9 else 'VIOLATION'}")
