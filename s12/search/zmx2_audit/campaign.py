"""Adversarial soundness campaign for zmx2: random adversarial mixed covers + targeted poses
(germs, walls at the admissibility boundary, points/segment ends/line crossings exactly on edges
and corners of Q, tiny and zero u, u on bin boundaries), nested boxes around each pose (probe.py),
check bound <= exact mu."""
import os
import random, sys, math
from fractions import Fraction as Fr
from aud import write, load, inside
import probe


def rand_cover(rng, path):
    s = rng.choice([Fr(2), Fr(3), Fr(5, 2), Fr(12, 5)])
    D = rng.choice([10, 100, 1000, 20, 7 * 8, 1 << 18])
    S = s * D
    if S.denominator != 1:
        D *= S.denominator
        S = s * D
    S = int(S)
    W = rng.choice([1000, 10 ** 6, 97, 10 ** 12])
    # lines: grid, boundary, off-grid, and partners at distance 1
    xs = set()
    for k in range(int(s) + 1):
        if rng.random() < 0.7:
            xs.add(k * D)
    for _ in range(3):
        a = rng.randint(0, S)
        xs.add(a)
        if a + D <= S and rng.random() < 0.7:
            xs.add(a + D)
    xs = sorted(xs)
    segs = []
    for vert in (True, False):
        for p in xs:
            if rng.random() < 0.25:
                continue
            # a few segments along the line, possibly overlapping
            for _ in range(rng.randint(1, 6)):
                L = rng.choice([x for x in (1, 2, 3, D // 10, D // 5, D // 2, D) if 0 < x <= S])
                a = rng.randint(0, S - L)
                b = a + L
                w = rng.randint(1, W)
                segs.append((p, a, p, b, w) if vert else (a, p, b, p, w))
    pts = []
    for _ in range(rng.randint(0, 40)):
        kind = rng.random()
        if kind < 0.4 and xs:
            x, y = rng.choice(xs), rng.randint(0, S)
        elif kind < 0.6 and xs:
            x, y = rng.choice(xs), rng.choice(xs)
        else:
            x, y = rng.randint(0, S), rng.randint(0, S)
        if rng.random() < 0.5:
            x, y = y, x
        pts.append((x, y, rng.randint(0, W)))
    cv = dict(s=s, D=D, W=W, pts=pts, segs=segs)
    write(path, cv)
    return load(path), xs


def rot(u, v):
    n = 1 + u * u
    C, S = (1 - u * u) / n, 2 * u / n
    return lambda a, b: (a * C - b * S, a * S + b * C)


def rand_u(rng):
    k = rng.random()
    if k < 0.2:
        return Fr(0)
    if k < 0.45:
        return Fr(1, rng.choice([10 ** 9, 10 ** 7, 10 ** 5, 10 ** 3, 2 ** 30, 3 * 10 ** 6]))
    if k < 0.6:
        return Fr(rng.randint(1, 4), 8)  # bin boundaries (1/8 .. 1/2)
    return Fr(rng.randint(1, 10 ** 6), 2 * 10 ** 6)


def poses(rng, cv, xs, n):
    s, D = cv['s'], cv['D']
    out = []
    tries = 0
    while len(out) < n and tries < 50 * n:
        tries += 1
        u = rand_u(rng)
        R = rot(u, None)
        k = rng.random()
        half = Fr(1, 2)
        if k < 0.25 and xs:
            # germ: an edge of Q (nearly) on a line, other coordinate random or on a line
            l = Fr(rng.choice(xs), D)
            d = rng.choice([Fr(0), Fr(0), Fr(rng.choice([1, -1]), rng.choice([10 ** 7, 10 ** 4, 10 ** 9]))])
            cx = l + rng.choice([half, -half]) + d
            cy = Fr(rng.choice(xs), D) + rng.choice([half, -half, Fr(rng.randint(-500, 500), 1000)]) if rng.random() < .5 \
                else Fr(rng.randint(0, 10 ** 6), 10 ** 6) * s
            c = (cx, cy)
        elif k < 0.4:
            # wall: admissibility boundary
            w = (1 - u * u + 2 * u) / (1 + u * u)
            cx = rng.choice([w / 2, s - w / 2])
            cy = rng.choice([w / 2, s - w / 2, Fr(rng.randint(0, 10 ** 6), 10 ** 6) * s,
                             (Fr(rng.choice(xs), D) + half) if xs else s / 2])
            c = (cx, cy)
        elif k < 0.7:
            # a point / segment endpoint / line crossing exactly on an edge or a corner of Q
            cand = [(Fr(p[0], D), Fr(p[1], D)) for p in cv['pts']]
            cand += [(Fr(g[0], D), Fr(g[1], D)) for g in cv['segs']] + [(Fr(g[2], D), Fr(g[3], D)) for g in cv['segs']]
            if xs:
                cand += [(Fr(rng.choice(xs), D), Fr(rng.choice(xs), D)) for _ in range(5)]
            if not cand:
                continue
            p = rng.choice(cand)
            e = rng.choice([half, -half])
            t = rng.choice([half, -half, Fr(rng.randint(-1000, 1000), 2000)])
            a, b = (e, t) if rng.random() < .5 else (t, e)
            ra, rb = R(a, b)
            c = (p[0] - ra, p[1] - rb)
        else:
            c = (Fr(rng.randint(0, 10 ** 6), 10 ** 6) * s, Fr(rng.randint(0, 10 ** 6), 10 ** 6) * s)
        P = (c[0], c[1], u)
        if inside(cv, *P):
            out.append(P)
    return out


if __name__ == '__main__':
    seed0 = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    ncov = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    npose = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    os.makedirs(os.environ.get('ZT', 'runs/zmx2_audit'), exist_ok=True)
    tot = 0
    for i in range(ncov):
        rng = random.Random(seed0 * 1000 + i)
        path = os.environ.get('ZT', 'runs/zmx2_audit') + '/rc_%d_%d.txt' % (seed0, i)
        cv, xs = rand_cover(rng, path)
        P = poses(rng, cv, xs, npose)
        f = probe.run(path, P, seed=i, kmax=rng.choice([22, 27]), refl=rng.random() < 0.3, extra=rng.choice([[], [], ['--pair-points'], ['--no-atoms']]))
        tot += f
    print('campaign seed %d: %d FAIL' % (seed0, tot))
