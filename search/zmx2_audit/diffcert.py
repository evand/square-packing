"""Differential test of the full `zmx2 cert` path (with float pre-filters, inherited point sets,
root enumeration, both passes): random mixed covers scaled so that the exact mu at a known pose is
1 - eps < 1.  zmx2 must answer NOT VERIFIED.  (Also runs the 1 + eps' version, informational.)"""
import os
import random, sys, math, subprocess
from fractions import Fraction as Fr
import numpy as np
from aud import write, load, mu, inside, FCover

Z = os.environ.get('ZMX2', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'verify2', 'target', 'release', 'zmx2'))
T = os.environ.get('ZT', 'runs/zmx2_audit') + '/'


def gen(rng):
    s = rng.choice([Fr(2), Fr(12, 5), Fr(3)])
    D = rng.choice([100, 1000, 60])
    S = int(s * D)
    W = 10 ** 9
    q = D // 10
    lines = set(range(0, S + 1, D))
    for _ in range(rng.randint(0, 3)):
        lines.add(rng.randrange(0, S + 1, q))
    segs, pts = [], []
    for vert in (True, False):
        for p in sorted(lines):
            base = rng.random()
            for i in range(S // q):
                w = int(W * q / D * max(0.0, base + rng.uniform(-.5, .5)))
                if w > 0:
                    segs.append((p, i * q, p, (i + 1) * q, w) if vert else (i * q, p, (i + 1) * q, p, w))
    for _ in range(rng.randint(0, 30)):
        x = rng.choice(sorted(lines)) if rng.random() < .5 else rng.randint(0, S)
        y = rng.randint(0, S)
        if rng.random() < .5:
            x, y = y, x
        pts.append((x, y, rng.randint(1, W // 5)))
    return dict(s=s, D=D, W=W, pts=pts, segs=segs)


def fmin(cv, rng):
    fc = FCover(cv)
    s = float(cv['s'])
    best = []
    ths = [0.0, 1e-6] + list(np.linspace(0.02, math.pi / 2 - 0.02, 24))
    for th in ths:
        w = math.cos(th) + math.sin(th)
        for x in np.linspace(w / 2, s - w / 2, 25):
            for y in np.linspace(w / 2, s - w / 2, 25):
                best.append((fc.mu(x, y, th), x, y, th))
    best.sort()
    cands = best[:6]
    out = []
    for (m, x, y, th) in cands:
        step = 0.02
        for it in range(300):
            nx, ny, nt = x + rng.gauss(0, step), y + rng.gauss(0, step), abs(th + rng.gauss(0, step))
            if nt >= math.pi / 2 or not fc.adm(nx, ny, nt):
                continue
            mm = fc.mu(nx, ny, nt)
            if mm < m:
                m, x, y, th = mm, nx, ny, nt
            if it % 60 == 59:
                step /= 3
        out.append((m, x, y, th))
    out.sort()
    m, x, y, th = out[0]
    # rational pose, exact mu; fold theta into [0, 90)
    u = Fr(math.tan(th / 2)).limit_denominator(10 ** 9)
    X, Y = Fr(x).limit_denominator(10 ** 9), Fr(y).limit_denominator(10 ** 9)
    if not inside(cv, X, Y, u):
        return None
    return mu(cv, X, Y, u), (X, Y, u)


def scaled(cv, f):
    c = dict(cv)
    c['pts'] = [(a, b, int(w * f)) for (a, b, w) in cv['pts']]
    c['segs'] = [(a, b, cc, d, int(w * f)) for (a, b, cc, d, w) in cv['segs']]
    return c


LAST = []


def verdict(path):
    r = subprocess.run([Z, 'cert', path, '--full', '--threads', '2', '--uncert-cap', '3', '--node-cap', '2000000'],
                       capture_output=True, text=True)
    LAST.clear()
    for l in r.stdout.split('\n'):
        if l.startswith('UNCERT'):
            f = l.split()
            pas = int(f[3])
            b = [Fr(t) for t in f[5].split(',')]
            LAST.append((pas, b))
    for l in r.stdout.split('\n'):
        if l.startswith(('VERIFIED', 'NOT VERIFIED', 'INCOMPLETE', 'ERROR')):
            return l.split(':')[0]
    return 'rc%d %s' % (r.returncode, (r.stdout + r.stderr)[-200:])


if __name__ == '__main__':
    seed0, n = int(sys.argv[1]), int(sys.argv[2])
    eps = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(1, 10 ** 4)
    os.makedirs(T, exist_ok=True)
    bad = 0
    for i in range(n):
        rng = random.Random(seed0 * 7919 + i)
        cv = gen(rng)
        r = fmin(cv, rng)
        if r is None or r[0] == 0:
            print(i, 'skip (no usable minimum)')
            continue
        m0, P = r
        f = (1 - eps) / m0
        cvr = scaled(cv, f)
        mr = mu(cvr, *P)
        assert mr < 1
        pr = T + 'dc_%d_%d_rej.txt' % (seed0, i)
        write(pr, cvr)
        v = verdict(pr)
        ok = v == 'NOT VERIFIED'
        # distance from the known violating pose to the nearest uncertified box (pass 1 = reflected frame)
        dist = 9
        for pas, b in LAST:
            X, Y, U = P
            if pas == 1:
                X, Y, U = X, cv['s'] - Y, (1 - U) / (1 + U)  # theta -> 90 - theta
            dd = max(max(b[0] - X, X - b[1], 0), max(b[2] - Y, Y - b[3], 0), max(b[4] - U, U - b[5], 0))
            dist = min(dist, float(dd))
        v += ' (nearest uncert box %.2g from P)' % dist
        bad += not ok
        cva = scaled(cv, Fr(1003, 1000) / m0)
        pa = T + 'dc_%d_%d_ok.txt' % (seed0, i)
        write(pa, cva)
        va = verdict(pa)
        print('%d s=%s D=%d pieces=%d pts=%d  min pose (%.5f,%.5f,th %.4f deg) exact mu*f = %.9f -> %s %s | x1.003/min -> %s' % (
            i, cv['s'], cv['D'], len(cv['segs']), len(cv['pts']), float(P[0]), float(P[1]),
            2 * math.degrees(math.atan(float(P[2]))), float(mr), v, '' if ok else '  <<<<< SOUNDNESS?', va), flush=True)
    print('diffcert seed %d: %d unexpected' % (seed0, bad))
