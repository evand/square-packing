#!/usr/bin/env python3
"""Constructive channel layouts (contact-preserving by construction) + exact axis fill.

Genome = (left, bottom): each a list of rows (w, theta_deg).
Left channel: row r = w unit squares along u = (cos t, sin t); its leftmost point touches x = 0; row 0 is placed as high as
possible (touching the top wall), each later row as high as possible without overlapping anything already placed
("slide down from the top until contact").  Bottom channel: the same construction in transposed coordinates
(x <-> y), with the left channel as an obstacle, so its columns stand on the bottom wall and stack leftward from the right wall.
A row that cannot be placed inside the box ends its channel.
Score = |tilted| + exact axis fill (gen.mis_fill2).

  layout.py --s S --n N --seed k --time T --out prefix [--init 3x6@27/3x6@62]
"""
import argparse, math, random, time
from gen import pen, wall_viol, mis_fill2, write, corners


def row_squares(w, t, k=0):
    """Row = w + k squares along u anchored at x = 0, minus its first k (a lattice segment set back from the wall)."""
    tr = math.radians(t)
    ki = int(math.floor(k)); kf = k - ki
    full = [(j * math.cos(tr), j * math.sin(tr), t % 90) for j in range(w + ki)]
    pts = [p for q in full for p in corners(*q)]
    mx = min(p[0] for p in pts)
    sq = [(x - mx + kf * math.cos(tr), y + kf * math.sin(tr), a) for x, y, a in full[ki:]]
    pts = [p for q in sq for p in corners(*q)]
    return sq, max(p[1] for p in pts) - min(p[1] for p in pts), min(p[1] for p in pts)


def place_channel(s, rows, obstacles):
    """Rows anchored at x = 0, slid down from the top until contact.  Returns placed squares."""
    placed = []
    for row in rows:
        w, t = row[0], row[1]
        k = row[2] if len(row) > 2 else 0
        base, height, ymin0 = row_squares(w, t, k)
        if max(x for x, _, _ in base) + 0.5 * (abs(math.cos(math.radians(t))) + abs(math.sin(math.radians(t)))) > s + 1e-9:
            break
        hi = s - (ymin0 + height)          # offset putting the row's top point on the top wall
        lo = -ymin0                         # offset putting its bottom point on the bottom wall
        if hi < lo:
            break
        obs = obstacles + placed

        def ok(dy):
            for x, y, a in base:
                q = (x, y + dy, a)
                if any(pen(q, p) > 1e-12 for p in obs if abs(p[0] - q[0]) < 1.5 and abs(p[1] - q[1]) < 1.5):
                    return False
            return True
        if ok(hi):
            dy = hi
        else:
            # find the highest free offset below hi: scan down, then bisect
            step = 0.05
            d = hi
            while d >= lo and not ok(d):
                d -= step
            if d < lo:
                break
            a, b = d, min(d + step, hi)
            for _ in range(45):
                m = 0.5 * (a + b)
                if ok(m):
                    a = m
                else:
                    b = m
            dy = a
        placed += [(x, y + dy, a_) for x, y, a_ in base]
    return placed


def transpose(sq):
    return [(y, x, (90 - a) % 90) for x, y, a in sq]


def build(s, left, bottom):
    L = place_channel(s, left, [])
    B = transpose(place_channel(s, bottom, transpose(L)))
    T = [q for q in L + B if wall_viol(q, s) <= 1e-9]
    return T


def score(s, g, bonus=0.25):
    T = build(s, *g)
    ax = mis_fill2(s, T, 1e-6)
    f = len(T) + len(ax)
    if bonus:
        f += bonus * max(0, len(mis_fill2(s, T, 0.03)) - len(ax))
    return f, T, ax


def mutate(g, rng):
    left, bottom = [[(r[0], r[1], r[2] if len(r) > 2 else 0) for r in c] for c in g]
    ch = left if rng.random() < 0.5 else bottom
    k = rng.randrange(6)
    if k <= 2 and ch:                                  # angle of one row, or of all rows
        i = rng.randrange(len(ch))
        d = rng.gauss(0, [0.3, 1.5, 4][rng.randrange(3)])
        if k == 2:
            for j in range(len(ch)):
                ch[j] = (ch[j][0], ch[j][1] + d, ch[j][2])
        else:
            ch[i] = (ch[i][0], ch[i][1] + d, ch[i][2])
    elif k == 3 and ch:                                # width or set-back of one row
        i = rng.randrange(len(ch))
        if rng.random() < 0.5:
            ch[i] = (max(1, min(6, ch[i][0] + rng.choice((-1, 1)))), ch[i][1], ch[i][2])
        else:
            if rng.random() < 0.5:
                ch[i] = (ch[i][0], ch[i][1], max(0, min(4, round(ch[i][2]) + rng.choice((-1, 1)))))
            else:
                ch[i] = (ch[i][0], ch[i][1], max(0.0, min(4.0, ch[i][2] + rng.gauss(0, 0.15))))
    elif k == 4:                                       # add a row (copy of the last)
        ch.append(ch[-1] if ch else (3, 27.0, 0))
    elif ch:                                           # remove the last row
        ch.pop()
    return (left, bottom)


def parse(spec):
    """'3x6@27/3x6@62' -> left 6 rows of 3 at 27 deg, bottom 6 rows of 3 at 62 deg (bottom angles given in the
    transposed frame: 62 there is 28 deg in the picture)."""
    out = []
    for part in spec.split('/'):
        w, rest = part.split('x'); L, t = rest.split('@')
        out.append([(int(w), float(t))] * int(L))
    return tuple(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--s', type=float, required=True); ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--time', type=float, default=600)
    ap.add_argument('--temp', type=float, default=0.5); ap.add_argument('--init', default='3x6@27/3x6@27')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    g = parse(a.init)
    f, T, ax = score(a.s, g)
    best = f
    t0 = time.time(); it = 0
    print(f'init f={f:.2f} |T|={len(T)} need {a.n}', flush=True)
    while time.time() - t0 < a.time:
        it += 1
        g2 = mutate(g, rng)
        f2, T2, ax2 = score(a.s, g2)
        if f2 >= f or rng.random() < math.exp((f2 - f) / a.temp):
            g, f, T, ax = g2, f2, T2, ax2
        if f > best:
            best = f
            write(f'{a.out}.best.txt', a.s, T + ax)
            print(f'it={it} f={f:.2f} |T|={len(T)} left={[(r[0], round(r[1], 2), r[2]) for r in g[0]]} '
                  f'bottom={[(r[0], round(r[1], 2), r[2]) for r in g[1]]} t={time.time() - t0:.0f}s', flush=True)
        if len(T) + len(ax) >= a.n:
            write(f'{a.out}.hit.txt', a.s, (T + ax)[:a.n])
            print(f'HIT it={it} count={len(T) + len(ax)}', flush=True)
            break
    print(f'done best={best:.2f} its={it}', flush=True)


if __name__ == '__main__' and __import__('sys').argv[1:2] != ['hyb']:
    main()


def count_at(s, g):
    T = build(s, *g)
    return len(T) + len(mis_fill2(s, T, 1e-6)), T


def s_star(g, n, lo, hi, tol=1e-5):
    """Smallest side (bisection, assuming count is monotone in s) at which genome g reaches n squares; hi if never."""
    if count_at(hi, g)[0] < n:
        return float('inf')
    while hi - lo > tol:
        m = 0.5 * (lo + hi)
        if count_at(m, g)[0] >= n:
            hi = m
        else:
            lo = m
    return hi


def main_sstar():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True); ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--time', type=float, default=600)
    ap.add_argument('--temp', type=float, default=0.003); ap.add_argument('--init', default='3x6@27/3x6@27')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    lo, hi = a.k - 0.3, a.k + 0.02
    g = parse(a.init)
    f = s_star(g, a.n, lo, hi)
    best = f
    t0 = time.time(); it = 0
    print(f'init s*={f:.6f}', flush=True)
    while time.time() - t0 < a.time:
        it += 1
        g2 = mutate(g, rng)
        f2 = s_star(g2, a.n, lo, hi)
        if f2 <= f or (f2 < float('inf') and rng.random() < math.exp(-(f2 - f) / a.temp)):
            g, f = g2, f2
        if f < best - 1e-7:
            best = f
            T = build(f, *g)
            write(f'{a.out}.best.txt', f, T + mis_fill2(f, T, 1e-6))
            print(f'it={it} s*={f:.6f} left={[(r[0], round(r[1], 2), r[2]) for r in g[0]]} '
                  f'bottom={[(r[0], round(r[1], 2), r[2]) for r in g[1]]} t={time.time() - t0:.0f}s', flush=True)
    print(f'done best s*={best:.6f} its={it}', flush=True)


def axis_chain(ax, c):
    """Longest chain of axis squares along coordinate c (0: x, 1: y) whose consecutive members overlap in the other
    projection (|d_other| < 1).  Two axis squares overlapping in one projection are >= 1 apart in the other, so a chain of
    length m forces side >= m.  A straight full line is the special case; staggered chains are the rest (10-06: run B
    jammed at exactly 11 on one with no full line)."""
    o = 1 - c
    P = sorted(ax, key=lambda q: q[c])
    best = [1] * len(P)
    for j in range(len(P)):
        for i in range(j):
            if P[j][c] - P[i][c] > 0.5 and abs(P[j][o] - P[i][o]) < 1 - 1e-9 and best[i] + 1 > best[j]:
                best[j] = best[i] + 1
    return max(best, default=0)


def full_lines(path, k):
    """Number of horizontal/vertical lines crossed by >= k near-axis squares (each such line forces side >= k); if there is
    none, 1 when a staggered axis chain (axis_chain) of length >= k exists (also forces side >= k), else 0."""
    L = [l.split() for l in open(path) if l.strip()]
    n = int(L[0][0])
    ax = [(float(x), float(y)) for x, y, a, *_ in L[1:n + 1] if min(float(a) % 90, 90 - float(a) % 90) < 1.0]
    cnt = _straight_lines(ax, k)
    if cnt == 0 and max(axis_chain(ax, 0), axis_chain(ax, 1)) >= k:
        return 1
    return cnt


def _straight_lines(ax, k):
    cnt = 0
    for c in (0, 1):
        seen = set()
        for p in ax:
            for off in (-0.45, 0.0, 0.45):
                y0 = p[c] + off
                m = frozenset(i for i, q in enumerate(ax) if abs(q[c] - y0) < 0.5 - 1e-9)
                if len(m) >= k and m not in seen:
                    seen.add(m); cnt += 1
    return cnt


def s_hyb(g, n, k, tmp, step=0.15, span=0.8):
    """Build at the first side (scanning up from k - 0.1) where the count reaches n, then squeeze with packer."""
    import subprocess, os
    s = k - 0.1
    while s < k + span:
        c, T = count_at(s, g)
        if c >= n:
            ax = mis_fill2(s, T, 1e-6)
            a = os.path.join(tmp, 'h.txt'); b = os.path.join(tmp, 'hs.txt')
            write(a, s, (T + ax)[:n])
            r = subprocess.run([os.path.join(os.path.dirname(os.path.abspath(__file__)), 'target/release/packer'),
                                'relax', '--in', a, '--squeeze-pen', '--out', b], capture_output=True, text=True)
            try:
                sq_s = float(r.stdout.split(' s=')[1].split()[0])
            except Exception:
                return float('inf'), None
            if full_lines(b, k) > 0:
                return float('inf'), None
            return sq_s, b
        s += step
    return float('inf'), None


def main_hyb():
    import tempfile, shutil
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True); ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--seed', type=int, default=1); ap.add_argument('--time', type=float, default=600)
    ap.add_argument('--temp', type=float, default=0.005); ap.add_argument('--init', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    tmp = tempfile.mkdtemp()
    g = eval(a.init) if a.init.startswith('(') else parse(a.init)
    f, p = s_hyb(g, a.n, a.k, tmp)
    best = f
    if p: shutil.copy(p, f'{a.out}.best.txt')
    t0 = time.time(); it = 0
    print(f'init s={f:.6f}', flush=True)
    while time.time() - t0 < a.time:
        it += 1
        g2 = mutate(g, rng)
        f2, p2 = s_hyb(g2, a.n, a.k, tmp)
        if f2 < best - 1e-7 and p2:
            best = f2
            shutil.copy(p2, f'{a.out}.best.txt')
            print(f'it={it} s={f2:.6f} left={[(r[0], round(r[1], 2), round(r[2], 2)) for r in g2[0]]} '
                  f'bottom={[(r[0], round(r[1], 2), round(r[2], 2)) for r in g2[1]]} t={time.time() - t0:.0f}s', flush=True)
        if f2 <= f or (f2 < float('inf') and rng.random() < math.exp(-(f2 - f) / a.temp)):
            g, f = g2, f2
    print(f'done best s={best:.6f} its={it}', flush=True)


if __name__ == '__main__' and len(__import__('sys').argv) > 1 and __import__('sys').argv[1] == 'hyb':
    __import__('sys').argv.pop(1)
    main_hyb()
