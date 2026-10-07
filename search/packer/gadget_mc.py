#!/usr/bin/env python3
"""Gadget Monte Carlo: search over the tilted squares only; axis squares are solved exactly.

State T = tilted squares in [0,s]^2 (kept feasible among themselves by `packer relax`).
Score f(T) = |T| + MIS(T)  (exact max axis fill incl. staircase candidates, gen.mis_fill2)
         + bonus * (MIS at tol 0.03 - MIS at tol 0)   (near-fits: a tie-break toward configurations about to gain one)
Goal: f >= n at side s.  Then the full configuration (T + fill) is relaxed and squeezed to confirm.

  gadget_mc.py --n N --s S --init file.txt --seed k --iters I --out prefix [--temp 0.3]
init: any packer-format file; its non-axis squares (centres inside the box) form the initial T.
"""
import argparse, math, os, random, subprocess, tempfile, time
from gen import mis_fill2, is_axis, write, wall_viol

HERE = os.path.dirname(os.path.abspath(__file__))
PACKER = os.path.join(HERE, 'target/release/packer')


def load(path):
    L = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    n, s = int(L[0][0]), float(L[0][1])
    return s, [(float(a), float(b), float(c)) for a, b, c, *_ in L[1:n + 1]]


def relax(T, s, tmp, extra=()):
    """Relax T (degrees) at side s with packer; returns (E, T')."""
    a, b = os.path.join(tmp, 'a.txt'), os.path.join(tmp, 'b.txt')
    write(a, s, list(T) + list(extra))
    r = subprocess.run([PACKER, 'relax', '--in', a, '--s', repr(s), '--out', b, '--maxit', '5000'],
                       capture_output=True, text=True)
    E = float(r.stdout.split('E=')[1].split()[0])
    return E, load(b)[1]


def score(T, s, bonus):
    a0 = mis_fill2(s, T, 1e-6)
    f = len(T) + len(a0)
    if bonus > 0:
        a1 = mis_fill2(s, T, 0.03)
        f += bonus * max(0, len(a1) - len(a0))
    return f, a0


def propose(T, s, rng):
    T = list(T)
    kind = rng.randrange(6)
    if not T:
        kind = 4
    c = T[rng.randrange(len(T))] if T else (s / 2, s / 2, 30.0)
    if kind == 0:      # disk jiggle
        r = 1 + 2 * rng.random(); sp = 0.03 + 0.15 * rng.random(); sa = 3 * rng.random()
        T = [(x + rng.gauss(0, sp), y + rng.gauss(0, sp), a + rng.gauss(0, sa))
             if math.hypot(x - c[0], y - c[1]) < r else (x, y, a) for x, y, a in T]
    elif kind == 1:    # cluster rotation about c
        r = 1 + 3 * rng.random(); ph = math.radians(rng.gauss(0, 3))
        cp, sp = math.cos(ph), math.sin(ph)
        T = [(c[0] + cp * (x - c[0]) - sp * (y - c[1]), c[1] + sp * (x - c[0]) + cp * (y - c[1]), a + math.degrees(ph))
             if math.hypot(x - c[0], y - c[1]) < r else (x, y, a) for x, y, a in T]
    elif kind == 2:    # cluster translation
        r = 1 + 3 * rng.random(); tx, ty = rng.gauss(0, 0.15), rng.gauss(0, 0.15)
        T = [(x + tx, y + ty, a) if math.hypot(x - c[0], y - c[1]) < r else (x, y, a) for x, y, a in T]
    elif kind == 3:    # delete
        T.remove(c)
    elif kind == 4:    # add a lattice neighbour of c (same angle)
        t = math.radians(c[2]); u = (math.cos(t), math.sin(t)); v = (-u[1], u[0])
        d = rng.choice([u, v, (-u[0], -u[1]), (-v[0], -v[1])])
        T.append((c[0] + d[0], c[1] + d[1], c[2]))
    else:              # global angle drift of one angle class
        a0 = c[2] % 90; da = rng.gauss(0, 1.0)
        T = [(x, y, a + da) if abs((a % 90) - a0) < 3 else (x, y, a) for x, y, a in T]
    T = [(min(max(x, 0.5), s - 0.5), min(max(y, 0.5), s - 0.5), a) for x, y, a in T]
    return T, kind


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True); ap.add_argument('--s', type=float, required=True)
    ap.add_argument('--init', required=True); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--iters', type=int, default=100000); ap.add_argument('--temp', type=float, default=0.3)
    ap.add_argument('--bonus', type=float, default=0.25); ap.add_argument('--out', required=True)
    ap.add_argument('--time', type=float, default=1e9)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    tmp = tempfile.mkdtemp()
    s0, sq = load(a.init)
    T = [q for q in sq if not is_axis(q, 1e-3) and 0.5 <= q[0] <= a.s - 0.5 and 0.5 <= q[1] <= a.s - 0.5]
    _, T = relax(T, a.s, tmp)
    f, ax = score(T, a.s, a.bonus)
    best = f
    t0 = time.time()
    print(f'init |T|={len(T)} f={f:.2f} need {a.n}', flush=True)
    for it in range(a.iters):
        if time.time() - t0 > a.time:
            break
        T2, kind = propose(T, a.s, rng)
        E, T2 = relax(T2, a.s, tmp)
        if E > 1e-20:
            continue
        f2, ax2 = score(T2, a.s, a.bonus)
        if f2 >= f or rng.random() < math.exp((f2 - f) / a.temp):
            T, f, ax = T2, f2, ax2
        if f > best:
            best = f
            print(f'it={it} |T|={len(T)} f={f:.2f} t={time.time() - t0:.0f}s', flush=True)
            write(f'{a.out}.best.txt', a.s, T + ax)
        if len(T) + len(ax) >= a.n:
            full = (T + ax)[:a.n]
            E, full2 = relax(full, a.s, tmp)
            write(f'{a.out}.hit.txt', a.s, full2)
            print(f'HIT it={it} count={len(T) + len(ax)} full relax E={E:.2e}', flush=True)
            if E < 1e-26:
                r = subprocess.run([PACKER, 'relax', '--in', f'{a.out}.hit.txt', '--squeeze', '--out', f'{a.out}.sq.txt'],
                                   capture_output=True, text=True)
                print('SQUEEZE', r.stdout.strip(), flush=True)
                break
    print(f'done best f={best:.2f} t={time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main()
