#!/usr/bin/env python3
"""Free-tilted Monte Carlo: search over the tilted squares as free pieces; axis squares re-solved exactly every step.

State: a squeezed packing C (n squares, side s).  T = its non-axis squares (> 0.5 deg from axis).
Step:  perturb T  ->  packer-relax T alone (looser box)  ->  exact axis fill (gen.mis_fill2) at the smallest scanned side
       where |T| + fill >= n  ->  penalty squeeze (packer relax --squeeze-pen)  ->  reject if a full line of k axis squares
       ->  Metropolis on the squeezed side.  The squeezed configuration becomes the new state.

  ftmc.py --n N --k K --init file.txt --seed s --time T --out prefix [--temp 0.003] [--mu0 1e3]
"""
import argparse, math, os, random, shutil, subprocess, tempfile, time
from gen import mis_fill2, write, wall_viol
from layout import full_lines

HERE = os.path.dirname(os.path.abspath(__file__))
PACKER = os.path.join(HERE, 'target/release/packer')
FTOL = 0.03   # allowed axis-tilted overlap in the fill (the squeeze resolves it); axis-axis stays exact


def load(path):
    L = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    n, s = int(L[0][0]), float(L[0][1])
    return s, [(float(a), float(b), float(c)) for a, b, c, *_ in L[1:n + 1]]


def tilted(C, eps=0.5):
    return [q for q in C if min(q[2] % 90, 90 - q[2] % 90) > eps]


def run_packer(args):
    r = subprocess.run([PACKER] + args, capture_output=True, text=True)
    out = r.stdout
    E = float(out.split('E=')[1].split()[0]) if 'E=' in out else float('inf')
    s = float(out.split(' s=')[1].split()[0]) if ' s=' in out else float('inf')
    return E, s


def perturb(T, s, rng):
    T = list(T)
    kind = rng.randrange(7)
    c = T[rng.randrange(len(T))] if T else (s / 2, s / 2, 27.0)
    if kind == 0:                                   # disk jiggle
        r = 0.8 + 2 * rng.random(); sp = 0.01 + 0.08 * rng.random(); sa = 2 * rng.random()
        T = [(x + rng.gauss(0, sp), y + rng.gauss(0, sp), a + rng.gauss(0, sa))
             if math.hypot(x - c[0], y - c[1]) < r else (x, y, a) for x, y, a in T]
    elif kind == 1:                                 # cluster rotation
        r = 1 + 2.5 * rng.random(); ph = math.radians(rng.gauss(0, 2))
        cp, sp = math.cos(ph), math.sin(ph)
        T = [(c[0] + cp * (x - c[0]) - sp * (y - c[1]), c[1] + sp * (x - c[0]) + cp * (y - c[1]), a + math.degrees(ph))
             if math.hypot(x - c[0], y - c[1]) < r else (x, y, a) for x, y, a in T]
    elif kind == 2 and T:                           # delete a tilted square
        T.remove(c)
    elif kind == 3:                                 # add a lattice neighbour of c
        t = math.radians(c[2]); u = (math.cos(t), math.sin(t)); v = (-u[1], u[0])
        d = rng.choice([u, v, (-u[0], -u[1]), (-v[0], -v[1])])
        T.append((c[0] + d[0], c[1] + d[1], c[2]))
    elif kind == 4:                                 # shift one angle group (channel) rigidly
        a0 = c[2] % 90; tx, ty = rng.gauss(0, 0.08), rng.gauss(0, 0.08)
        T = [(x + tx, y + ty, a) if abs(((a % 90) - a0 + 45) % 90 - 45) < 6 else (x, y, a) for x, y, a in T]
    elif kind == 5:                                 # re-angle one square
        T = [(x, y, a + rng.gauss(0, 3)) if (x, y, a) == c else (x, y, a) for x, y, a in T]
    else:                                           # delete one, add one elsewhere next to another
        if len(T) > 1:
            T.remove(c)
            d0 = T[rng.randrange(len(T))]
            t = math.radians(d0[2]); u = (math.cos(t), math.sin(t)); v = (-u[1], u[0])
            d = rng.choice([u, v, (-u[0], -u[1]), (-v[0], -v[1])])
            T.append((d0[0] + d[0], d0[1] + d[1], d0[2]))
    return [(min(max(x, 0.5), s - 0.5), min(max(y, 0.5), s - 0.5), a) for x, y, a in T], kind


def evaluate(T, n, k, s_cur, tmp, mu0):
    """Relax T, fill exactly, squeeze.  Returns (side, path) or (inf, None)."""
    a, b, c = (os.path.join(tmp, x) for x in ('t.txt', 'tr.txt', 'f.txt'))
    sb = s_cur + 0.3
    write(a, sb, T)
    E, _ = run_packer(['relax', '--in', a, '--s', repr(sb), '--out', b, '--maxit', '5000'])
    if E > 1e-20:
        return float('inf'), None
    _, T = load(b)
    s = s_cur - 0.03
    while s < s_cur + 0.6:
        if all(wall_viol(q, s) <= 1e-9 for q in T):
            ax = mis_fill2(s, T, FTOL)
            if len(T) + len(ax) >= n:
                write(c, s, (T + ax)[:n])
                _, sq = run_packer(['relax', '--in', c, '--squeeze-pen', '--mu0', repr(mu0), '--out', b])
                if sq == float('inf') or full_lines(b, k) > 0:
                    return float('inf'), None
                return sq, b
        s += 0.04
    return float('inf'), None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True); ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--init', required=True); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--time', type=float, default=600); ap.add_argument('--temp', type=float, default=0.003)
    ap.add_argument('--mu0', type=float, default=1e3); ap.add_argument('--out', required=True)
    ap.add_argument('--ftol', type=float, default=0.03)
    a = ap.parse_args()
    global FTOL
    FTOL = a.ftol
    rng = random.Random(a.seed)
    tmp = tempfile.mkdtemp()
    s_cur, C = load(a.init)
    best = s_cur
    shutil.copy(a.init, f'{a.out}.best.txt')
    t0 = time.time(); it = 0; acc = [0] * 7; tried = [0] * 7
    print(f'init s={s_cur:.6f} |T|={len(tilted(C))}', flush=True)
    while time.time() - t0 < a.time:
        it += 1
        T2, kind = perturb(tilted(C), s_cur, rng)
        tried[kind] += 1
        s2, p = evaluate(T2, a.n, a.k, s_cur, tmp, a.mu0)
        if p is None:
            continue
        if s2 <= s_cur or rng.random() < math.exp(-(s2 - s_cur) / a.temp):
            _, C = load(p); s_cur = s2; acc[kind] += 1
            if s2 < best - 1e-7:
                best = s2
                shutil.copy(p, f'{a.out}.best.txt')
                print(f'it={it} best s={best:.6f} |T|={len(tilted(C))} kind={kind} t={time.time() - t0:.0f}s', flush=True)
    print(f'done best={best:.6f} its={it} acc={acc} tried={tried}', flush=True)


if __name__ == '__main__':
    main()
