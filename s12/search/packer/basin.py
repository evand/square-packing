#!/usr/bin/env python3
"""Basin-entry measurement: perturb a packing by Gaussian noise sigma (positions; angles sigma*20 deg), squeeze with
packer --squeeze-pen, and record whether it returns to the reference side (within tol).  Prints return fraction per sigma
and the distribution of sides reached.

  basin.py ref.txt --sigmas 1e-3,3e-3,1e-2,3e-2,0.1 --trials 20 [--mu0 1e3] [--loosen 0.02]
"""
import argparse, os, random, subprocess, tempfile, statistics
from gen import write

PACKER = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'target/release/packer')


def load(path):
    L = [l.split() for l in open(path) if l.strip()]
    n, s = int(L[0][0]), float(L[0][1])
    return s, [(float(a), float(b), float(c)) for a, b, c, *_ in L[1:n + 1]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('ref'); ap.add_argument('--sigmas', default='1e-3,3e-3,1e-2,3e-2,0.1')
    ap.add_argument('--trials', type=int, default=20); ap.add_argument('--mu0', type=float, default=1e3)
    ap.add_argument('--loosen', type=float, default=0.02); ap.add_argument('--k', type=float, default=0); ap.add_argument('--tol', type=float, default=1e-5)
    a = ap.parse_args()
    s0, C = load(a.ref)
    tmp = tempfile.mkdtemp(); rng = random.Random(1)
    # reference side under this descent map (squeeze the unperturbed loosened configuration)
    for sig in [0.0] + [float(x) for x in a.sigmas.split(',')]:
        sides = []
        for t in range(a.trials if sig > 0 else 1):
            f = 1 + a.loosen
            P = [(x * f + rng.gauss(0, sig), y * f + rng.gauss(0, sig), ang + rng.gauss(0, 20 * sig)) for x, y, ang in C]
            p = os.path.join(tmp, 'p.txt'); q = os.path.join(tmp, 'q.txt')
            write(p, s0 * f + 4 * sig, P)
            r = subprocess.run([PACKER, 'relax', '--in', p, '--squeeze-pen', '--mu0', repr(a.mu0), '--out', q],
                               capture_output=True, text=True).stdout
            sides.append(float(r.split(' s=')[1].split()[0]) if ' s=' in r else float('inf'))
        if sig == 0.0:
            sref = sides[0]
            print(f'ref {a.ref}: file s={s0:.9f}, loosened+squeezed s={sref:.9f}')
            continue
        ret = sum(abs(x - min(sref, s0)) < a.tol or x < min(sref, s0) for x in sides)
        from collections import Counter
        keys = Counter(round(x, 8) for x in sides if x < float('inf'))
        f1 = sum(1 for v in keys.values() if v == 1); f2 = sum(1 for v in keys.values() if v == 2)
        chao = len(keys) + (f1 * f1 / (2 * f2) if f2 else f1 * (f1 - 1) / 2)
        below = sum(x < a.k for x in sides) if a.k else 0
        fin = sorted(x for x in sides if x < float('inf'))
        print(f'sigma={sig:g}: distinct minima {len(keys)}/{len(sides)} (seen once {f1}, twice {f2}, Chao1 ~{chao:.0f}); below k {below}; '
              f'returned {ret}; sides min {fin[0]:.6f} median {statistics.median(fin):.6f} max {fin[-1]:.6f}',
              flush=True)


if __name__ == '__main__':
    main()
