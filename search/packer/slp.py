#!/usr/bin/env python3
"""Sequential-LP squeeze: repeat (shrink LP in a trust region -> apply -> packer relax at the new side) until the
LP finds no shrink.  Reports the side reached and the final first-order status.

  slp.py in.txt out.txt [--R 1e-3] [--iters 200]
"""
import math, subprocess, sys, os, tempfile
import numpy as np
from rigid import load, shrink_lp

PACKER = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'target/release/packer')


def save(path, s, sq):
    with open(path, 'w') as f:
        f.write(f'{len(sq)} {s!r}\n')
        for x, y, t in sq:
            f.write(f'{x!r} {y!r} {math.degrees(t) % 90!r}\n')


def relax(path, s, out):
    r = subprocess.run([PACKER, 'relax', '--in', path, '--s', repr(s), '--out', out], capture_output=True, text=True)
    if "E=" not in r.stdout:
        raise RuntimeError(f"packer failed: {r.stdout!r} {r.stderr[-500:]!r}")
    return float(r.stdout.split('E=')[1].split()[0])


def slp(s, sq, R=1e-3, iters=200, verbose=False):
    tmp = tempfile.mkdtemp()
    a, b = os.path.join(tmp, 'a.txt'), os.path.join(tmp, 'b.txt')
    for it in range(iters):
        ds, z, nc = shrink_lp(s, sq, delta=max(10 * R, 1e-7), R=R)
        if ds is None or ds > -1e-13:
            if R > 1e-9:
                R *= 0.3
                continue
            break
        n = len(sq)
        new = [(float(x + z[3 * i]), float(y + z[3 * i + 1]), float(t + z[3 * i + 2])) for i, (x, y, t) in enumerate(sq)]
        s_new = float(s + ds)
        save(a, s_new, new)
        e = relax(a, s_new, b)
        if e < 1e-26:
            s, sq = load(b)
            R = min(R * 1.5, 1e-2)
            if verbose:
                print(f'it={it} s={s:.12f} contacts={nc} ds={ds:.2e}')
        else:
            R *= 0.3
            if R < 1e-10:
                break
    ds0, _, nc = shrink_lp(s, sq, delta=1e-7, exact=True)
    return s, sq, ds0


if __name__ == '__main__':
    s, sq = load(sys.argv[1])
    R = float(sys.argv[sys.argv.index('--R') + 1]) if '--R' in sys.argv else 1e-3
    s2, sq2, ds0 = slp(s, sq, R, verbose='-v' in sys.argv)
    save(sys.argv[2], s2, sq2)
    print(f'{sys.argv[1]}: s {s:.12f} -> {s2:.12f}; final exact-contact ds* = {ds0:.2e}')
