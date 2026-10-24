#!/usr/bin/env python3
"""Kaggle submission CSV -> the store (10-23; Evan: "import and use _all_ known packings, regardless of source ... use the
store db as a db": holdouts are select statements).  CSV rows (id, n, square_id, x, y, angle[rad]) -> per-n files in our
frame (box [0, s], degrees), every n (grid ones too: the store keeps them), then
`pk.py import --source kaggle:<name>`.

  kaggle_import.py runs/kaggle/pub/lightningv08_v12_2026-10-22.csv [--name lightningv08_v12] [--finder F]
"""
import argparse, csv, math, os, subprocess, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))


def half_ext(a_deg):
    """half the axis-aligned extent of a unit square at a_deg degrees"""
    r = np.radians(a_deg)
    return 0.5 * (np.abs(np.cos(r)) + np.abs(np.sin(r)))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('csv'); ap.add_argument('--name'); ap.add_argument('--finder')
    a = ap.parse_args()
    name = a.name or os.path.basename(a.csv).rsplit('.', 1)[0]
    out = os.path.join(os.path.dirname(os.path.abspath(a.csv)), name)
    os.makedirs(out, exist_ok=True)
    by = {}
    with open(a.csv) as fh:
        for r in csv.DictReader(fh):
            by.setdefault(int(r['n']), []).append((float(r['x']), float(r['y']), float(r['angle'] or 0.0)))
    files = []
    for n, rows in sorted(by.items()):
        P = np.array(rows)
        ang = np.degrees(P[:, 2]) % 90.0
        h = half_ext(ang)
        x, y = P[:, 0] - (P[:, 0] - h).min(), P[:, 1] - (P[:, 1] - h).min()
        s = max((x + h).max(), (y + h).max())
        f = os.path.join(out, f'n{n}.txt')
        with open(f, 'w') as fh:
            fh.write(f'{n} {float(s)!r}\n' + ''.join(f'{float(u)!r} {float(v)!r} {float(t)!r}\n' for u, v, t in zip(x, y, ang)))
        files.append(f)
    print(f'{len(files)} n written to {out}', flush=True)
    cmd = [sys.executable, os.path.join(HERE, 'pk.py'), 'import', *files, '--source', f'kaggle:{name}', '--ref', os.path.abspath(a.csv)]
    if a.finder:
        cmd += ['--finder', a.finder]
    subprocess.run(cmd, check=True)


if __name__ == '__main__':
    main()
