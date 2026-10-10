#!/usr/bin/env python3
"""Import third-party packings (data only) into our text format, check them, and pool them per n as search seeds.

  import_ext.py scan DIR... --tag TAG --out POOL        read every packing file under DIR (formats below), write
                                                        POOL/n<NNN>/<tag>_<side>.txt for each feasible one (dedup)
  import_ext.py best POOL --out DIR                     DIR/n-<n>.txt = smallest-side packing per n in POOL

Formats: Ellsworth text (`s: ...` + `Square i: x=, y=, deg=`, box centred at 0; SQUISH, ry-xu, Mishapolk `.cert.txt`);
Couzo text (`# n =`, `# s =`, then `x y theta(rad)`, lower-left origin); Witness/v2 yaml (jlevy register, ry-xu), gz ok;
our text (`n s` then `x y deg`).  Feasibility: f64 SAT, max penetration <= --tol (seeds only; nothing here certifies).
"""
import argparse, glob, gzip, math, os, re, sys
import numpy as np


def _open(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def side(x):
    from fractions import Fraction
    return float(Fraction(x)) if '/' in x else float(x)


def read_any(p):
    """-> (s, [(x, y, deg)], s_str) in lower-left coordinates, or None if the format is not recognised."""
    txt = _open(p).read()
    if 'Square 1:' in txt and re.search(r'^(?:Final )?s:\s', txt, re.M):
        s_str = re.search(r'^(?:Final )?s:\s*([0-9.eE+-]+)', txt, re.M).group(1); s = float(s_str)
        sq = [(float(a) + s / 2, float(b) + s / 2, float(c)) for a, b, c in
              re.findall(r'Square \d+:\s*x=([^,\s]+),\s*y=([^,\s]+),\s*deg=([^,\s]+)', txt)]
        return s, sq, s_str
    if txt.startswith('# n ='):
        s_str = re.search(r'# s =\s*(\S+)', txt).group(1)
        sq = [tuple(map(float, l.split()[:3])) for l in txt.splitlines() if l.strip() and not l.startswith('#')]
        return float(s_str), [(x, y, math.degrees(t)) for x, y, t in sq], s_str
    if 'witness:' in txt and 'representation' in txt:
        import yaml
        w = yaml.safe_load(txt)['witness']
        if w.get('representation') != 'center-angle' or w['coordinates'].get('origin') != 'lower-left':
            return None
        f = 1.0 if w['coordinates']['angle_unit'] == 'degrees' else 180 / math.pi
        sq = [(float(q['center'][0]), float(q['center'][1]), f * float(q['angle'])) for q in w['squares']]
        return side(str(w['side'])), sq, str(w['side'])
    first = txt.split('\n', 1)[0].split()
    if len(first) == 2 and first[0].isdigit():
        sq = [tuple(map(float, l.split()[:3])) for l in txt.splitlines()[1:] if l.strip()]
        if len(sq) == int(first[0]):
            return side(first[1]), sq, first[1]
    return None


def max_pen(s, sq):
    """max over walls and pairs of the SAT penetration depth (0 = feasible), f64."""
    P = np.array(sq, float); th = np.radians(P[:, 2]); c, si = np.cos(th), np.sin(th)
    off = np.array([[.5, .5], [-.5, .5], [-.5, -.5], [.5, -.5]])
    C = P[:, None, :2] + np.stack([off[:, 0] * c[:, None] - off[:, 1] * si[:, None],
                                   off[:, 0] * si[:, None] + off[:, 1] * c[:, None]], -1)
    worst = max(0.0, -C.min(), C.max() - s)
    ctr = P[:, :2]
    d = np.linalg.norm(ctr[:, None] - ctr[None], axis=-1)
    I, J = np.nonzero(np.triu(d < math.sqrt(2) + 1e-9, 1))
    for i, j in zip(I, J):
        pen = math.inf
        for k in (i, j):
            for ax in ((c[k], si[k]), (-si[k], c[k])):
                a, b = C[i] @ ax, C[j] @ ax
                pen = min(pen, min(a.max() - b.min(), b.max() - a.min()))
        worst = max(worst, pen)
    return worst


def scan(a):
    n_ok = n_bad = 0
    for root in a.dirs:
        files = [root] if os.path.isfile(root) else sorted(glob.glob(f'{root}/**/*', recursive=True))
        for p in files:
            if not os.path.isfile(p) or not re.search(r'\.(txt|yaml|yml|gz)$', p) or re.search(r'check|output|cert50', p):
                continue
            try:
                r = read_any(p)
            except Exception as e:
                print(f'  unreadable {p}: {e!r}'[:160]); continue
            if r is None or not r[1]:
                continue
            s, sq, s_str = r; n = len(sq)
            if a.ns and n not in a.ns:
                continue
            pen = max_pen(s, sq)
            if pen > a.tol:
                n_bad += 1; print(f'  REJECT n={n} s={s:.12f} pen={pen:.2e} {p}'); continue
            d = f'{a.out}/n{n:03d}'; os.makedirs(d, exist_ok=True)
            if any(abs(side(open(q).readline().split()[1]) - s) < 1e-11 and a.tag in os.path.basename(q)
                   for q in glob.glob(f'{d}/*.txt') if not q.endswith('SOURCES.txt')):
                continue
            q = f'{d}/{a.tag}_{s:.12f}.txt'
            with open(q, 'w') as f:
                f.write(f'{n} {s_str}\n' + ''.join(f'{x!r} {y!r} {t!r}\n' for x, y, t in sq))
            with open(f'{d}/SOURCES.txt', 'a') as f:
                f.write(f'{os.path.basename(q)}\t{p}\tpen={pen:.1e}\n')
            n_ok += 1
    print(f'{a.tag}: {n_ok} imported, {n_bad} rejected')


def best(a):
    os.makedirs(a.out, exist_ok=True)
    for d in sorted(glob.glob(f'{a.pool}/n[0-9]*')):
        fs = [(side(open(q).readline().split()[1]), q) for q in glob.glob(f'{d}/*.txt') if not q.endswith('SOURCES.txt')]
        if fs:
            s, q = min(fs); n = int(os.path.basename(d)[1:])
            os.system(f'cp {q} {a.out}/n-{n}.txt')
            print(f'{n:4d} {s:.12f} {os.path.basename(q)}  ({len(fs)} in pool)')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('scan'); p.add_argument('dirs', nargs='+'); p.add_argument('--tag', required=True)
    p.add_argument('--out', required=True); p.add_argument('--tol', type=float, default=1e-9)
    p.add_argument('--ns', type=lambda x: {int(v) for v in x.split(',')}, default=None)
    p = sp.add_parser('best'); p.add_argument('pool'); p.add_argument('--out', required=True)
    a = ap.parse_args()
    {'scan': scan, 'best': best}[a.cmd](a)
