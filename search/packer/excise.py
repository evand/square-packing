#!/usr/bin/env python3
"""Strip excision: carry a packing family one triangle row down (10-14 night; Evan's family thesis).

Friedman's L-step gives s(n + 2k - 1) <= s(n) + 1 for floor(s(n)) = k - 1.  Run backwards: take best(n') with n' = n + 2k
- 1 (same column m = n - (k-1)^2, one row up), delete every square meeting a vertical unit strip x in (x0, x0 + 1) and a
horizontal one y in (y0, y0 + 1), shift what lies beyond each strip back by 1: side s(n') - 1.  If the cut removes at
most 2k - 1 squares, the excised layout holds >= n squares in side s(n') - 1 (seams may overlap a little: squares that
touched the strip edges now touch each other, tilted ones can overlap).  Extras are dropped (fewest neighbours first).
The seams are then healed by `harden` (box held, release, settle) + fq polish (harden.py).

  python3 excise.py scan N'|ID --n N [--step 0.05]          best cuts: removed count vs the 2k - 1 allowance
  python3 excise.py build N'|ID --n N --out DIR [--top 8]   write the top cuts as harden inputs (DIR/cut*.txt)
"""
from __future__ import annotations
import argparse, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def half_extent(t_deg):
    t = math.radians(t_deg % 90.0)
    return 0.5 * (abs(math.cos(t)) + abs(math.sin(t)))


def cut(sq, s, x0, y0):
    """Squares kept (shifted) after excising (x0, x0+1) and (y0, y0+1); a square is removed when its axis-aligned
    extent meets an open strip (conservative for tilted squares)."""
    keep = []
    for x, y, t in sq:
        w = half_extent(t)
        if x + w > x0 + 1e-9 and x - w < x0 + 1 - 1e-9:
            continue
        if y + w > y0 + 1e-9 and y - w < y0 + 1 - 1e-9:
            continue
        keep.append((x - (1.0 if x - w >= x0 + 1 - 1e-9 else 0.0), y - (1.0 if y - w >= y0 + 1 - 1e-9 else 0.0), t))
    return np.array(keep)


def load(ref):
    from pk.store import Store
    st = Store()
    if ref.isdigit() and int(ref) > 400:
        p = st.get(int(ref))
    else:
        r = st.db.execute("SELECT id FROM packing WHERE n=? AND status!='infeasible' AND (pen IS NULL OR pen<=1e-9) "
                          "ORDER BY s LIMIT 1", (int(ref),)).fetchone()
        p = st.get(r['id'])
    return p.n, float(p.s), np.asarray(p.sq, float)


def scan(n1, s1, sq, n, step):
    res = []
    xs = np.arange(0.0, s1 - 1 + 1e-9, step)
    for x0 in xs:
        for y0 in xs:
            k = cut(sq, s1, x0, y0)
            res.append((len(k), x0, y0))
    res.sort(key=lambda r: -r[0])
    return res


def drop_extras(K, n):
    """Drop the squares with the fewest close neighbours until n remain."""
    while len(K) > n:
        d = np.linalg.norm(K[:, None, :2] - K[None, :, :2], axis=2)
        deg = (d < 1.2).sum(1)
        K = np.delete(K, int(np.argmin(deg)), axis=0)
    return K


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    for c in ('scan', 'build'):
        p = sp.add_parser(c); p.add_argument('src'); p.add_argument('--n', type=int, required=True)
        p.add_argument('--step', type=float, default=0.05)
        if c == 'build':
            p.add_argument('--out', required=True); p.add_argument('--top', type=int, default=8)
    a = ap.parse_args()
    n1, s1, sq = load(a.src)
    allow = n1 - a.n
    res = scan(n1, s1, sq, a.n, a.step)
    best = res[0][0]
    print(f"n' {n1} s' {s1:.10f} -> n {a.n} at side {s1 - 1:.10f}; allowance {allow} squares; best cut keeps {best} "
          f"(need {a.n}: {'OK' if best >= a.n else f'short {a.n - best}'}); cuts keeping >= n: {sum(r[0] >= a.n for r in res)}")
    if a.cmd == 'build':
        os.makedirs(a.out, exist_ok=True)
        seen, k = set(), 0
        for c, x0, y0 in res:
            if k >= a.top:
                break
            key = (round(x0, 1), round(y0, 1))
            if key in seen:
                continue
            seen.add(key)
            K = cut(sq, s1, x0, y0)
            if len(K) < a.n:
                break
            K = drop_extras(K, a.n)
            with open(os.path.join(a.out, f'cut{k}.txt'), 'w') as f:
                f.write(f'{a.n} {s1 - 1!r}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(t)!r}\n' for x, y, t in K))
            print(f'  cut{k}: x0 {x0:.2f} y0 {y0:.2f} kept {c}')
            k += 1


if __name__ == '__main__':
    main()
