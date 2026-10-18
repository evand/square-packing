"""Data for the k^2 - c explorer page (site/www/k2c/explore.html).

* site/www/k2c/data/fang-pieces.json: Fang's 18 gluing pieces (certificates/fang-gluing/inputs/*.cert, exact
  rationals) as display floats: {name: {n, W, H, sq: [[x, y, t], ...]}}, centre (x, y), t = tan(theta/2).
* site/www/k2c/data/best-known.json: per k, the largest n <= k^2 with a known packing in side < k, from the packing
  store's frontier table (register and issue-reported entries), with finder and source.  Shown on the page only as
  a reference mark: these packings are not drawn there.

  python3 search/ryu/explorer_data.py
"""
import json, math, os, sqlite3, sys
from fractions import Fraction as F

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, 'site/www/k2c/data')
CERTS = os.path.join(ROOT, 'certificates/fang-gluing/inputs')
STORE = os.environ.get('PK_STORE', os.path.join(ROOT, 'search/packer/runs/store.sqlite'))


def rd(x, nd=9):
    return float(round(float(x), nd))


def pieces():
    out = {}
    for fn in sorted(os.listdir(CERTS)):
        if not fn.endswith('.cert'):
            continue
        rows = [l.split() for l in open(os.path.join(CERTS, fn)) if l.strip() and not l.startswith('#')]
        head, body = rows[0], rows[1:]
        n = int(head[0])
        W = F(head[1])
        H = F(head[2]) if len(head) > 2 else W
        assert len(body) == n, fn
        sq = [[rd(F(x)), rd(F(y)), rd(F(t))] for x, y, t in body]
        out[fn[:-5]] = dict(n=n, W=rd(W, 12), H=rd(H, 12), sq=sq)
    return out


def best_known(kmax=100):
    db = sqlite3.connect(f'file:{STORE}?mode=ro', uri=True)
    best = {}
    for n, s, src, pend, finder, ref in db.execute('select n, s, source, pending, finder, ref from frontier'):
        if n not in best or s < best[n][0]:
            best[n] = (s, src, pend, finder, ref)
    out = {}
    for k in range(2, kmax + 1):
        ns = [n for n in best if (k - 1) ** 2 < n <= k * k and best[n][0] < k - 1e-9]
        if ns:
            n = max(ns)
            s, src, pend, finder, ref = best[n]
            out[k] = dict(n=n, c=k * k - n, s=s, finder=finder, source=src, registered=src == 'register')
    return out


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    P = pieces()
    with open(os.path.join(OUT, 'fang-pieces.json'), 'w') as f:
        json.dump(P, f, separators=(',', ':'))
    B = best_known()
    with open(os.path.join(OUT, 'best-known.json'), 'w') as f:
        json.dump(B, f, indent=0, ensure_ascii=False)
    print(f'{len(P)} pieces, {sum(p["n"] for p in P.values())} squares; best known for k = {min(B)}..{max(B)}')
