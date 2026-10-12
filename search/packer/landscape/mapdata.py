#!/usr/bin/env python3
"""Connectivity-map data for the s(110) explainer (10-11): certified minima below 11 from the store (`minimum` table;
weak / unresolved left out until re-certified), junction-aware role counts, all-pairs fingerprint distance.

  python3 landscape/mapdata.py [--out runs/landscape110/map]      (run from packer/; ~15 s)

Writes minima.json (key, S, family, display orientation g = (transpose, flip y, flip x) applied to the stored
packing, packing in 1e-4 / 1e-2 deg ints) and fpD.npy (n x n, pk.fingerprint distance)."""
import argparse, collections, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import numpy as np
from pk.store import Store, RANK
from pk import fingerprint as fp
from explore import roles_j


SYMS = [(t, fy, fx) for t in (0, 1) for fy in (0, 1) for fx in (0, 1)]


def sym(sq, s, t, fy, fx):
    """Apply a box symmetry to (x, y, theta deg): transpose, then flip y, then flip x."""
    sq = np.asarray(sq, float); x, y, th = sq[:, 0].copy(), sq[:, 1].copy(), sq[:, 2].copy()
    if t: x, y, th = y, x, 90 - th
    if fy: y, th = s - y, -th
    if fx: x, th = s - x, -th
    return np.c_[x, y, th % 90]


def orient_tree(M, D, root):
    """Display orientation: root as stored; then along the minimum spanning tree of D (Prim from root), each minimum gets
    the box symmetry whose fingerprint image is closest to its already-oriented tree parent.  Neighbours agree; tilt
    counts play no part (so squares flipping grid <-> band cannot flip the orientation).  Returns list of syms."""
    n = len(M); g = [None] * n; g[root] = (0, 0, 0)
    ims = [[fp.image(m['s'], sym(m['sq'], m['s'], *e)) for e in SYMS] for m in M]
    done = np.zeros(n, bool); done[root] = True
    best = D[root].copy(); par = np.full(n, root)
    for _ in range(n - 1):
        j = int(np.argmin(np.where(done, np.inf, best)))
        P = ims[par[j]][SYMS.index(g[par[j]])]
        g[j] = SYMS[int(np.argmin([np.abs(im - P).sum() for im in ims[j]]))]
        done[j] = True
        upd = (~done) & (D[j] < best); best[upd] = D[j][upd]; par[upd] = j
    return g


def family(sq):
    L, B, A, J, O = roles_j(sq)
    return (max(L, B), min(L, B), A, J, O)          # bands swap under the diagonal reflection


def load(n=110, k=11, cls=('certified',)):
    st = Store()
    best = {}
    for pid, s, status, pen in st.db.execute('SELECT id, s, status, pen FROM packing WHERE n=? AND s<?', (n, k)):
        key = f'{s:.9f}'; c = (-RANK.get(status, 0), pen or 0.0, pid)
        if key not in best or c < best[key]:
            best[key] = c
    cp = collections.defaultdict(list)
    for sk, mk in st.db.execute('SELECT side_key, min_key FROM copy WHERE n=?', (n,)):
        cp[mk].append(sk)
    out = []
    q = f"SELECT key, S, cls, pid FROM minimum WHERE n=? AND cls IN ({','.join('?' * len(cls))})"
    for r in st.db.execute(q, (n, *cls)):
        if float(r['S']) >= k:
            continue
        pid = r['pid'] or min((best[sk] for sk in cp[r['key']] if sk in best), default=(0, 0, None))[2]
        if pid is None:
            print('no packing for', r['key'], file=sys.stderr); continue
        s, coords = st.db.execute('SELECT s, coords FROM packing WHERE id=?', (pid,)).fetchone()
        sq = np.frombuffer(coords, dtype=np.float64).reshape(-1, 3)
        out.append(dict(key=r['key'], S=r['S'], cls=r['cls'], pid=pid, s=s, sq=sq))
    out.sort(key=lambda m: float(m['S']))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=os.path.join(os.path.dirname(HERE), 'runs/landscape110/map'))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    M = load()
    n = len(M); t = time.time()
    I = np.stack([fp.image(m['s'], m['sq']) for m in M]).reshape(n, -1)
    SY = np.stack([fp.syms(fp.image(m['s'], m['sq'])).reshape(8, -1) for m in M])
    D = np.empty((n, n))
    for i in range(n):
        D[i] = np.abs(SY[None, i] - I[:, None]).sum(-1).min(-1) / 110
    D = np.minimum(D, D.T); np.fill_diagonal(D, 0)
    np.save(f'{a.out}/fpD.npy', D)
    G = orient_tree(M, D, int(np.argmin([float(m['S']) for m in M])))
    for m, e in zip(M, G):
        m['sq'] = sym(m['sq'], m['s'], *e); m['g'] = e
    json.dump([dict(key=m['key'], S=m['S'][:16], pid=m['pid'], fam=family(m['sq']), g=m['g'],
                    sq=[v for x, y, th in m['sq'] for v in (round(x * 1e4), round(y * 1e4), round(th * 100))])
               for m in M], open(f'{a.out}/minima.json', 'w'), separators=(',', ':'))
    print(f'{n} certified minima, distances in {time.time() - t:.0f} s -> {a.out}')


if __name__ == '__main__':
    main()
