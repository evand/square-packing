"""Parent / start selection from the store (10-14): one pool query and named policies, shared by every driver.

  pool = Pool(n, ks=(0,), asof=None, above_k=0.05)        rows at n (and n + k for k in ks), feasible, s < k(n) + above_k
  pick(pool, m, policy='cluster', rng=random.Random(0))    -> [store ids]
  starts(n, m, out_dir, ...)                               -> files out_dir/lib_<i>.txt for legacy drivers (chain)

Policies
  best      lowest side first (the old default everywhere)
  spread    rank-uniform over the side range (pk/corpus 'spread')
  banded    a gap decade (log10 gap to the frontier) uniformly, then `cluster` inside it: spread in depth and in region
            (C33: deep parents pay most, far ones some; at most n clusters are near-singletons, so `cluster` alone
            would draw mostly far stepping stones).  Chain default (10-14).
  cluster   uniform over region clusters, then uniform over minima (side keys) inside, then the least-expanded copy
            (Evan 10-14: diversity = spread across distinct regions); the default.  Minima weighted down by ledger CPU
            already spent on them (1 / (1 + cpu / CPU0)).
Clusters: pk/cluster.py (leader balls in matching distance, r = 0.42 / sqrt n; analysis/similarity110.md), computed on
the pool itself (so an as-of pool never sees later packings), cached in runs/clusters/.  cluster=False: a coarse
provisional key (tilted-count bin, log10 gap bin).
As-of: rows whose `added` (store-local time) is before the cutoff only (time travel; pk/eventbench cutoff names or UTC ISO).
"""
from __future__ import annotations
import math, os, random, sqlite3
import numpy as np
from collections import defaultdict
from .store import Store, DEFAULT

CPU0 = 600.0                 # ledger weight scale: a packing that already got CPU0 s of expansion counts half


def asof_local(T):
    from . import eventbench as EB
    return EB.cutoff(T).astimezone().strftime('%Y-%m-%dT%H:%M:%S')


def kgrid(n):
    return math.ceil(math.sqrt(n) - 1e-12)


def gap_bin(gap):
    """log10 bins of the gap above the frontier: '<=0', '-9' .. '-1', '>=0' (gap in side units)"""
    if gap <= 0:
        return '<=0'
    return str(max(-9, min(0, math.floor(math.log10(gap)))))


def cap_pool(rows, cap, seed=0):
    """At most `cap` distinct side keys (rows of one key kept together), sampled evenly over gap decades (gap_bin):
    pk/cluster.py is O(keys x leaders); 110 alone has ~150k rows.  Deterministic per seed, so cached labels stay valid."""
    keys = defaultdict(list)
    for d in rows:
        keys[round(d['s'] / 2e-9)].append(d)
    if len(keys) <= cap:
        return rows
    band = defaultdict(list)
    for k, ds in keys.items():
        band[gap_bin(ds[0]['gap'])].append(k)
    rng = random.Random(seed)
    for b in band.values():
        b.sort(); rng.shuffle(b)
    keep, bs = [], sorted(band, key=lambda b: len(band[b]))
    left = cap
    for i, b in enumerate(bs):                       # equal share per band; small bands give their unused share on
        take = band[b][:left // (len(bs) - i)]
        keep += take; left -= len(take)
    return [d for k in keep for d in keys[k]]


class Pool:
    def __init__(self, n, ks=(0,), asof=None, above_k=0.05, frontier=None, statuses=None, path=DEFAULT, ledger=True,
                 cluster=True, cap=6000):
        self.n, self.ks, self.cluster, self.cap = n, tuple(ks), cluster, cap
        db = sqlite3.connect(f'file:{path}?mode=ro', uri=True)
        db.row_factory = sqlite3.Row
        self.rows = []
        for k in self.ks:
            m = n + k
            if frontier is None:
                from . import eventbench as EB
                fr = EB.frontier_at(m, asof) if asof else EB.frontier_now(m)
            else:
                fr = frontier(m)
            q = ("SELECT id, n, s, ntilted, grid, status, cluster, coords FROM packing WHERE n=? AND s<? AND status != 'infeasible'"
                 " AND (grid IS NULL OR grid=0)")
            a = [m, kgrid(m) + above_k]
            if asof:
                q += ' AND added < ?'; a.append(asof_local(asof))
            if statuses:
                q += f' AND status IN ({",".join("?" * len(statuses))})'; a += list(statuses)
            got = []
            for r in db.execute(q, a):
                d = dict(r); d['k'] = k; d['gap'] = d['s'] - fr
                d['sq'] = np.frombuffer(d.pop('coords'), '<f8').reshape(-1, 3) if self.cluster else None
                got.append(d)
            got = cap_pool(got, self.cap, seed=m)
            if self.cluster and got:                   # region clusters (pk/cluster.py), computed on this pool only
                from . import cluster as CL
                tag = f"asof{asof_local(asof).replace(':', '')}" if asof else 'live'
                L = CL.labels(m, [(d['id'], d['s'], d['sq']) for d in got], tag=tag)
                for d in got:
                    d['cl'] = L[d['id']]; d['sq'] = None
            for d in got:
                d.setdefault('cl', d['cluster'] or f"p:{m}:t{(d['ntilted'] or 0) // 4}:g{gap_bin(d['gap'])}")
                d['key'] = f"{m}:{d['s']:.9f}"
                self.rows.append(d)
        self.cpu = defaultdict(float)
        if ledger and self.rows:
            q = 'SELECT packing, sum(cpu_s) c FROM expansion WHERE n IN (%s)' % ','.join('?' * len(self.ks))
            a = [n + k for k in self.ks]
            if asof:
                q += ' AND added < ?'; a.append(asof_local(asof))
            for r in db.execute(q + ' GROUP BY packing', a):
                self.cpu[r['packing']] = r['c']
        db.close()

    def clusters(self):
        C = defaultdict(list)
        for r in self.rows:
            C[r['cl']].append(r)
        return C

    def summary(self):
        C = self.clusters()
        return dict(n=self.n, rows=len(self.rows), clusters=len(C), largest=max((len(v) for v in C.values()), default=0),
                    expanded=sum(1 for r in self.rows if self.cpu.get(r['id'])))


def _w(pool, r):
    return 1.0 / (1.0 + pool.cpu.get(r['id'], 0.0) / CPU0)


def pick(pool: Pool, m, policy='cluster', rng=None, exclude=()):
    """m distinct store ids by `policy` (see module doc)."""
    rng = rng or random.Random(0)
    rows = [r for r in pool.rows if r['id'] not in set(exclude)]
    if not rows:
        return []
    if policy == 'best':
        return [r['id'] for r in sorted(rows, key=lambda r: r['s'])[:m]]
    if policy == 'spread':
        rows.sort(key=lambda r: r['s'])
        if len(rows) <= m:
            return [r['id'] for r in rows]
        return [rows[round(i * (len(rows) - 1) / (m - 1 if m > 1 else 1))]['id'] for i in range(m)]
    if policy == 'banded':                           # gap decade uniform, then cluster, then minimum (default 10-14)
        B = defaultdict(list)
        for r in rows:
            B[gap_bin(r['gap'])].append(r)
        out, bands = [], list(B)
        while len(out) < m and bands:
            b = rng.choice(bands)
            sub = Pool.__new__(Pool); sub.rows = B[b]; sub.cpu = pool.cpu
            got = pick(sub, 1, 'cluster', rng, exclude=out)
            if not got:
                bands.remove(b); continue
            out += got
        return out
    if policy == 'cluster':
        C = defaultdict(lambda: defaultdict(list))   # cluster -> side key (one minimum) -> rows
        for r in rows:
            C[r['cl']][r['key']].append(r)
        def kw(rs):                                  # a minimum's weight: down by ledger CPU on any of its copies
            return 1.0 / (1.0 + sum(pool.cpu.get(x['id'], 0.0) for x in rs) / CPU0)
        cl = list(C)
        rng.shuffle(cl)
        out = []
        while len(out) < m and cl:                   # uniform over clusters (without replacement), then over minima
            for c in cl[:m - len(out)]:
                ks = list(C[c])
                key = rng.choices(ks, weights=[kw(C[c][x]) for x in ks])[0]
                rs = C[c].pop(key)
                out.append(min(rs, key=lambda x: (pool.cpu.get(x['id'], 0.0), x['id']))['id'])
            cl = [c for c in cl if C[c]]
            rng.shuffle(cl)
        return out
    raise ValueError(f'unknown policy {policy}')


def starts(n, m, out_dir, policy='cluster', ks=(0,), asof=None, seed=0, exclude=(), prefix='lib'):
    """Write m picked packings to out_dir/<prefix>_<i>.txt (ours format); returns [(path, id, s, cluster)]."""
    pool = Pool(n, ks=ks, asof=asof)
    ids = pick(pool, m, policy, random.Random(seed), exclude)
    st = Store()
    by = {r['id']: r for r in pool.rows}
    out = []
    try:
        for i, pid in enumerate(ids):
            p = st.get(pid)
            path = os.path.join(out_dir, f'{prefix}_{i}.txt')
            p.write(path)
            out.append((path, pid, by[pid]['s'], by[pid]['cl']))
    finally:
        st.close()
    return out


if __name__ == '__main__':
    import argparse, json
    ap = argparse.ArgumentParser(description='parent pools / picks from the store')
    ap.add_argument('n', type=int); ap.add_argument('--m', type=int, default=8); ap.add_argument('--policy', default='cluster')
    ap.add_argument('--ks', default='0'); ap.add_argument('--asof'); ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()
    P = Pool(a.n, ks=[int(x) for x in a.ks.split(',')], asof=a.asof)
    print(json.dumps(P.summary()))
    by = {r['id']: r for r in P.rows}
    for pid in pick(P, a.m, a.policy, random.Random(a.seed)):
        r = by[pid]
        print(f"{pid:8d} n={r['n']} s={r['s']:.10f} gap={r['gap']:+.2e} {r['cl']} cpu={P.cpu.get(pid, 0):.0f}")
