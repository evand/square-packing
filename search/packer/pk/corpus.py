"""Replay bench: frozen proposal corpora, replayed under pipeline variants; certification shared through the store's basin
table; paired metrics in minutes.

  corpus.py build NAME --n 110 --parents 'query' --arms kick lkick ... --per 4 [--seed 1]
      parents: store query (e.g. "n=110 below=11 status=certified limit=40" or ids "ids=12,15,...");
      writes runs/corpus/NAME/{items.jsonl, props.npz}: item = (parent id, arm, move params, seed), proposal coords frozen.
  corpus.py replay NAME VARIANT [--procs 16] [--known ref|none]    -> runs/corpus/NAME/out_VARIANT.jsonl (+ .npz finals)
  corpus.py certify NAME [VARIANT...] [--procs 16]   exactsolve every distinct new? side below k (cache: store basin table)
  corpus.py report NAME [VARIANT...] [--by arm|param:sigma|loosen]

Variants are named Policy settings (VARIANTS below; add more freely).  The known-basin set for a corpus is frozen at build
time (`known.json`: certified sides at that n in the store, minus nothing), so every variant sees the same reference.
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):    # before numpy: one BLAS thread per worker
    os.environ.setdefault(_v, '1')
import argparse, collections, json, math, os, random, statistics as stt, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pk.store import Store
from pk.packing import Packing
from pk import moves
from pk.engine import QOpt, SCREEN, POLISH
from pk.pipeline import Policy, Known, evaluate

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs', 'corpus')

VARIANTS = {
    'base': Policy(),                                             # explore.py 10-08 recipe
    'staged20': Policy(stage_pit=20),                             # chunked polish with late-return
    'staged50': Policy(stage_pit=50),
    'same1e5': Policy(same=1e-5),                                 # looser screen-return test
    'nolines': Policy(grid_check=False),
    'loosen1': Policy(loosen=(1.0,)), 'loosen102': Policy(loosen=(1.02,)), 'loosen105': Policy(loosen=(1.05,)),
    'loosen110': Policy(loosen=(1.10,)),
    'mu1': Policy(screen=SCREEN.but(mu0=1.0)), 'mu100': Policy(screen=SCREEN.but(mu0=100.0)),
    'pit0': Policy(screen=SCREEN.but(pit=0)), 'pit30': Policy(screen=SCREEN.but(pit=30)),
    'full': Policy(same=0.0),                                     # polish everything (reference for recall)
}


def parse_query(q):
    kw = {}
    for tok in q.split():
        k, v = tok.split('=', 1)
        kw[k] = v
    return kw


def select_parents(st, n, q):
    kw = parse_query(q)
    if 'ids' in kw:
        return [int(x) for x in kw['ids'].split(',')]
    rows = st.query(n=n, below=float(kw['below']) if 'below' in kw else None,
                    above=float(kw['above']) if 'above' in kw else None, status=kw.get('status'),
                    source=kw.get('source'), tag=kw.get('tag'), nongrid=kw.get('nongrid') == '1', limit=None)
    rows = list(rows)
    lim = int(kw.get('limit', 40))
    if kw.get('pick', 'spread') == 'best':
        rows = rows[:lim]
    else:                                     # spread over the side range (rank-uniform)
        idx = np.linspace(0, len(rows) - 1, min(lim, len(rows))).round().astype(int)
        rows = [rows[i] for i in sorted(set(idx.tolist()))]
    return [r['id'] for r in rows]


def build(a):
    st = Store()
    d = os.path.join(ROOT, a.name)
    os.makedirs(d, exist_ok=True)
    pids = select_parents(st, a.n, a.parents)
    rng = random.Random(a.seed)
    items, props = [], []
    for pid in pids:
        p = st.get(pid)
        for arm in a.arms:
            for j in range(a.per):
                seed = rng.randrange(1 << 30)
                q, info = moves.apply(p, arm, seed)
                items.append(dict(i=len(items), parent=pid, parent_s=p.s, arm=arm, move=info['move'], params=info['params'],
                                  seed=seed, s0=q.s, moved=info.get('moved')))
                props.append(q.sq)
    np.savez_compressed(f'{d}/props.npz', *props)
    with open(f'{d}/items.jsonl', 'w') as f:
        for it in items:
            f.write(json.dumps(it) + '\n')
    known = sorted({float(r['key']) for r in st.db.execute("SELECT key FROM basin WHERE n=? AND cls='certified'", (a.n,))})
    k = math.ceil(math.sqrt(a.n) - 1e-12)
    best = st.frontier().get(a.n)
    json.dump(dict(n=a.n, k=k, known=known, best=min(([best['s']] if best else []) + known[:1]), parents=pids,
                   arms=a.arms, per=a.per, built=time.strftime('%F %T')), open(f'{d}/meta.json', 'w'))
    print(f'{a.name}: {len(items)} proposals from {len(pids)} parents x {len(a.arms)} arms x {a.per}; '
          f'{len(known)} known certified basins at n={a.n}')


_CTX = {}


def _init(d, variant, known_mode):
    meta = json.load(open(f'{d}/meta.json'))
    Z = np.load(f'{d}/props.npz')
    _CTX.update(meta=meta, Z=Z, pol=VARIANTS[variant],
                known=Known(meta['known'] if known_mode == 'ref' else []))


def _work(item):
    meta, Z = _CTX['meta'], _CTX['Z']
    sq = Z[f'arr_{item["i"]}']
    prop = Packing(item['s0'], sq)
    r = evaluate(prop, _CTX['pol'], _CTX['known'], meta['best'], item['seed'] ^ 0x5eed, meta['k'])
    p = r.pop('p', None)
    return item['i'], r, (None if p is None else (p.s, p.sq))


def replay(a):
    d = os.path.join(ROOT, a.name)
    items = [json.loads(l) for l in open(f'{d}/items.jsonl')]
    if a.limit:
        items = items[:a.limit]
    t0 = time.time()
    out = {}
    finals = {}
    with ProcessPoolExecutor(a.procs, initializer=_init, initargs=(d, a.variant, a.known)) as ex:
        for i, r, fin in ex.map(_work, items, chunksize=2):
            out[i] = r
            if fin is not None and r['st'] == 'new?':
                finals[str(i)] = np.concatenate([[fin[0], 0, 0], fin[1].ravel()])
    tag = a.variant + ('' if a.known == 'ref' else f'_{a.known}')
    with open(f'{d}/out_{tag}.jsonl', 'w') as f:
        for i in sorted(out):
            f.write(json.dumps(dict(i=i, **out[i])) + '\n')
    np.savez_compressed(f'{d}/fin_{tag}.npz', **finals)
    cpu = sum(r['sec'] for r in out.values())
    print(f'{a.name}/{tag}: {len(out)} proposals, wall {time.time() - t0:.0f}s, CPU {cpu:.0f}s; '
          + ' '.join(f'{k}={v}' for k, v in collections.Counter(r['st'] for r in out.values()).most_common()))


# ---------- certification (shared, cached in the store's basin table) ----------
def _solve(args):
    key, s, sq = args
    import tempfile, contextlib, io
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), '..', 'exact'))
    import census_exact
    with tempfile.TemporaryDirectory() as tmp:
        Packing(s, sq).write(f'{tmp}/in.txt')
        try:
            r = census_exact.solve((f'{tmp}/in.txt', tmp))
        except Exception as e:
            r = dict(cls='unresolved', status=repr(e)[:60])
    return key, r


def certify(a):
    st = Store()
    d = os.path.join(ROOT, a.name)
    meta = json.load(open(f'{d}/meta.json'))
    n, k = meta['n'], meta['k']
    variants = a.variants or [f[4:-6] for f in os.listdir(d) if f.startswith('out_')]
    todo = {}
    for v in variants:
        F = np.load(f'{d}/fin_{v}.npz')
        for l in open(f'{d}/out_{v}.jsonl'):
            r = json.loads(l)
            if r['st'] != 'new?' or r.get('unpolished') or r.get('lines') or r['s'] >= k:
                continue
            key = basin_key(r['s'])
            if key in todo or st.basin(n, key) is not None:
                continue
            x = F[str(r['i'])]
            todo[key] = (key, float(x[0]), x[3:].reshape(-1, 3))
    print(f'certify {a.name}: {len(todo)} distinct uncached sides below {k}')
    t0 = time.time()
    with ProcessPoolExecutor(a.procs) as ex:
        for key, r in ex.map(_solve, todo.values()):
            st.set_basin(n, key, r.get('cls', 'unresolved'), r.get('S'), r.get('status'))
    print(f'  done in {time.time() - t0:.0f}s')


def basin_key(s):
    return f'{s:.9f}'


# ---------- report ----------
def report(a):
    st = Store()
    d = os.path.join(ROOT, a.name)
    meta = json.load(open(f'{d}/meta.json'))
    n, k, best = meta['n'], meta['k'], meta['best']
    known = set(basin_key(s) for s in meta['known'])
    items = {it['i']: it for it in map(json.loads, open(f'{d}/items.jsonl'))}
    variants = a.variants or sorted(f[4:-6] for f in os.listdir(d) if f.startswith('out_'))
    B = {r['key']: dict(r) for r in st.db.execute('SELECT * FROM basin WHERE n=?', (n,))}
    rows = []
    union = set()
    per = {}
    for v in variants:
        R = [json.loads(l) for l in open(f'{d}/out_{v}.jsonl')]
        per[v] = R
        for r in R:
            if r['st'] == 'new?' and not r.get('unpolished') and r['s'] < k and B.get(basin_key(r['s']), {}).get('cls') == 'certified':
                union.add(basin_key(r['s']))
    print(f'{a.name}: n={n} k={k} best={best:.10f}; {len(items)} proposals; certified basins found by any variant: '
          f'{len(union)} ({len(union - known)} not in the frozen known set)')
    hdr = f'{"variant":12s} {"CPU s":>7s} {"/prop":>6s} {"screen":>6s} {"polish":>6s} {"grid":>5s} {"disc":>5s} {"ret":>5s} {"late":>5s} {"new?":>5s} {"cert":>5s} {"newcert":>7s} {"recall":>6s} {"<1e-3":>5s} {"notmin":>6s}'
    print(hdr)
    for v in variants:
        R = per[v]
        c = collections.Counter(r['st'] for r in R)
        cpu = sum(r['sec'] for r in R)
        ts = sum(r.get('t_screen', 0) for r in R); tp = sum(r.get('t_polish', 0) for r in R)
        found = set(); notmin = 0
        for r in R:
            if r['st'] == 'new?' and not r.get('unpolished') and r['s'] < k:
                b = B.get(basin_key(r['s']))
                if b and b['cls'] == 'certified':
                    found.add(basin_key(r['s']))
                elif b and b['cls'] == 'not-min':
                    notmin += 1
        near = {x for x in found if float(x) < best + 1e-3}
        print(f'{v:12s} {cpu:7.0f} {cpu / len(R):6.2f} {ts:6.0f} {tp:6.0f} {c["screen-grid"]:5d} {c["screen-discard"]:5d} '
              f'{c["return"]:5d} {c["late-return"]:5d} {c["new?"]:5d} {len(found):5d} {len(found - known):7d} '
              f'{len(found) / max(1, len(union)):6.2f} {len(near):5d} {notmin:6d}')
    if a.by:
        v = variants[0]
        print(f'\nby {a.by} ({v}): proposals, CPU s, grid share, certified distinct, new certified, per CPU-h')
        G = collections.defaultdict(list)
        for r in per[v]:
            it = items[r['i']]
            if a.by == 'arm':
                g = it['arm']
            elif a.by.startswith('param:'):
                x = it['params'].get(a.by[6:])
                g = 'na' if not isinstance(x, (int, float)) else f'{10 ** math.floor(math.log10(x) * 4) / 4 if x > 0 else 0:.3g}'
            elif a.by == 'loosen':
                g = r.get('loosen')
            elif a.by == 'parent_gap':
                gap = it['parent_s'] - best
                g = '<1e-4' if gap < 1e-4 else '<1e-3' if gap < 1e-3 else '<1e-2' if gap < 1e-2 else '>=1e-2'
            else:
                g = str(it.get(a.by))
            G[g].append(r)
        for g, R in sorted(G.items(), key=lambda t: str(t[0])):
            cpu = sum(r['sec'] for r in R)
            found = {basin_key(r['s']) for r in R if r['st'] == 'new?' and not r.get('unpolished') and r['s'] < k
                     and B.get(basin_key(r['s']), {}).get('cls') == 'certified'}
            grid = sum(r['st'] == 'screen-grid' for r in R) / len(R)
            print(f'  {str(g):10s} {len(R):5d} {cpu:7.0f} {grid:5.2f} {len(found):5d} {len(found - known):5d} '
                  f'{len(found - known) / max(cpu, 1) * 3600:7.1f}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('build'); p.add_argument('name'); p.add_argument('--n', type=int, required=True)
    p.add_argument('--parents', required=True); p.add_argument('--arms', nargs='+', default=list(moves.WEIGHTS))
    p.add_argument('--per', type=int, default=4); p.add_argument('--seed', type=int, default=1)
    p = sp.add_parser('replay'); p.add_argument('name'); p.add_argument('variant'); p.add_argument('--procs', type=int, default=16)
    p.add_argument('--known', default='ref', choices=['ref', 'none']); p.add_argument('--limit', type=int)
    p = sp.add_parser('certify'); p.add_argument('name'); p.add_argument('variants', nargs='*'); p.add_argument('--procs', type=int, default=16)
    p = sp.add_parser('report'); p.add_argument('name'); p.add_argument('variants', nargs='*'); p.add_argument('--by')
    a = ap.parse_args()
    dict(build=build, replay=replay, certify=certify, report=report)[a.cmd](a)


if __name__ == '__main__':
    main()
