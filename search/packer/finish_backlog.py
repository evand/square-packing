#!/usr/bin/env python3
"""Cheap finisher on the census's not-min backlog (10-11; PACKER.md "Finisher on the backlog").

Points (n = 110 side keys below k): (1) basin cls not-min never descended (no row in runs/census/n<n>/descend.jsonl);
(2) descended but still not-min, started from descend's output point; (3) basin cls unresolved or descend unresolved /
timeout, from the stored packing.  Each: explore.kick_retry(alm=False) = kicked (sigma 1e-4) uniform repair + SLP polish,
up to --tries, keep strictly lower; then exactsolve (census_exact.solve).  Certified / weak results go into the store's
`minimum` (outcome new / known by the table before the write) and `copy` (input side key -> minimum).
Primary metric: new certified minima per measured CPU-hour (manifest cpu_s), vs the explorer's ~1 / worker-h at 110.

  finish_backlog.py --n 110 --procs 4 [--tries 4] [--limit N] [--minutes M]
"""
import argparse, collections, contextlib, io, json, math, os, random, sys, tempfile, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from jobpool import run_jobs


def todo(st, n, k):
    D = {}
    p = f'{HERE}/runs/census/n{n}/descend.jsonl'
    if os.path.exists(p):
        for l in open(p):
            r = json.loads(l); D[r['key']] = r
    T = []
    def stored(key):
        s = float(key)
        r = st.db.execute('SELECT id FROM packing WHERE n=? AND s BETWEEN ? AND ? ORDER BY pen LIMIT 1',
                          (n, s - 5e-10, s + 5e-10)).fetchone()
        return r[0] if r else None
    for key, cls in st.db.execute("SELECT key, cls FROM basin WHERE n=? AND cls IN ('not-min','unresolved')", (n,)):
        if float(key) >= k:
            continue
        d = D.get(key)
        if d is None:
            pid = stored(key)
            if pid: T.append((key, 'never-descended' if cls == 'not-min' else 'unresolved', pid, None))
        elif d.get('cls') == 'not-min' and os.path.exists(f'{HERE}/runs/census/n{n}/descend/{key}/in.txt'):
            T.append((key, 'descend-still-not-min', None, f'{HERE}/runs/census/n{n}/descend/{key}/in.txt'))
        elif d.get('cls') not in ('certified', 'weak', 'not-min'):
            pid = stored(key)
            if pid: T.append((key, 'descend-unresolved', pid, None))
    return T


def job(arg):
    key, group, pid, path, n, tries, sigma, seed = arg
    import explore, census_exact, mcmin
    from pk.store import Store
    from pk.manifest import cpu_seconds
    c0 = cpu_seconds()
    if path:
        s, sq = mcmin.load_deg(path)
    else:
        p = Store().get(pid); s, sq = float(p.s), [tuple(map(float, r)) for r in p.sq]
    k = math.ceil(math.sqrt(n) - 1e-12)
    with tempfile.TemporaryDirectory() as tmp:
        s2, sq2, info = explore.kick_retry(s, sq, random.Random(seed), tmp, k, tries, sigma, alm=False)
        mcmin.write_deg(f'{tmp}/in.txt', s2, sq2)
        t1 = cpu_seconds()
        try:
            r = census_exact.solve((f'{tmp}/in.txt', tmp))
        except Exception as e:
            r = dict(cls='unresolved', status=repr(e)[:80])
    return dict(key=key, group=group, pid=pid, s_in=s, s_out=s2, drops=info['drops'], traj=info.get('traj', []),
                cls=r.get('cls'), S=r.get('S'), status=r.get('status'), cpu_finish=round(t1 - c0, 2),
                cpu_total=round(cpu_seconds() - c0, 2), sq=[list(x) for x in sq2] if r.get('cls') in ('certified', 'weak') else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=110); ap.add_argument('--procs', type=int, default=4)
    ap.add_argument('--tries', type=int, default=4); ap.add_argument('--sigma', type=float, default=1e-4)
    ap.add_argument('--limit', type=int); ap.add_argument('--timeout', type=float, default=300)
    ap.add_argument('--seed', type=int, default=1011); ap.add_argument('--out', default=None)
    a = ap.parse_args()
    import scratch; scratch.scope('finish_backlog')
    from pk.store import Store
    from pk import manifest
    from pk.packing import Packing
    import numpy as np
    out = a.out or f'{HERE}/runs/census/n{a.n}/finish_backlog'
    st = Store()
    k = math.ceil(math.sqrt(a.n) - 1e-12)
    T = todo(st, a.n, k)
    if os.path.exists(f'{out}/results.jsonl'):                  # resume: skip keys already finished
        done = {json.loads(l)['key'] for l in open(f'{out}/results.jsonl')}
        T = [t for t in T if t[0] not in done]
    rng = random.Random(a.seed); rng.shuffle(T)                 # random order: a partial run is a fair sample
    T = T[:a.limit or None]
    man = manifest.start(out, 'finish_backlog', seed=a.seed, args=vars(a), procs=a.procs,
                         extra=dict(todo=len(T), groups=dict(collections.Counter(t[1] for t in T))))
    print(f'n={a.n}: {len(T)} backlog points {dict(collections.Counter(t[1] for t in T))}', flush=True)
    jobs = [(key, g, pid, path, a.n, a.tries, a.sigma, rng.randrange(1 << 30)) for key, g, pid, path in T]
    res = run_jobs(job, jobs, procs=a.procs, timeout=a.timeout)
    C = collections.Counter(); G = collections.defaultdict(collections.Counter)
    with open(f'{out}/results.jsonl', 'a') as fo:
        for (key, g, pid, path), r in zip(T, res):
            r = r if r else dict(key=key, group=g, cls='unresolved', status='timeout')
            oc = r.get('cls') or 'unresolved'
            if oc in ('certified', 'weak') and r.get('S'):
                known = st.find_minimum(a.n, r['S']) is not None
                pid2, _ = st.add(Packing(r['s_out'], np.array(r['sq'])), source='finish:backlog', finder='ours',
                                 ref=f'from {key}', status=oc)
                mk, _ = st.add_minimum(a.n, r['S'], oc, r.get('status'), src='finish:backlog', pid=pid2)
                st.add_copy(a.n, key, mk, via='finish:backlog')
                oc = f'{oc}-{"known" if known else "new"}'; r['min_key'] = mk
            r.pop('sq', None); r['outcome'] = oc
            C[oc] += 1; G[g][oc] += 1
            fo.write(json.dumps(r) + '\n')
    man.finish(summary=dict(outcomes=dict(C), by_group={g: dict(c) for g, c in G.items()}))
    cpu_h = man.d['cpu_s'] / 3600
    newc = C['certified-new']
    print(f'outcomes {dict(C)}')
    for g, c in G.items(): print(f'  {g}: {dict(c)}')
    print(f'new certified minima {newc} in {cpu_h:.2f} measured CPU-h = {newc / max(cpu_h, 1e-9):.2f} per CPU-h '
          f'(+ {C["weak-new"]} new weak)')


if __name__ == '__main__':
    main()
