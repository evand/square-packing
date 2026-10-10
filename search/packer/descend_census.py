#!/usr/bin/env python3
"""Where do the census's not-min points go when descended properly?  (10-10)

For every `basin` row of class not-min (side below k) that still has a packing in the store: exact/descend.py's
descend() (corner-corner MILP kick + slp2, until the MILP first-order dS >= -thr), then exactsolve (census_exact.solve)
on the result.  Outcome per point: known (lands on an exact side already certified in `basin` before this run), new
(certified, side not certified before), not-min (still not a local minimum), unresolved / failed.  Results are appended
to runs/census/n<n>/descend.jsonl (re-runs skip keys already there); landed minima go into `basin` and the store
(source descend:census).  The point -> minimum map also measures duplication: how many census points collapse onto one
minimum.

  descend_census.py run    --n 110 [--procs 14] [--limit N] [--thr 1e-9] [--kick 1e-5] [--rounds 20] [--timeout 900]
  descend_census.py report --n 110
"""
import argparse, collections, contextlib, io, json, math, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'exact'))
from jobpool import run_jobs
from pk.store import Store, RANK
from census_certify import kgrid


def todo(st, n, k, done):
    rows = st.db.execute("SELECT key FROM basin WHERE n=? AND cls='not-min'", (n,)).fetchall()
    T = []
    for (key,) in rows:
        if key in done or float(key) >= k:
            continue
        s = float(key)
        r = st.db.execute('SELECT id FROM packing WHERE n=? AND s BETWEEN ? AND ? ORDER BY pen LIMIT 1',
                          (n, s - 5e-10, s + 5e-10)).fetchone()
        if r:
            T.append((key, r[0]))
    return T


def job(arg):
    key, pid, n, root, kick, rounds, thr = arg
    import descend as D, census_exact
    from pk.packing import Packing
    st = Store()
    p = st.get(pid)
    out = f'{root}/{key}'
    os.makedirs(out, exist_ok=True)
    t0 = time.time()
    log = io.StringIO()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            s2, sq2 = D.descend(p.s, [(x, y, math.radians(d)) for x, y, d in p.sq], kick, rounds, thr=thr,
                                 log=lambda *a: print(*a, file=log))
        sq2 = [(x, y, math.degrees(t) % 90.0) for x, y, t in sq2]
    except Exception as e:
        return dict(key=key, pid=pid, s_in=p.s, cls='failed', status=repr(e)[:80], t=time.time() - t0)
    Packing(s2, sq2).write(f'{out}/in.txt')
    t1 = time.time()
    try:
        r = census_exact.solve((f'{out}/in.txt', out))
    except Exception as e:
        r = dict(cls='unresolved', status=repr(e)[:60])
    open(f'{out}/descend.log', 'w').write(log.getvalue())
    return dict(key=key, pid=pid, s_in=p.s, s_out=s2, cls=r.get('cls', 'unresolved'), S=r.get('S'), status=r.get('status'),
                rounds=log.getvalue().count('round '), t_descend=t1 - t0, t_solve=time.time() - t1)


def run(a):
    st = Store()
    k = a.k or kgrid(a.n)
    root = f'{HERE}/runs/census/n{a.n}/descend'
    os.makedirs(root, exist_ok=True)
    outp = f'{HERE}/runs/census/n{a.n}/descend.jsonl'
    done = set(json.loads(l)['key'] for l in open(outp)) if os.path.exists(outp) else set()
    if a.skip_started:                                      # concurrent workers: skip keys another run has started
        done |= set(os.listdir(root))
    known = set((S or '')[:23] for (S,) in st.db.execute("SELECT S FROM basin WHERE n=? AND cls='certified'", (a.n,)) if S)
    T = todo(st, a.n, k, done)
    import random; random.Random(a.seed).shuffle(T)        # random order: a partial run is a fair sample
    T = T[:a.limit or None]
    print(f'n={a.n}: {len(T)} not-min points to descend ({len(done)} done before); {len(known)} certified exact sides known', flush=True)
    t0 = time.time()
    res = run_jobs(job, [(key, pid, a.n, root, a.kick, a.rounds, a.thr) for key, pid in T], procs=a.procs, timeout=a.timeout)
    C = collections.Counter()
    with open(outp, 'a') as fo:
        for (key, pid), r in zip(T, res):
            r = r or dict(key=key, pid=pid, cls='unresolved', status='timeout')
            if r['cls'] == 'certified':
                r['outcome'] = 'known' if r['S'][:23] in known else 'new'
                k2 = f"{float(r['S'][:20]):.9f}"
                if not st.basin(a.n, k2) or st.basin(a.n, k2)['cls'] != 'certified':
                    st.set_basin(a.n, k2, 'certified', r['S'], r['status'])
                from pk.packing import read
                pid2, _ = st.add(read(f"{root}/{key}/in.txt"), source='descend:census', finder='ours', ref=f'from {pid}',
                                 status='certified')
                r['pid_out'] = pid2
            else:
                r['outcome'] = r['cls']
            C[r['outcome']] += 1
            fo.write(json.dumps(r) + '\n')
    print(f'  done in {time.time() - t0:.0f}s: {dict(C)}')


def report(a):
    outp = f'{HERE}/runs/census/n{a.n}/descend.jsonl'
    R = [json.loads(l) for l in open(outp)]
    C = collections.Counter(r['outcome'] for r in R)
    print(f'n={a.n}: {len(R)} not-min points descended: {dict(C)}')
    ok = [r for r in R if 's_out' in r]
    if ok:
        dsv = sorted(r['s_in'] - r['s_out'] for r in ok)
        td = sorted(r['t_descend'] + r.get('t_solve', 0) for r in ok)
        print(f'  side drop: median {dsv[len(dsv) // 2]:.1e}, p90 {dsv[int(.9 * len(dsv))]:.1e}, max {dsv[-1]:.1e};'
              f' time per point: median {td[len(td) // 2]:.0f}s, p90 {td[int(.9 * len(td))]:.0f}s')
    land = collections.Counter(r['S'][:23] for r in R if r.get('S'))
    if land:
        print(f'  certified landings: {sum(land.values())} points -> {len(land)} distinct minima;'
              f' points per minimum: max {max(land.values())}, top {land.most_common(5)}')
        new = sorted(set(S for S, c in land.items() if any((r.get('S') or '')[:23] == S and r['outcome'] == 'new' for r in R)))
        print(f'  new minima: {len(new)}; lowest {new[:5]}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run'); r.add_argument('--n', type=int, required=True); r.add_argument('--k', type=float)
    r.add_argument('--procs', type=int, default=14); r.add_argument('--limit', type=int)
    r.add_argument('--thr', type=float, default=1e-9); r.add_argument('--kick', type=float, default=1e-5)
    r.add_argument('--rounds', type=int, default=20); r.add_argument('--seed', type=int, default=0); r.add_argument('--skip-started', action='store_true'); r.add_argument('--timeout', type=float, default=900)
    p = sp.add_parser('report'); p.add_argument('--n', type=int, required=True)
    a = ap.parse_args()
    {'run': run, 'report': report}[a.cmd](a)
