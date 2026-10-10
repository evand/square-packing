#!/usr/bin/env python3
"""Bring a per-n basin census up to date (10-10).

The store keeps two layers: `packing` rows (deduplicated by exact coordinates) and the `basin` table (one row per polished
side to 1e-9, with exactsolve's class: certified / not-min / unresolved).  `certify` solves every side below k that has
packing rows but no basin row: one representative each (best status, then lowest penalty), exactsolve via
census_exact.solve, result into `basin` and the representative's status.  Unlike corpus.py certify, the exactsolve
output directory is kept (runs/census/n<n>/<key>/: .exact.txt, .cert, contacts.json), so contact graphs can be compared
later ("same basin" is a contact-graph question; matching distance is the parent-usefulness proxy).

  census_certify.py certify --n 110 [--k 11] [--procs 14] [--timeout 900] [--limit N]
  census_certify.py report  --n 110 [--k 11]
"""
import argparse, collections, json, math, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from jobpool import run_jobs
from pk.store import Store, RANK


def kgrid(n):
    r = math.isqrt(n)
    return r if r * r == n else r + 1


def todo(st, n, k):
    have = set(r[0] for r in st.db.execute('SELECT key FROM basin WHERE n=?', (n,)))
    rep = {}
    for pid, s, status, pen in st.db.execute('SELECT id, s, status, pen FROM packing WHERE n=? AND s<?', (n, k)):
        key = f'{s:.9f}'
        if key in have:
            continue
        cand = (-RANK.get(status, 0), pen if pen is not None else 0.0, pid)
        if key not in rep or cand < rep[key][0]:
            rep[key] = (cand, pid)
    return sorted((key, pid) for key, (_, pid) in rep.items())


def job(arg):
    key, pid, n, root = arg
    import census_exact
    st = Store()
    out = f'{root}/{key}'
    os.makedirs(out, exist_ok=True)
    st.get(pid).write(f'{out}/in.txt')
    try:
        r = census_exact.solve((f'{out}/in.txt', out))
    except Exception as e:
        r = dict(cls='unresolved', status=repr(e)[:60])
    json.dump(dict(r, pid=pid, key=key), open(f'{out}/class.json', 'w'))
    return key, pid, r


def certify(a):
    st = Store()
    k = a.k or kgrid(a.n)
    T = todo(st, a.n, k)[:a.limit or None]
    root = f'{HERE}/runs/census/n{a.n}'
    print(f'n={a.n}: {len(T)} sides below {k} without a basin row; solving on {a.procs} procs', flush=True)
    t0 = time.time()
    res = run_jobs(job, [(key, pid, a.n, root) for key, pid in T], procs=a.procs, timeout=a.timeout)
    C = collections.Counter()
    for (key, pid), r in zip(T, res):
        r = r[2] if r else dict(cls='unresolved', status='timeout')
        cls = r.get('cls', 'unresolved')
        st.set_basin(a.n, key, cls, r.get('S'), r.get('status'))
        if RANK.get(cls, 0) > RANK.get(st.db.execute('SELECT status FROM packing WHERE id=?', (pid,)).fetchone()[0], 0):
            st.set(pid, status=cls)
        C[cls] += 1
    print(f'  done in {time.time() - t0:.0f}s: {dict(C)}')


def report(a):
    st = Store()
    k = a.k or kgrid(a.n)
    B = st.db.execute('SELECT key, cls, S FROM basin WHERE n=?', (a.n,)).fetchall()
    sides = set(f'{s:.9f}' for (s,) in st.db.execute('SELECT s FROM packing WHERE n=? AND s<?', (a.n, k)))
    sub = [b for b in B if float(b[0]) < k]
    cls = collections.Counter(b[1] for b in sub)
    cert = [b for b in sub if b[1] == 'certified']
    byS = collections.defaultdict(list)
    noS = [key for key, _, S in cert if not S]
    for key, _, S in cert:
        if S:
            byS[S[:23]].append(key)
    print(f'n={a.n}, below {k}: {len(sub)} side keys in `basin` ({dict(cls)}); {len(sides - set(b[0] for b in B))} sides still unsolved')
    print(f'  certified local minima: {len(cert)} side keys -> {len(byS)} distinct exact sides (S to 1e-20)'
          f'{f"; {len(noS)} certified rows without a stored S" if noS else ""}')
    S = sorted(float(s[:20]) for s in byS)
    for g in (1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2):
        print(f'    within {g:g} of the best certified: {sum(x - S[0] <= g for x in S)}')
    multi = [v for v in byS.values() if len(v) > 1]
    print(f'  exact sides reached under more than one 1e-9 side key: {len(multi)} (largest {max(map(len, multi)) if multi else 0})')
    kept = sum(os.path.exists(f'{HERE}/runs/census/n{a.n}/{b[0]}/class.json') for b in sub)
    print(f'  with kept exactsolve output (contacts): {kept}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    for c in ('certify', 'report'):
        p = sp.add_parser(c); p.add_argument('--n', type=int, required=True); p.add_argument('--k', type=float)
        if c == 'certify':
            p.add_argument('--procs', type=int, default=14); p.add_argument('--timeout', type=float, default=900)
            p.add_argument('--limit', type=int)
    a = ap.parse_args()
    {'certify': certify, 'report': report}[a.cmd](a)
