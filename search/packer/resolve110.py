#!/usr/bin/env python3
"""Re-solve the s(110) weak minima and unresolved sides with exactsolve output kept (10-11), for diagnosis.
Does not write the store.  One representative packing per target (best status, lowest penalty), as census_certify.

  python3 resolve110.py [--procs 15] [--timeout 600]      -> runs/resolve110/<key>/ (+ summary.json)"""
import argparse, collections, contextlib, io, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'exact'))
from jobpool import run_jobs
from pk.store import Store, RANK
OUT = os.path.join(HERE, 'runs', 'resolve110')


def job(arg):
    key, pid, kind = arg
    import exactsolve as E, census_exact
    st = Store(); d = f'{OUT}/{key}'; os.makedirs(d, exist_ok=True)
    st.get(pid).write(f'{d}/in.txt'); t = time.time()
    try:
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            rep = E.run(f'{d}/in.txt', outdir=d, quiet=False)[0]
        open(f'{d}/log.txt', 'w').write(buf.getvalue())
        so = rep.get('second_order') or {}
        r = dict(cls=census_exact.classify(rep)['cls'], status=so.get('status') if isinstance(so, dict) else None,
                 n_neg=so.get('n_neg'), n_zero=so.get('n_zero'), realised=so.get('zero_modes_realised'),
                 probe=so.get('zero_mode_probe'), S=rep.get('S_exact'))
    except Exception as e:
        r = dict(cls='error', status=repr(e)[:200])
    r.update(key=key, pid=pid, kind=kind, secs=round(time.time() - t, 1))
    json.dump(r, open(f'{d}/result.json', 'w'), default=str)
    return r


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--procs', type=int, default=15); ap.add_argument('--timeout', type=int, default=600)
    a = ap.parse_args(); os.makedirs(OUT, exist_ok=True); st = Store()
    best = {}
    for pid, s, status, pen in st.db.execute('SELECT id, s, status, pen FROM packing WHERE n=110 AND s<11'):
        k = f'{s:.9f}'; c = (-RANK.get(status, 0), pen or 0.0, pid)
        if k not in best or c < best[k]: best[k] = c
    cp = collections.defaultdict(list)
    for sk, mk in st.db.execute('SELECT side_key, min_key FROM copy WHERE n=110'): cp[mk].append(sk)
    T = []
    for (mk,) in st.db.execute("SELECT key FROM minimum WHERE n=110 AND cls='weak' AND CAST(S AS REAL) < 11"):
        c = min((best[sk] for sk in cp[mk] if sk in best), default=None)
        if c: T.append((mk[:20], c[2], 'weak'))
    for (k,) in st.db.execute("SELECT key FROM basin WHERE n=110 AND cls='unresolved'"):
        if k in best: T.append((k, best[k][2], 'unresolved'))
    print(f'{len(T)} targets ({collections.Counter(t[2] for t in T)}) on {a.procs} procs', flush=True)
    t0 = time.time()
    res = run_jobs(job, T, procs=a.procs, timeout=a.timeout, on_timeout=lambda j: dict(key=j[0], kind=j[2], cls='timeout'),
                   on_error=lambda j, x: dict(key=j[0], kind=j[2], cls='error', status=str(x)[:200]))
    json.dump(res, open(f'{OUT}/summary.json', 'w'), indent=0, default=str)
    print(f'done in {time.time() - t0:.0f} s')
    for kind in ('weak', 'unresolved'):
        R = [r for r in res if r and r.get('kind') == kind]
        print(kind, collections.Counter(r['cls'] for r in R))
        print('  statuses:', collections.Counter((r.get('status') or '')[:60] for r in R).most_common(6))


if __name__ == '__main__':
    main()
