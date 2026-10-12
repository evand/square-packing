#!/usr/bin/env python3
"""Signed side change along the unrealised zero modes of the PSD-weak s(110) minima (10-11, diagnosis).
exactsolve.probe_flat records |dS|; this wraps it to log the signed dS per mode and direction (+t, -t), t = 1e-3.
dS > 0 both ways for every mode: the side grows along the mode at fixed contacts (higher-order minimum, not a proof);
dS < 0 anywhere: the contacts allow a smaller side (a descent candidate: check feasibility before calling it).

  python3 resolve110_sign.py [--procs 15]      (reads runs/resolve110/summary.json; writes runs/resolve110/<key>/sign.json)"""
import argparse, contextlib, inspect, io, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'exact'))
from jobpool import run_jobs
OUT = os.path.join(HERE, 'runs', 'resolve110')


def job(key):
    import exactsolve as E
    src = inspect.getsource(E.probe_flat).replace('dS = float(abs(Q.S - P.S))',
                                                  'dS = float(abs(Q.S - P.S)); SIGNED.append((k, sg, float(Q.S - P.S), float(nr)))')
    ns = vars(E); SIGNED = []; ns['SIGNED'] = SIGNED; exec(src, ns)
    d = f'{OUT}/{key}'
    with contextlib.redirect_stdout(io.StringIO()):
        E.run(f'{d}/in.txt', outdir=f'{d}/sign_run', quiet=True)
    json.dump(SIGNED, open(f'{d}/sign.json', 'w'))
    return key, SIGNED


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--procs', type=int, default=15); a = ap.parse_args()
    R = json.load(open(f'{OUT}/summary.json'))
    keys = [r['key'] for r in R if r and r.get('cls') == 'weak' and (r.get('status') or '').startswith('PSD')]
    for k in keys: os.makedirs(f'{OUT}/{k}/sign_run', exist_ok=True)
    t0 = time.time()
    res = run_jobs(job, keys, procs=a.procs, timeout=600, on_timeout=lambda j: None, on_error=lambda j, x: None)
    out = {}
    for r in res:
        if not r: continue
        key, S = r
        modes = {}
        for k, sg, ds, nr in S: modes.setdefault(k, {})[sg] = ds
        allpos = all(v > 0 for m in modes.values() for v in m.values())
        anyneg = any(v < 0 for m in modes.values() for v in m.values())
        out[key] = dict(modes=len(modes), all_positive=allpos, any_negative=anyneg,
                        min_dS=min((v for m in modes.values() for v in m.values()), default=None))
    json.dump(out, open(f'{OUT}/sign_summary.json', 'w'), indent=0)
    print(f'{len(out)} of {len(keys)} done in {time.time() - t0:.0f} s')
    print('dS > 0 both ways on every mode:', sum(v['all_positive'] for v in out.values()),
          '| some dS < 0:', sum(v['any_negative'] for v in out.values()))
    for k, v in sorted(out.items(), key=lambda kv: kv[1]['min_dS'] or 0)[:8]:
        print(' ', k[:16], v)


if __name__ == '__main__':
    main()
