#!/usr/bin/env python3
"""Exact dedupe + coupon statistics for census3 runs.

  census_exact.py runs/X/rec.json runs/X/cand.json [--procs 14] [--k 11]

For every trial output below k without full lines: exactsolve (../exact) -> class
  certified   KKT point with a valid certificate and a local-minimum second-order status; key = S_exact to 1e-20
              (also: second order inconclusive, flagged 'inconclusive' in census_exact.json)
  not-min     exactsolve says not a local minimum (not jammed in the corner-corner MILP, or a second-order descent)
  unresolved  anything else (refused, inconclusive, solve failed, killed)
Trials at or above k (or with a full line) are keyed in f64: side to 1e-7 + tilted-count fingerprint.
Per reference and sigma, and pooled: N, hits per class, distinct D, singletons f1, doubletons f2,
Good-Turing unseen mass f1/N (P(next trial lands in an unseen class)), Chao1 = D + f1^2/(2 f2).
Exact results are cached in <run dir>/exact/<name>.json (re-runs only solve what is missing).
"""
import argparse, collections, contextlib, io, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
EXACT = os.path.join(HERE, '..', 'exact')
sys.path.insert(0, EXACT)
from jobpool import run_jobs

MIN_OK = ('PD', 'local minimum', 'PSD', 'rigid')


def solve(args):
    path, outdir = args
    import exactsolve as E
    with contextlib.redirect_stdout(io.StringIO()):
        rep = E.run(path, outdir=outdir, quiet=True)[0]
    return classify(rep)


def classify(rep):
    so = rep.get('second_order') or {}
    st = so.get('status', '') if isinstance(so, dict) else ''
    if rep.get('S_exact') and rep.get('cert_valid') and st.startswith(MIN_OK):
        return dict(cls='certified', key=rep['S_exact'][:23], S=rep['S_exact'][:32], status=st[:40])
    if rep.get('S_exact') and rep.get('cert_valid') and st.startswith('inconclusive'):
        # certified KKT point, second-order test blocked by weakly active contacts: keyed by S, flagged
        return dict(cls='certified', key=rep['S_exact'][:23], S=rep['S_exact'][:32], status=st[:40], inconclusive=True)
    if rep.get('status') == 'not jammed' or st.startswith('NOT'):
        return dict(cls='not-min', status=(rep.get('status') or st)[:60])
    return dict(cls='unresolved', status=(rep.get('status') or st or 'no status')[:60])


def load_run(p):
    d = json.load(open(p))
    if isinstance(d, list):                       # pre-10-05 census3 format: failed trials were dropped
        d = dict(ref=None, trials_out=d)
    pref = p[:-5]
    for r in d['trials_out']:
        r['path'] = r.get('path') or f"{pref}_{r['sig']:g}_{r['t']}.txt"
        r['run'] = os.path.basename(pref)
    return d['trials_out']


def stats(keys):
    c = collections.Counter(keys)
    N, D = len(keys), len(c)
    f1 = sum(1 for v in c.values() if v == 1)
    f2 = sum(1 for v in c.values() if v == 2)
    chao = D + (f1 * f1 / (2 * f2) if f2 else f1 * (f1 - 1) / 2)
    return N, D, f1, f2, (f1 / N if N else float('nan')), chao, c


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('runs', nargs='+'); ap.add_argument('--procs', type=int, default=14)
    ap.add_argument('--k', type=float, default=11); ap.add_argument('--timeout', type=float, default=900)
    a = ap.parse_args()
    T = [r for p in a.runs for r in load_run(p)]
    for r in T:
        r['outdir'] = os.path.join(os.path.dirname(r['path']), 'exact')
        r['name'] = os.path.splitext(os.path.basename(r['path']))[0]
        r['sub'] = 's' in r and r['s'] < a.k and r.get('lines', 1) == 0
    todo = [r for r in T if r['sub'] and not os.path.exists(os.path.join(r['outdir'], r['name'] + '.json'))]
    print(f'{len(T)} trials, {sum(r["sub"] for r in T)} below {a.k:g} without lines; {len(todo)} to solve', flush=True)
    res = run_jobs(solve, [(r['path'], r['outdir']) for r in todo], procs=a.procs, timeout=a.timeout,
                   on_timeout=lambda j: dict(cls='unresolved', status='killed'),
                   on_error=lambda j, e: dict(cls='unresolved', status=f'error {e}'[:60]))
    for r, x in zip(todo, res):
        r['ex'] = x or dict(cls='unresolved', status='lost')
        if r['ex']['status'] in ('killed', 'lost') or r['ex']['status'].startswith('error'):
            json.dump(dict(status=r['ex']['status']), open(os.path.join(r['outdir'], r['name'] + '.json'), 'w'))
    for r in T:
        if r['sub'] and 'ex' not in r:
            r['ex'] = classify(json.load(open(os.path.join(r['outdir'], r['name'] + '.json'))))
        if 's' not in r:
            r['key'] = 'trial ' + r.get('status', 'failed')
        elif r['sub']:
            e = r['ex']
            r['key'] = e['key'] if e['cls'] == 'certified' else f"{e['cls']} ~{r['s']:.7f}"
        else:
            r['key'] = f">=k {r['s']:.7f} t{r['fp'][0]}{' L' if r.get('lines') else ''}"
    json.dump([{k: v for k, v in r.items() if k not in ('outdir',)} for r in T],
              open(os.path.join(os.path.dirname(a.runs[0]), 'census_exact.json'), 'w'), indent=0)

    def report(label, R):
        cert = [r['key'] for r in R if r['sub'] and r['ex']['cls'] == 'certified']
        n_inc = sum(1 for r in R if r['sub'] and r['ex'].get('inconclusive'))
        n_nm = sum(1 for r in R if r['sub'] and r['ex']['cls'] == 'not-min')
        n_un = sum(1 for r in R if r['sub'] and r['ex']['cls'] == 'unresolved')
        n_up = sum(1 for r in R if 's' in r and not r['sub'])
        n_fail = sum(1 for r in R if 's' not in r)
        N, D, f1, f2, gt, chao, _ = stats(cert)
        print(f'{label:>18}: N={len(R):4d}  sub-k certified {N:4d} (distinct {D:3d}, f1 {f1:3d}, f2 {f2:3d}, '
              f'f1/N {gt:.3f}, Chao1 {chao:6.1f}; 2nd-order inconclusive {n_inc})  not-min {n_nm}  unresolved {n_un}  >=k/lines {n_up}  failed {n_fail}')

    print('\nper run and sigma (sub-k classes keyed by exact S to 1e-20):')
    for run in sorted({r['run'] for r in T}):
        for sg in sorted({r['sig'] for r in T if r['run'] == run}):
            report(f'{run} s={sg:g}', [r for r in T if r['run'] == run and r['sig'] == sg])
        report(f'{run} all', [r for r in T if r['run'] == run])
    report('pooled', T)
    _, _, _, _, _, _, c = stats([r['key'] for r in T if r['sub'] and r['ex']['cls'] == 'certified'])
    print('\ncertified sub-k classes (S_exact, hits by run):')
    for key, v in sorted(c.items()):
        by = collections.Counter(r['run'] for r in T if r['sub'] and r['key'] == key)
        print(f'  {key}  {v:4d}  {dict(by)}')
    other = collections.Counter(r['key'] for r in T if not (r['sub'] and r['ex']['cls'] == 'certified'))
    print('\nother outcomes:')
    for key, v in sorted(other.items()):
        print(f'  {key}  {v}')
