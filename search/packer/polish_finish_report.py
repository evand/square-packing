#!/usr/bin/env python3
"""Report for `polish_replay.py finish`: does an exactsolve (KKT Newton) finish from the polish state at iteration k
reproduce the full polish's side?  Per k: success (Newton converged, jammed, rational certificate valid) and side vs the
full polish (same within 1e-9 / lower / higher), failures by status.  Policy simulation: try the finish at the listed ks
in turn (each try costs the exactsolve time), else fall back to the full polish.

  polish_finish_report.py bench/polish110 [--tag fin] [--base base]
"""
import argparse, collections, json, statistics as st


def ok(r):
    return r.get('status') == 'ok' and r.get('newton') and r.get('cert') and r.get('S') is not None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('set'); ap.add_argument('--tag', default='fin'); ap.add_argument('--base', default='base')
    ap.add_argument('--policy', type=int, nargs='+', default=None)
    a = ap.parse_args()
    S = json.load(open(f'{a.set}/set.json'))
    base = json.load(open(f'{a.set}/out_{a.base}/res.json'))['res']
    F = json.load(open(f'{a.set}/out_{a.tag}/finish.json'))
    ks = sorted({r['k'] for v in F.values() for r in v})
    print(f'== {a.set}: {len(F)} states, finish tried at k = {ks}; base polish median {st.median(b["slp_it"] for b in base if b)} it')
    for k in ks:
        c = collections.Counter(); dsl = []; es = []
        for it, b in zip(S['items'], base):
            r = next((x for x in F.get(f'p{it["i"]:04d}.txt', []) if x['k'] == k), None)
            if not r or not b:
                continue
            es.append(r.get('es_sec', 0))
            if ok(r):
                d = r['S'] - b['s']; dsl.append(d)
                c['same' if abs(d) < 1e-9 else 'lower' if d < 0 else 'higher'] += 1
            else:
                c['fail:' + str(r.get('status'))] += 1
        lo = [d for d in dsl if d < -1e-9]; hi = [d for d in dsl if d > 1e-9]
        print(f'  k={k:3d}: {dict(c)}; lower by median {st.median(lo) if lo else 0:.1e}, higher by median '
              f'{st.median(hi) if hi else 0:.1e} (max {max(hi) if hi else 0:.1e}); exactsolve median {st.median(es):.2f} s')
    pol = a.policy or [k for k in ks if k >= 8]
    tb = tf = tpy = 0.0; c = collections.Counter(); dd = []
    for it, b in zip(S['items'], base):
        R = {x['k']: x for x in F.get(f'p{it["i"]:04d}.txt', [])}
        if not b:
            continue
        tb += b['wall']
        done = False; py = 0.0
        for k in pol:
            r = R.get(k)
            if not r:
                continue
            py += r.get('es_sec', 0)
            if ok(r):
                tf += r['wall_k']; d = r['S'] - b['s']; dd.append(d)
                c['same' if abs(d) < 1e-9 else 'lower' if d < 0 else 'higher'] += 1
                done = True; break
        if not done:
            tf += b['wall']; c['fallback'] += 1
        tpy += py
    print(f'  policy {pol}: {dict(c)}; polish time {tf:.0f} s vs base {tb:.0f} s ({tf / tb:.2f}x), plus exactsolve (Python) '
          f'{tpy:.0f} s; worst higher {max(dd, default=0):.1e}')


if __name__ == '__main__':
    main()
