#!/usr/bin/env python3
"""Lineage of a chain.py result: which rounds, counts and start packings it descends from (for credit; 10-12, after
the s(154) registration had to be traced by hand).

Walks back from a packing (default: the chain's best at n) through the round archives:
  round r of n: the archive entry at side s -> its root start ('start:<file>', explore.py)
    best_<n>.txt / seed_c* (carry) / seed_l* (library): an earlier result at n with the root's side (latest earlier
                 round of n whose archive holds that side to --tol), or n's chain start if no round does
    seed_k<k>.txt: built from best(n + k) at the round's start (k > 0: squares removed, k < 0: inserted): continue at
                 n + k from its chain best at that time (state.json log), recorded side if logged (seeds_src, 10-12+)
Carry / library / seed files stay in the round dirs, so a start is matched by its file's header side (exact copy of an
earlier archive entry); best_<n>.txt is overwritten, so it is matched by the polished side to --tol.  A step that does
not match is reported and the walk stops there.

  chain_lineage.py RUNDIR N [--side S] [--json]
"""
import argparse, glob, json, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))


def load_rounds(out):
    R = {}
    for d in glob.glob(f'{out}/r[0-9]*_n*'):
        m = re.match(r'r(\d+)_n(\d+)$', os.path.basename(d))
        if m and os.path.exists(f'{d}/archive.jsonl'):
            R[int(m[1])] = dict(n=int(m[2]), dir=d, E=None)
    return R


def archive(r):
    if r['E'] is None:
        r['E'] = {}
        for l in open(f'{r["dir"]}/archive.jsonl'):
            e = json.loads(l); r['E'][e['i']] = e
    return r['E']


def trace(out, n, side=None, tol=1e-8):
    S = json.load(open(f'{out}/state.json'))
    R = load_rounds(out)
    starts = {e['round']: e for e in S['log'] if e['ev'] == 'start'}
    done = [e for e in S['log'] if e['ev'] == 'done']
    if side is None:
        side = S['n'][str(n)]['best']
    steps, before = [], max(R) + 1 if R else 0
    while True:
        # latest round of n before `before` holding `side`: an entry produced in the round first, else one of its
        # starts (an unchanged carry / library pick / best, or a seed_k origin)
        hit = None
        for rn in sorted((rn for rn, r in R.items() if r['n'] == n and rn < before), reverse=True):
            E = archive(R[rn])
            es = [e for e in E.values() if abs(e['s'] - side) < tol]
            if es:
                hit = (rn, min(es, key=lambda e: (e['parent'] < 0, abs(e['s'] - side), e['i']))); break
        if hit is None:
            v = S['n'][str(n)]
            ok = abs(v.get('start', v['reg']) - side) < tol
            steps.append(dict(n=n, side=side, origin='chain start' if ok else 'UNMATCHED',
                              start_side=v.get('start', v['reg'])))
            return steps
        rn, e = hit
        E = archive(R[rn])
        root = E[e.get('root', e['i'])]
        path, hops = [], e
        while hops['parent'] >= 0:
            path.append(hops['kind']); hops = E[hops['parent']]
        kind = root['kind'].split(':', 1)[1] if root['kind'].startswith('start:') else root['kind']
        steps.append(dict(round=rn, n=n, side=side, root_side=root['s'], root=kind, moves=list(reversed(path))))
        m = re.match(r'seed_k(-?\d+)\.txt$', kind)
        if m:
            k = int(m[1])
            t0 = starts.get(rn, {}).get('t', 1e18)
            src = (starts.get(rn, {}).get('seeds_src') or {}).get(str(k))
            if src is None:      # n + k's chain best at the round's start: last 'done' of n + k before t0, else its start
                prev = [d for d in done if d['n'] == n + k and d['t'] <= t0]
                v = S['n'].get(str(n + k))
                src = min([d['best'] for d in prev], default=v.get('start', v['reg']) if v else None)
            if src is None:
                steps.append(dict(n=n + k, origin='UNMATCHED (not in the chain)')); return steps
            n, side, before = n + k, src, rn
        else:       # the start file as copied in (header side, before explore's start polish), else the polished side
            f = os.path.join(R[rn]['dir'], kind)
            side = float(open(f).readline().split()[1]) if os.path.exists(f) else root['s']
            before = rn


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('out'); ap.add_argument('n', type=int)
    ap.add_argument('--side', type=float); ap.add_argument('--tol', type=float, default=1e-8)
    ap.add_argument('--json', action='store_true')
    a = ap.parse_args()
    steps = trace(a.out, a.n, a.side, a.tol)
    if a.json:
        print(json.dumps(steps, indent=1)); return
    for s in steps:
        if 'round' in s and not s['moves'] and not s['root'].startswith('seed_k'):
            continue            # pass-through (the packing was a start of this round, unchanged); kept in --json
        if 'round' in s:
            mv = ' '.join(s['moves'])
            print(f'r{s["round"]:04d} n={s["n"]}  {s["root_side"]:.10f} ({s["root"]}) -> {s["side"]:.10f}  [{len(s["moves"])} moves: {mv[:80]}]')
        else:
            print(f'      n={s["n"]}  {s["origin"]}' + (f' {s["start_side"]:.10f}' if 'start_side' in s else ''))


if __name__ == '__main__':
    main()
