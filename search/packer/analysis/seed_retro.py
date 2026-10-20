#!/usr/bin/env python3
"""SR (10-19; Evan: "records make the best seeds, but play out relatively soon, and there are far more near-record seeds
available").  Retrospective over chain.py rounds (explore.py archives): per (round, start) the start kind, its gap to the
round's best start, effort (archive entries rooted at it) and whether any descendant beat the round's best side.

  kinds: best (best_<n>.txt = the chain best at n), carry (seed_c*), lib (seed_l*), kseed (seed_k*), other
  record age: earlier rounds of the same chain at the same n whose best start had the same side (to 1e-9)

  seed_retro.py [CHAINDIR ...]   (default: runs/tt2/S runs/tt2/C runs/tt3/A runs/tt3/G runs/tt3/L runs/tt3/M)
"""
import collections, glob, json, os, re, sys
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEF = [f'{HERE}/runs/tt2/{a}' for a in 'SC'] + [f'{HERE}/runs/tt3/{a}' for a in 'AGLM']
EPS = 1e-9


def kind_of(k):
    f = k.split(':', 1)[1] if ':' in k else k
    for pre, name in (('best_', 'best'), ('seed_c', 'carry'), ('seed_l', 'lib'), ('seed_k', 'kseed')):
        if f.startswith(pre):
            return name
    return 'other'


def rounds(chain):
    R = []
    for d in glob.glob(f'{chain}/r[0-9]*_n*'):
        m = re.match(r'r(\d+)_n(\d+)$', os.path.basename(d))
        if m and os.path.exists(f'{d}/archive.jsonl'):
            R.append((int(m[1]), int(m[2]), d))
    return sorted(R)


def rows(chains):
    out = []
    for ch in chains:
        age = collections.defaultdict(list)            # n -> best sides of earlier rounds
        for r, n, d in rounds(ch):
            try:
                A = [json.loads(l) for l in open(f'{d}/archive.jsonl')]
            except Exception:
                continue
            S = {a['i']: a for a in A if a.get('parent', 0) < 0}
            if not S:
                continue
            bs = [a['s'] for a in S.values() if kind_of(a['kind']) == 'best']
            best = min(bs) if bs else min(a['s'] for a in S.values())
            ag = sum(abs(x - best) < EPS for x in age[n])
            age[n].append(best)
            eff, lo = collections.Counter(), {}
            for a in A:
                rt = a.get('root', a['i'])
                eff[rt] += 1
                lo[rt] = min(lo.get(rt, 9e9), a['s'])
            for i, a in S.items():
                out.append(dict(chain=os.path.relpath(ch, HERE), round=r, n=n, kind=kind_of(a['kind']), gap=a['s'] - best,
                                eff=eff[i], hit=lo[i] < best - EPS, gain=max(0.0, best - lo[i]), age=ag,
                                has_best=bool(bs)))
    return out


def gbin(g):
    return '0' if g < EPS else '<1e-3' if g < 1e-3 else '1e-3..1e-2' if g < 1e-2 else '>1e-2'


def table(R, key, title):
    print(f'\n{title}')
    print(f'  {"group":24s} {"starts":>7s} {"hits":>5s} {"hit/start":>9s} {"desc":>8s} {"hit/1k desc":>11s}')
    G = collections.defaultdict(list)
    for x in R:
        G[key(x)].append(x)
    for k in sorted(G, key=str):
        X = G[k]; h = sum(x['hit'] for x in X); e = sum(x['eff'] for x in X)
        print(f'  {str(k):24s} {len(X):7d} {h:5d} {h / len(X):9.4f} {e:8d} {1000 * h / max(e, 1):11.3f}')


def main():
    chains = sys.argv[1:] or DEF
    R = [x for x in rows(chains) if x['has_best']]
    rounds_ = {(x['chain'], x['round']) for x in R}
    print(f'{len(rounds_)} rounds with a best start, {len(R)} starts, chains {sorted({x["chain"] for x in R})}')
    table(R, lambda x: x['kind'], 'H1 by kind')
    table(R, lambda x: (x['kind'] == 'best', gbin(x['gap'])), 'H1 by gap (is_best, gap bin)')
    table([x for x in R if x['kind'] == 'best'], lambda x: min(x['age'], 3), 'H2 best starts by record age (3 = 3+)')
    table([x for x in R if x['kind'] != 'best'], lambda x: min(x['age'], 3), 'H2b non-best starts by record age')
    H = [x for x in R if x['hit']]
    for lo_age in (0, 1, 2):
        Y = [x for x in H if x['age'] >= lo_age]
        nb = sum(x['kind'] != 'best' for x in Y)
        print(f'H3 hits in rounds with age >= {lo_age}: {len(Y)}, rooted at non-best {nb} ({nb / max(len(Y), 1):.0%})')
    # gains: who made the big moves
    gb = sum(x['gain'] for x in H if x['kind'] == 'best'); gn = sum(x['gain'] for x in H if x['kind'] != 'best')
    print(f'total gain (sum over hits of best side - lowest descendant): best {gb:.4g}, non-best {gn:.4g}')
    json.dump(R, open(f'{HERE}/runs/seed_retro.json', 'w'))


if __name__ == '__main__':
    main()
