#!/usr/bin/env python3
"""Bandit reward replay (10-21; Evan: "should we be trying to figure out better reward criteria in cases where we're
chasing a new below-grid packing?" -> "Yes, let's try the replay!").  Logs only, no search CPU.

Unit: an archive entry created by a move (not a start) in an explorer run, created in the first half of the run's
archive time (later entries have too little future to be judged).  Candidate rewards, computable at creation:
  legacy   s < k                                                         (explore --reward legacy, the default)
  value    min(1, 1e-3 / gap) x parent robustness, gap = s - min(public frontier(n), run best so far)
                                                                         (explore --reward value)
  closek   min(1, 1e-3 / max(s - k, 1e-12))  (closeness to the grid side itself)
  improve  s < parent s - 2e-9
  nnd      fingerprint distance to the run's fp bank (structural novelty; runs with --fp-bank only)
  blend    value x (1 + nnd rank within run)
Outcomes (the future of the run, through the entry's descendants, parent links):
  prog     some descendant (not the entry) sets a new run best (< running best at its time - 2e-9)
  subk     some descendant is below k                                    (runs that ever go below k)
  near     some descendant is within 1e-4 of the run's final best and distinct (> 2e-9) from it
Selection bias: the explorer expands low-side entries more, so every low-side reward "predicts" descendants.  Each
outcome is therefore also scored on expanded entries only (children > 0).
Metric: AUC of each reward for each outcome within a run (ties 0.5), averaged over runs with both classes, weighted by
the minority count; by n bucket.  Kind level: per run, Spearman between kinds' reward per CPU-s and outcome per CPU-s.

  python3 analysis/reward_replay.py --procs 2 [--json out.json] [--limit N]
"""
import argparse, collections, glob, json, math, os, sys
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
REW = ('legacy', 'value', 'closek', 'improve', 'nnd', 'blend')
OUT = ('prog', 'subk', 'near')
FRONT = {}


def auc(x, y):
    """Mann-Whitney AUC of score x for binary y (average ranks for ties)"""
    pos = sum(y); neg = len(y) - pos
    if not pos or not neg:
        return None
    o = sorted(range(len(x)), key=lambda i: x[i])
    r = [0.0] * len(x); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and x[o[j + 1]] == x[o[i]]:
            j += 1
        for t in range(i, j + 1):
            r[o[t]] = (i + j) / 2 + 1
        i = j + 1
    return (sum(r[i] for i in range(len(y)) if y[i]) - pos * (pos + 1) / 2) / (pos * neg)


def spearman(a, b):
    if len(a) < 3:
        return None
    def rk(v):
        o = sorted(range(len(v)), key=lambda i: v[i]); r = [0] * len(v)
        for p, i in enumerate(o): r[i] = p
        return r
    ra, rb = rk(a), rk(b); m = (len(a) - 1) / 2
    num = sum((x - m) * (y - m) for x, y in zip(ra, rb))
    den = math.sqrt(sum((x - m) ** 2 for x in ra) * sum((y - m) ** 2 for y in rb))
    return num / den if den else None


def one(run):
    try:
        E = []
        for l in open(os.path.join(run, 'archive.jsonl')):
            try:
                E.append(json.loads(l))
            except Exception:
                pass
        if len(E) < 30:
            return None
        n = None
        try:
            m = json.load(open(os.path.join(run, 'manifest.json'))); a = m.get('argv', [])
            n = int(a[a.index('--n') + 1]) if '--n' in a else None
        except Exception:
            pass
        if n is None:
            from mcmin import load_deg
            p = E[0].get('path')
            if not p or not os.path.exists(p):
                return None
            n = len(load_deg(p)[1])
        k = math.ceil(math.sqrt(n) - 1e-12)
        F = FRONT.get(n, float('inf'))
        byi = {e['i']: e for e in E}
        # robustness of each parent at each child's creation, from proposals in time order
        rob = {}
        pf = os.path.join(run, 'proposals.jsonl'); cpu_kind = collections.Counter()
        if os.path.exists(pf):
            nk, ng = collections.Counter(), collections.Counter()
            for l in open(pf):
                try:
                    q = json.loads(l)
                except Exception:
                    continue
                cpu_kind[q.get('kind')] += q.get('cpu') or q.get('sec') or 0
                p = q.get('par')
                if q.get('new') is not None:
                    rob[q['new']] = (1 + nk[p] - ng[p]) / (2 + nk[p])
                if q.get('st') != 'fail':
                    nk[p] += 1; ng[p] += (q.get('st') == 'screen-grid')
        order = sorted(E, key=lambda e: (e.get('t', 0), e['i']))
        tmax = max(e.get('t', 0) for e in E)
        best = float('inf'); prog = {}; runbest_at = {}
        for e in order:
            runbest_at[e['i']] = best
            prog[e['i']] = e['s'] < best - 2e-9 and not str(e.get('kind', '')).startswith('start')
            best = min(best, e['s'])
        final = best
        kids = collections.defaultdict(list)
        for e in E:
            if e.get('parent', -1) in byi:
                kids[e['parent']].append(e['i'])
        sub = {}                                          # (prog, subk, near) in the subtree incl. self
        for e in sorted(E, key=lambda e: -e['i']):
            i = e['i']
            me = (prog[i], e['s'] < k - 1e-9, abs(e['s'] - final) > 2e-9 and e['s'] < final + 1e-4)
            agg = list(me)
            for c in kids[i]:
                agg = [x or y for x, y in zip(agg, sub[c])]
            sub[i] = tuple(agg)
        everk = final < k - 1e-9
        rows = []
        nn = [e.get('nnd') for e in E if e.get('nnd') is not None]
        nnr = sorted(nn)
        for e in E:
            if str(e.get('kind', '')).startswith('start') or e.get('parent', -1) not in byi:
                continue
            if e.get('t', 0) > tmax / 2:
                continue
            i = e['i']; par = byi[e['parent']]
            y = [any(sub[c][j] for c in kids[i]) for j in range(3)]
            gap = max(e['s'] - min(F, runbest_at[i]), 1e-12)
            r = rob.get(i, 0.5)
            v = min(1.0, 1e-3 / gap) * r
            nd = e.get('nnd')
            nr = (sum(1 for x in nnr if x < nd) / len(nnr)) if nd is not None else None
            rows.append(dict(kind=e.get('kind'), exp=bool(kids[i]),
                             r=dict(legacy=float(e['s'] < k - 1e-9), value=v, closek=min(1.0, 1e-3 / max(e['s'] - k, 1e-12)),
                                    improve=float(e['s'] < par['s'] - 2e-9), nnd=nd, blend=v * (1 + nr) if nr is not None else None),
                             y=dict(prog=y[0], subk=y[1] if everk else None, near=y[2])))
        if len(rows) < 20:
            return None
        res = dict(run=run, n=n, k=k, rows=len(rows), everk=everk, auc={}, kinds={})
        for sel in ('all', 'exp'):
            R = rows if sel == 'all' else [x for x in rows if x['exp']]
            for o in OUT:
                for w in REW:
                    Z = [(x['r'][w], x['y'][o]) for x in R if x['r'][w] is not None and x['y'][o] is not None]
                    if len(Z) < 20:
                        continue
                    a = auc([z[0] for z in Z], [int(z[1]) for z in Z])
                    if a is not None:
                        res['auc'][f'{sel}/{o}/{w}'] = (a, min(sum(z[1] for z in Z), len(Z) - sum(z[1] for z in Z)))
        # kind level
        K = collections.defaultdict(lambda: collections.defaultdict(float))
        for x in rows:
            for w in REW:
                if x['r'][w] is not None:
                    K[x['kind']][w] += x['r'][w]
            for o in OUT:
                if x['y'][o] is not None:
                    K[x['kind']]['y_' + o] += x['y'][o]
        kinds = [kd for kd in K if cpu_kind.get(kd, 0) > 0]
        for o in OUT:
            for w in REW:
                a_ = [K[kd][w] / cpu_kind[kd] for kd in kinds]; b_ = [K[kd]['y_' + o] / cpu_kind[kd] for kd in kinds]
                if sum(b_) > 0 and sum(a_) > 0:
                    s_ = spearman(a_, b_)
                    if s_ is not None:
                        res['kinds'][f'{o}/{w}'] = s_
        return res
    except Exception as ex:
        return dict(run=run, error=repr(ex)[:200])


def init(F):
    FRONT.update(F)


def bucket(n):
    return '<70' if n < 70 else '70-129' if n < 130 else '130-199' if n < 200 else '200-259' if n < 260 else '>=260'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--procs', type=int, default=2); ap.add_argument('--json'); ap.add_argument('--limit', type=int)
    a = ap.parse_args()
    from pk.store import Store
    F = {n: r['s'] for n, r in Store().frontier().items()}
    runs = sorted(os.path.dirname(f) for f in glob.glob(os.path.join(HERE, 'runs', '**', 'archive.jsonl'), recursive=True)
                  if '/ext/' not in f)
    if a.limit:
        runs = runs[:a.limit]
    with ProcessPoolExecutor(a.procs, initializer=init, initargs=(F,)) as ex:
        R = [r for r in ex.map(one, runs, chunksize=8) if r]
    err = [r for r in R if 'error' in r]; R = [r for r in R if 'error' not in r]
    print(f'{len(runs)} runs, {len(R)} usable, {len(err)} errors, {sum(r["rows"] for r in R)} entries')
    for sel in ('all', 'exp'):
        for o in OUT:
            print(f'\n[{sel} entries] outcome {o}: weighted mean AUC (runs)  by n bucket')
            B = sorted({bucket(r['n']) for r in R}, key=lambda b: ['<70', '70-129', '130-199', '200-259', '>=260'].index(b))
            print('  reward   ' + ''.join(f'{b:>18s}' for b in B + ['all']))
            for w in REW:
                cells = []
                for b in B + ['all']:
                    X = [r['auc'][f'{sel}/{o}/{w}'] for r in R if f'{sel}/{o}/{w}' in r['auc'] and (b == 'all' or bucket(r['n']) == b)]
                    W = sum(m for _, m in X)
                    cells.append(f'{sum(a_ * m for a_, m in X) / W:.3f} ({len(X):4d})' if W else '        -')
                print(f'  {w:8s} ' + ''.join(f'{c:>18s}' for c in cells))
    print('\nkind level: mean Spearman (reward per CPU-s vs outcome per CPU-s across kinds), runs')
    for o in OUT:
        for w in REW:
            X = [r['kinds'][f'{o}/{w}'] for r in R if f'{o}/{w}' in r['kinds']]
            if X:
                print(f'  {o:5s} {w:8s} {sum(X) / len(X):+.3f} ({len(X)})')
    if a.json:
        json.dump(R, open(a.json, 'w'))
    if err:
        print('errors e.g.', err[:3])


if __name__ == '__main__':
    main()
