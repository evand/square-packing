#!/usr/bin/env python3
"""Offline parent / move value model (10-14, Evan: is parent choice a modeling problem that generalizes?).

Dataset: every proposal of the TT3 replay arms (39 n <= 240, same schedule), in time order per arm.  Features of the
parent at pick time (gap to frontier_at(T), polished or screened stepping stone, start seed or found, lineage depth,
times expanded so far, archive-new children so far), the move kind, n.  Labels:
  y   = the proposal is a distinct new near-frontier minimum: polished, within 1e-3 of frontier_at(T), > 2e-9 from every
        as-of-T reference side and from every earlier result of the same arm at that n (first sighting in the arm)
  yl  = lineage credit: y + GAMMA x (sum of y over the descendants of the archive entry this proposal created)
  cpu = measured CPU s of the proposal.
Model: Poisson GLM (log link, offset log cpu) by IRLS on one-hot / binned features; validation by 5-fold CV over n.
Report: yield per CPU-h by feature bins; held-out deviance vs a gap-only model and a constant; the yield the logged
CPU would have earned if reallocated to the top predicted share (an upper bound: ignores saturation beyond the
expansion features).

  python3 analysis/parent_model.py build [--arms G M L A]   -> runs/tt3/parent_ds.json
  python3 analysis/parent_model.py report
"""
import argparse, bisect, collections, glob, json, math, os, re, sqlite3, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE)
sys.path.insert(0, PK); sys.path.insert(0, HERE)
import tt2_census as TC                                                      # noqa: E402
from pk import eventbench as EB, parents as PA                               # noqa: E402

NEAR, SAME, GAMMA = 1e-3, 2e-9, 0.5
DS = os.path.join(PK, 'runs', 'tt3', 'parent_ds.json')
GAPB = (-1, 0, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 1e9)          # parent gap bins (left edges; first = below frontier_T)


def build(arms):
    os.chdir(PK)
    sched = [tuple(r) for r in json.load(open('runs/tt3/schedule.json'))]
    ns = sorted({n for _, n in sched})
    db = sqlite3.connect(f'file:{TC.DB}?mode=ro', uri=True)
    asof = PA.asof_local(TC.T)
    fT = {n: EB.frontier_at(n, TC.T) for n in ns}
    ref = {n: sorted(s for (s,) in db.execute("SELECT s FROM packing WHERE n=? AND s<? AND status!='infeasible' AND added<? "
                                               "AND (grid IS NULL OR grid=0)", (n, math.ceil(math.sqrt(n) - 1e-12), asof)))
           for n in ns}

    def near_key(L, s):
        i = bisect.bisect_left(L, s - SAME)
        return i < len(L) and abs(L[i] - s) < SAME
    rows = []
    for arm in arms:
        seen = {n: [] for n in ns}
        for r, n in sorted(sched):
            rd = f'runs/tt3/{arm}/r{r:04d}_n{n}'
            if not os.path.exists(f'{rd}/proposals.jsonl'):
                continue
            A = {e['i']: e for e in map(json.loads, open(f'{rd}/archive.jsonl'))}
            depth = {}
            for i in sorted(A):
                p = A[i]['parent']; depth[i] = depth.get(p, -1) + 1 if p >= 0 else 0
            P = [json.loads(l) for l in open(f'{rd}/proposals.jsonl')]
            P.sort(key=lambda p: p['t'])
            exp_, kids = collections.Counter(), collections.Counter()
            made = {}                                  # archive index created -> row index
            base = len(rows)
            for p in P:
                par = A.get(p.get('par'))
                if par is None:
                    continue
                y = 0
                s = p.get('s')
                if s is not None and p['st'] == 'new?' and not p.get('unp') and s - fT[n] < NEAR and not near_key(ref[n], s) \
                        and not near_key(seen[n], s):
                    y = 1
                if s is not None and p['st'] == 'new?' and not p.get('unp'):
                    bisect.insort(seen[n], s)
                rows.append(dict(arm=arm, n=n, r=r, kind=p['kind'], cpu=max(p['cpu'], 1e-3), y=y,
                                 gap=par['s'] - fT[n], unp=int(bool(par.get('unpolished'))),
                                 start=int(par['kind'].startswith('start')), depth=depth.get(par['i'], 0),
                                 exp=exp_[par['i']], kids=kids[par['i']], par=par['i'], new=p.get('new')))
                exp_[par['i']] += 1
                if p.get('new') is not None:
                    kids[par['i']] += 1; made[p['new']] = len(rows) - 1
            # lineage credit within the round: yl(row) = y + GAMMA * sum over rows whose parent is the entry this row made
            R = rows[base:]
            by_par = collections.defaultdict(list)
            for j, x in enumerate(R):
                by_par[x['par']].append(j)
            yl = [None] * len(R)
            for j in range(len(R) - 1, -1, -1):        # children come later in time than the row that made their parent
                x = R[j]; c = 0.0
                if x['new'] is not None:
                    c = sum(yl[k] for k in by_par.get(x['new'], ()) if k > j and yl[k] is not None)
                yl[j] = x['y'] + GAMMA * c
            for j, x in enumerate(R):
                x['yl'] = yl[j]
        print(arm, sum(1 for x in rows if x['arm'] == arm), 'rows', sum(x['y'] for x in rows if x['arm'] == arm), 'y', flush=True)
    json.dump(dict(T=TC.T, gamma=GAMMA, rows=rows), open(DS, 'w'))


def gbin(g):
    return max(0, bisect.bisect_right(GAPB, g) - 1)


def design(R, kinds, feats):
    cols = ['1']
    X = []
    for x in R:
        v = [1.0]
        if 'gap' in feats:
            b = gbin(x['gap']); v += [float(b == k) for k in range(1, len(GAPB) - 1)]
        if 'unp' in feats:
            v += [float(x['unp']), float(x['start'])]
        if 'exp' in feats:
            v += [math.log1p(x['exp']), math.log1p(x['kids'])]
        if 'depth' in feats:
            v += [math.log1p(x['depth'])]
        if 'kind' in feats:
            v += [float(x['kind'] == k) for k in kinds[1:]]
        if 'n' in feats:
            v += [float(x['n'] >= 72), float(x['n'] >= 151)]
        if 'gapxkind' in feats:                      # near-parent indicator x kind
            near = x['gap'] < 1e-3
            v += [float(near and x['kind'] == k) for k in kinds[1:]]
        X.append(v)
    return np.array(X)


def poisson_fit(X, y, off, l2=1e-2, it=50):
    b = np.zeros(X.shape[1]); b[0] = math.log(max(y.sum(), 1) / np.exp(off).sum())
    for _ in range(it):
        eta = X @ b + off; mu = np.exp(np.clip(eta, -30, 30))
        g = X.T @ (y - mu) - l2 * b; H = (X * mu[:, None]).T @ X + l2 * np.eye(len(b))
        st = np.linalg.solve(H, g); b += st
        if np.abs(st).max() < 1e-8:
            break
    return b


def dev(y, mu):
    t = np.where(y > 0, y * np.log(np.maximum(y, 1e-300) / mu), 0.0)
    return 2 * float(np.sum(t - (y - mu)))


def report(label='y', arms=('G', 'M', 'L')):
    D = json.load(open(DS)); R = [x for x in D['rows'] if x['arm'] in arms]
    kinds = sorted({x['kind'] for x in R})
    y = np.array([x[label] for x in R], float); cpu = np.array([x['cpu'] for x in R]); off = np.log(cpu)
    print(f'{len(R)} proposals, arms {arms}, label {label}: total {y.sum():.0f}, CPU-h {cpu.sum() / 3600:.2f}, '
          f'{y.sum() / (cpu.sum() / 3600):.1f} per CPU-h')
    # tables
    def tab(key, name):
        T = collections.defaultdict(lambda: [0.0, 0.0, 0])
        for x, yy in zip(R, y):
            k = key(x); T[k][0] += yy; T[k][1] += x['cpu']; T[k][2] += 1
        print(f'  by {name}:')
        tot = cpu.sum()
        for k in sorted(T):
            a, c, m = T[k]
            print(f'    {str(k):>22s}  CPU {c / tot:5.1%}  props {m:6d}  {label} {a:6.0f}  per CPU-h {a / (c / 3600):8.1f}')
    tab(lambda x: f'{GAPB[gbin(x["gap"])]:g}', 'parent gap bin (left edge)')
    tab(lambda x: ('start' if x['start'] else 'screened' if x['unp'] else 'polished'), 'parent type')
    tab(lambda x: f"{x['exp']:02d}" if x['exp'] < 4 else '04-15' if x['exp'] < 16 else '16-63' if x['exp'] < 64 else '64+', 'times parent expanded before')
    tab(lambda x: ('near' if x['gap'] < 1e-3 else 'far') + ' ' + x['kind'], 'near(<1e-3)/far parent x kind')
    # CV over n
    ns = sorted({x['n'] for x in R}); fold = {n: i % 5 for i, n in enumerate(ns)}
    F = np.array([fold[x['n']] for x in R])
    models = {'const': [], 'gap': ['gap'], 'gap+type': ['gap', 'unp'], 'gap+type+sat': ['gap', 'unp', 'exp', 'depth'],
              'full': ['gap', 'unp', 'exp', 'depth', 'kind', 'n'], 'full+gapxkind': ['gap', 'unp', 'exp', 'depth', 'kind', 'n', 'gapxkind']}
    print('  held-out (5-fold over n) Poisson deviance; lower is better; and top-share reallocation bound:')
    preds = {}
    for name, feats in models.items():
        X = design(R, kinds, feats); mu = np.zeros(len(R))
        for f in range(5):
            tr, te = F != f, F == f
            b = poisson_fit(X[tr], y[tr], off[tr]); mu[te] = np.exp(X[te] @ b + off[te])
        preds[name] = mu
        rate = mu / cpu
        o = np.argsort(-rate); cc = np.cumsum(cpu[o]); yy = np.cumsum(y[o])
        bd = []
        for q in (0.1, 0.25, 0.5):
            j = np.searchsorted(cc, q * cc[-1]); bd.append((yy[j] / cc[j]) / (yy[-1] / cc[-1]))
        print(f'    {name:16s} dev {dev(y, mu):9.1f}   yield per CPU in the top 10 / 25 / 50 % of CPU by predicted rate: '
              + ' / '.join(f'{v:.2f}x' for v in bd))
    b = poisson_fit(design(R, kinds, models['full']), y, off)
    names = ['1'] + [f'gap>={GAPB[k]:g}' for k in range(1, len(GAPB) - 1)] + ['par screened', 'par start', 'log1p exp',
             'log1p kids', 'log1p depth'] + [f'kind {k}' for k in kinds[1:]] + ['n>=72', 'n>=151']
    print('  full-model coefficients (log rate per CPU-s; ref: gap < 0, polished found parent, kind ' + kinds[0] + ', n < 72):')
    for nm, v in zip(names, b):
        print(f'    {nm:22s} {v:+.2f}  (x{math.exp(v):.2f})')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('cmd', choices=['build', 'report'])
    ap.add_argument('--arms', nargs='+', default=['G', 'M', 'L', 'A']); ap.add_argument('--label', default='y')
    a = ap.parse_args()
    if a.cmd == 'build':
        build(a.arms)
    else:
        report(a.label, tuple(x for x in a.arms if x != 'A') or ('G',))
