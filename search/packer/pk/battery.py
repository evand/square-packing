"""Idea battery: one command from a new move / schedule to a verdict against a frozen baseline, with uncertainty.

  battery.py freeze [--n 110]                       write runs/battery/suite_<n>.json (frozen parents, tasks, withheld set)
  battery.py run NAME ARM [--minutes 10] [--procs 16] [--n 110]   ARM: pk.moves arm spec ('ashallow', 'kick:sigma=0.02',
                                                    or 'mix' = the explorer's weighted arm mix = the baseline)
  battery.py report NAME [--vs baseline]

Suites (all n = 110 for now; frozen at `freeze` time):
  bench   certified sub-k parents spread over the side range: new certified basins / CPU-h, near-frontier new basins
          (< best + 1e-3) / CPU-h, distinct polished sub-k sides / CPU-h (screening proxy, many more counts)
  reach   parents at matching distance >= 0.02 from the record, the record's 0.02-neighbourhood withheld: withheld basins
          rediscovered / CPU-h; continuous: dd = d(child, record) - d(parent, record) (fraction with dd < -0.005)
  cross   clean above-grid starts (no sub-k ancestry, >= 0.1 from sub-k): grid-collapse fraction, crossings; positive
          control (sub-k pushed to k + 0.01): return-below-k fraction
Rounds of ~1 minute (all suites interleaved); new sub-k sides are certified after each round (store basin cache).
Statistics: counts per CPU -> Gamma(x + 0.5, CPU) posterior; proportions -> Beta(x + 0.5, N - x + 0.5);
P(candidate > baseline) by sampling.  Early stop when every headline metric has P < 0.05 or > 0.95 (after >= 3 rounds).
"""
from __future__ import annotations
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import argparse, collections, json, math, random, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pk.store import Store
from pk.packing import Packing
from pk import moves
from pk.pipeline import Policy, Known, evaluate

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs', 'battery')
POLICY = Policy(same=0.0)                 # polish everything (10-09 bench)
key = lambda s: f'{s:.9f}'


def kof(n):
    return math.ceil(math.sqrt(n) - 1e-12)


# ---------- freeze ----------
def freeze(a):
    from pk.crossing import clean_above, match_distance
    st = Store(); n = a.n; k = kof(n)
    rows = list(st.query(n=n, below=k, status='certified', limit=None))
    reps = {}
    for r in rows:
        reps.setdefault(key(r['s']), r['id'])
    rec_id = rows[0]['id']; rec = st.get(rec_id)
    dist = json.load(open(a.dist)) if a.dist and os.path.exists(a.dist) else None
    if dist is None:
        dist = dict(basins=[dict(id=i, s=st.get(i).s, d=match_distance(st.get(i), rec)) for i in reps.values()])
    D = dist['basins']
    withheld = sorted(key(x['s']) for x in D if x['d'] < 0.02)
    rng = random.Random(7)
    far = [x for x in D if x['d'] >= 0.02]
    reach = [x['id'] for x in rng.sample(far, min(30, len(far)))]
    allb = sorted(reps.items())
    bench = [allb[i][1] for i in np.linspace(0, len(allb) - 1, 40).round().astype(int)]
    above = [i for i, _ in clean_above(st, n, 0.1, cap=60, seed=7)]
    pos = [allb[i][1] for i in np.linspace(0, len(allb) - 1, 20).round().astype(int)]
    known = sorted(float(r['key']) for r in st.db.execute("SELECT key FROM basin WHERE n=? AND cls='certified'", (n,)))
    fr = st.frontier().get(n)
    suite = dict(n=n, k=k, best=min([fr['s']] if fr else [] + known[:1]), record=rec_id, rec_s=rec.s, withheld=withheld,
                 known=[s for s in known if key(s) not in set(withheld)], bench=bench, reach=reach,
                 reach_d={str(x['id']): x['d'] for x in D}, above=above, pos=pos, frozen=time.strftime('%F %T'))
    os.makedirs(ROOT, exist_ok=True)
    json.dump(suite, open(f'{ROOT}/suite_{n}.json', 'w'))
    print(f'suite n={n}: bench {len(bench)}, reach {len(reach)} (withheld {len(withheld)}), above {len(above)}, pos {len(pos)}')


# ---------- one proposal ----------
_C = {}


def _init(suite_path):
    S = json.load(open(suite_path))
    st = Store()
    P = {i: st.get(i) for i in set(S['bench'] + S['reach'] + S['above'] + S['pos'] + [S['record']])}
    _C.update(S=S, P=P, rec=P[S['record']], known=Known(S['known']))


def _arm(spec, rng):
    if spec == 'mix':
        kinds, w = zip(*[(kk, ww) for kk, ww in moves.WEIGHTS.items() if ww > 0])
        return rng.choices(kinds, weights=w)[0]
    return spec


def _job(j):
    suite, pid, spec, seed = j
    S, P = _C['S'], _C['P']
    rng = random.Random(seed)
    k = S['k']
    p = P[pid]
    if suite == 'pos':
        f = (k + 0.01) / p.s
        sq = p.sq.copy(); sq[:, :2] = sq[:, :2] * f + np.array([[rng.gauss(0, .01), rng.gauss(0, .01)] for _ in range(p.n)])
        p = Packing(k + 0.01, sq)
    arm = _arm(spec, rng)
    t0 = time.time()
    extra = {}
    if arm.split(':')[0] == 'cross':                    # mate: another parent of the same suite
        pool = {'bench': S['bench'], 'reach': S['reach'], 'cross': S['above'], 'pos': S['pos']}[suite]
        mid = rng.choice([i for i in pool if i != pid])
        extra['mate'] = P[mid]
    try:
        q, info = moves.apply(p, arm, rng.randrange(1 << 30), **extra)
    except Exception as e:
        return dict(suite=suite, pid=pid, arm=arm, st='move-fail', sec=time.time() - t0, err=repr(e)[:100])
    r = evaluate(q, POLICY, _C['known'], S['best'], seed, k)
    fin = r.pop('p', None)
    out = dict(suite=suite, pid=pid, arm=arm, st=r['st'], s=r.get('s'), sec=time.time() - t0, params=info.get('params'))
    if 'mate_dist' in info:
        out['mate'] = extra['mate'].meta.get('id'); out['mate_dist'] = info['mate_dist']
    if fin is not None:
        out['lines'] = bool(fin.grid_lines(k))
        if suite == 'reach' and r['st'] == 'new?':
            from pk.crossing import match_distance
            out['dd'] = match_distance(fin, _C['rec']) - S['reach_d'][str(pid)]
        if r['st'] == 'new?' and fin.s < k and not out['lines']:
            out['sq'] = fin.sq.tolist()
    return out


def run(a):
    sp = f'{ROOT}/suite_{a.n}.json'
    S = json.load(open(sp))
    d = f'{ROOT}/{a.name}'; os.makedirs(d, exist_ok=True)
    json.dump(dict(arm=a.arm, n=a.n, minutes=a.minutes), open(f'{d}/meta.json', 'w'))
    rng = random.Random(hash(a.name) & 0xffffffff)
    t_end = time.time() + 60 * a.minutes
    rnd = 0
    with ProcessPoolExecutor(a.procs, initializer=_init, initargs=(sp,)) as ex:
        while time.time() < t_end:
            rnd += 1
            jobs = []
            for suite, ids, m in (('bench', S['bench'], 2), ('reach', S['reach'], 2), ('cross', S['above'], 1), ('pos', S['pos'], 1)):
                for pid in ids:
                    for _ in range(m):
                        jobs.append((suite, pid, a.arm, rng.randrange(1 << 30)))
            rng.shuffle(jobs)
            with open(f'{d}/trials.jsonl', 'a') as f:
                for r in ex.map(_job, jobs, chunksize=2):
                    sq = r.pop('sq', None)
                    if sq is not None:
                        r['fin'] = sq
                    f.write(json.dumps(dict(r, round=rnd)) + '\n')
            certify(d, S)
            if rnd >= 3 and a.vs and decisive(a.name, a.vs):
                print(f'round {rnd}: decisive, stopping'); break
            print(f'round {rnd} done', flush=True)
    report(argparse.Namespace(name=a.name, vs=a.vs, n=a.n))


def _cert1(args):
    s, sq = args
    import tempfile
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    import census_exact
    with tempfile.TemporaryDirectory() as tmp:
        Packing(s, np.array(sq)).write(f'{tmp}/in.txt')
        try:
            r = census_exact.solve((f'{tmp}/in.txt', tmp))
        except Exception as e:
            r = dict(cls='unresolved', status=repr(e)[:60])
    return key(s), r


def certify(d, S, procs=16):
    st = Store()
    todo = {}
    lines = []
    for l in open(f'{d}/trials.jsonl'):
        r = json.loads(l)
        if 'fin' in r and r['s'] is not None:
            kk = key(r['s'])
            if kk not in todo and st.basin(S['n'], kk) is None:
                todo[kk] = (r['s'], r['fin'])
    if todo:
        with ProcessPoolExecutor(procs) as ex:
            for kk, r in ex.map(_cert1, todo.values()):
                st.set_basin(S['n'], kk, r.get('cls', 'unresolved'), r.get('S'), r.get('status'))


# ---------- metrics ----------
def metrics(name, n=110):
    S = json.load(open(f'{ROOT}/suite_{n}.json'))
    st = Store()
    B = {r['key']: r['cls'] for r in st.db.execute('SELECT key, cls FROM basin WHERE n=?', (S['n'],))}
    known = {key(s) for s in S['known']}; wh = set(S['withheld']); k = S['k']; best = S['best']
    R = [json.loads(l) for l in open(f'{ROOT}/{name}/trials.jsonl')]
    M = {}
    for suite in ('bench', 'reach'):
        rs = [r for r in R if r['suite'] == suite]
        cpu = sum(r['sec'] for r in rs) / 3600
        cert = {key(r['s']) for r in rs if 'fin' in r and B.get(key(r['s'])) == 'certified'}
        scr = {key(r['s']) for r in rs if 'fin' in r}
        M[f'{suite}.cpu_h'] = cpu
        M[f'{suite}.new_cert'] = ('rate', len(cert - known - wh), cpu)
        M[f'{suite}.near_new'] = ('rate', len({x for x in cert - known - wh if float(x) < best + 1e-3}), cpu)
        M[f'{suite}.sub_k_sides'] = ('rate', len(scr), cpu)
        if suite == 'reach':
            M['reach.withheld'] = ('rate', len(cert & wh), cpu)
            dd = [r['dd'] for r in rs if 'dd' in r]
            M['reach.toward'] = ('prop', sum(x < -0.005 for x in dd), len(dd))
    rs = [r for r in R if r['suite'] == 'cross']
    M['cross.grid'] = ('prop', sum(1 for r in rs if r.get('lines', True)), len(rs))
    M['cross.crossed'] = ('prop', sum(1 for r in rs if 'fin' in r), len(rs))
    rs = [r for r in R if r['suite'] == 'pos']
    M['pos.return'] = ('prop', sum(1 for r in rs if 'fin' in r), len(rs))
    M['cpu_per_prop'] = sum(r['sec'] for r in R) / max(1, len(R))
    M['proposals'] = len(R)
    return M


def _draw(m, size=20000, rng=np.random.default_rng(0)):
    kind, x, y = m
    if kind == 'rate':
        return rng.gamma(x + 0.5, 1.0 / max(y, 1e-9), size)
    return rng.beta(x + 0.5, y - x + 0.5, size)


HEADLINE = ['bench.new_cert', 'bench.near_new', 'reach.withheld', 'reach.toward', 'cross.grid', 'pos.return']
LOWER_BETTER = {'cross.grid'}


def compare(Mc, Mb):
    out = {}
    for kk, v in Mc.items():
        if isinstance(v, tuple) and kk in Mb:
            pc, pb = _draw(v), _draw(Mb[kk])
            p = float(np.mean(pc > pb))
            out[kk] = 1 - p if kk in LOWER_BETTER else p
    return out


def decisive(name, vs):
    try:
        P = compare(metrics(name), metrics(vs))
    except FileNotFoundError:
        return False
    return all(P.get(h, 0.5) < 0.05 or P.get(h, 0.5) > 0.95 for h in HEADLINE)


def fmt(m):
    if not isinstance(m, tuple):
        return f'{m:.3g}'
    kind, x, y = m
    if kind == 'rate':
        lo, hi = np.quantile(_draw(m), [0.05, 0.95])
        return f'{x:4d} = {x / max(y, 1e-9):6.1f}/CPUh [{lo:.1f}, {hi:.1f}]'
    lo, hi = np.quantile(_draw(m), [0.05, 0.95])
    return f'{x:4d}/{y:<4d} = {x / max(y, 1):.2f} [{lo:.2f}, {hi:.2f}]'


def report(a):
    Mc = metrics(a.name, a.n)
    Mb = metrics(a.vs, a.n) if a.vs and os.path.exists(f'{ROOT}/{a.vs}/trials.jsonl') else None
    P = compare(Mc, Mb) if Mb else {}
    arm = json.load(open(f'{ROOT}/{a.name}/meta.json'))['arm']
    print(f'battery {a.name} (arm {arm}): {Mc["proposals"]} proposals, {Mc["cpu_per_prop"]:.2f} CPU s each'
          + (f'; vs {a.vs} ({Mb["proposals"]} proposals)' if Mb else ''))
    print(f'{"metric":18s} {"candidate (90% interval)":42s} ' + (f'{"baseline":42s} P(cand better)' if Mb else ''))
    for kk, v in Mc.items():
        if not isinstance(v, tuple):
            continue
        line = f'{kk:18s} {fmt(v):42s} '
        if Mb and kk in Mb:
            line += f'{fmt(Mb[kk]):42s} {P[kk]:.2f}' + ('  *' if kk in HEADLINE and (P[kk] > 0.95 or P[kk] < 0.05) else '')
        print(line)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('freeze'); p.add_argument('--n', type=int, default=110); p.add_argument('--dist', default='runs/dist110.json')
    p = sp.add_parser('run'); p.add_argument('name'); p.add_argument('arm'); p.add_argument('--minutes', type=float, default=10)
    p.add_argument('--procs', type=int, default=16); p.add_argument('--n', type=int, default=110); p.add_argument('--vs', default='baseline')
    p = sp.add_parser('report'); p.add_argument('name'); p.add_argument('--vs', default='baseline'); p.add_argument('--n', type=int, default=110)
    a = ap.parse_args()
    {'freeze': freeze, 'run': run, 'report': report}[a.cmd](a)


if __name__ == '__main__':
    main()
