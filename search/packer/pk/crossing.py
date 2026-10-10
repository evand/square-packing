"""Grid-crossing tests: from states with side >= k and no sub-k ancestry, how often does a method reach a sub-k minimum?

  crossing.py run NAME --n 110 --classes grid above scratch pos neg --methods quench kick:0.05 anneal:diagonal:0.3:20000
              [--trials 20] [--procs 16] [--snapq 0]
  crossing.py report NAME [--by method,cls]

Start classes (per n; k = ceil sqrt n):
  grid     k x k lattice at side k with k^2 - n random holes (a fresh hole pattern per trial)
  above    store packings with k < s < k + 0.05, no full lines, no sub-k packing in their parent chain, and matching
           distance >= --dmin from the reference sub-k set (best certified basins); drawn per trial
  scratch  no start (constructor methods only; others skip)
  pos      positive control: certified sub-k packing scaled to side k + 0.01, kicked 0.01 (should cross back)
  neg      negative control at n' = k^2 - 4 (k >= 5, s(n') = k proved): crossing there is a bug
Methods (spec strings):
  quench                         screen + polish the start (no move)
  kick:SIGMA                     gaussian kick (angles 20 sigma deg), then quench
  anneal:PATH:RMAX:SWEEPS[:BP1]  schedule-driven anneal (pk.anneal paths), then quench; from scratch for class scratch
  melt:RAD:RMAX:SWEEPS           regional anneal around a random square (path hold), then quench
Outcome per trial: final polished side, full lines, crossed = side < k and no full line; matched certified basin if any.
Trajectory (anneal methods): line occupancy and side per snapshot; with --snapq K also the quenched side of the last K
snapshots (where the trajectory locks into the grid).
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
from pk import anneal as AN
from pk.pipeline import Policy, Known, evaluate

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs', 'crossing')
POLICY = Policy(same=0.0, smax_above_k=10.0, grid_check=False)   # polish everything; grid tested on the final state


def kof(n):
    return math.ceil(math.sqrt(n) - 1e-12)


def match_distance(p: Packing, q: Packing):
    """min over the 8 box symmetries of the optimal-matching RMS corner-displacement distance (rediscover.distance)."""
    import rediscover
    return rediscover.distance(p.tuples(), p.s, q.tuples(), q.s)


def grid_start(n, rng):
    k = kof(n)
    cells = [(i + 0.5, j + 0.5, 0.0) for i in range(k) for j in range(k)]
    rng.shuffle(cells)
    return Packing(k + 1e-9, cells[:n])


# ---------- task pools ----------
def clean_above(st, n, dmin, nref=30, cap=400, seed=1):
    """Above-grid library entries without sub-k ancestry, >= dmin from the reference sub-k set."""
    k = kof(n)
    rows = st.query(n=n, above=k, below=k + 0.05, nongrid=True, limit=None)
    rng = random.Random(seed)
    rows = list(rows); rng.shuffle(rows)
    refs = [st.get(r['id']) for r in st.query(n=n, below=k, status='certified', limit=nref)]
    out = []
    for r in rows:
        if len(out) >= cap:
            break
        q, leak = r['parent'], False
        for _ in range(500):
            if q is None:
                break
            a = st.db.execute('SELECT s, parent FROM packing WHERE id=?', (q,)).fetchone()
            if a['s'] < k:
                leak = True; break
            q = a['parent']
        if leak:
            continue
        p = st.get(r['id'])
        d = min((match_distance(p, f) for f in refs), default=9.0)
        if d >= dmin:
            out.append((r['id'], d))
    return out


# ---------- one trial ----------
def method_apply(spec, start, n, rng, snaps):
    """-> (proposal Packing, info dict, trajectory)"""
    parts = spec.split(':')
    kind = parts[0]
    if kind == 'quench':
        return start, {}, []
    if kind == 'kick':
        sig = float(parts[1])
        sq = start.sq.copy()
        sq[:, :2] += np.array([[rng.gauss(0, sig), rng.gauss(0, sig)] for _ in range(start.n)])
        sq[:, 2] += np.array([rng.gauss(0, 20 * sig) for _ in range(start.n)])
        return Packing(start.s, sq), dict(sigma=sig), []
    if kind == 'anneal':
        # anneal:PATH:RMAX:SWEEPS[:BP1][:k=v,...]  extras: bias=lin|harm, lam=, q0=, release= (bias held to `release`,
        # off by release + 0.05), bp0= (from a packing: starting pressure; default bp1 / 10 = no expansion)
        path, rmax, sweeps = parts[1], float(parts[2]), int(parts[3])
        bp1 = float(parts[4]) if len(parts) > 4 and '=' not in parts[4] else 3000.0
        kv = dict(x.split('=') for p_ in parts[4:] if '=' in p_ for x in p_.split(','))
        sch = AN.Schedule.from_path(path, sweeps=sweeps, rmax=rmax, bp1=bp1)
        if start is not None:
            sch = sch.but(bp=AN.curve([(0, float(kv.get('bp0', bp1 / 10))), (1, bp1)]))
        if 'bias' in kv:
            rel = float(kv.get('release', 0.9)); lam = float(kv['lam'])
            sch = sch.but(bias=kv['bias'], bias_lam=AN.curve([(0, lam), (rel, lam), (min(1, rel + 0.05), 0), (1, 0)]),
                          bias_q0=kv.get('q0', '0'))
        fin, info, traj = AN.run(sch, start=start, n=n, seed=rng.randrange(1 << 30), snaps=snaps)
        return fin, info, traj
    if kind == 'melt':
        rad, rmax, sweeps = float(parts[1]), float(parts[2]), int(parts[3])
        i = rng.randrange(start.n)
        sch = AN.Schedule.from_path('hold', sweeps=sweeps, rmax=rmax).but(
            bp='2000', box_every=4, region=(start.sq[i, 0], start.sq[i, 1], rad),
            r=AN.curve([(0, 0), (0.25, rmax), (0.55, rmax), (0.85, 0), (1, 0)]))
        fin, info, traj = AN.run(sch, start=start, seed=rng.randrange(1 << 30), snaps=snaps)
        return fin, info, traj
    raise ValueError(spec)


def trial(job):
    task, spec, seed, snaps, snapq, known_keys = job
    rng = random.Random(seed)
    n, k = task['n'], kof(task['n'])
    if task['cls'] == 'grid' or task['cls'] == 'neg':
        start = grid_start(n, rng)
    elif task['cls'] == 'scratch':
        start = None
    else:
        start = Packing(task['s'], np.array(task['sq']))
        if task['cls'] == 'pos':
            f = (k + 0.01) / start.s
            sq = start.sq.copy(); sq[:, :2] *= f
            sq[:, :2] += np.array([[rng.gauss(0, 0.01), rng.gauss(0, 0.01)] for _ in range(start.n)])
            start = Packing(k + 0.01, sq)
    if start is None and not spec.startswith('anneal'):
        return None
    t0 = time.time()
    prop, info, traj = method_apply(spec, start, n, rng, snaps)
    t_move = time.time() - t0
    if prop is None:
        return dict(task=task['name'], cls=task['cls'], n=n, method=spec, seed=seed, st='move-fail', err=str(info)[:200])
    r = evaluate(prop, POLICY, Known(), best=0.0, seed=seed, k=k)
    p = r.pop('p', None)
    s = r.get('s')
    lines = bool(p.grid_lines(k)) if p is not None else True
    crossed = r['st'] == 'new?' and s is not None and s < k and not lines
    out = dict(task=task['name'], cls=task['cls'], n=n, k=k, method=spec, seed=seed, st=r['st'], s=s, lines=lines,
               crossed=crossed, s_prop=prop.s, t_move=t_move, sec=time.time() - t0,
               basin=(f'{s:.9f}' in known_keys) if crossed else None, info={k_: v for k_, v in info.items() if k_ != 'error'})
    if traj:
        out['traj'] = [(round(x['t'], 3), round(x['s'], 6), AN.line_occupancy(x['p']), x['p'].ntilted()) for x in traj]
        if snapq:
            qs = []
            for x in traj[-snapq:]:
                rq = evaluate(x['p'], POLICY, Known(), best=0.0, seed=seed, k=k)
                pq = rq.get('p')
                qs.append((round(x['t'], 3), rq.get('s'), bool(pq.grid_lines(k)) if pq is not None else None))
            out['snapq'] = qs
    if crossed and p is not None:
        out['sq'] = p.sq.tolist(); out['s_full'] = p.s
    return out


# ---------- driver ----------
def run(a):
    st = Store()
    d = os.path.join(ROOT, a.name); os.makedirs(d, exist_ok=True)
    rng = random.Random(a.seed)
    tasks = []
    for n in a.n:
        k = kof(n)
        keys = {r['key'] for r in st.db.execute("SELECT key FROM basin WHERE n=? AND cls='certified'", (n,))}
        for cls in a.classes:
            if cls in ('grid', 'scratch'):
                for t in range(a.trials):
                    tasks.append((dict(name=f'{cls}{n}', cls=cls, n=n), t, keys))
            elif cls == 'neg':
                m = k * k - 4 if k >= 5 else None
                if m:
                    mk = {r['key'] for r in st.db.execute("SELECT key FROM basin WHERE n=? AND cls='certified'", (m,))}
                    for t in range(a.trials):
                        tasks.append((dict(name=f'neg{m}', cls='neg', n=m), t, mk))
            elif cls == 'above':
                pool = clean_above(st, n, a.dmin, seed=a.seed)
                print(f'n={n}: {len(pool)} clean above-grid starts (dmin {a.dmin})', flush=True)
                for t in range(a.trials):
                    if not pool:
                        break
                    pid, dist = pool[t % len(pool)]
                    p = st.get(pid)
                    tasks.append((dict(name=f'above{n}#{pid}', cls='above', n=n, s=p.s, sq=p.sq.tolist(), dist=dist, id=pid), t, keys))
            elif cls == 'pos':
                subs = list(st.query(n=n, below=k, status='certified', limit=None))
                for t in range(a.trials):
                    r = subs[rng.randrange(len(subs))]
                    p = st.get(r['id'])
                    tasks.append((dict(name=f'pos{n}#{r["id"]}', cls='pos', n=n, s=p.s, sq=p.sq.tolist(), id=r['id']), t, keys))
    jobs = []
    for task, t, keys in tasks:
        for spec in a.methods:
            if task['cls'] == 'scratch' and not spec.startswith('anneal'):
                continue
            jobs.append((task, spec, rng.randrange(1 << 30), a.snaps if spec.split(':')[0] in ('anneal', 'melt') else 0,
                         a.snapq, keys))
    print(f'{a.name}: {len(jobs)} trials', flush=True)
    t0 = time.time()
    with open(f'{d}/trials.jsonl', 'a') as f, ProcessPoolExecutor(a.procs) as ex:
        for r in ex.map(trial, jobs, chunksize=1):
            if r is not None:
                f.write(json.dumps(r) + '\n'); f.flush()
    print(f'done in {time.time() - t0:.0f}s', flush=True)
    report(a)


def wilson_hi(x, n, z=1.96):
    if n == 0:
        return 1.0
    p = x / n
    return min(1.0, (p + z * z / (2 * n) + z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / (1 + z * z / n))


def report(a):
    d = os.path.join(ROOT, a.name)
    R = [json.loads(l) for l in open(f'{d}/trials.jsonl')]
    G = collections.defaultdict(list)
    for r in R:
        G[(r['n'], r['cls'], r['method'])].append(r)
    print(f'{"n":>4} {"class":8s} {"method":32s} {"N":>4} {"cross":>5} {"95%hi":>6} {"known":>5} {"grid":>5} {"best s":>12} {"median s":>12} {"CPU s/trial":>11}')
    for (n, cls, m), rs in sorted(G.items()):
        ok = [r for r in rs if r.get('s') is not None]
        x = sum(r['crossed'] for r in rs)
        kn = sum(1 for r in rs if r['crossed'] and r.get('basin'))
        grid = sum(1 for r in ok if r['lines'])
        ss = sorted(r['s'] for r in ok)
        print(f'{n:4d} {cls:8s} {m:32s} {len(rs):4d} {x:5d} {wilson_hi(x, len(rs)):6.3f} {kn:5d} {grid:5d} '
              f'{(ss[0] if ss else float("nan")):12.7f} {(ss[len(ss) // 2] if ss else float("nan")):12.7f} '
              f'{sum(r["sec"] for r in rs) / len(rs):11.1f}')


def traj_report(a):
    """Per (class, method) with trajectories: mean side / line occupancy / tilted count at each snapshot time, lock-in time
    (first snapshot from which occupancy >= k for the rest of the run), and the quenched-snapshot sides (--snapq)."""
    d = os.path.join(ROOT, a.name)
    R = [json.loads(l) for l in open(f'{d}/trials.jsonl')]
    G = collections.defaultdict(list)
    for r in R:
        if r.get('traj'):
            G[(r['n'], r['cls'], r['method'])].append(r)
    for (n, cls, m), rs in sorted(G.items()):
        k = kof(n)
        T = collections.defaultdict(list)
        locks = []
        for r in rs:
            tr = r['traj']
            for t_, s_, occ, til in tr:
                T[t_].append((s_, occ, til))
            lk = None
            for i in range(len(tr)):
                if all(x[2] >= k for x in tr[i:]):
                    lk = tr[i][0]; break
            locks.append(lk)
        nl = [x for x in locks if x is not None]
        print(f'{n} {cls} {m}: {len(rs)} runs; locked into a full line during the anneal: {len(nl)}'
              + (f', median lock-in t {sorted(nl)[len(nl) // 2]:.2f}' if nl else ''))
        print('   t     s       occ/k  tilted')
        for t_ in sorted(T)[::max(1, len(T) // 10)]:
            v = T[t_]
            print(f'   {t_:.2f} {sum(x[0] for x in v) / len(v):8.4f} {sum(x[1] for x in v) / len(v) / k:6.2f} '
                  f'{sum(x[2] for x in v) / len(v):6.1f}')
        sq = [q for r in rs for q in r.get('snapq', [])]
        if sq:
            by = collections.defaultdict(list)
            for t_, s_, g in sq:
                if s_ is not None:
                    by[t_].append((s_, g))
            print('   quenched snapshots: ' + '  '.join(f't {t_:.2f}: grid {sum(g for _, g in v)}/{len(v)} best {min(x for x, _ in v):.5f}'
                                                     for t_, v in sorted(by.items())))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('run'); p.add_argument('name'); p.add_argument('--n', type=int, nargs='+', required=True)
    p.add_argument('--classes', nargs='+', default=['grid', 'above', 'scratch', 'pos', 'neg'])
    p.add_argument('--methods', nargs='+', required=True); p.add_argument('--trials', type=int, default=20)
    p.add_argument('--procs', type=int, default=16); p.add_argument('--seed', type=int, default=1)
    p.add_argument('--dmin', type=float, default=0.1); p.add_argument('--snaps', type=int, default=20)
    p.add_argument('--snapq', type=int, default=0)
    p = sp.add_parser('report'); p.add_argument('name')
    p = sp.add_parser('traj'); p.add_argument('name')
    a = ap.parse_args()
    {'run': run, 'report': report, 'traj': traj_report}[a.cmd](a)


if __name__ == '__main__':
    main()
