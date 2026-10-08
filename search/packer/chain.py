#!/usr/bin/env python3
"""Chaining across sizes (10-08; SQUISH-style lineage): a window of n, each starting from its frozen 10-05 register record
(../exact/batch/inputs: pre-SQUISH, so SQUISH's sides are ground truth for "how much does the chain recover").

Round for one n: an explorer run (P procs x T min) from best(n) plus removal seeds best(n + k) - k, k = 1, 2 (best of 4
random removals each, fq-quenched).  A new best at n marks n - 1 and n - 2 dirty (fresh seeds).  Scheduling: dirty n
first (most recently dirtied), otherwise the n with the best recent gain per round (ties: fewest rounds).  The top of the
window (n_hi + 1, n_hi + 2) stays at its register record (seed source only).

  chain.py --lo 123 --hi 131 --minutes 4 --procs-per 2 --slots 8 --hours 1 --out runs/ch1
  chain.py --report runs/ch1
"""
import argparse, json, math, os, random, shutil, signal, subprocess, sys, tempfile, time
from jobpool import run_jobs
import mcmin, hop

HERE = os.path.dirname(os.path.abspath(__file__))
TRIES, SEED_FIRST = 12, False                # removal tries per seed; v2 (10-08): non-shrinking removals rejected


def seed_job(j):
    src, k, n, out, seed = j
    s1, sq1 = mcmin.load_deg(src)
    rng = random.Random(seed)
    tmp = tempfile.mkdtemp()
    hop.EXTRA[:] = ['--loosen', '1.0']
    best = None
    for _ in range(TRIES):
        drop = set(rng.sample(range(len(sq1)), k))
        r = hop.quench(s1, [q for i, q in enumerate(sq1) if i not in drop], tmp)
        if not r or r[0] > s1 - 1e-7:          # did not shrink (removed a rattler / polished back): a copy, not a seed
            continue
        if best is None or r[0] < best[0]:
            best = (r[0], r[1])
        if SEED_FIRST and best:
            break
    if best:
        mcmin.write_deg(out, *best)
        return best[0]
    return None


def archive_best(d):
    best = None
    if os.path.exists(f'{d}/archive.jsonl'):
        for l in open(f'{d}/archive.jsonl'):
            e = json.loads(l)
            if best is None or e['s'] < best[0]:
                best = (e['s'], e['path'])
    return best


def main(a):
    os.makedirs(a.out, exist_ok=True)
    stp = f'{a.out}/state.json'
    T = hop.squish_targets()
    if os.path.exists(stp):
        S = json.load(open(stp))
    else:
        S = {'n': {}, 'log': [], 't0': time.time(), 'round': 0}
        src = a.src or hop.BATCH
        ns = json.load(open(a.ns_file)) if a.ns_file else list(range(a.lo, a.hi + 1))
        ours = dict(x.split(':') for x in (a.ours or []))
        for n in sorted(set(ns) | {m + k for m in ns for k in (1, 2)}):
            p = f'{a.out}/best_{n}.txt'
            sp = f'{src}/n-{n}.txt' if os.path.exists(f'{src}/n-{n}.txt') else f'{hop.BATCH}/n-{n}.txt'
            if not os.path.exists(sp):
                continue
            shutil.copy(sp, p)
            reg, _ = mcmin.load_deg(p)
            if str(n) in ours:
                so, _ = mcmin.load_deg(ours[str(n)])
                if so < reg - 1e-12:
                    shutil.copy(ours[str(n)], p)
            s, _ = mcmin.load_deg(p)
            S['n'][str(n)] = dict(best=s, reg=reg, start=s, squish=None if a.src else T.get(n), path=p, rounds=0,
                                  dirty=(n in ns), gains=[], dirtied=0.0, active=(n in ns))
    live = {}
    t_end = time.time() + 3600 * a.hours
    def pick():
        cand = [int(n) for n, v in S['n'].items() if v.get('active', int(n) <= (a.hi or 0)) and int(n) not in live.values()]
        dirty = [n for n in cand if S['n'][str(n)]['dirty']]
        if dirty:
            return max(dirty, key=lambda n: (S['n'][str(n)]['dirtied'], -n))
        if not cand:
            return None
        if a.first_pass and any(S['n'][str(n)]['rounds'] == 0 for n in cand):
            cand = [n for n in cand if S['n'][str(n)]['rounds'] == 0]
        return max(cand, key=lambda n: (sum(S['n'][str(n)]['gains'][-2:]), -S['n'][str(n)]['rounds']))
    while time.time() < t_end or live:
        while time.time() < t_end and len(live) < a.slots:
            n = pick()
            if n is None:
                break
            v = S['n'][str(n)]
            v['dirty'] = False
            S['round'] += 1
            rd = f'{a.out}/r{S["round"]:04d}_n{n}'
            os.makedirs(rd)
            seeds = [(S['n'][str(n + k)]['path'], k, n, f'{rd}/seed_k{k}.txt', S['round'] * 10 + k) for k in (1, 2) if str(n + k) in S['n']]
            res = run_jobs(seed_job, seeds, procs=2, timeout=1800)
            carried = []
            for q, cp in enumerate(v.get('carry', [])):         # lineage carry from the previous round(s)
                dst = f'{rd}/seed_c{q}.txt'; shutil.copy(cp, dst); carried.append(dst)
            starts = [v['path']] + carried + [sd[3] for sd, r in zip(seeds, res) if isinstance(r, float)]
            cmd = [sys.executable, os.path.join(HERE, 'explore.py'), '--n', str(n), '--starts', *starts, '--procs', str(a.procs_per),
                   '--minutes', str(a.minutes), '--seed', str(S['round']), '--out', rd, '--adapt', '--elite-share', '0.3',
                   '--frontier', f'{rd}/seed', '--frontier-share', str(a.frontier_share), '--bandit-state', f'{a.out}/bandit_{n}.json'] + (a.explore_extra.split() if a.explore_extra else [])
            p = subprocess.Popen(cmd, stdout=open(f'{rd}/out', 'w'), stderr=subprocess.STDOUT, start_new_session=True)
            live[p] = n
            S['log'].append(dict(t=round(time.time() - S['t0']), ev='start', n=n, round=S['round'], seeds=[r for r in res]))
        time.sleep(3)
        for p, n in list(live.items()):
            if p.poll() is None:
                continue
            del live[p]
            v = S['n'][str(n)]
            rd = p.args[p.args.index('--out') + 1]
            b = archive_best(rd)
            v['rounds'] += 1
            gain = 0.0
            if b and b[0] < v['best'] - 1e-9:
                gain = v['best'] - b[0]
                shutil.copy(b[1], v['path'])
                v['best'] = b[0]
                if b[0] < min(v['reg'], v.get('start', v['reg'])) - 1e-9:      # beyond the live register and our own start
                    cp = f'{a.out}/cand_{n}_{b[0]:.10f}.txt'; shutil.copy(b[1], cp)
                    print(f'CANDIDATE n={n} {b[0]:.12f} (register {v["reg"]:.12f}, start {v.get("start", v["reg"]):.12f}) -> {cp}', flush=True)
                for m in (n - 1, n - 2):
                    if str(m) in S['n'] and S['n'][str(m)].get('active', m >= (a.lo or 0)):
                        S['n'][str(m)]['dirty'] = True
                        S['n'][str(m)]['dirtied'] = time.time()
            # carry: best 3 distinct basins of the non-record lineages (roots = seed / carried starts), below the ceiling
            E = [json.loads(l) for l in open(f'{rd}/archive.jsonl')] if os.path.exists(f'{rd}/archive.jsonl') else []
            roots = {e['i'] for e in E if e['parent'] < 0 and e['kind'].startswith('start:seed')}
            k_n = math.ceil(math.sqrt(n) - 1e-12)
            lin = sorted((e for e in E if e.get('root', e['i']) in roots and e['s'] < k_n + 0.05), key=lambda e: e['s'])
            car, seen = [], []
            for e in lin:
                if all(abs(e['s'] - x) > 1e-7 for x in seen):
                    seen.append(e['s']); dst = f'{a.out}/carry_{n}_{len(car)}.txt'; shutil.copy(e['path'], dst + '.tmp'); car.append(dst)
                if len(car) == 3:
                    break
            for c in car: os.replace(c + '.tmp', c)
            v['carry'] = car
            lb = seen[0] if seen else None
            lgain = (v.get('lin_best', float('inf')) - lb) if lb is not None and v.get('lin_best') is not None else 0.0
            if lb is not None and (v.get('lin_best') is None or lb < v['lin_best']): v['lin_best'] = lb
            gain = gain + max(0.0, lgain) * 0.1        # scheduling credit for lineage progress too (record gains dominate)
            v['gains'].append(gain)
            S['log'].append(dict(t=round(time.time() - S['t0']), ev='done', n=n, best=v['best'], gain=gain,
                                 vs_reg=v['best'] - v['reg'], vs_squish=(v['best'] - v['squish']) if v['squish'] else None))
            print(f'{time.time() - S["t0"]:7.0f}s n={n} lineage {v.get("lin_best") or 0:.10f} best {v["best"]:.10f} gain {gain:.2e} vs reg {v["best"] - v["reg"]:+.2e}'
                  + (f' vs SQUISH {v["best"] - v["squish"]:+.2e}' if v['squish'] else ''), flush=True)
            json.dump(S, open(stp + '.tmp', 'w')); os.replace(stp + '.tmp', stp)
    json.dump(S, open(stp + '.tmp', 'w')); os.replace(stp + '.tmp', stp)
    report(a.out)


def report(out):
    S = json.load(open(f'{out}/state.json'))
    print(f'== {out}: {S["round"]} rounds')
    cands = [(n, v) for n, v in S['n'].items() if v['best'] < min(v['reg'], v.get('start', v['reg'])) - 1e-9]
    print(f'  candidates below register and start: {[(n, round(v["best"], 10), round(v["reg"] - v["best"], 10)) for n, v in cands]}')
    print(f'  {"n":>4} {"rounds":>6} {"register":>15} {"best":>15} {"gain":>9} {"SQUISH":>15} {"frac of SQUISH gain":>19} {"lineage best":>15}')
    for n, v in sorted(S['n'].items(), key=lambda q: int(q[0])):
        sq = v['squish']
        fr = (v['reg'] - v['best']) / (v['reg'] - sq) if sq and v['reg'] - sq > 1e-12 else None
        print(f'  {n:>4} {v["rounds"]:6d} {v["reg"]:15.10f} {v["best"]:15.10f} {v["reg"] - v["best"]:9.2e} '
              f'{(f"{sq:15.10f}") if sq else " " * 15} {(f"{fr:19.2f}") if fr is not None else " " * 19} {v.get("lin_best") or 0:15.10f}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--lo', type=int); ap.add_argument('--hi', type=int); ap.add_argument('--minutes', type=float, default=4)
    ap.add_argument('--procs-per', type=int, default=2); ap.add_argument('--slots', type=int, default=8)
    ap.add_argument('--hours', type=float, default=1); ap.add_argument('--frontier-share', type=float, default=0.7);
    ap.add_argument('--explore-extra', default='', help='extra explore.py args (e.g. "--melt 2")'); ap.add_argument('--first-pass', action='store_true', help='every n once before any repeat'); ap.add_argument('--ns-file', help='json list of n to work on (instead of --lo/--hi)'); ap.add_argument('--src', help='start packings dir (n-<n>.txt), e.g. ../exact/batch/inputs_live')
    ap.add_argument('--ours', nargs='*', help='n:path of our own better packings (used as start if better than src)'); ap.add_argument('--out'); ap.add_argument('--report')
    a = ap.parse_args()
    if a.report:
        report(a.report)
    else:
        main(a)
