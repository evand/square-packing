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


def seed_job(j):
    src, k, n, out, seed = j
    s1, sq1 = mcmin.load_deg(src)
    rng = random.Random(seed)
    tmp = tempfile.mkdtemp()
    hop.EXTRA[:] = ['--loosen', '1.0']
    best = None
    for _ in range(4):
        drop = set(rng.sample(range(len(sq1)), k))
        r = hop.quench(s1, [q for i, q in enumerate(sq1) if i not in drop], tmp)
        if r and (best is None or r[0] < best[0]):
            best = (r[0], r[1])
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
        for n in range(a.lo, a.hi + 3):
            p = f'{a.out}/best_{n}.txt'
            shutil.copy(f'{hop.BATCH}/n-{n}.txt', p)
            s, _ = mcmin.load_deg(p)
            S['n'][str(n)] = dict(best=s, reg=s, squish=T.get(n), path=p, rounds=0, dirty=(n <= a.hi), gains=[], dirtied=0.0)
    live = {}
    t_end = time.time() + 3600 * a.hours
    def pick():
        cand = [int(n) for n, v in S['n'].items() if int(n) <= a.hi and int(n) not in live.values()]
        dirty = [n for n in cand if S['n'][str(n)]['dirty']]
        if dirty:
            return max(dirty, key=lambda n: (S['n'][str(n)]['dirtied'], -n))
        if not cand:
            return None
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
            seeds = [(S['n'][str(n + k)]['path'], k, n, f'{rd}/seed_k{k}.txt', S['round'] * 10 + k) for k in (1, 2)]
            res = run_jobs(seed_job, seeds, procs=2, timeout=1800)
            starts = [v['path']] + [sd[3] for sd, r in zip(seeds, res) if isinstance(r, float)]
            cmd = [sys.executable, os.path.join(HERE, 'explore.py'), '--n', str(n), '--starts', *starts, '--procs', str(a.procs_per),
                   '--minutes', str(a.minutes), '--seed', str(S['round']), '--out', rd, '--adapt', '--elite-share', '0.3']
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
                for m in (n - 1, n - 2):
                    if str(m) in S['n'] and m >= a.lo:
                        S['n'][str(m)]['dirty'] = True
                        S['n'][str(m)]['dirtied'] = time.time()
            v['gains'].append(gain)
            S['log'].append(dict(t=round(time.time() - S['t0']), ev='done', n=n, best=v['best'], gain=gain,
                                 vs_reg=v['best'] - v['reg'], vs_squish=(v['best'] - v['squish']) if v['squish'] else None))
            print(f'{time.time() - S["t0"]:7.0f}s n={n} best {v["best"]:.10f} gain {gain:.2e} vs reg {v["best"] - v["reg"]:+.2e}'
                  + (f' vs SQUISH {v["best"] - v["squish"]:+.2e}' if v['squish'] else ''), flush=True)
            json.dump(S, open(stp + '.tmp', 'w')); os.replace(stp + '.tmp', stp)
    json.dump(S, open(stp + '.tmp', 'w')); os.replace(stp + '.tmp', stp)
    report(a.out)


def report(out):
    S = json.load(open(f'{out}/state.json'))
    print(f'== {out}: {S["round"]} rounds')
    print(f'  {"n":>4} {"rounds":>6} {"register":>15} {"best":>15} {"gain":>9} {"SQUISH":>15} {"frac of SQUISH gain":>19}')
    for n, v in sorted(S['n'].items(), key=lambda q: int(q[0])):
        sq = v['squish']
        fr = (v['reg'] - v['best']) / (v['reg'] - sq) if sq and v['reg'] - sq > 1e-12 else None
        print(f'  {n:>4} {v["rounds"]:6d} {v["reg"]:15.10f} {v["best"]:15.10f} {v["reg"] - v["best"]:9.2e} '
              f'{(f"{sq:15.10f}") if sq else " " * 15} {(f"{fr:19.2f}") if fr is not None else ""}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--lo', type=int); ap.add_argument('--hi', type=int); ap.add_argument('--minutes', type=float, default=4)
    ap.add_argument('--procs-per', type=int, default=2); ap.add_argument('--slots', type=int, default=8)
    ap.add_argument('--hours', type=float, default=1); ap.add_argument('--out'); ap.add_argument('--report')
    a = ap.parse_args()
    if a.report:
        report(a.report)
    else:
        main(a)
