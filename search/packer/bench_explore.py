#!/usr/bin/env python3
"""Explorer benchmark suite (10-08): fixed tasks with known better answers, R replicates of small explorer runs per task,
paired comparisons between variants.  Built because single runs could not rank variants (replicate spread ~10x).

Task = (start packing, target side).  Two sources, both frozen under bench/:
  * SQUISH seeds (bench/seeds/, bench_search.py): s108, s179, r126, r237; target = SQUISH's side.
  * rd1 pairs (rediscover.py, runs/rd1): older record (09-30 site mirror) -> newer register record; the start is the
    mirror packing polished by fq (--no-alm --loosen 1.0), frozen in bench/tasks/rd<n>.txt; target = register side.
    Picked for a spread of difficulty under rd1's kick protocol (hits / best fraction of the gap there):
    106 1/38, 177 4/40, 236 3/40, 130 0/40 (0.73), 263 (0.97), 307 (0.54), 237 (0.43), 210 (0.29), 271 (0.21),
    131, 180, 305 (~0: far pairs).
Replicate = one `explore.py` run (P procs x T min, its own seed), all variants interleaved so they share machine load.
Score per replicate: frac = (start - best) / (start - target) clipped to [0, 1] at checkpoints (fractions of T); score = mean
of the checkpoint fracs (area under the curve: speed counts); hit = best <= target + 1e-7.
Summary: per task median frac per variant; paired difference (B - A) per task, mean over tasks with a bootstrap CI over
replicates.  Tiers: quick (8 tasks), full (16 tasks).

  bench_explore.py --make-tasks
  bench_explore.py --variants base elite0 --tier quick --reps 2 --procs-per 2 --minutes 8 --total-procs 16 --out runs/be_x
  bench_explore.py --summary runs/be_x [--ref base]
"""
import argparse, json, math, os, random, signal, statistics as st, subprocess, sys, tempfile, time
import mcmin, hop

HERE = os.path.dirname(os.path.abspath(__file__))
TASKS = os.path.join(HERE, 'bench/tasks')
SITE = os.path.join(HERE, '../../site/www/data/p')
RD_N = (106, 177, 236, 130, 263, 307, 237, 210, 271, 131, 180, 305)
SQUISH = {'s108': (108, 10.909940073445), 's179': (179, 13.883795490512), 'r126': (126, 11.773303606607),
          'r237': (237, 15.903676235191)}
TIERS = {'quick': ['s108', 's179', 'rd106', 'rd130', 'rd177', 'rd210', 'rd236', 'rd263'],
         'full': ['s108', 's179', 'r126', 'r237'] + [f'rd{n}' for n in RD_N]}
CHECK = (0.125, 0.25, 0.5, 1.0)

# Variants: extra explorer arguments.  base = the record-hunting recipe of the tx runs (adaptive kinds, elite 0.3).
VARIANTS = {
    'base': ['--adapt', '--elite-share', '0.3'],
    'elite0': ['--adapt'],
    'fixed': ['--elite-share', '0.3'],
    'pm3': ['--adapt', '--elite-share', '0.3', '--polish-margin', '3e-3'],
}


def make_tasks():
    os.makedirs(TASKS, exist_ok=True)
    rd = {p['n']: p for p in json.load(open(os.path.join(HERE, 'runs/rd1/rd.json')))['pairs'] if 'mirror' in p['label']}
    man = {}
    for name, (n, tgt) in SQUISH.items():
        s, _ = mcmin.load_deg(f'{HERE}/bench/seeds/{name}.txt')
        man[name] = dict(n=n, path=f'bench/seeds/{name}.txt', start=s, target=tgt, source='SQUISH seed (bench_search.py)')
    hop.EXTRA[:] = ['--loosen', '1.0']
    tmp = tempfile.mkdtemp()
    for n in RD_N:
        d = json.load(open(f'{SITE}/square-{n}.json'))
        s0, sq0 = float(d['s']), [(float(x), float(y), float(a)) for x, y, a in d['squares']]
        r = hop.quench(s0, sq0, tmp, extra=('--no-alm',))
        assert r, n
        s, sq = r[0], r[1]
        p = f'bench/tasks/rd{n}.txt'
        mcmin.write_deg(os.path.join(HERE, p), s, sq)
        man[f'rd{n}'] = dict(n=n, path=p, start=s, target=rd[n]['s_new'], mirror=s0, source='rd1: 09-30 mirror -> register 10-05')
        print(f'rd{n}: mirror {s0:.10f} polished {s:.10f} target {rd[n]["s_new"]:.10f} gap {s - rd[n]["s_new"]:.2e}', flush=True)
    json.dump(man, open(f'{TASKS}/tasks.json', 'w'), indent=1)


def run(a):
    man = json.load(open(f'{TASKS}/tasks.json'))
    tasks = a.tasks or TIERS[a.tier]
    os.makedirs(a.out, exist_ok=True)
    json.dump(dict(argv=sys.argv, variants={v: VARIANTS[v] for v in a.variants}, tasks=tasks, reps=a.reps,
                   procs_per=a.procs_per, minutes=a.minutes), open(f'{a.out}/meta.json', 'w'), indent=1)
    # big n first (longest polishes), variants interleaved innermost so A and B run side by side
    jobs = [(t, r, v) for r in range(a.reps) for t in sorted(tasks, key=lambda t: -man[t]['n']) for v in a.variants]
    jobs = [j for j in jobs if not os.path.exists(f'{a.out}/{j[2]}/{j[0]}_r{j[1]}/done')]
    slots = max(1, a.total_procs // a.procs_per)
    live = {}
    print(f'{len(jobs)} replicates, {slots} at a time, ~{math.ceil(len(jobs) / slots) * a.minutes:.0f} min', flush=True)
    while jobs or live:
        while jobs and len(live) < slots:
            t, r, v = jobs.pop(0)
            d = f'{a.out}/{v}/{t}_r{r}'
            os.makedirs(d, exist_ok=True)
            for f in ('archive.jsonl', 'proposals.jsonl'):
                if os.path.exists(f'{d}/{f}'): os.remove(f'{d}/{f}')
            seed = 1000 * r + sum(map(ord, t))          # same seed for every variant: paired replicates
            cmd = [sys.executable, os.path.join(HERE, 'explore.py'), '--n', str(man[t]['n']), '--starts',
                   os.path.join(HERE, man[t]['path']), '--procs', str(a.procs_per), '--minutes', str(a.minutes),
                   '--seed', str(seed), '--out', d] + VARIANTS[v]
            p = subprocess.Popen(cmd, stdout=open(f'{d}/out', 'w'), stderr=subprocess.STDOUT, start_new_session=True)
            live[p] = (d, time.time())
        time.sleep(2)
        for p, (d, t0) in list(live.items()):
            if p.poll() is not None:
                open(f'{d}/done', 'w').write(str(p.returncode)); del live[p]
            elif time.time() - t0 > 60 * a.minutes + 900:           # hung: kill the whole session (pool workers, fq)
                os.killpg(p.pid, signal.SIGKILL); p.wait()
                open(f'{d}/done', 'w').write('killed'); del live[p]
    summary([a.out], a.ref)


def curve(d, minutes):
    E = [json.loads(l) for l in open(f'{d}/archive.jsonl')] if os.path.exists(f'{d}/archive.jsonl') else []
    st0 = [e['s'] for e in E if e['parent'] < 0]
    if not st0:
        return None
    return [min(e['s'] for e in E if e['t'] <= c * 60 * minutes or e['parent'] < 0) for c in CHECK], len(E)


def load(out):
    meta = json.load(open(f'{out}/meta.json'))
    man = json.load(open(f'{TASKS}/tasks.json'))
    R = {}                                                  # R[v][t][r] = (fracs at CHECK, hit, nbasins)
    for v in meta['variants']:
        for t in meta['tasks']:
            for r in range(meta['reps']):
                d = f'{out}/{v}/{t}_r{r}'
                c = curve(d, meta['minutes']) if os.path.exists(f'{d}/done') else None
                if c is None:
                    continue
                m = man[t]; gap = m['start'] - m['target']
                R.setdefault(v, {}).setdefault(t, {})[r] = ([min(1.0, max(0.0, (m['start'] - b) / gap)) for b in c[0]],
                                                            c[0][-1] <= m['target'] + 1e-7, c[1], c[0][-1])
    return meta, R


def summary(outs, ref=None):
    for out in outs:
        meta, R = load(out)
        vs = list(R)
        print(f'== {out}: {meta["procs_per"]} procs x {meta["minutes"]} min, reps {meta["reps"]}; frac = gain / gap '
              f'(median over reps) at {", ".join(f"{c:g}" for c in CHECK)} of T; hits; basins')
        for t in meta['tasks']:
            row = []
            for v in vs:
                X = list(R.get(v, {}).get(t, {}).values())
                if not X:
                    row.append(f'{v}: -'); continue
                f = [st.median(x[0][i] for x in X) for i in range(len(CHECK))]
                row.append(f'{v}: ' + '/'.join(f'{y:.2f}' for y in f) + f' h{sum(x[1] for x in X)}/{len(X)} '
                           f'b{st.median(x[2] for x in X):.0f}')
            print(f'  {t:6} ' + ' | '.join(row))
        rf = ref or vs[0]
        for v in vs:
            if v == rf:
                continue
            ds = paired(R[rf], R[v])
            if ds is None:
                continue
            m, lo, hi, nt = ds
            print(f'  {v} - {rf}: mean over {nt} tasks of (median score) difference, score = mean clipped frac over checkpoints, {m:+.3f}  [95% bootstrap {lo:+.3f}, {hi:+.3f}]')


def paired(A, B, nboot=4000):
    ts = [t for t in A if t in B and A[t] and B[t]]
    if not ts:
        return None
    rng = random.Random(0)
    def stat(pick):
        return st.mean(st.median(pick(B[t])) - st.median(pick(A[t])) for t in ts)
    auc = lambda x: st.mean(x[0])                         # mean clipped frac over the checkpoints: rewards speed too
    m = stat(lambda X: [auc(x) for x in X.values()])
    bs = sorted(stat(lambda X: [auc(rng.choice(list(X.values()))) for _ in X]) for _ in range(nboot))
    return m, bs[int(0.025 * nboot)], bs[int(0.975 * nboot)], len(ts)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--make-tasks', action='store_true'); ap.add_argument('--summary', nargs='+')
    ap.add_argument('--variants', nargs='+', default=['base']); ap.add_argument('--tier', default='quick', choices=TIERS)
    ap.add_argument('--tasks', nargs='+'); ap.add_argument('--reps', type=int, default=2)
    ap.add_argument('--procs-per', type=int, default=2); ap.add_argument('--minutes', type=float, default=8)
    ap.add_argument('--total-procs', type=int, default=16); ap.add_argument('--out'); ap.add_argument('--ref')
    a = ap.parse_args()
    if a.make_tasks:
        make_tasks()
    elif a.summary:
        summary(a.summary, a.ref)
    else:
        run(a)
