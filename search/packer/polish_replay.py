#!/usr/bin/env python3
"""Polish replay (10-08): a frozen set of screened states (what the explorer hands to the full polish), replayed through
polish variants.  Deterministic and paired: the tier-0 check for polish changes (same final basin? side? time?),
before any search-level benchmark.

  polish_replay.py make --archive runs/ex10 --n 110 --count 200 --out bench/polish110
  polish_replay.py run bench/polish110 --tag base [--fq target/release/fq] [-- extra fq args]
  polish_replay.py finish bench/polish110 --ks 2 4 8 16 32 --tag fin    # prototype: stop at k, exactsolve (Newton) finish
  polish_replay.py compare bench/polish110 base other
"""
import argparse, json, math, os, random, statistics as st, subprocess, sys, tempfile, time
from jobpool import run_jobs
import mcmin, explore, polish_trace

HERE = os.path.dirname(os.path.abspath(__file__))
EXACT = os.path.join(HERE, '../exact/exactsolve.py')
FQ = os.path.join(HERE, 'target/release/fq')


def make_job(j):
    path, kind, seed, k, smax = j
    from layout import full_lines
    rng = random.Random(seed)
    s, sq = mcmin.load_deg(path)
    prop, _ = explore.MOVES[kind](s, sq, rng, explore.A)
    tmp = tempfile.mkdtemp()
    loosen = rng.choice(('1.0', '1.02', '1.05'))
    r1, d1 = polish_trace.fq(s, prop, tmp, ['--loosen', loosen, '--pit', '8', '--flip-top', '0'])
    if r1 is None or r1[0] > smax + 1e-3:
        return None
    p1 = os.path.join(tmp, 's1.txt'); mcmin.write_deg(p1, *r1)
    if full_lines(p1, k):
        return None
    return dict(parent=path, kind=kind, loosen=loosen, s=r1[0], sq=r1[1])


def make(a):
    k = math.ceil(math.sqrt(a.n) - 1e-12); smax = k + 0.05
    E = [json.loads(l) for l in open(f'{a.archive}/archive.jsonl')]
    E = [e for e in E if e['s'] < smax and os.path.exists(e['path'])]
    rng = random.Random(a.seed)
    kinds, wts = zip(*explore.WEIGHTS.items())
    os.makedirs(a.out, exist_ok=True)
    got = []
    while len(got) < a.count:
        jobs = [(rng.choice(E)['path'], rng.choices(kinds, weights=wts)[0], rng.randrange(1 << 30), k, smax)
                for _ in range(int(1.6 * (a.count - len(got))) + 4)]
        got += [r for r in run_jobs(make_job, jobs, procs=a.procs, timeout=900) if isinstance(r, dict)]
    meta = []
    for i, r in enumerate(got[:a.count]):
        mcmin.write_deg(f'{a.out}/p{i:04d}.txt', r['s'], r['sq'])
        meta.append(dict(i=i, parent=r['parent'], kind=r['kind'], loosen=r['loosen'], s_screen=r['s']))
    json.dump(dict(n=a.n, archive=a.archive, items=meta), open(f'{a.out}/set.json', 'w'), indent=0)
    print(a.out, len(meta), 'states')


def fq_run(fq, path, args, outp):
    t0 = time.time()
    r = subprocess.run([fq, 'quench', '--in', path, '--out', outp, *args], capture_output=True, text=True, timeout=1800)
    d = json.loads(r.stdout)
    d['wall'] = time.time() - t0
    d.pop('trace', None); d.pop('fps', None)
    return d


def run_job(j):
    fq, path, args, outp = j
    return fq_run(fq, path, args, outp)


def run(a, extra):
    S = json.load(open(f'{a.set}/set.json'))
    od = f'{a.set}/out_{a.tag}'; os.makedirs(od, exist_ok=True)
    args = ['--loosen', '1.0', '--no-alm'] + extra
    jobs = [(a.fq, f'{a.set}/p{it["i"]:04d}.txt', args, f'{od}/p{it["i"]:04d}.txt') for it in S['items']]
    res = run_jobs(run_job, jobs, procs=a.procs, timeout=2400)
    json.dump(dict(args=args, fq=a.fq, res=[r if isinstance(r, dict) else None for r in res]), open(f'{od}/res.json', 'w'))
    print(od, 'done')


def finish_job(j):
    fq, path, k, od, xa = j
    tmp = tempfile.mkdtemp()
    sk = f'{tmp}/k.txt'
    d = fq_run(fq, path, ['--loosen', '1.0', '--no-alm', '--pit', str(k)], sk)
    t0 = time.time()
    r = subprocess.run([sys.executable, EXACT, sk, '--dps', '30', '--eps', '1e-12', '--out', tmp, '-q'] + list(xa), capture_output=True,
                       text=True, timeout=1200)
    try:
        e = json.load(open(f'{tmp}/k.json'))
    except Exception:
        e = {'status': 'crash', 'err': r.stderr[-300:]}
    so = e.get('second_order')
    return dict(k=k, s_k=d['s'], it=d['slp_it'], wall_k=d['wall'], es_sec=time.time() - t0, status=e.get('status', 'ok'),
                S=float(e['S_exact']) if e.get('S_exact') else None, newton=e.get('newton_converged'),
                cert=e.get('cert_valid'), lam_min=(e.get('lambda_A_maxmin') or [None, None])[-1] if isinstance(e.get('lambda_A_maxmin'), list) else e.get('lambda_A_maxmin'),
                so=so if isinstance(so, (str, int, float, bool)) or so is None else str(so)[:120])


def finish(a):
    S = json.load(open(f'{a.set}/set.json'))
    jobs = [(a.fq, f'{a.set}/p{it["i"]:04d}.txt', k, None, ['--force-kkt'] if a.force_kkt else []) for it in S['items'] for k in a.ks]
    res = run_jobs(finish_job, jobs, procs=a.procs, timeout=2400)
    out = {}
    for (fq, p, k, _, _), r in zip(jobs, res):
        out.setdefault(os.path.basename(p), []).append(r if isinstance(r, dict) else {'k': k, 'status': 'timeout'})
    od = f'{a.set}/out_{a.tag}'; os.makedirs(od, exist_ok=True)
    json.dump(out, open(f'{od}/finish.json', 'w'))
    print(od, 'done')


def compare(a):
    S = json.load(open(f'{a.set}/set.json'))
    A = json.load(open(f'{a.set}/out_{a.a}/res.json'))['res']
    B = json.load(open(f'{a.set}/out_{a.b}/res.json'))['res']
    ok = [(x, y) for x, y in zip(A, B) if x and y]
    ds = [y['s'] - x['s'] for x, y in ok]
    same = sum(abs(d) < 1e-9 for d in ds); lower = sum(d < -1e-9 for d in ds); higher = sum(d > 1e-9 for d in ds)
    ta = sum(x['wall'] for x, _ in ok); tb = sum(y['wall'] for _, y in ok)
    print(f'{a.set}: {a.b} vs {a.a}: {len(ok)} states; same side (1e-9) {same}, lower {lower}, higher {higher}; '
          f'ds p10/50/90 {sorted(ds)[len(ds) // 10]:.1e}/{st.median(ds):.1e}/{sorted(ds)[9 * len(ds) // 10]:.1e}; '
          f'time {tb:.0f} s vs {ta:.0f} s ({tb / ta:.2f}x); iterations median {st.median(y["slp_it"] for _, y in ok)} vs '
          f'{st.median(x["slp_it"] for x, _ in ok)}')


if __name__ == '__main__':
    argv = sys.argv[1:]
    extra = argv[argv.index('--') + 1:] if '--' in argv else []
    argv = argv[:argv.index('--')] if '--' in argv else argv
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    m = sp.add_parser('make'); m.add_argument('--archive'); m.add_argument('--n', type=int); m.add_argument('--count', type=int, default=200)
    m.add_argument('--out'); m.add_argument('--seed', type=int, default=11); m.add_argument('--procs', type=int, default=16)
    r = sp.add_parser('run'); r.add_argument('set'); r.add_argument('--tag'); r.add_argument('--fq', default=FQ); r.add_argument('--procs', type=int, default=16)
    f = sp.add_parser('finish'); f.add_argument('set'); f.add_argument('--tag', default='fin'); f.add_argument('--fq', default=FQ)
    f.add_argument('--ks', type=int, nargs='+', default=[2, 4, 8, 16, 32]); f.add_argument('--force-kkt', action='store_true'); f.add_argument('--procs', type=int, default=16)
    c = sp.add_parser('compare'); c.add_argument('set'); c.add_argument('a'); c.add_argument('b')
    a = ap.parse_args(argv)
    {'make': make, 'run': lambda a: run(a, extra), 'finish': finish, 'compare': compare}[a.cmd](a)
