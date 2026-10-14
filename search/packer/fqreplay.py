#!/usr/bin/env python3
"""fq replay equivalence benchmark (10-13 performance session; PACKER "Q0").

An fq-local change (build flags, LP settings, a faster kernel) is accepted by replaying a corpus of real explorer fq calls
under each arm and comparing per-call outputs, instead of a whole-search race.  Null for "same search, re-rolled": the
base binary on inputs jittered by ~1e-13 (arms `jitter=SEED`).  A change whose outputs differ from base no more than in
an unbiased way (no class-frequency shift) is a re-roll of the same search, so its speedup carries over to the search.

  fqreplay.py record --out runs/fqcorpus/c237 -- <explore.py args>      (runs explore.py with every fq call logged)
  fqreplay.py run --corpus runs/fqcorpus/c110 --sample 4000 --out runs/fqreplay/q0 \
        --arm base:BIN --arm j1:BIN:jitter=1 --arm noir:BIN:FQ_LP_IR=0
  fqreplay.py compare runs/fqreplay/q0 [--base base]
  fqreplay.py certify runs/fqreplay/q0 --arms base noir [--calls 60]   (exactsolve on polish outputs; diverging calls first)

Arm = name:binary[:KEY=VAL,...]; KEY jitter = jitter seed, any other KEY is an environment variable for fq.
Results: OUT/<arm>.jsonl, one line per call (id, cat, s, feasible, cpu, class, fp).  Calls of all arms are interleaved in
one pool (same load conditions for the CPU ratio); per-call CPU is the fq child's rusage.
"""
from __future__ import annotations
import argparse, json, math, os, random, subprocess, sys, tempfile, textwrap, zlib
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
BEST_VISIBLE = 10.9968869868          # bench110: best visible side of the h110_basinA holdout
SAME_SIDE = 2e-9                      # pk.identity.SAME_SIDE


def cat_of(args):
    s = ' '.join(args)
    if '--pit 8' in s:
        return 'screen'
    return 'polish+kick' if '--kick' in s else ('polish' if '--no-alm' in s else 'other')


# ---------------------------------------------------------------- record

WRAP = """#!/usr/bin/env bash
in=""; prev=""; for a in "$@"; do [ "$prev" = "--in" ] && in="$a"; prev="$a"; done
id=$$_$RANDOM$RANDOM
if [ "$1" = quench ] && [ -n "$in" ]; then cp "$in" {out}/in/$id.txt; fi
{fq} "$@"; rc=$?
if [ "$1" = quench ] && [ -n "$in" ]; then
  printf '%s\\n' "$id $*" >> {out}/calls.log
fi
exit $rc
"""


def record_cmd(a, rest):
    out = os.path.abspath(a.out)
    os.makedirs(f'{out}/in', exist_ok=True)
    wrap = f'{out}/fq-wrap.sh'
    open(wrap, 'w').write(WRAP.replace('{out}', out).replace('{fq}', os.path.join(HERE, 'target/release/fq')))
    os.chmod(wrap, 0o755)
    drive = (f"import sys; sys.path.insert(0, {HERE!r}); import hop, explore; hop.FQ = {wrap!r}\n"
             f"sys.argv = ['explore.py'] + {rest!r}\n"
             f"import textwrap; exec(textwrap.dedent(open({os.path.join(HERE, 'explore.py')!r}).read()"
             f".split(\"if __name__ == '__main__':\", 1)[1]), vars(explore))\n")
    subprocess.run([sys.executable, '-c', drive], check=True)
    with open(f'{out}/calls.jsonl', 'w') as f:
        for l in open(f'{out}/calls.log'):
            id_, *args = l.split()
            i = args.index('--out')
            f.write(json.dumps(dict(id=id_, args=args[i + 2:])) + '\n')
    print(sum(1 for _ in open(f'{out}/calls.jsonl')), 'calls recorded in', out)


# ---------------------------------------------------------------- run

def jitter_file(src, dst, seed):
    rng = random.Random(seed)
    L = open(src).read().split('\n')
    out = [L[0]]
    for l in L[1:]:
        f = l.split()
        if len(f) >= 3:
            x, y, t = (float(v) for v in f[:3])
            l = f'{x + rng.uniform(-1e-13, 1e-13)!r} {y + rng.uniform(-1e-13, 1e-13)!r} {t + rng.uniform(-1e-11, 1e-11)!r}'
        out.append(l)
    open(dst, 'w').write('\n'.join(out))


def classify(n, s, path, cat):
    from layout import full_lines
    k = math.ceil(math.sqrt(n) - 1e-12)
    if cat == 'screen' and s > k + 0.051:
        return 'discard'
    if full_lines(path, k):
        return 'grid'
    if s < k and n == 110 and s < BEST_VISIBLE - 1e-9:
        return 'below-best'
    return 'sub-k' if s < k else 'above-k'


def one(job):
    arm, binary, env, jit, corpus, c = job
    src = os.path.join(corpus, 'in', c['id'] + '.txt')
    with tempfile.TemporaryDirectory(prefix='fqr-') as tmp:
        inp = src
        if jit is not None:
            inp = os.path.join(tmp, 'in.txt')
            # per-call seed: the same jitter seed gives different (but reproducible) noise on each call
            jitter_file(src, inp, int(jit) * 1_000_003 + zlib.crc32(c['id'].encode()))
        outp = os.path.join(tmp, 'out.txt')
        p = subprocess.Popen([binary, 'quench', '--in', inp, '--out', outp, *c['args']], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, text=True, env={**os.environ, **env})
        so = p.stdout.read()
        _, status, ru = os.wait4(p.pid, 0)
        p.returncode = status
        r = dict(arm=arm, id=c['id'], cat=cat_of(c['args']), cpu=round(ru.ru_utime + ru.ru_stime, 4))
        try:
            d = json.loads(so)
            n = int(open(src).readline().split()[0])
            r.update(s=d['s'], feasible=d['min_gap'] >= 0 and d['min_wall'] >= 0, fp=d.get('fp'),
                     slp_it=d.get('slp_it'), cls=classify(n, d['s'], outp, r['cat']))
        except Exception as e:
            r.update(s=None, feasible=False, cls='fail', err=str(e)[:200])
    return r


def run_cmd(a):
    from pk import manifest
    calls = [json.loads(l) for l in open(os.path.join(a.corpus, 'calls.jsonl'))]
    calls = [c for c in calls if cat_of(c['args']) in ('screen', 'polish')]
    random.Random(a.seed).shuffle(calls)
    if a.stage != 'all':                         # --stage screen / polish: one stage only (polish = --no-alm calls)
        calls = [c for c in calls if ('--no-alm' in c['args']) == (a.stage == 'polish')]
    calls = calls[a.skip:a.skip + a.sample]     # --skip: disjoint slices of one shuffle (dev vs held-out samples)
    arms = []
    for spec in a.arm:
        name, binary, *kv = spec.split(':')
        env, jit = {}, None
        for item in (kv[0].split(',') if kv else []):
            k, v = item.split('=', 1)
            if k == 'jitter':
                jit = v
            else:
                env[k] = v
        arms.append((name, os.path.abspath(binary), env, jit))
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'fqreplay', args=vars(a), procs=a.procs)
    jobs = [(name, b, env, jit, a.corpus, c) for c in calls for (name, b, env, jit) in arms]   # interleaved
    fh = {name: open(os.path.join(a.out, f'{name}.jsonl'), 'w') for name, *_ in arms}
    with ProcessPoolExecutor(a.procs) as ex:
        for i, r in enumerate(ex.map(one, jobs, chunksize=4)):
            fh[r['arm']].write(json.dumps(r) + '\n')
            if (i + 1) % 2000 == 0:
                print(f'{i + 1} / {len(jobs)}', flush=True)
    for f in fh.values():
        f.close()
    json.dump(dict(corpus=a.corpus, sample=a.sample, skip=a.skip, stage=a.stage, seed=a.seed, arms=[dict(name=x[0], binary=x[1], env=x[2], jitter=x[3])
              for x in arms]), open(os.path.join(a.out, 'arms.json'), 'w'), indent=1)
    man.finish()


# ---------------------------------------------------------------- compare

def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


def sign_p(lo, hi):
    """Two-sided exact binomial sign test."""
    n, k = lo + hi, min(lo, hi)
    if n == 0:
        return 1.0
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n)


def compare_cmd(a):
    arms = json.load(open(os.path.join(a.dir, 'arms.json')))['arms']
    R = {x['name']: {r['id']: r for r in map(json.loads, open(os.path.join(a.dir, x['name'] + '.jsonl')))} for x in arms}
    B = R[a.base]
    for cat in ('screen', 'polish'):
        ids = [i for i, r in B.items() if r['cat'] == cat]
        if not ids:
            continue
        print(f'\n== {cat}: {len(ids)} calls')
        print(f"{'arm':8s} {'cpu/call':>8s} {'speed':>6s} {'fail':>5s} {'diverge':>8s} {'lower':>6s} {'higher':>6s} {'sign p':>7s}  class frequencies vs base (diff, 95% paired CI)")
        cb = sum(B[i]['cpu'] for i in ids)
        classes = sorted({r['cls'] for X in R.values() for i, r in X.items() if r['cat'] == cat})
        for name, X in R.items():
            got = [i for i in ids if i in X]
            cx = sum(X[i]['cpu'] for i in got)
            fail = sum(not X[i].get('feasible') for i in got)
            dv = [(X[i]['s'] - B[i]['s']) for i in got if X[i]['s'] is not None and B[i]['s'] is not None]
            div = [d for d in dv if abs(d) > SAME_SIDE]
            lo, hi = sum(d < 0 for d in div), sum(d > 0 for d in div)
            parts = []
            for c in classes:
                b_ = sum(B[i]['cls'] == c for i in got); x_ = sum(X[i]['cls'] == c for i in got)
                # paired difference CI: discordant pairs (McNemar-style normal approximation)
                n10 = sum(X[i]['cls'] == c and B[i]['cls'] != c for i in got); n01 = sum(X[i]['cls'] != c and B[i]['cls'] == c for i in got)
                N = len(got); diff = (n10 - n01) / N
                se = math.sqrt(max(n10 + n01 - (n10 - n01) ** 2 / N, 0)) / N
                parts.append(f'{c} {x_ / N:.3f} ({diff:+.4f} [{diff - 1.96 * se:+.4f},{diff + 1.96 * se:+.4f}])')
            print(f'{name:8s} {cx / len(got):8.3f} {cb / cx:6.2f} {fail:5d} {len(div) / len(dv):8.3f} {lo:6d} {hi:6d} {sign_p(lo, hi):7.3f}  ' + '; '.join(parts))
        # magnitude of divergence
        for name, X in R.items():
            if name == a.base:
                continue
            dv = sorted(abs(X[i]['s'] - B[i]['s']) for i in ids if X[i].get('s') is not None and B[i].get('s') is not None)
            q = lambda f: dv[min(len(dv) - 1, int(f * len(dv)))]
            print(f'   |ds| {name:8s} q50 {q(.5):.1e} q90 {q(.9):.1e} q99 {q(.99):.1e} max {dv[-1]:.1e}')


# ---------------------------------------------------------------- certify (M4: does the change keep outputs certifiable?)

EXACT = os.path.join(HERE, '..', 'exact')


def cert_one(job):
    arm, binary, env, corpus, c, dst = job
    dst = os.path.abspath(dst); os.makedirs(dst, exist_ok=True)
    outp = os.path.join(dst, 'fq.txt')
    d = json.loads(subprocess.run([binary, 'quench', '--in', os.path.join(corpus, 'in', c['id'] + '.txt'), '--out', outp,
                                   *c['args']], capture_output=True, text=True, env={**os.environ, **env}).stdout)
    subprocess.run([os.path.join(EXACT, 'env.sh'), 'python3', os.path.join(EXACT, 'exactsolve.py'), outp, '--out', dst, '-q'],
                   cwd=EXACT, capture_output=True, text=True, timeout=3600)
    J = os.path.join(dst, 'fq.json')
    J = json.load(open(J)) if os.path.exists(J) else {}
    S = J.get('S_exact')
    return dict(arm=arm, id=c['id'], s=d['s'], valid=bool(J.get('cert_valid')),
                prec=abs(float(S[:40]) - d['s']) if S else None)


def certify_cmd(a):
    meta = json.load(open(os.path.join(a.dir, 'arms.json')))
    arms = [x for x in meta['arms'] if x['name'] in a.arms]
    R = {x['name']: {r['id']: r for r in map(json.loads, open(os.path.join(a.dir, x['name'] + '.jsonl')))} for x in arms}
    ok = lambda r: r['cat'] == 'polish' and r['cls'] in ('sub-k', 'below-best')
    ids = [i for i in R[a.arms[0]] if all(i in R[x] and ok(R[x][i]) for x in a.arms)]
    rng = random.Random(a.seed)
    div = [i for i in ids if max(R[x][i]['s'] for x in a.arms) - min(R[x][i]['s'] for x in a.arms) > SAME_SIDE]
    rest = [i for i in ids if i not in set(div)]
    rng.shuffle(div); rng.shuffle(rest)
    pick = (div + rest)[:a.calls]               # diverging calls first: they are where a change could hurt
    calls = {c['id']: c for c in map(json.loads, open(os.path.join(meta['corpus'], 'calls.jsonl')))}
    jobs = [(x['name'], x['binary'], x['env'], meta['corpus'], calls[i], os.path.join(a.dir, 'cert', x['name'], i))
            for i in pick for x in arms]
    with ProcessPoolExecutor(a.procs) as ex:
        res = list(ex.map(cert_one, jobs))
    json.dump(res, open(os.path.join(a.dir, 'cert.json'), 'w'))
    print(f'{len(pick)} calls ({min(len(div), a.calls)} diverging between arms)')
    for x in a.arms:
        r = [q for q in res if q['arm'] == x]
        p = sorted(q['prec'] for q in r if q['prec'] is not None)
        print(f"{x:8s} certified {sum(q['valid'] for q in r)}/{len(r)}  |S_exact - s_fq| median {p[len(p) // 2]:.1e} "
              f"q90 {p[int(.9 * len(p))]:.1e} max {p[-1]:.1e}" if p else f'{x}: no exact results')


if __name__ == '__main__':
    argv = sys.argv[1:]
    rest = []
    if '--' in argv:
        i = argv.index('--'); argv, rest = argv[:i], argv[i + 1:]
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('record'); p.add_argument('--out', required=True)
    p = sub.add_parser('run'); p.add_argument('--corpus', required=True); p.add_argument('--out', required=True)
    p.add_argument('--arm', action='append', required=True); p.add_argument('--sample', type=int, default=4000)
    p.add_argument('--skip', type=int, default=0)
    p.add_argument('--stage', choices=['all', 'screen', 'polish'], default='all')
    p.add_argument('--seed', type=int, default=1); p.add_argument('--procs', type=int, default=16)
    p = sub.add_parser('compare'); p.add_argument('dir'); p.add_argument('--base', default='base')
    p = sub.add_parser('certify'); p.add_argument('dir'); p.add_argument('--arms', nargs='+', default=['base', 'noir'])
    p.add_argument('--calls', type=int, default=60); p.add_argument('--seed', type=int, default=1); p.add_argument('--procs', type=int, default=16)
    a = ap.parse_args(argv)
    {'record': lambda: record_cmd(a, rest), 'run': lambda: run_cmd(a), 'compare': lambda: compare_cmd(a),
     'certify': lambda: certify_cmd(a)}[a.cmd]()
