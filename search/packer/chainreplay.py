#!/usr/bin/env python3
# Screen -> polish chain replay (perf session 10-13, PACKER "perfbench noise"): run each recorded screen under arms
# off (FQ_FASTSC=0) / on (defaults) / jit (off + 1e-13 input jitter), then polish every output with deployed defaults;
# reports paired polish and whole-chain CPU.  Usage: SCR=<tmp dir> [PROCS=16] chainreplay.py runs/fqcorpus/c237c N [SKIP].
# Edit ARMS for other knobs.  Needs >= ~1,000 chains to resolve ~5 % (polish cost per chain is heavy-tailed).
"""screen -> polish chain replay: does the screen variant change the cost of the polish that follows it?"""
import json, os, random, subprocess, sys, resource, concurrent.futures as cf, statistics as st, math
B = 'target/release/fq'; corpus = sys.argv[1]; N = int(sys.argv[2]); skip = int(sys.argv[3]) if len(sys.argv) > 3 else 0
S = os.environ['SCR']
calls = [json.loads(l) for l in open(f'{corpus}/calls.jsonl')]
random.Random(11).shuffle(calls)
scr = [c for c in calls if '--no-alm' not in c['args']][skip:skip + N]
ARMS = {'off': {'FQ_FASTSC': '0'}, 'on': {}, 'jit': {'FQ_FASTSC': '0', 'JIT': '1'}}
def jitter(src, dst, seed):
    rng = random.Random(seed); L = open(src).read().split('\n'); out = [L[0]]
    for l in L[1:]:
        f = l.split()
        if len(f) >= 3: x, y, t = map(float, f[:3]); l = f'{x + rng.uniform(-1e-13, 1e-13)!r} {y + rng.uniform(-1e-13, 1e-13)!r} {t + rng.uniform(-1e-11, 1e-11)!r}'
        out.append(l)
    open(dst, 'w').write('\n'.join(out))
def run(args, env):
    r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    p = subprocess.run([B, 'quench', *args], env={**os.environ, **env}, capture_output=True, text=True)
    r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    return r1.ru_utime - r0.ru_utime + r1.ru_stime - r0.ru_stime, p.stdout
def one(job):
    c, arm = job
    env = {k: v for k, v in ARMS[arm].items() if k != 'JIT'}
    inp = f"{corpus}/in/{c['id']}.txt"
    if 'JIT' in ARMS[arm]: j = f"{S}/ch_{c['id']}_jin.txt"; jitter(inp, j, hash(c['id']) & 0xffff); inp = j
    so = f"{S}/ch_{c['id']}_{arm}_s.txt"; po = f"{S}/ch_{c['id']}_{arm}_p.txt"
    t1, out1 = run(['--in', inp, '--out', so, *c['args']], env)
    try: s1 = json.loads(out1)['s']
    except Exception: s1 = None
    if not os.path.exists(so): return (c['id'], arm, t1, None, s1, None)
    t2, out2 = run(['--in', so, '--out', po, '--loosen', '1.0', '--no-alm'], {})     # polish always with deployed defaults
    try: s2 = json.loads(out2)['s']
    except Exception: s2 = None
    return (c['id'], arm, t1, t2, s1, s2)
# each job runs in its own thread: RUSAGE_CHILDREN deltas are per process, so use a process pool
with cf.ProcessPoolExecutor(int(os.environ.get('PROCS', '16'))) as ex:
    R = list(ex.map(one, [(c, a) for c in scr for a in ARMS]))
by = {}
for i, a, t1, t2, s1, s2 in R: by.setdefault(i, {})[a] = (t1, t2, s1, s2)
ok = [i for i, v in by.items() if all(v[a][1] is not None for a in ARMS)]
print(f'{corpus}: {len(ok)} chains'); json.dump(by, open(f"{S}/chain_{os.path.basename(corpus)}_{skip}.json", 'w'))
for a in ARMS:
    print(f"  {a:4s} screen cpu {sum(by[i][a][0] for i in ok):7.1f}  polish cpu {sum(by[i][a][1] for i in ok):7.1f}")
tot = {a: [by[i][a][0] + by[i][a][1] for i in ok] for a in ARMS}
for a in ('on', 'jit'):
    d = [tot[a][k] - tot['off'][k] for k in range(len(ok))]
    m = st.mean(d); se = st.stdev(d) / math.sqrt(len(d))
    print(f"  chain (screen + polish) cpu {a} - off: {m:+.3f} s/chain [{m-1.96*se:+.3f}, {m+1.96*se:+.3f}]  (off mean {st.mean(tot['off']):.3f}); ratio of sums {sum(tot[a])/sum(tot['off']):.3f}")
for a in ('on', 'jit'):
    lr = [math.log(by[i][a][1] / by[i]['off'][1]) for i in ok if by[i]['off'][1] > 0.02 and by[i][a][1] > 0.02]
    m = st.mean(lr); se = st.stdev(lr) / math.sqrt(len(lr))
    print(f"  polish cpu {a}/off: geo-mean ratio {math.exp(m):.3f} [{math.exp(m-1.96*se):.3f}, {math.exp(m+1.96*se):.3f}] over {len(lr)} chains")
