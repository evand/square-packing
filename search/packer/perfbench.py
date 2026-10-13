#!/usr/bin/env python3
"""End-to-end explorer throughput benchmark (10-13 performance session; PACKER "perf score").

Perf score = explorer proposals per measured CPU-hour (Python + fq + anneal, the explorer's own manifest cpu_s), with
fixed seeds and fixed start sets, at n = 110 (bench110 recipe) and n = 237; the aggregate is the geometric mean of the
two speedups against a reference run.  It counts only together with an equivalence check (fqreplay.py): a change that
is a re-roll of the same search turns more proposals per CPU-h into more events per CPU-h.

  perfbench.py run --bin DIR --out runs/perfbench/NAME [--reps 8 --minutes 6]
      DIR holds fq, anneal, packer (e.g. target/release); 2 x reps single-proc explorers run concurrently
  perfbench.py report runs/perfbench/A runs/perfbench/B ...       (first = reference)

Starts: runs/perfbench/starts{110,237}/rep<r>.lst, created on first use (110: bench110 'rand' pool, 40 per replicate;
237: 4 of the 200 best store packings per replicate) and then frozen, so every run sees the same starts.
"""
from __future__ import annotations
import argparse, glob, json, math, os, random, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, 'runs', 'perfbench')
RECIPE = ['--adapt', '--reward', 'legacy', '--sparse-share', '0.5']


def _write(path, n, p):
    with open(path, 'w') as f:
        f.write(f'{n} {float(p.s)!r}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(t)!r}\n' for x, y, t in p.sq))


def starts(n, reps):
    d = os.path.join(ROOT, f'starts{n}')
    if not os.path.exists(os.path.join(d, f'rep{reps - 1}.lst')):
        os.makedirs(d, exist_ok=True)
        if n == 110:
            import bench110
            P = bench110.Pools()
            for r in range(reps):
                sd = os.path.join(d, f'rep{r}'); os.makedirs(sd, exist_ok=True)
                open(os.path.join(d, f'rep{r}.lst'), 'w').write(' '.join(P.starts('rand')(1000 + r, sd)))
        else:
            from pk.store import Store
            st = Store()
            ids = [r[0] for r in st.db.execute("SELECT id FROM packing WHERE n=? AND status!='infeasible' ORDER BY s LIMIT 200", (n,))]
            rng = random.Random(n)
            for r in range(reps):
                fs = []
                for i in rng.sample(ids, 4):
                    f = os.path.join(d, f'p{i}.txt'); _write(f, n, st.get(i)); fs.append(f)
                open(os.path.join(d, f'rep{r}.lst'), 'w').write(' '.join(fs))
    return [open(os.path.join(d, f'rep{r}.lst')).read().split() for r in range(reps)]


DRIVE = """import sys, textwrap; sys.path.insert(0, {here!r})
import hop, explore, mcmin
from pk import moves as pkm
B = {bin!r}
hop.FQ = B + '/fq'; explore.ANNEAL = B + '/anneal'; pkm.ANNEAL = B + '/anneal'; mcmin.PACKER = B + '/packer'
sys.argv = ['explore.py'] + {argv!r}
exec(textwrap.dedent(open({explore!r}).read().split("if __name__ == '__main__':", 1)[1]), vars(explore))
"""


def run_cmd(a):
    bin_ = os.path.abspath(a.bin)
    os.makedirs(a.out, exist_ok=True)
    procs = []
    for n in (110, 237):
        for r, st in enumerate(starts(n, a.reps)):
            out = os.path.join(a.out, f'n{n}_r{r}')
            argv = ['--n', str(n), *RECIPE, '--starts', *st, '--minutes', str(a.minutes), '--procs', '1',
                    '--seed', str(1000 + r), '--out', out]
            code = DRIVE.format(here=HERE, bin=bin_, argv=argv, explore=os.path.join(HERE, 'explore.py'))
            procs.append(subprocess.Popen([sys.executable, '-c', code], stdout=open(out + '.log', 'w'), stderr=subprocess.STDOUT))
    for p in procs:
        p.wait()
    json.dump(dict(bin=bin_, reps=a.reps, minutes=a.minutes, t=time.strftime('%F %T')), open(os.path.join(a.out, 'run.json'), 'w'))
    report([a.out])


def stats(d):
    out = {}
    for n in (110, 237):
        P, C = 0, 0.0
        for m in glob.glob(os.path.join(d, f'n{n}_r*', 'manifest.json')):
            C += json.load(open(m)).get('cpu_s') or 0
            P += sum(1 for _ in open(os.path.join(os.path.dirname(m), 'proposals.jsonl')))
        out[n] = (P, C, 3600 * P / C if C else float('nan'))
    return out


def report(dirs):
    S = [stats(d) for d in dirs]
    print(f"{'run':40s} {'n=110 prop/CPU-h':>18s} {'n=237 prop/CPU-h':>18s} {'speedup 110':>12s} {'speedup 237':>12s} {'perf score':>10s}")
    for d, s in zip(dirs, S):
        r0 = S[0]
        x110, x237 = s[110][2] / r0[110][2], s[237][2] / r0[237][2]
        print(f'{os.path.basename(d.rstrip("/")):40s} {s[110][2]:18.0f} {s[237][2]:18.0f} {x110:12.2f} {x237:12.2f} {math.sqrt(x110 * x237):10.2f}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('run'); p.add_argument('--bin', required=True); p.add_argument('--out', required=True)
    p.add_argument('--reps', type=int, default=8); p.add_argument('--minutes', type=float, default=6)
    p = sub.add_parser('report'); p.add_argument('dirs', nargs='+')
    a = ap.parse_args()
    run_cmd(a) if a.cmd == 'run' else report(a.dirs)
