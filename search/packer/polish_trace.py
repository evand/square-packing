#!/usr/bin/env python3
"""Polish-length census (10-08, Evan: polish until the contact graph is identified, dedupe on graphs, finish later).
Replays the explorer's proposal pipeline on parents drawn from a real archive: move -> screen (<= 8 SLP iterations, no
flips) -> full polish (--no-alm), both with `fq --ident` (per accepted step: load-bearing network size, change count,
fingerprint).  Every proposal that the explorer would polish is polished (returns included, flagged), so the log answers:
how many iterations does a polish take, when does the load network stop changing, and how often do proposals that end
in the same basin already share the network earlier.

  polish_trace.py --archive runs/ex10 --n 110 --count 400 --procs 14 --out runs/pt_110.jsonl
"""
import argparse, json, math, os, random, tempfile, time
from jobpool import run_jobs
import mcmin, hop, explore

HERE = os.path.dirname(os.path.abspath(__file__))
hop.FQ = os.path.join(HERE, 'target-dev/release/fq')       # build with the per-step fingerprints (--ident "fps")


def fq(s, sq, tmp, args):
    import subprocess
    a, b = os.path.join(tmp, 'a.txt'), os.path.join(tmp, 'b.txt')
    mcmin.write_deg(a, s, sq)
    r = subprocess.run([hop.FQ, 'quench', '--in', a, '--out', b, *args], capture_output=True, text=True, timeout=1200)
    d = json.loads(r.stdout)
    if not (d['min_gap'] >= 0 and d['min_wall'] >= 0):
        return None, d
    return mcmin.load_deg(b), d


def job(j):
    path, kind, seed, k, smax, sides = j
    from layout import full_lines
    rng = random.Random(seed)
    s, sq = mcmin.load_deg(path)
    prop, desc = explore.MOVES[kind](s, sq, rng, explore.A)
    tmp = tempfile.mkdtemp()
    loosen = rng.choice(('1.0', '1.02', '1.05'))
    t0 = time.time()
    r1, d1 = fq(s, prop, tmp, ['--loosen', loosen, '--pit', '8', '--flip-top', '0', '--ident'])
    out = dict(parent=path, kind=kind, loosen=loosen, scr=dict(s=d1['s'], it=d1['slp_it'], sec=d1['sec'], fp=d1['fp'], fps=d1.get('fps', [])))
    if r1 is None:
        return dict(out, status='fail')
    s1, sq1 = r1
    if s1 > smax + 1e-3:
        return dict(out, status='screen-discard')
    p1 = os.path.join(tmp, 's1.txt'); mcmin.write_deg(p1, s1, sq1)
    if full_lines(p1, k):
        return dict(out, status='screen-grid')
    out['return'] = any(abs(x - s1) < 1e-6 for x in sides)
    r2, d2 = fq(s1, sq1, tmp, ['--loosen', '1.0', '--no-alm', '--ident'])
    if r2 is None:
        return dict(out, status='fail2')
    out['full'] = dict(s=d2['s'], it=d2['slp_it'], sec=d2['sec'], flips=d2['flips'], fp=d2['fp'], fps=d2.get('fps', []),
                       trace=d2['trace'])
    out['status'] = 'polished'
    out['sec'] = time.time() - t0
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--archive'); ap.add_argument('--n', type=int); ap.add_argument('--count', type=int, default=300)
    ap.add_argument('--procs', type=int, default=14); ap.add_argument('--out'); ap.add_argument('--seed', type=int, default=7)
    a = ap.parse_args()
    k = math.ceil(math.sqrt(a.n) - 1e-12); smax = k + 0.05
    E = [json.loads(l) for l in open(f'{a.archive}/archive.jsonl')]
    E = [e for e in E if e['s'] < smax and os.path.exists(e['path'])]
    sides = [e['s'] for e in E]
    rng = random.Random(a.seed)
    kinds, wts = zip(*explore.WEIGHTS.items())
    jobs = [(rng.choice(E)['path'], rng.choices(kinds, weights=wts)[0], rng.randrange(1 << 30), k, smax, sides)
            for _ in range(a.count)]
    res = run_jobs(job, jobs, procs=a.procs, timeout=1800)
    with open(a.out, 'w') as f:
        for r in res:
            if isinstance(r, dict):
                f.write(json.dumps(r) + '\n')
    print(a.out, sum(isinstance(r, dict) and r.get('status') == 'polished' for r in res), 'polished of', len(res))
