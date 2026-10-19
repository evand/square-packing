#!/usr/bin/env python3
"""Stage 0 of the symmetric-subspace test (10-18): free vs Md-symmetric random starts, the SAME soft-to-rigid run
(`harden reharden`, box free under pressure while soft = Nakajima morph) + fq polish.  Arm sym: symmetric start AND
dynamics projected onto the Md subspace every step (`harden --sym-md 1`); the fq polish after is unconstrained (= the
release).  (A first smoke with symmetric starts but free dynamics lost the symmetry in every run, noise or not.)

  symharden.py run --ns 52,82,125 --k 32 --out runs/symsub/s0 [--cells "r0=0.2,0.35;hold=30000"] [--procs 15]
  symharden.py report runs/symsub/s0
Outcome per run: polished side s; hit = |s - record| < 1e-7; subgrid = s < ceil(record) - 1e-6; sym_out = output
still Md-invariant (1e-4).
"""
import sys, os, json, argparse, itertools, subprocess, tempfile, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from symstart import free_start, sym_start, ondiag
from symscan import invariant

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HARDEN = os.path.join(HERE, 'target', 'release', 'harden')


def write_cfg(path, s, sq):
    with open(path, 'w') as f:
        f.write(f'{len(sq)} {float(s)!r}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(t)!r}\n' for x, y, t in sq))


def job(t):
    from pk.engine import fq_quench, POLISH
    from pk.packing import read
    arm, n, rec, m, seed, cell = t
    rng = np.random.default_rng(seed)
    S = rec * 1.03
    sq = free_start(n, S, rng) if arm == 'free' else sym_start(n, S, None, rng)
    c0 = time.process_time(); w0 = time.time()
    with tempfile.TemporaryDirectory() as tmp:
        a, b = os.path.join(tmp, 'a.txt'), os.path.join(tmp, 'b.txt')
        write_cfg(a, S, sq)
        args = [HARDEN, 'reharden', '--in', a, '--out', b, '--seed', str(seed), '--hold-move', '1']
        if arm == 'sym':
            args += ['--sym-md', '1']             # dynamics projected onto the Md subspace every step
        for k, v in cell.items():
            args += [f'--{k.replace("_", "-")}', str(v)]
        r = subprocess.run(args, capture_output=True, text=True, timeout=900)
        if r.returncode or not os.path.exists(b):
            return dict(arm=arm, n=n, seed=seed, cell=cell, s=None, err=r.stderr[-200:])
        sym_h = None
        try:
            ph = read(b); sym_h = bool(invariant(ph.sq, ph.s, "Md", 1e-4)[0])
        except Exception:
            pass
        run = fq_quench(b, POLISH, tmp=tmp, timeout=900)
        p = read(run.out) if run.out else None
    if p is None:
        return dict(arm=arm, n=n, seed=seed, cell=cell, s=None)
    s = float(p.s)
    return dict(arm=arm, n=n, seed=seed, cell=cell, s=s, rec=rec, hit=bool(abs(s - rec) < 1e-7),
                subgrid=bool(s < np.ceil(rec - 1e-9) - 1e-6), pen=float(p.max_pen()),
                sym_hard=sym_h, sym_out=bool(invariant(p.sq, p.s, "Md", 1e-4)[0]), wall=time.time() - w0,
                sq=p.sq.tolist())


def parse_cells(g):
    keys, vals = [], []
    for part in filter(None, g.split(';')):
        k, v = part.split('='); keys.append(k); vals.append(v.split(','))
    return [dict(zip(keys, c)) for c in itertools.product(*vals)]


def run(a):
    from multiprocessing import Pool
    from pk.store import Store
    from pk import manifest
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'symsub-s0', args=vars(a), procs=a.procs)
    sc = json.load(open('runs/symsub/symscan.json'))
    st = Store()
    T = []
    for n in map(int, a.ns.split(',')):
        rec = st.get(sc[str(n)]['pid']); m = ondiag(rec.sq)
        for cell in parse_cells(a.cells):
            for i in range(a.k):
                for arm in ('free', 'sym'):
                    T.append((arm, n, float(rec.s), m, 10000 * n + i, cell))
    print(len(T), 'runs', flush=True)
    with Pool(a.procs) as pool, open(os.path.join(a.out, 'runs.jsonl'), 'w') as fo:
        for r in pool.imap_unordered(job, T):
            fo.write(json.dumps(r) + '\n'); fo.flush()
    man.finish()
    report(a)


def report(a):
    import collections
    R = [json.loads(l) for l in open(os.path.join(a.out, 'runs.jsonl'))]
    g = collections.defaultdict(list)
    for r in R:
        if r.get('s') is not None and r.get('pen', 1) <= 1e-9:
            g[(r['n'], json.dumps(r['cell'], sort_keys=True), r['arm'])].append(r)
    print('n cell arm: runs hits subgrid best gap median gap sym_hard sym_out  wall')
    for k in sorted(g):
        v = g[k]; gaps = sorted(r['s'] - r['rec'] for r in v)
        print(k[0], k[1], k[2], len(v), sum(r['hit'] for r in v), sum(r['subgrid'] for r in v),
              f'{gaps[0]:+.2e} {gaps[len(gaps)//2]:+.2e}', sum(bool(r.get('sym_hard')) for r in v),
              sum(r['sym_out'] for r in v), f"{np.mean([r['wall'] for r in v]):.1f}s")
    for arm in ('free', 'sym'):
        v = [r for kk, vv in g.items() if kk[2] == arm for r in vv]
        print('TOTAL', arm, len(v), 'hits', sum(r['hit'] for r in v), 'subgrid', sum(r['subgrid'] for r in v))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest='cmd')
    r = sub.add_parser('run'); r.add_argument('--ns', required=True); r.add_argument('--k', type=int, default=32)
    r.add_argument('--out', required=True); r.add_argument('--cells', default='r0=0.2,0.35;hold=30000')
    r.add_argument('--procs', type=int, default=15)
    q = sub.add_parser('report'); q.add_argument('out')
    a = ap.parse_args()
    run(a) if a.cmd == 'run' else report(a)
