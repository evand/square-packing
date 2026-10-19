#!/usr/bin/env python3
"""DUP3 (10-18): duplicate-half as a START DISTRIBUTION (Evan: "a new idea never produces one packing; read it as
starts for refinement").  Small: 6 n, ~96 variants each, filter, then equal explorer refinement vs a baseline.

  gen   : per n, variants = mirror {Md, Ma} x half {lo, hi} x seam band w {0, 0.7, 1.5} x r0 {0.05, 0.15} x box {loose
          (x1.03, box free while soft), tight (held at the record side, overlaps persist)} (+ a second seed at w = 0).
          Seam band w > 0: squares within w of the diagonal deleted after the mirror, the count refilled with random
          poses inside the band, band left free during the symmetric soft stage (harden --sym-band w).  Free polish.
  filter: feasible, within 5e-3 of the record, new (Hungarian matching distance >= 0.005 to the record, best of 8
          symmetries) and distinct (sides 1e-7 apart); up to 6 per n, best first.
  refine: explore.py 5 min / 1 proc per kept start (harden_explore's flags); baseline arm: explore.py from the record
          itself, 6 seeds per n, same 5 min.
  dup3.py gen|refine|report --out runs/symsub/dup3
"""
import sys, os, json, argparse, tempfile, subprocess, collections, time
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(HERE, 'core'))
import numpy as np
import symtt
from symharden import write_cfg, HARDEN

NS = (51, 110, 153, 172, 207, 297)


def variant(t):
    from pk.engine import fq_quench, POLISH
    from pk.packing import read
    n, s, sq, el, half, w, r0, box, seed = t
    rng = np.random.default_rng(seed + 31 * n)
    x = symtt.to_md(np.array(sq), s, el) if el == 'Md' else symtt.to_md(np.array(sq), s, 'Ma')
    y = symtt.duplicate_half(x, s, n, half, rng)
    if w > 0:
        d = np.abs(y[:, 0] - y[:, 1]) / np.sqrt(2)
        y = y[d >= w]
        k = n - len(y)
        if k > 0:                                                    # refill the band: random poses near the diagonal
            u = rng.uniform(0.5, s - 0.5, k); v = rng.uniform(-w, w, k) / np.sqrt(2)
            y = np.r_[y, np.c_[np.clip(u + v, 0.5, s - 0.5), np.clip(u - v, 0.5, s - 0.5), rng.uniform(0, 90, k)]]
        y = y[:n]
    L = 1.03 if box == 'loose' else 1.0
    with tempfile.TemporaryDirectory() as tmp:
        a, b = os.path.join(tmp, 'a'), os.path.join(tmp, 'b')
        c = s / 2; z = y.copy(); z[:, :2] = c * L + (z[:, :2] - c) * (1 if L == 1.0 else L)
        write_cfg(a, s * L, z)
        args = [HARDEN, 'reharden', '--in', a, '--out', b, '--seed', str(seed), '--r0', str(r0), '--sym-md', '1']
        args += ['--hold-move', '1'] if box == 'loose' else ['--hold', '20000']
        if w > 0:
            args += ['--sym-band', str(w)]
        subprocess.run(args, capture_output=True, text=True, timeout=900)
        if not os.path.exists(b):
            return None
        run = fq_quench(b, POLISH, tmp=tmp, timeout=900)
        p = read(run.out) if run.out else None
    if p is None or p.n != n or p.max_pen() > 1e-9:
        return None
    return dict(n=n, el=el, half=half, w=w, r0=r0, box=box, seed=seed, s=float(p.s), sq=p.sq.tolist())


def gen(a):
    from multiprocessing import Pool
    from pk.store import Store
    from coremap import align
    st = Store(); os.makedirs(a.out, exist_ok=True)
    recs = {n: st.best(n) for n in NS}
    T = []
    for n, p in recs.items():
        for el in ('Md', 'Ma'):
            for half in ('lo', 'hi'):
                for w in (0.0, 0.7, 1.5):
                    for r0 in (0.05, 0.15):
                        for box in ('loose', 'tight'):
                            for seed in ((1, 2) if w == 0 else (1,)):
                                T.append((n, p.s, p.sq.tolist(), el, half, w, r0, box, seed))
    print(len(T), 'variants', flush=True)
    t0 = time.time()
    with Pool(a.procs) as pool:
        R = [r for r in pool.map(variant, T, chunksize=1) if r]
    print(f'gen {time.time()-t0:.0f}s wall, {len(R)} feasible', flush=True)
    keep = collections.defaultdict(list)
    for r in sorted(R, key=lambda r: r['s']):
        n = r['n']; rec = recs[n]; gap = r['s'] - rec.s
        r['gap'] = gap
        if gap > 5e-3 or len(keep[n]) >= 6:
            continue
        if any(abs(r['s'] - k['s']) < 1e-7 for k in keep[n]):
            continue
        md = align(rec.sq, np.array(r['sq']), rec.s)[0]
        r['md'] = float(md)
        if md < 0.005:
            continue
        keep[n].append(r)
    json.dump(dict(all=[{k: v for k, v in r.items() if k != 'sq'} for r in R], keep=keep), open(os.path.join(a.out, 'gen.json'), 'w'))
    for n in NS:
        g = sorted(r['gap'] for r in R if r['n'] == n)
        print(n, 'variants', len(g), 'within 5e-3:', sum(x <= 5e-3 for x in g), 'within 1e-3:', sum(x <= 1e-3 for x in g),
              'kept new:', len(keep[n]), [f"{k['gap']:+.1e}/md{k['md']:.3f}/{k['el']}{k['half']}w{k['w']}{k['box'][0]}" for k in keep[n]])
        d = os.path.join(a.out, 'starts', f'n{n}'); os.makedirs(d, exist_ok=True)
        for i, k in enumerate(keep[n]):
            write_cfg(os.path.join(d, f'dup_{i}.txt'), k['s'], np.array(k['sq']))
        write_cfg(os.path.join(d, 'record.txt'), recs[n].s, recs[n].sq)


def explore_job(t):
    n, start, out, seed, minutes, frontier = t
    cmd = [sys.executable, os.path.join(HERE, 'explore.py'), '--n', str(n), '--adapt', '--reward', 'value', '--elite-share', '0.3',
           '--polish-gate', '2e-3', '--gate-ref', 'run', '--frontier-s', repr(frontier), '--starts', start,
           '--minutes', str(minutes), '--procs', '1', '--seed', str(seed), '--out', out]
    r = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))
    return out, r.returncode


def refine(a):
    from multiprocessing import Pool
    from pk.store import Store
    st = Store()
    G = json.load(open(os.path.join(a.out, 'gen.json')))
    T = []
    for n in NS:
        rec = st.best(n); d = os.path.join(a.out, 'starts', f'n{n}')
        k = len(G['keep'].get(str(n), []))
        for i in range(k):
            T.append((n, os.path.join(d, f'dup_{i}.txt'), os.path.join(a.out, 'ref', f'n{n}', f'dup_{i}'), 7000 + 10 * n + i, a.minutes, rec.s))
        for i in range(k):                                             # baseline: same number of runs from the record
            T.append((n, os.path.join(d, 'record.txt'), os.path.join(a.out, 'ref', f'n{n}', f'base_{i}'), 9000 + 10 * n + i, a.minutes, rec.s))
    print(len(T), 'explorer runs', flush=True)
    with Pool(a.procs) as pool:
        for out, rc in pool.imap_unordered(explore_job, T):
            print(out, rc, flush=True)


def report(a):
    import glob
    from pk.store import Store
    st = Store()
    print('n arm runs | best refined - record | distinct minima within 1e-4 of the record (sides 1e-9 apart) | <= record')
    tot = collections.Counter()
    for n in NS:
        rec = st.best(n).s
        for arm in ('dup', 'base'):
            sides, best = set(), []
            runs = sorted(glob.glob(os.path.join(a.out, 'ref', f'n{n}', f'{arm}_*')))
            for d in runs:
                b = None
                for l in open(os.path.join(d, 'archive.jsonl')) if os.path.exists(os.path.join(d, 'archive.jsonl')) else []:
                    e = json.loads(l)
                    if e.get('unpolished'):
                        continue
                    b = e['s'] if b is None else min(b, e['s'])
                    if e['s'] - rec < 1e-4:
                        sides.add(round(e['s'], 9))
                if b is not None:
                    best.append(b - rec)
            le = sum(1 for x in sides if x <= rec + 1e-9)
            tot[(arm, 'near')] += len(sides); tot[(arm, 'le')] += le
            print(n, arm, len(runs), '|', ' '.join(f'{x:+.1e}' for x in sorted(best)), '|', len(sides), '|', le)
    print('TOTAL distinct near-record minima: dup', tot[('dup', 'near')], 'base', tot[('base', 'near')],
          '; at/below record: dup', tot[('dup', 'le')], 'base', tot[('base', 'le')])


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('cmd'); ap.add_argument('--out', default='runs/symsub/dup3')
    ap.add_argument('--procs', type=int, default=15); ap.add_argument('--minutes', type=float, default=5)
    a = ap.parse_args()
    {'gen': gen, 'refine': refine, 'report': report}[a.cmd](a)
