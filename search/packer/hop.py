#!/usr/bin/env python3
"""Basin hopping on fq quenches (10-07).  Benchmark use: start from the pre-SQUISH register packing
(../exact/batch/inputs/n-N.txt, register of 10-05) and see how close a fixed CPU budget gets to SQUISH's side.

Moves (per proposal, from the current state): global kick (sigma 0.003/0.01/0.03, angles 20 sigma deg) or local kick
(squares within radius r of a random square get sigma 0.05-0.2, the rest 0.002).  Quench = `fq quench` (loosen 1.02).
Metropolis at temperature T (side units); reset to the best state after `--patience` proposals without a new best.

  hop.py --n 108 130 --minutes 10 --chains 1 --procs 15 --out runs/hop0
  hop.py --summary runs/hop0
"""
import argparse, json, math, os, random, subprocess, tempfile, time
from jobpool import run_jobs
import mcmin

HERE = os.path.dirname(os.path.abspath(__file__))
FQ = os.path.join(HERE, 'target/release/fq')
BATCH = os.path.join(HERE, '../exact/batch/inputs')
EXTRA = []
SQUISH = os.path.expanduser('~/math/_untrusted-third-party/itsnaka-squish-certs/README.md')


def squish_targets():
    out = {}
    for line in open(SQUISH):
        f = [c.strip() for c in line.split('|')]
        if len(f) > 3 and f[1].isdigit() and f[2].replace('.', '').isdigit():
            out.setdefault(int(f[1]), float(f[2]))       # first table = current packings
    return out


def quench(s, sq, tmp, extra=()):
    a, b = os.path.join(tmp, 'a.txt'), os.path.join(tmp, 'b.txt')
    mcmin.write_deg(a, s, sq)
    r = subprocess.run([FQ, 'quench', '--in', a, '--out', b, *EXTRA, *extra], capture_output=True, text=True, timeout=900)
    d = json.loads(r.stdout)
    if not (d['min_gap'] >= 0 and d['min_wall'] >= 0):
        return None
    s2, sq2 = mcmin.load_deg(b)
    return s2, sq2, d


def propose(sq, rng):
    if rng.random() < 0.5:
        sig = rng.choice((0.003, 0.01, 0.03))
        return [(x + rng.gauss(0, sig), y + rng.gauss(0, sig), a + rng.gauss(0, 20 * sig)) for x, y, a in sq], f'kick{sig}'
    c = rng.randrange(len(sq))
    r = rng.choice((1.5, 2.5, 3.5))
    sig = rng.choice((0.05, 0.1, 0.2))
    cx, cy = sq[c][0], sq[c][1]
    out = []
    for x, y, a in sq:
        sg = sig if (x - cx) ** 2 + (y - cy) ** 2 < r * r else 0.002
        out.append((x + rng.gauss(0, sg), y + rng.gauss(0, sg), a + rng.gauss(0, 20 * sg)))
    return out, f'local{r}/{sig}'


def seed_remove(n, k, m, rng, tmp):
    """Best of m quenches of record(n + k) with k random squares removed."""
    s1, sq1 = mcmin.load_deg(f'{BATCH}/n-{n + k}.txt')
    best = None
    for _ in range(m):
        drop = set(rng.sample(range(n + k), k))
        r = quench(s1, [q for i, q in enumerate(sq1) if i not in drop], tmp)
        if r and (best is None or r[0] < best[0]):
            best = (r[0], r[1])
    return best


def chain(j):
    n, cid, minutes, T, patience, out, sk, sm = j
    rng = random.Random(n * 1000 + cid)
    tmp = tempfile.mkdtemp()
    s0, sq0 = mcmin.load_deg(f'{BATCH}/n-{n}.txt')
    if sk:
        cur = seed_remove(n, sk, sm, rng, tmp)
    else:
        r = quench(s0, sq0, tmp)
        cur = (r[0], r[1]) if r else (s0, sq0)
    best = cur
    mcmin.write_deg(f'{out}/n{n}_c{cid}_best.txt', *best)
    log = open(f'{out}/n{n}_c{cid}.jsonl', 'w')
    log.write(json.dumps(dict(start=s0, polished=cur[0], seed=f'{n + sk}-{sk}' if sk else 'register')) + '\n')
    t_end = time.time() + 60 * minutes
    k = since = 0
    while time.time() < t_end:
        k += 1
        prop, kind = propose(cur[1], rng)
        t0 = time.time()
        r = quench(cur[0], prop, tmp)
        if r is None:
            log.write(json.dumps(dict(k=k, kind=kind, s=None)) + '\n'); continue
        s, sq, d = r
        acc = s < cur[0] or rng.random() < math.exp(-(s - cur[0]) / T)
        if acc:
            cur = (s, sq)
        if s < best[0] - 1e-11:
            best = (s, sq); since = 0
            mcmin.write_deg(f'{out}/n{n}_c{cid}_best.txt', s, sq)
        else:
            since += 1
            if since >= patience:
                cur = best; since = 0
        log.write(json.dumps(dict(k=k, kind=kind, s=s, acc=acc, best=best[0], sec=round(time.time() - t0, 2))) + '\n')
        log.flush()
    return dict(n=n, cid=cid, start=s0, best=best[0], proposals=k)


def summary(out):
    tg = squish_targets()
    rows = {}
    for f in sorted(os.listdir(out)):
        if f.endswith('.jsonl'):
            L = [json.loads(l) for l in open(os.path.join(out, f))]
            if not L: continue
            n = int(f[1:].split('_')[0])
            b = min([L[0]['polished']] + [l['best'] for l in L[1:] if l.get('best')])
            rows.setdefault(n, []).append((L[0]['start'], b, len(L) - 1))
    print(f"{'n':>4} {'start':>16} {'ours best':>16} {'SQUISH':>16} {'our gain':>9} {'their gain':>10} {'frac':>5} props")
    for n in sorted(rows):
        st = rows[n][0][0]; b = min(r[1] for r in rows[n]); P = sum(r[2] for r in rows[n])
        t = tg.get(n, float('nan'))
        print(f"{n:4d} {st:16.10f} {b:16.10f} {t:16.10f} {st - b:9.2e} {st - t:10.2e} {(st - b) / (st - t):5.2f} {P}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, nargs='+')
    ap.add_argument('--minutes', type=float, default=10)
    ap.add_argument('--chains', type=int, default=1)
    ap.add_argument('--T', type=float, default=1e-5)
    ap.add_argument('--patience', type=int, default=40)
    ap.add_argument('--procs', type=int, default=15)
    ap.add_argument('--out', default='runs/hop0')
    ap.add_argument('--summary')
    ap.add_argument('--seed-k', type=int, default=0, help='start from record(n+k) minus k random squares')
    ap.add_argument('--seed-m', type=int, default=8, help='removal trials for the seed (best kept)')
    ap.add_argument('--loosen', type=float, default=1.02)
    a = ap.parse_args()
    if a.summary:
        summary(a.summary); raise SystemExit
    EXTRA[:] = ['--loosen', str(a.loosen)]
    os.makedirs(a.out, exist_ok=True)
    jobs = [(n, c, a.minutes, a.T, a.patience, a.out, a.seed_k, a.seed_m) for n in a.n for c in range(a.chains)]
    run_jobs(chain, jobs, procs=a.procs, timeout=60 * a.minutes + 1800)
    summary(a.out)
