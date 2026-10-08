#!/usr/bin/env python3
"""Tiling starts for "tilings are never optimal" (10-07; WISHLIST N17/P6, TODO "Tiling starts").

From the register record for n (../exact/batch/inputs/n-N.txt), build a packing of 4n squares in side 2·s(n):
  T  plain 2×2 copies
  M  mirrored copies (copy (i, j) reflected in x if i = 1, in y if j = 1; the seams meet their mirror images)
  R  pinwheel (copy q rotated by 90°·q about its own centre)
  S  scale ×2 and split every square into a 2×2 block (tilted blocks, the other construction)
then quench (fq, loosen 1.0 and 1.02) and hop (hop.search) for a fixed budget.  The conjecture predicts every start can
be improved (side < 2·s(n)).

  tile_hop.py --n 241 273 307 65 --variants T M R S --minutes 20 --procs 8 --out runs/tile0
  tile_hop.py --summary runs/tile0
"""
import argparse, json, math, os, random, tempfile, time
import hop, mcmin
from jobpool import run_jobs


def tiles(s, sq, v):
    out = []
    if v == 'S':
        for x, y, a in sq:
            t = math.radians(a)
            c, si = math.cos(t), math.sin(t)
            for u in (-0.5, 0.5):
                for w in (-0.5, 0.5):
                    out.append((2 * x + c * u - si * w, 2 * y + si * u + c * w, a))
        return 2 * s, out
    for q, (i, j) in enumerate(((0, 0), (1, 0), (1, 1), (0, 1))):
        for x, y, a in sq:
            if v == 'M':
                if i: x, a = s - x, -a
                if j: y, a = s - y, -a
            elif v == 'R':
                for _ in range(q):
                    x, y = s - y, x            # rotate by 90° about (s/2, s/2); angle mod 90 unchanged
            out.append((x + i * s, y + j * s, a))
    return 2 * s, out


def job(j):
    n, v, minutes, T, patience, out, seed = j
    rng = random.Random(seed)
    tmp = tempfile.mkdtemp()
    s0, sq0 = mcmin.load_deg(f'{hop.BATCH}/n-{n}.txt')
    S, sq = tiles(s0, sq0, v)
    tag = f'n{4 * n}_{v}'
    mcmin.write_deg(f'{out}/{tag}_seed.txt', S, sq)
    log = open(f'{out}/{tag}.jsonl', 'w')
    rec = dict(n=4 * n, base=n, variant=v, base_s=s0, tiled=S)
    t0 = time.time()
    hop.EXTRA[:] = ['--loosen', '1.0']
    r10 = hop.quench(S, sq, tmp)
    rec['q10'] = r10[0] if r10 else None; rec['q10_sec'] = round(time.time() - t0, 1)
    t0 = time.time()
    hop.EXTRA[:] = ['--loosen', '1.02']
    r102 = hop.quench(S, sq, tmp)
    rec['q102'] = r102[0] if r102 else None; rec['q102_sec'] = round(time.time() - t0, 1)
    log.write(json.dumps(rec) + '\n'); log.flush()
    cands = [r for r in (r10, r102) if r]
    if not cands:
        return rec
    s1, sq1, _ = min(cands, key=lambda r: r[0])
    mcmin.write_deg(f'{out}/{tag}_best.txt', s1, sq1)
    if minutes > 0:
        best = hop.search(s1, sq1, 60 * minutes, rng, tmp, log, T=T, patience=patience)
        mcmin.write_deg(f'{out}/{tag}_best.txt', *best)
        rec['hop'] = best[0]
    rec['best'] = min(x for x in (rec.get('q10'), rec.get('q102'), rec.get('hop')) if x is not None)
    rec['gain'] = S - rec['best']
    log.write(json.dumps(dict(final=rec)) + '\n')
    return rec


def summary(out):
    rows = []
    for f in sorted(os.listdir(out)):
        if f.endswith('.jsonl'):
            L = [json.loads(l) for l in open(os.path.join(out, f))]
            if not L: continue
            r = L[0]
            hb = [l['best'] for l in L[1:] if 'best' in l]
            r['hop'] = min(hb) if hb else None
            r['props'] = len(hb)
            rows.append(r)
    print(f"{'4n':>5} {'v':>2} {'tiled 2s(n)':>14} {'q loosen1.0':>14} {'q 1.02':>14} {'hop best':>14} {'gain':>9} props  qsec")
    for r in rows:
        b = min(x for x in (r.get('q10'), r.get('q102'), r.get('hop')) if x is not None) if any(
            r.get(k) is not None for k in ('q10', 'q102', 'hop')) else float('nan')
        f = lambda x: f'{x:14.9f}' if x is not None else f"{'-':>14}"
        print(f"{r['n']:5d} {r['variant']:>2} {r['tiled']:14.9f} {f(r.get('q10'))} {f(r.get('q102'))} {f(r.get('hop'))} "
              f"{r['tiled'] - b:9.2e} {r['props']:5d}  {r.get('q10_sec')}/{r.get('q102_sec')}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, nargs='+')
    ap.add_argument('--variants', nargs='+', default=['T', 'M', 'R', 'S'])
    ap.add_argument('--minutes', type=float, default=20)
    ap.add_argument('--T', type=float, default=1e-5)
    ap.add_argument('--patience', type=int, default=40)
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--out', default='runs/tile0')
    ap.add_argument('--summary')
    a = ap.parse_args()
    if a.summary:
        summary(a.summary); raise SystemExit
    os.makedirs(a.out, exist_ok=True)
    jobs = [(n, v, a.minutes, a.T, a.patience, a.out, 1000 * n + k) for n in a.n for k, v in enumerate(a.variants)]
    run_jobs(job, jobs, procs=a.procs, timeout=60 * a.minutes + 3600)
    summary(a.out)
