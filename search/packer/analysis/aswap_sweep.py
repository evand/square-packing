#!/usr/bin/env python3
"""AS1 (10-16): exhaustive greedy angle-swap descent.

The explorer's `aswap` move (mcmin.mv_aswap: swap the tilts of two neighbours, centre distance < 1.6, tilts differing
by > 3 deg mod 90) found s271_ex1.  Here every such pair is tried each round, each through the explorer's two fq
stages (SCREEN at loosen 1.0 / 1.02 / 1.05, then POLISH); the best side below the current one is accepted and the
sweep repeats from it until a round finds none.

  python3 analysis/aswap_sweep.py START.txt --out runs/as1 [--procs 4 --rounds 10 --loosen 1.0,1.02,1.05]
"""
import argparse, json, math, os, sys, tempfile, time
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, HERE)
import hop, mcmin
from pk.engine import SCREEN, POLISH


def adiff(u, v):
    return min((u - v) % 90, (v - u) % 90)


def pairs(sq, r=1.6, dmin=3.0):
    out = []
    for i in range(len(sq)):
        for j in range(i + 1, len(sq)):
            if math.hypot(sq[i][0] - sq[j][0], sq[i][1] - sq[j][1]) < r and adiff(sq[i][2], sq[j][2]) > dmin:
                out.append((i, j))
    return out


def trial(args):
    s, sq, i, j, loosen = args
    prop = list(sq)
    prop[i] = (sq[i][0], sq[i][1], sq[j][2]); prop[j] = (sq[j][0], sq[j][1], sq[i][2])
    with tempfile.TemporaryDirectory(prefix='as1_') as tmp:
        r1 = hop.quench(s, prop, tmp, opt=SCREEN.but(loosen=loosen))
        if r1 is None:
            return dict(i=i, j=j, loosen=loosen, s=None)
        s1, sq1, _ = r1
        if s1 > s + 2e-3:                                    # the explorer's gate: far above, no polish
            return dict(i=i, j=j, loosen=loosen, s=None, s_screen=s1)
        r2 = hop.quench(s1, sq1, tmp, opt=POLISH)
        if r2 is None:
            return dict(i=i, j=j, loosen=loosen, s=None, s_screen=s1)
        return dict(i=i, j=j, loosen=loosen, s=r2[0], s_screen=s1, sq=r2[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('start'); ap.add_argument('--out', required=True)
    ap.add_argument('--procs', type=int, default=4); ap.add_argument('--rounds', type=int, default=10)
    ap.add_argument('--loosen', default='1.0,1.02,1.05')
    a = ap.parse_args()
    from pk import manifest
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'aswap_sweep', args=vars(a), procs=a.procs)
    s, sq = mcmin.load_deg(a.start)
    L = [float(x) for x in a.loosen.split(',')]
    log = open(os.path.join(a.out, 'rounds.jsonl'), 'a')
    print(f'start {s:.12f}  n {len(sq)}', flush=True)
    hist = []
    with ProcessPoolExecutor(a.procs) as ex:
        for rd in range(a.rounds):
            t0 = time.time()
            P = pairs(sq)
            R = list(ex.map(trial, [(s, sq, i, j, l) for i, j in P for l in L], chunksize=4))
            ok = sorted((r for r in R if r['s'] is not None), key=lambda r: r['s'])
            best = ok[0] if ok else None
            gains = [r for r in ok if r['s'] < s - 1e-10]
            row = dict(round=rd, s_in=s, pairs=len(P), trials=len(R), feasible=len(ok), below=len(gains),
                       best=best['s'] if best else None, best_pair=(best['i'], best['j'], best['loosen']) if best else None,
                       top=[(r['s'], r['i'], r['j'], r['loosen']) for r in ok[:10]], sec=round(time.time() - t0))
            log.write(json.dumps(row) + '\n'); log.flush(); hist.append(row)
            print(f"round {rd}: {len(P)} pairs x {len(L)}, {len(ok)} feasible, {len(gains)} below; best "
                  f"{row['best']} ({row['best']-s if best else float('nan'):+.3e}) pair {row['best_pair']}  {row['sec']} s", flush=True)
            if not gains:
                break
            s, sq = gains[0]['s'], gains[0]['sq']
            mcmin.write_deg(os.path.join(a.out, f'best_r{rd:02d}.txt'), s, sq)
    mcmin.write_deg(os.path.join(a.out, 'final.txt'), s, sq)
    print(f'final {s:.12f}', flush=True)
    man.finish(summary=dict(final=s, rounds=len(hist)))


if __name__ == '__main__':
    main()
