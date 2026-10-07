#!/usr/bin/env python3
"""Diagnostics for an mcmin.py run (log.jsonl): proposals and acceptance per move kind and per temperature, swap acceptance per
pair and per walker, distinct minima, excursions below thresholds with N_eff = N / (1 + CV^2) .

  ptstat.py runs/pt1 [--thr 11,10.9972,10.9968]
"""
import argparse, collections, json, math, statistics


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('--thr', default='11,10.9972,10.9968')
    a = ap.parse_args()
    L = [json.loads(l) for l in open(f'{a.run}/log.jsonl')]
    meta = next(x['meta'] for x in L if 'meta' in x)
    Ts = meta['Ts']; R = len(Ts)
    P = [x for x in L if 'step' in x]
    S = [x for x in L if 'swaps' in x]
    print(f'{a.run}: {len(S)} epochs, {len(P)} proposals, {sum(x["dt"] for x in P) / 3600:.1f} CPU-h in proposals, '
          f'median dt {statistics.median(x["dt"] for x in P):.1f}s')
    # per kind (first move's kind; single-move proposals only for clean attribution)
    kd = collections.defaultdict(lambda: [0, 0, 0, 0])      # n, accepted, lines, sub11
    for x in P:
        k = x['moves'][0].split()[0] if len(x['moves']) == 1 else 'multi'
        c = kd[k]; c[0] += 1; c[1] += x.get('acc', False); c[2] += bool(x.get('lines')); c[3] += x.get('s', 99) < 11 and not x.get('lines')
    print('kind: proposals, acc, full-line, sub-11 results')
    for k, (n, ac, li, sb) in sorted(kd.items()):
        print(f'  {k:9s} {n:5d}  acc {ac / n:.2f}  lines {li / n:.2f}  sub11 {sb / n:.2f}')
    print('per T: proposals, acc, mean E, best s')
    for r in range(R):
        X = [x for x in P if x['rid'] == r]
        if X:
            print(f'  T={Ts[r]:.1e} {len(X):4d} acc {sum(x["acc"] for x in X) / len(X):.2f}  '
                  f'best {min((x["s"] for x in X if "s" in x and not x["lines"]), default=float("nan")):.7f}')
    # swaps per pair
    pc = collections.defaultdict(lambda: [0, 0])
    for x in S:
        for r, ok in x['swaps']:
            pc[r][0] += 1; pc[r][1] += ok
    print('swap acc per pair: ' + ' '.join(f'{pc[r][1] / max(pc[r][0], 1):.2f}' for r in range(R - 1)))
    # walker occupancy of the coldest slot
    cold = collections.Counter(x['slot'][0] for x in S)
    print(f'walkers at the coldest T (epochs): {dict(cold)}')
    # distinct minima (f64 side to 1e-8, line-free)
    sides = collections.Counter(round(x['s'], 8) for x in P if 's' in x and not x['lines'])
    print(f'distinct line-free minima visited {len(sides)}; below 11: {sum(1 for s in sides if s < 11)}; best {min(sides):.10f}')
    # excursions: per walker, runs of consecutive accepted states below threshold (in that walker's own proposal sequence)
    for th in map(float, a.thr.split(',')):
        durs = []
        for w in range(R):
            seq = [x for x in P if x['w'] == w]
            cur, ins = 0, False
            st = None
            for x in seq:
                if x['acc']:
                    st = x['s'] if not x['lines'] else 99
                if st is None:
                    continue
                if st < th:
                    cur += 1; ins = True
                elif ins:
                    durs.append(cur); cur, ins = 0, False
            if ins:
                durs.append(cur)
        if durs:
            m = statistics.mean(durs); cv = statistics.pstdev(durs) / m if m else 0
            print(f'excursions below {th}: N {len(durs)}, mean length {m:.1f} proposals, N_eff {len(durs) / (1 + cv * cv):.1f}')
        else:
            print(f'excursions below {th}: none')


if __name__ == '__main__':
    main()
