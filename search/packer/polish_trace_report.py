#!/usr/bin/env python3
"""Report for polish_trace.py logs: iteration histogram, settle point of the load network, work after settling,
and how well a fingerprint -> final-basin table (built online from earlier proposals) predicts the basin early.

  polish_trace_report.py runs/pt_110.jsonl [...]
"""
import collections, json, statistics as st, sys


def q(xs, ps=(0.1, 0.5, 0.9)):
    xs = sorted(xs)
    return '/'.join(f'{xs[min(len(xs) - 1, int(p * len(xs)))]:.3g}' for p in ps) if xs else '-'


def hist(xs, edges):
    c = collections.Counter()
    for x in xs:
        for lo, hi in zip(edges, edges[1:]):
            if lo <= x < hi:
                c[(lo, hi)] += 1; break
    return '  '.join(f'[{lo},{hi}):{c[(lo, hi)]}' for lo, hi in zip(edges, edges[1:]))


def basin_ids(P, tol=2e-9):
    ss = sorted(set(p['full']['s'] for p in P)); ids = {}; cur = -1; last = None
    for s in ss:
        if last is None or s - last > tol: cur += 1
        ids[s] = cur; last = s
    return {id(p): ids[p['full']['s']] for p in P}


def report(path):
    R = [json.loads(l) for l in open(path)]
    st_c = collections.Counter(r['status'] for r in R)
    P = [r for r in R if r['status'] == 'polished']
    print(f'== {path}: {len(R)} proposals {dict(st_c)}; returns (screen within 1e-6 of an archive side) {sum(p["return"] for p in P)}')
    it_s = [p['scr']['it'] for p in P]; it_f = [p['full']['it'] for p in P]
    print(f'  screen iterations p10/50/90 {q(it_s)}; full polish iterations {q(it_f)}; full polish sec {q([p["full"]["sec"] for p in P])}')
    print('  full-polish iteration histogram:', hist(it_f, [0, 3, 6, 10, 20, 40, 80, 150, 299, 301]))
    tot = sum(p['full']['sec'] for p in P)
    big = [p for p in P if p['full']['it'] >= 80]
    print(f'  share of full-polish time in polishes >= 80 it: {sum(p["full"]["sec"] for p in big) / tot:.2f} ({len(big)} polishes); '
          f'>= 299 it: {sum(p["full"]["sec"] for p in P if p["full"]["it"] >= 299) / tot:.2f}')
    # settle point: first accepted step j with fps[j:] constant and equal to the final network fingerprint
    after_t, after_ds, after_dm, never, settle_it, frac_it = [], [], [], 0, [], []
    for p in P:
        f = p['full']; F = f['fps']; T = f['trace']
        if not F or F[-1] != f['fp']:
            never += 1
            continue
        j = len(F) - 1
        while j > 0 and F[j - 1] == F[-1]: j -= 1
        settle_it.append(T[j][0]); frac_it.append(T[j][0] / max(1, f['it']))
        after_t.append((f['sec'] - T[j][1]) / max(1e-9, f['sec']))
        after_ds.append(T[j][2] - f['s'])
        after_dm.append(sum(t[4] for t in T[j + 1:]))
    print(f'  load network settles (last step = final fp): {len(settle_it)}/{len(P)} (final network differs from last step: {never})')
    print(f'    settle iteration p10/50/90 {q(settle_it)}; as fraction of the polish\'s iterations {q(frac_it)}')
    print(f'    share of polish time after settling {q(after_t)} (time-weighted overall '
          f'{sum(a * p["full"]["sec"] for a, p in zip(after_t, P)) / tot:.2f})')
    print(f'    side error at settle (s_settle - s_final) {q(after_ds)}; summed max-corner displacement after settle {q(after_dm)}')
    # online fingerprint table: fp -> final basins seen through that fp on earlier trajectories (screen + full)
    bid = basin_ids(P)
    table = collections.defaultdict(set)
    correct = wrong = none = 0; saved = 0.0
    for p in P:
        traj = p['scr']['fps'] + p['full']['fps']
        b = bid[id(p)]; T = p['full']['trace']; ns = len(p['scr']['fps'])
        hit = None
        for j, h in enumerate(traj):
            if h in table and len(table[h]) == 1:
                hit = j; break
        if hit is None:
            none += 1
        elif next(iter(table[traj[hit]])) == b:
            correct += 1
            jt = hit - ns                    # index into the full polish (negative = already at the screen)
            saved += p['full']['sec'] - (T[jt][1] if jt >= 0 else 0.0)
        else:
            wrong += 1
        for h in traj:
            table[h].add(b)
    amb = sum(len(v) > 1 for v in table.values())
    print(f'  online fp->basin table: early prediction correct {correct}, wrong {wrong}, no prediction {none}; '
          f'polish time saved if stopped at the first correct match {saved / tot:.2f}; '
          f'{len(table)} fps, {amb} ambiguous; distinct final basins {len(set(bid.values()))}')
    fin = collections.Counter(bid.values())
    print(f'  final-basin multiplicity: {sum(c > 1 for c in fin.values())} basins reached > once, covering '
          f'{sum(c for c in fin.values() if c > 1)}/{len(P)} polishes; final fp per basin unique for '
          f'{sum(len({p["full"]["fp"] for p in P if bid[id(p)] == b}) == 1 for b in fin)}/{len(fin)}')


if __name__ == '__main__':
    for f in sys.argv[1:]:
        report(f)
