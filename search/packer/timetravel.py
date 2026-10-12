#!/usr/bin/env python3
"""Time-travel record test report (10-12; pre-registration: PACKER.md "Night plan 10-11").

Chains from the frozen 10-05 register (runs/timetravel/{plain,sparse}, chain.py) scored against runs/timetravel/truth.json
(per n: s1005 = the 10-05 register side, s_now = the store's best when the test was frozen).  Per n and arm:
recovered gap fraction f = (s1005 - best) / (s1005 - s_now), capped at 1; full recovery = best <= s_now + 1e-9;
below = best < s_now - 1e-9 (check against the live register before calling anything a record).
Paired sparse vs plain: Wilcoxon signed-rank on f, median paired difference.

  timetravel.py [--root runs/timetravel]
"""
import argparse, json, os, statistics, sys
HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', default=f'{HERE}/runs/timetravel')
    a = ap.parse_args()
    T = json.load(open(f'{a.root}/truth.json'))
    arms = [d for d in ('plain', 'sparse') if os.path.exists(f'{a.root}/{d}/state.json')]
    S = {d: json.load(open(f'{a.root}/{d}/state.json'))['n'] for d in arms}
    rows = []
    print(f'{"n":>4} {"s1005":>14} {"s_now":>14} {"gap":>8} ' + ' '.join(f'{d + " best":>15} {"f":>5} {"rnd":>3}' for d in arms))
    for n in T['ns']:
        t = T['truth'][str(n)]
        r = dict(n=n, gap=t['gap'])
        line = f'{n:4d} {t["s1005"]:14.10f} {t["s_now"]:14.10f} {t["gap"]:8.1e} '
        for d in arms:
            v = S[d].get(str(n))
            best = v['best'] if v else t['s1005']
            f = min(1.0, max(0.0, (t['s1005'] - best) / t['gap']))
            r[d] = dict(best=best, f=f, rounds=v['rounds'] if v else 0, full=best <= t['s_now'] + 1e-9, below=best < t['s_now'] - 1e-9)
            line += f'{best:15.10f} {f:5.2f} {r[d]["rounds"]:3d}' + (' FULL' if r[d]['full'] else '') + (' BELOW' if r[d]['below'] else '') + ' '
        rows.append(r)
        print(line)
    for d in arms:
        F = [r[d]['f'] for r in rows]
        print(f'{d}: median f {statistics.median(F):.3f}, mean {statistics.mean(F):.3f}; full recovery {sum(r[d]["full"] for r in rows)} / '
              f'{len(rows)}; below s_now {[r["n"] for r in rows if r[d]["below"]]}; rounds {sum(r[d]["rounds"] for r in rows)}')
    if len(arms) == 2:
        from scipy.stats import wilcoxon
        D = [r['sparse']['f'] - r['plain']['f'] for r in rows]
        nz = [x for x in D if abs(x) > 1e-12]
        p = wilcoxon(nz).pvalue if len(nz) >= 5 else float('nan')
        print(f'sparse - plain: median paired difference {statistics.median(D):+.3f}, sparse better in {sum(x > 1e-12 for x in D)}, '
              f'worse in {sum(x < -1e-12 for x in D)}, tied {len(D) - len(nz)}; Wilcoxon p = {p:.3g}')
    json.dump(rows, open(f'{a.root}/report.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
