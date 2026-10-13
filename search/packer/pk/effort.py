"""Effort map (10-13): measured CPU per n over a time window, from run manifests, against where records happened.

Every launcher writes runs/**/manifest.json (pk.manifest: start, cpu_s, args).  This module attributes each run's CPU to
the n it worked on and plots it next to other people's record events (bench/events.json + the store's pending
sightings) and our certified candidates.  Evan 10-13: "every run should have a target compute budget and evaluate
after the fact what it got done" -- this is the after-the-fact view.

Attribution, per run kind (manifest name):
  compose2         fills.jsonl / results.jsonl 'sec' per n (else equal shares over results' n)
  relax            results.jsonl rows per n
  explore, candidate  args.n (chain: its explore sub-runs carry their own manifests)
  fresh            its compose sub-runs carry their own manifests; the crop step is counted as pieces
  rectlib, rectsearch (pieces, not n-specific), others: 'unattributed' by name
Running jobs without cpu_s count as 0 (finish() writes cpu_s).

  python3 -m pk.effort --since 2026-10-12T21:00 [--until ISO] [--plot runs/plots/effort.png] [--nmax 480]
"""
from __future__ import annotations
import argparse, collections, glob, json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, '..', 'runs')
PIECE_KINDS = ('rectlib', 'rectsearch', 'rectsearch-bandit')


def manifests(since, until=None):
    for m in glob.glob(os.path.join(RUNS, '**', 'manifest.json'), recursive=True):
        if '/ext/' in m:
            continue
        try:
            d = json.load(open(m))
        except Exception:
            continue
        if not isinstance(d, dict) or not d.get('start'):
            continue
        if d['start'] < since or (until and d['start'] >= until):
            continue
        yield os.path.dirname(m), d


def _weights(run):
    w = collections.Counter()
    for f in ('fills.jsonl', 'results.jsonl'):
        p = os.path.join(run, f)
        if os.path.exists(p):
            for l in open(p):
                try:
                    x = json.loads(l)
                except Exception:
                    continue
                if isinstance(x, dict) and x.get('n'):
                    w[int(x['n'])] += x.get('sec', 1) or 1
            if w:
                break
    return w


def attribute(since, until=None):
    """{n: {kind: cpu_h}}, {kind: cpu_h} unattributed, [(run, kind, cpu_h)]"""
    per_n = collections.defaultdict(collections.Counter)
    other, runs = collections.Counter(), []
    for run, d in manifests(since, until):
        kind, cpu = d.get('name') or '?', (d.get('cpu_s') or 0) / 3600
        args = d.get('args') or {}
        if kind == 'chain' or kind == 'fresh':
            # sub-runs (explore, compose2, rectlib) have their own manifests; count the parent only for its own CPU
            # beyond them is not separable here, so skip it (cpu_s of the parent includes children)
            runs.append((run, kind, 0.0))
            continue
        runs.append((run, kind, cpu))
        if kind in PIECE_KINDS:
            other['pieces'] += cpu; continue
        if isinstance(args.get('n'), int):            # explore, candidate: one n per run
            per_n[int(args['n'])][kind] += cpu; continue
        w = _weights(run)
        if w:
            tot = sum(w.values())
            for n, v in w.items():
                per_n[n][kind] += cpu * v / tot
        else:
            other[kind] += cpu
    return per_n, other, runs


def record_events(since_events='2026-10-05'):
    """n -> count of others' record events (events table; posted only)."""
    p = os.path.join(HERE, '..', 'bench', 'events.json')
    c = collections.Counter()
    if os.path.exists(p):
        for e in json.load(open(p))['events']:
            if not e.get('ours') and e['date'] >= since_events:
                c[e['n']] += 1
    return c


def ours_certified(since):
    from .store import Store
    db = Store().db
    return sorted({r[0] for r in db.execute(
        "SELECT p.n FROM packing p JOIN sighting g ON g.packing = p.id WHERE g.source LIKE 'candidate:%' AND g.added >= ?",
        (since,))})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--since', required=True); ap.add_argument('--until')
    ap.add_argument('--plot'); ap.add_argument('--nmax', type=int, default=480)
    a = ap.parse_args()
    per_n, other, runs = attribute(a.since, a.until)
    tot = collections.Counter()
    for v in per_n.values():
        tot.update(v)
    lo = sum(sum(v.values()) for n, v in per_n.items() if n <= 324)
    hi = sum(sum(v.values()) for n, v in per_n.items() if n > 324)
    print(f'window {a.since} .. {a.until or "now"}: {len(runs)} runs')
    print('by kind (n-attributed): ' + ', '.join(f'{k} {v:.1f}' for k, v in tot.most_common()))
    print('not n-specific: ' + ', '.join(f'{k} {v:.1f}' for k, v in other.most_common()))
    print(f'n <= 324: {lo:.1f} CPU-h over {sum(1 for n in per_n if n <= 324)} n;  n > 324: {hi:.1f} CPU-h')
    ev = record_events()
    hot = [n for n in range(1, 325) if ev[n] >= 2]
    print(f"others' record events since 10-05 at n <= 324: {sum(ev[n] for n in range(325))}; CPU on n with >= 2 events: "
          f"{sum(sum(per_n[n].values()) for n in hot if n in per_n):.1f} of {lo:.1f}")
    if a.plot:
        plot(per_n, other, ev, ours_certified(a.since), a)


def plot(per_n, other, ev, ours, a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    kinds = sorted({k for v in per_n.values() for k in v})
    N = range(1, a.nmax + 1)
    fig, ax = plt.subplots(2, 1, figsize=(15, 7), sharex=True, gridspec_kw={'height_ratios': [3, 2]})
    bot = [0.0] * len(N)
    for k in kinds:
        v = [per_n.get(n, {}).get(k, 0.0) for n in N]
        ax[0].bar(N, v, bottom=bot, width=1.0, label=k)
        bot = [x + y for x, y in zip(bot, v)]
    for n in ours:
        if n <= a.nmax:
            ax[0].annotate('*', (n, bot[n - 1]), ha='center', va='bottom', color='red', fontsize=16)
    tot = sum(bot)
    ax[0].set_title(f'CPU-h by n, {a.since} .. {a.until or "now"}: {tot:.1f} by n + '
                    f'{sum(other.values()):.1f} not n-specific (red * = our certified candidates)')
    ax[0].set_ylabel('CPU-h'); ax[0].legend(ncol=6, loc='upper right')
    ax[1].bar(N, [ev[n] for n in N], width=1.0, color='tab:red')
    ax[1].set_ylabel("others' record events\nsince 10-05")
    for x in ax:
        x.axvline(324.5, color='k', ls=':')
    ax[1].set_xlabel('n')
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(a.plot)), exist_ok=True)
    plt.savefig(a.plot, dpi=110)
    print(a.plot)


if __name__ == '__main__':
    main()
