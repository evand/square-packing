#!/usr/bin/env python3
"""Explorer from harden outputs (10-15, BM3): compose -> harden -> explorer, the "nearby search" after the squeeze.

For each event (event.json with n, ante, new / S) in the given harden output dirs (reach_harden.py layout; several dirs
per event are pooled by basename), take the --keep best distinct harden outputs (results.jsonl 'out', side gap > 2e-9),
then run explore.py (current config: gate 2e-3 gate-ref run, value reward, elite 0.3, frontier_s = ante) from them,
1 proc each, all events at once, spare cores to extra replicates.  Tiers as removal_bench / comp_l1.

  python3 harden_explore.py run --src runs/bm1 runs/bm2 --out runs/bm3 [--minutes 20 --cores 16 --keep 4]
  python3 harden_explore.py score --out runs/bm3
"""
import argparse, glob, json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)


def gather(srcs, starts_glob=None):
    E = {}
    for s in srcs:
        if starts_glob:                                  # control arm: rank plain start files by their header side
            for ej in sorted(glob.glob(f'{s}/*/event.json')):
                d = os.path.dirname(ej); b = os.path.basename(d)
                ev = json.load(open(ej)); ev.setdefault('new', ev.get('S'))
                e = E.setdefault(b, dict(ev=ev, res=[]))
                for g in starts_glob.split(','):
                    for f in glob.glob(f'{d}/{g}'):
                        e['res'].append((float(open(f).readline().split()[1]), f))
            continue
        for ej in sorted(glob.glob(f'{s}/*/event.json')):
            d = os.path.dirname(ej); b = os.path.basename(d)
            ev = json.load(open(ej)); ev.setdefault('new', ev.get('S'))
            e = E.setdefault(b, dict(ev=ev, res=[]))
            if os.path.exists(f'{d}/results.jsonl'):
                for l in open(f'{d}/results.jsonl'):
                    try:
                        r = json.loads(l)
                    except Exception:
                        continue
                    if r.get('ok') and isinstance(r.get('s'), float) and os.path.exists(r.get('out', '')):
                        e['res'].append((r['s'], r['out']))
    return E


def run(a):
    from pk import manifest
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'harden_explore', args=vars(a), procs=a.cores)
    runs = []
    for b, e in sorted(gather(a.src, a.starts_glob).items()):
        if a.only and int(e['ev']['n']) not in {int(x) for x in a.only.split(',')}:
            continue
        ev = e['ev']; d = f'{a.out}/{b}'; os.makedirs(d, exist_ok=True)
        keep, seen = [], []
        for s, f in sorted(e['res']):
            if any(abs(s - t) < 2e-9 for t in seen):
                continue
            seen.append(s); dst = f'{d}/start_{len(keep):02d}.txt'
            open(dst, 'w').write(open(f).read()); keep.append(dst)
            if len(keep) >= a.keep:
                break
        ev['harden_best'] = seen[0] if seen else None
        json.dump(ev, open(f'{d}/event.json', 'w'), indent=1)
        if keep:
            runs.append((ev, d, keep))
    slots, r = [(ev, d, k, 0) for ev, d, k in runs], 1
    while len(slots) < a.cores and runs:
        slots += [(ev, d, k, r) for ev, d, k in runs][:a.cores - len(slots)]; r += 1
    if a.max_need is not None:                           # live hunts: only events whose best start is near
        slots = [x for x in slots if x[0].get('harden_best') is not None and x[0]['harden_best'] - x[0]['ante'] < a.max_need]
    P, queue = [], list(slots)
    import time
    while queue or P:                                    # at most --cores explorers at a time (events beyond: queued)
        P = [p for p in P if p.poll() is None]
        while queue and len(P) < a.cores:
            ev, d, keep, rep = queue.pop(0)
            cmd = [sys.executable, f'{HERE}/explore.py', '--n', str(ev['n']), '--adapt', '--reward', 'value', '--elite-share', '0.3',
                   '--polish-gate', '2e-3', '--gate-ref', 'run', '--frontier-s', repr(ev['ante']), '--starts', *keep,
                   '--minutes', str(a.minutes), '--procs', '1', '--seed', str(ev['n'] * 10 + rep), '--out', f'{d}/explore{rep}']
            P.append(subprocess.Popen(cmd, stdout=open(f'{d}/explore{rep}.log', 'w'), stderr=subprocess.STDOUT, cwd=HERE))
        time.sleep(5)
    man.finish()
    score(a)


def score(a):
    from collections import Counter
    C = Counter()
    for ej in sorted(glob.glob(f'{a.out}/*/event.json')):
        ev = json.load(open(ej)); d = os.path.dirname(ej); best = ev.get('harden_best')
        for ad in glob.glob(f'{d}/explore*/archive.jsonl'):
            for l in open(ad):
                e = json.loads(l)
                if not e.get('unpolished') and (best is None or e['s'] < best):
                    best = e['s']
        tier = 'REPLICATED' if best is not None and best <= ev['new'] + 1e-9 else \
            'BEAT-ANTE' if best is not None and best < ev['ante'] - 1e-9 else 'no'
        C[tier] += 1
        hb = ev.get('harden_best')
        print(f"n={ev['n']:3d} harden {hb if hb is not None else float('nan'):.10f} final {best if best is not None else float('nan'):.10f} "
              f"need {(best - ev['ante']) if best is not None else float('nan'):+.2e}  posted {ev['new']:.10f} ante {ev['ante']:.10f} {tier}")
    print(dict(C))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run'); r.add_argument('--src', nargs='+', required=True); r.add_argument('--out', required=True)
    r.add_argument('--minutes', type=float, default=20); r.add_argument('--cores', type=int, default=16)
    r.add_argument('--keep', type=int, default=4)
    r.add_argument('--starts-glob', help='control: use these start files in each src event dir (no harden results)')
    r.add_argument('--only', help='comma list of n')
    r.add_argument('--max-need', type=float, help='skip events whose best start is more than this above ante')
    s = sp.add_parser('score'); s.add_argument('--out', required=True)
    a = ap.parse_args()
    {'run': run, 'score': score}[a.cmd](a)
