#!/usr/bin/env python3
"""Composition level 2 (blind time travel; 10-15): can the composer find an event's composition without being told the
pieces, and the explorer then squeeze it?  For an event set at cutoff T:

  1. compose2 run --blind T_local --blind-t T_utc --target cutoff:T_utc --probe P --only <ns>  (pieces: from-scratch
     rectsearch minima, Kingbird 09-30 crops, crops of store packings added before T_local or at the register frontier
     as of T_utc; every arm that fits n at frontier(T) + P is bisected to its smallest side S_min);
  2. per n, the --per-n tightest distinct near arms (layout at S_min, from near.jsonl) -> start files;
  3. explorer from them (harden_explore.py control arm: no harden), 1 proc x --minutes, queued --cores at a time.
Tiers vs the event's posted value / ante as usual.

  python3 analysis/l2_compose.py run --set 481 --out runs/l2/c481 [--probe 0.3 --budget 120 --minutes 20 --cores 16]
  python3 analysis/l2_compose.py score --out runs/l2/c481
"""
import argparse, json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE)
sys.path.insert(0, PK)


def event_set(name):
    sys.path.insert(0, PK)
    import comp_l1 as C
    import yaml
    E = yaml.safe_load(open(os.path.join(PK, '..', '..', 'tasks', 'composition-baseline', 'events.yaml')))
    E = E['events'] if isinstance(E, dict) else E
    evs = C.events(name)
    for ev in evs:
        e = next(x for x in E if x['n'] == ev['n'] and str(x['issue']).startswith(name))
        ev['T_utc'] = e['date'][:16]; ev['T_local'] = e['cutoff_local']; ev['new'] = ev['S']
    return evs


def run(a):
    from pk import manifest
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'l2_compose', args=vars(a), procs=a.cores)
    evs = event_set(a.set)
    if a.only:
        keep = {int(x) for x in a.only.split(',')}; evs = [e for e in evs if e['n'] in keep]
    T_utc, T_local = evs[0]['T_utc'], evs[0]['T_local']
    assert all(e['T_utc'] == T_utc for e in evs), 'one cutoff per run'
    ns = ','.join(str(e['n']) for e in evs)
    cdir = f'{a.out}/compose'
    if not os.path.exists(f'{cdir}/near.jsonl'):
        subprocess.run([sys.executable, f'{PK}/compose2.py', 'run', '--ns', '1-400', '--only', ns, '--blind', T_local,
                        '--blind-t', T_utc, '--target', f'cutoff:{T_utc}', '--probe', str(a.probe), '--budget', str(a.budget),
                        '--procs', str(a.cores), '--out', cdir], cwd=PK, check=True)
    near = {}
    for l in open(f'{cdir}/near.jsonl'):
        r = json.loads(l)
        if r.get('sq'):
            near.setdefault(r['n'], []).append(r)
    sdir = f'{a.out}/starts'
    for ev in evs:
        d = f"{sdir}/n{ev['n']}"; os.makedirs(d, exist_ok=True)
        R = sorted(near.get(ev['n'], []), key=lambda r: r['S_min'])
        k, seen = 0, []
        for r in R:
            if any(abs(r['S_min'] - t) < 1e-7 for t in seen):
                continue
            seen.append(r['S_min'])
            with open(f'{d}/start{k}.txt', 'w') as f:
                f.write(f"{ev['n']} {float(r['S_min'])!r}\n" + ''.join(f'{x!r} {y!r} {t!r}\n' for x, y, t in r['sq']))
            k += 1
            if k >= a.per_n:
                break
        ev['near_best'] = R[0]['S_min'] if R else None
        json.dump(ev, open(f'{d}/event.json', 'w'), indent=1)
        print(f"n={ev['n']}: {len(R)} near arms, tightest {ev['near_best']}  ante {ev['ante']:.10f} posted {ev['S']:.10f}", flush=True)
    subprocess.run([sys.executable, f'{PK}/harden_explore.py', 'run', '--src', sdir, '--starts-glob', 'start[0-9]*.txt',
                    '--out', f'{a.out}/explore', '--minutes', str(a.minutes), '--cores', str(a.cores)], cwd=PK, check=True)
    man.finish()


def score(a):
    subprocess.run([sys.executable, f'{PK}/harden_explore.py', 'score', '--out', f'{a.out}/explore'], cwd=PK)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run'); r.add_argument('--set', default='481'); r.add_argument('--out', required=True)
    r.add_argument('--only'); r.add_argument('--probe', type=float, default=0.3); r.add_argument('--budget', type=float, default=120)
    r.add_argument('--per-n', type=int, default=4); r.add_argument('--minutes', type=float, default=20)
    r.add_argument('--cores', type=int, default=16)
    s = sp.add_parser('score'); s.add_argument('--out', required=True)
    a = ap.parse_args()
    {'run': run, 'score': score}[a.cmd](a)
