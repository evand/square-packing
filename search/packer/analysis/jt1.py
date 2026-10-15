#!/usr/bin/env python3
"""JT1 (10-15): time travel on family J (Xu #432, LLM-evolved from-scratch programs, 10-08).  Can our explorer reach
Xu's values from the frontier packing before #432 (the ante)?  Top-10 J events by gain, n <= 324.

Per event: start = the store row at frontier_at(n, T) (the register packing at T), explorer (current config: gate 2e-3,
gate-ref run, value reward, elite 0.3, frontier_s = ante), --reps x 1 proc, all events at once, --minutes.
Tiers: REPLICATED (<= Xu + 1e-9), BEAT-ANTE (< ante - 1e-9), no.

  python3 analysis/jt1.py run [--minutes 20 --reps 1]   ->  runs/jt1
  python3 analysis/jt1.py score
"""
import argparse, glob, json, os, sqlite3, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE)
sys.path.insert(0, PK)
OUT = os.path.join(PK, 'runs', 'jt1')


def events(k=10):
    from pk import eventbench as EB
    E = json.load(open(os.path.join(PK, 'bench', 'events.json'))); E = E['events'] if isinstance(E, dict) else E
    J = sorted([e for e in E if e['family'] == 'J' and e['n'] <= 324], key=lambda e: -e['gain'])[:k]
    db = sqlite3.connect(f"file:{os.path.join(PK, 'runs', 'store.sqlite')}?mode=ro", uri=True)
    out = []
    for e in J:
        a = EB.frontier_at(e['n'], e['date'])
        r = db.execute("SELECT id FROM packing WHERE n=? AND abs(s-?)<2e-9 AND status!='infeasible' ORDER BY added LIMIT 1",
                       (e['n'], a)).fetchone()
        out.append(dict(n=e['n'], T=e['date'], ante=a, new=e['new'], ante_id=r and r[0]))
    return out


def run(a):
    from pk import manifest
    from pk.store import Store
    import mcmin
    os.makedirs(OUT, exist_ok=True)
    evs = [e for e in events(a.k) if e['ante_id']]
    man = manifest.start(OUT, 'jt1', args=vars(a), procs=len(evs) * a.reps)
    st = Store(); P = []
    for ev in evs:
        d = f"{OUT}/n{ev['n']}"; os.makedirs(d, exist_ok=True)
        p = st.get(ev['ante_id']); f = f'{d}/start.txt'
        mcmin.write_deg(f, float(p.s), [tuple(map(float, q)) for q in p.sq])
        json.dump(ev, open(f'{d}/event.json', 'w'), indent=1)
        for rep in range(a.reps):
            cmd = [sys.executable, f'{PK}/explore.py', '--n', str(ev['n']), '--adapt', '--reward', 'value', '--elite-share', '0.3',
                   '--polish-gate', '2e-3', '--gate-ref', 'run', '--frontier-s', repr(ev['ante']), '--starts', f,
                   '--minutes', str(a.minutes), '--procs', '1', '--seed', str(ev['n'] * 10 + rep), '--out', f'{d}/explore{rep}']
            P.append(subprocess.Popen(cmd, stdout=open(f'{d}/explore{rep}.log', 'w'), stderr=subprocess.STDOUT, cwd=PK))
    for p in P:
        p.wait()
    man.finish()
    score(a)


def score(a):
    from collections import Counter
    C = Counter()
    for ej in sorted(glob.glob(f'{OUT}/n*/event.json')):
        ev = json.load(open(ej)); d = os.path.dirname(ej); best = None
        for ad in glob.glob(f'{d}/explore*/archive.jsonl'):
            for l in open(ad):
                e = json.loads(l)
                if not e.get('unpolished') and (best is None or e['s'] < best):
                    best = e['s']
        tier = 'REPLICATED' if best is not None and best <= ev['new'] + 1e-9 else \
            'BEAT-ANTE' if best is not None and best < ev['ante'] - 1e-9 else 'no'
        C[tier] += 1
        print(f"n={ev['n']:3d} ante {ev['ante']:.10f} xu {ev['new']:.10f} ours {best if best is not None else float('nan'):.10f} "
              f"frac of Xu's gain {((ev['ante'] - best) / (ev['ante'] - ev['new'])) if best is not None else float('nan'):+.2f} {tier}")
    print(dict(C))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('cmd', choices=['run', 'score'])
    ap.add_argument('--minutes', type=float, default=20); ap.add_argument('--reps', type=int, default=1)
    ap.add_argument('--k', type=int, default=10)
    a = ap.parse_args()
    {'run': run, 'score': score}[a.cmd](a)
