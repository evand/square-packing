#!/usr/bin/env python3
"""Evaluator CLI for Claude Code subagents writing square-packing constructors (EV2, prompt evolution; 10-15).

  python3 analysis/evolve_eval.py eval PROG.py --n 126,127 --out DIR [--target live|cutoff:ISO] [--seeds 2] [--tag P003]
      runs PROG's construct(n, rng) (sandboxed: analysis/evolve_runner.py), optimises every candidate (fq squeeze +
      SLP polish), prints per n: best side, target, gap, normalised gap (0 = ties target, 1 = grid), a description, and
      writes DIR/n<N>.png (rendering of the best optimised layout) and DIR/n<N>.txt (its coordinates).  Every call is
      logged to runs/evolve2/evals.jsonl (with --tag).  A result below the target prints CANDIDATE RECORD.
  python3 analysis/evolve_eval.py ref --n 127 --out DIR [--target live|cutoff:ISO]
      renders the best known packing at the target time (store / register) to DIR/ref_n<N>.png for reference.
"""
import argparse, json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE); sys.path.insert(0, PK); sys.path.insert(0, HERE)
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
LOG = os.path.join(PK, 'runs', 'evolve2', 'evals.jsonl')


def target_of(n, spec):
    from pk import eventbench as EB
    if spec == 'live':
        return EB.frontier_at(n, None)
    return EB.frontier_at(n, spec.split(':', 1)[1])


def render(s, sq, path, title=''):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    fig, ax = plt.subplots(figsize=(6, 6))
    for x, y, a in sq:
        a = a % 90.0
        th = math.radians(a)
        c, si = math.cos(th), math.sin(th)
        pts = [(x + dx * c - dy * si, y + dx * si + dy * c) for dx, dy in ((.5, .5), (-.5, .5), (-.5, -.5), (.5, -.5))]
        tilt = min(a, 90 - a)
        col = '#d0d0d0' if tilt < 1 else plt.cm.viridis(min(tilt, 45) / 45)
        ax.add_patch(Polygon(pts, closed=True, fc=col, ec='k', lw=0.4))
    ax.plot([0, s, s, 0, 0], [0, 0, s, s, 0], 'r-', lw=1)
    ax.set_xlim(-0.2, s + 0.2); ax.set_ylim(-0.2, s + 0.2); ax.set_aspect('equal')
    ax.set_xticks(range(int(s) + 1)); ax.set_yticks(range(int(s) + 1)); ax.tick_params(labelsize=6)
    ax.grid(True, lw=0.2)
    ax.set_title(title + '  (grey: axis-aligned; colour: tilt 0-45 deg, yellow = 45)', fontsize=7)
    fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)


def cmd_eval(a):
    import evolve as E
    from concurrent.futures import ProcessPoolExecutor
    ns = [int(x) for x in a.n.split(',')]
    os.makedirs(a.out, exist_ok=True)
    prog = os.path.abspath(a.prog)
    bad = E.static_check(open(prog).read())
    if bad:
        print(f'REJECTED by the static check: {bad}'); return
    t0 = time.time()
    jobs = [(prog, n, 1000 * sd + n, a.maxc, (1.0, 1.03), a.fill) for n in ns for sd in range(a.seeds)]
    with ProcessPoolExecutor(min(a.procs, len(jobs))) as ex:
        R = list(ex.map(E.eval_job, jobs))
    best = {}
    errs = {}
    mins = {}
    for (_, n, *_), r in zip(jobs, R):
        for m in r.get('minima', []):
            if all(abs(m[0] - x[0]) > 1e-9 for x in mins.setdefault(n, [])):
                mins[n].append(m)
        if 's' in r:
            if n not in best or r['s'] < best[n]['s']:
                best[n] = r
        else:
            errs.setdefault(n, r['err'])
    rows = []
    for n in ns:
        T = target_of(n, a.target); g = math.ceil(math.sqrt(n) - 1e-12)
        if n not in best:
            print(f'n={n}: ERROR {errs.get(n)}'); rows.append(dict(n=n, err=errs.get(n))); continue
        r = best[n]; ng = (r['s'] - T) / max(g - T, 1e-9)
        flag = '  *** CANDIDATE RECORD (below target) ***' if r['s'] < T - 1e-9 else ''
        print(f"n={n}: side {r['s']:.9f}  target {T:.9f}  gap {r['s'] - T:+.2e}  normalised {ng:+.4f}  grid {g}  "
              f"({r['cands']} candidates)  {r['desc']}{flag}")
        with open(f'{a.out}/n{n}.txt', 'w') as f:
            f.write(f"{n} {r['s']!r}\n" + ''.join(f'{x!r} {y!r} {t!r}\n' for x, y, t in r['sq']))
        render(r['s'], r['sq'], f'{a.out}/n{n}.png', f"n={n} s={r['s']:.6f} target {T:.6f}")
        rows.append(dict(n=n, s=r['s'], target=T, ng=ng))
    ok = [r['ng'] for r in rows if 'ng' in r]
    if ok:
        T = 0.02                                          # soft-min: one great n counts (Evan 10-15), the mean for reference
        sm = -T * math.log(np.mean(np.exp(-np.array(ok) / T)))
        print(f'SCORE  softmin(T={T}) {sm:+.4f}   best {min(ok):+.4f}   mean {np.mean(ok):+.4f}   '
              f'(normalised gaps; failed n count as missing)')
    archive_minima(a, prog, mins)
    print(f'renderings: {a.out}/n<N>.png   coordinates: {a.out}/n<N>.txt   ({time.time() - t0:.0f} s)')
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with open(LOG, 'a') as f:
        f.write(json.dumps(dict(t=time.strftime('%Y-%m-%dT%H:%M:%S'), tag=a.tag, prog=prog, target=a.target, rows=rows,
                                sec=round(time.time() - t0, 1))) + '\n')


def archive_minima(a, prog, mins):
    """every distinct non-grid minimum (>= 5 tilted, feasible) -> <out>/n<N>/archive.jsonl in explore.py's format, so
    census_collect.py can certify / catalogue them like explorer runs (Evan 10-15: keep novel slightly-worse minima,
    as in the s(110) census); plus one line each in runs/evolve2/minima.jsonl (the cross-run index)."""
    import mcmin
    idx = os.path.join(PK, 'runs', 'evolve2', 'minima.jsonl')
    with open(idx, 'a') as fi:
        for n, L in mins.items():
            d = os.path.abspath(f'{a.out}/n{n}'); os.makedirs(d, exist_ok=True)
            T = target_of(n, a.target); g = math.ceil(math.sqrt(n) - 1e-12)
            with open(f'{d}/archive.jsonl', 'w') as fa:
                for i, (s, sq) in enumerate(sorted(L)):
                    path = f'{d}/b{i:05d}.txt'
                    mcmin.write_deg(path, s, [tuple(q) for q in sq])
                    a_ = np.asarray(sq)[:, 2] % 90.0
                    nt = int((np.minimum(a_, 90 - a_) > 1.0).sum())
                    fa.write(json.dumps(dict(i=i, s=s, path=path, kind=f'evolve:{a.tag}', t=0, known=False, parent=-1)) + '\n')
                    fi.write(json.dumps(dict(n=n, s=s, gap=s - T, ng=(s - T) / max(g - T, 1e-9), tilted=nt, path=path,
                                             tag=a.tag, prog=prog, t=time.strftime('%Y-%m-%dT%H:%M:%S'))) + '\n')


def cmd_near(a):
    """holdout references (Evan 10-15): k valid packings of n that are WORSE than the known record (gap to the target in
    [--lo, --hi]: same neighbourhood of structures, not the record's own basin), chosen for diversity (store cluster,
    then tilted count); rendered to DIR/near_n<N>_<i>.png.  The record itself is never shown.  The ids shown are logged
    to DIR/near_manifest.json for the leakage audit."""
    import sqlite3
    from pk.store import Store
    os.makedirs(a.out, exist_ok=True)
    db = sqlite3.connect(f"file:{os.path.join(PK, 'runs', 'store.sqlite')}?mode=ro", uri=True)
    st, man = Store(), {}
    for n in [int(x) for x in a.n.split(',')]:
        T = target_of(n, a.target)
        rows = db.execute("SELECT id, s, cluster, ntilted FROM packing WHERE n=? AND status!='infeasible' AND s>=? AND s<=? "
                          "ORDER BY s", (n, T + a.lo, T + a.hi)).fetchall()
        pick, seen_c, seen_s = [], set(), []
        for rnd in (0, 1):                               # first one per cluster / tilted count, then fill up
            for pid, s_, cl, nt in rows:
                key = (cl, nt)
                if len(pick) >= a.k or any(abs(s_ - t) < 1e-6 for t in seen_s):
                    continue
                if rnd == 0 and (key in seen_c or (cl is not None and any(c == cl for c, _ in seen_c))):
                    continue
                pick.append((pid, s_)); seen_c.add(key); seen_s.append(s_)
        rng = np.random.default_rng(n)
        for i, (pid, s_) in enumerate(pick):
            p = st.get(pid)
            render(float(p.s), np.asarray(p.sq, float), f'{a.out}/near_n{n}_{i}.png',
                   f'n={n}: a valid packing, s={s_:.6f} (the known record is {s_ - T:.4f} smaller)')
        man[n] = dict(target=T, shown=[dict(id=pid, s=s_) for pid, s_ in pick])
        print(f'n={n}: {len(pick)} near-record references ({len(rows)} in the gap band) -> {a.out}/near_n{n}_<i>.png')
    json.dump(man, open(f'{a.out}/near_manifest.json', 'w'), indent=1)


def cmd_ref(a):
    if not a.allow_best:
        print('ref (best known) is disabled for holdout runs; use `near` (pass --allow-best only when told to)'); return
    from pk.store import Store
    import sqlite3
    os.makedirs(a.out, exist_ok=True)
    db = sqlite3.connect(f"file:{os.path.join(PK, 'runs', 'store.sqlite')}?mode=ro", uri=True)
    st = Store()
    for n in [int(x) for x in a.n.split(',')]:
        T = target_of(n, a.target)
        r = db.execute("SELECT id FROM packing WHERE n=? AND abs(s-?)<2e-9 AND status!='infeasible' ORDER BY added LIMIT 1",
                       (n, T)).fetchone()
        if not r:
            print(f'n={n}: no stored packing at {T}'); continue
        p = st.get(r[0])
        render(float(p.s), np.asarray(p.sq, float), f'{a.out}/ref_n{n}.png', f'best known n={n} s={T:.6f}')
        print(f'n={n}: {a.out}/ref_n{n}.png  (s = {T:.9f})')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    e = sp.add_parser('eval'); e.add_argument('prog'); e.add_argument('--n', required=True); e.add_argument('--out', required=True)
    e.add_argument('--target', default='live'); e.add_argument('--seeds', type=int, default=2)
    e.add_argument('--maxc', type=int, default=24); e.add_argument('--procs', type=int, default=3); e.add_argument('--tag', default='')
    e.add_argument('--fill', action='store_true', help='candidates may hold only the tilted structure (< n squares): the exact MILP upright fill (mondrian.exact_fill) adds axis-aligned squares at the given s (retrying at s x 1.01 / 1.02 / 1.04 if short)')
    r = sp.add_parser('ref'); r.add_argument('--n', required=True); r.add_argument('--out', required=True)
    r.add_argument('--target', default='live'); r.add_argument('--allow-best', action='store_true')
    q = sp.add_parser('near'); q.add_argument('--n', required=True); q.add_argument('--out', required=True)
    q.add_argument('--target', default='live'); q.add_argument('--k', type=int, default=4)
    q.add_argument('--lo', type=float, default=5e-3); q.add_argument('--hi', type=float, default=5e-2)
    a = ap.parse_args()
    {'eval': cmd_eval, 'ref': cmd_ref, 'near': cmd_near}[a.cmd](a)
