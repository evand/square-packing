#!/usr/bin/env python3
"""Composition level 1 on SQUISH #500 (Mondrian v2, narrowed pieces) and the #481 events (10-15).

Given the named sources and piece sizes (tasks/composition-baseline/events.yaml), per event:
  1. crops = mondrian.crop_pieces(source, m) (8 orientations x up to R dropped rows / columns, tight corner cut);
  2. narrowing (#500: "pieces narrowed in their own rectangle"): each crop with all its squares (uprights included) is
     re-optimised in its own (h + d) x h rectangle, d = W - H fixed (rectlib's fq rectangle quench); the narrowed piece
     = its tilted squares in the smaller box.  Arms: plain crops, narrowed crops, both.
  3. mondrian.replicate (exact MILP upright fill + fq) at the composition sides near the posted S;
  4. explorer from the best distinct compositions, all events at once (gate-ref run), spare cores to replicates.
Tiers: REPLICATED (<= posted + 1e-9), BEAT-ANTE (< frontier_at(n, T) - 1e-9), no.

  python3 comp_l1.py run --set 500 --out runs/comp_l1/c500 [--minutes 15 --cores 16]
  python3 comp_l1.py score runs/comp_l1/c500
"""
import argparse, glob, json, math, os, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np                                                     # noqa: E402

# n: (T, posted S, (src A, m A), (src B, m B))   from events.yaml (#500, issue 2026-10-11T05:30Z)
SET500 = {
    181: ('store:212417', 37, 'store:212415', 42),
    238: ('store:212418', 68, 'store:212415', 42),
    239: ('store:92', 57, 'store:212416', 52),
    271: ('store:212418', 68, 'store:399', 47),
    272: ('store:399', 63, 'store:211856', 58),
    305: ('store:212414', 63, 'store:212414', 63),
    270: ('store:406', 209, 'store:211869', 3),
    304: ('store:406', 209, 'store:211783', 8),
}
T500 = '2026-10-11T05:30'


def events(which):
    from pk import eventbench as EB
    import yaml
    E = yaml.safe_load(open(os.path.join(HERE, '..', '..', 'tasks', 'composition-baseline', 'events.yaml')))
    E = E['events'] if isinstance(E, dict) else E
    out = []
    if which == '500':
        for n, (sa, ma, sb, mb) in SET500.items():
            e = next(x for x in E if x['n'] == n and str(x['issue']).startswith('500'))
            out.append(dict(n=n, T=T500, S=float(e['new']), ante=EB.frontier_at(n, T500), A=(sa, ma), B=(sb, mb)))
    elif which == '481':
        import mondrian as M
        for n, (S, (sa, ma), (sb, mb)) in M.LINEAGE_481.items():
            out.append(dict(n=n, T='2026-10-09T21:59', S=float(S), ante=float(M.ANTE_481[n]), A=(sa, ma), B=(sb, mb)))
    return out


def full_crop(src, desc):
    import mondrian as M
    sq, s = src
    o, i, j, us, vs, _ = desc
    q = M.orient(sq, s, o)
    h = M.half_ext(q[:, 2])
    x, y = q[:, 0] - i, q[:, 1] - j
    ok = (x - h >= -1e-9) & (y - h >= -1e-9) & (x + h <= us + 1e-7) & (y + h <= vs + 1e-7)
    return np.c_[x[ok], y[ok], q[ok, 2]], float(us), float(vs)


def narrow_job(args):
    """re-optimise one crop in its own rectangle at aspect offsets d + delta (delta < 0 narrows the width relative to the
    height, > 0 the height); -> (key, [(T', desc'), ...]) for every box no larger than the crop in both sides and
    smaller in one (the Pareto narrowing; #500: "searched toward a smaller rectangle of its own, by 7e-06 to 5e-02 in
    one side").  delta = 0 runs ALM + polish, other deltas polish only (stay in the crop's basin)."""
    key, P, W, H, desc, deltas = args
    import rectlib as RL
    import mondrian as M
    d0 = W - H
    out = []
    o, i, j, _, _, m = desc
    with tempfile.TemporaryDirectory() as td:
        inp, outp = os.path.join(td, 'in.txt'), os.path.join(td, 'out.txt')
        for dl in deltas:
            d = d0 + dl
            s0 = max(H, W - d)
            with open(inp, 'w') as f:
                f.write(f'{len(P)} {float(s0)!r}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(a)!r}\n' for x, y, a in P))
            best = None
            for extra in (([], ['--no-alm']) if dl == 0 else (['--no-alm'],)):
                r, _ = RL._fq(inp, outp, d, extra)
                if r is not None and (best is None or r[0] < best[0]):
                    best = r
            if best is None:
                continue
            h = best[0]; w = h + d
            if h <= H + 1e-9 and w <= W + 1e-9 and (h < H - 1e-9 or w < W - 1e-9):
                sq = best[2]
                out.append((sq[M.is_tilted(sq[:, 2])], (o, i, j, float(w), float(h), int(m))))
    return key, out


def run(a):
    import mondrian as M
    from pk import manifest
    from jobpool import run_jobs
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'comp_l1', args=vars(a), procs=a.cores)
    evs = events(a.set)
    if a.only:
        keep = {int(x) for x in a.only.split(',')}; evs = [e for e in evs if e['n'] in keep]
    # 1-2: crops and narrowing, every piece of every event in one queue
    jobs, pieces = [], {}
    deltas = [float(x) for x in a.narrow_deltas.split(',')]
    for ev in evs:
        for side in ('A', 'B'):
            spec, m = ev[side]
            src = M.load_src(spec)
            P = M.crop_pieces(src, m, a.R, 0)[:a.max_pieces]
            pieces[(ev['n'], side)] = P
            if a.arm in ('narrow', 'both'):
                for k, (T, desc) in enumerate(P):
                    Pf, W, H = full_crop(src, desc)
                    jobs.append(((ev['n'], side, k), Pf, W, H, desc, deltas))
    print(f'{len(evs)} events; {sum(len(v) for v in pieces.values())} crops; {len(jobs)} narrowings', flush=True)
    nar = {}
    for r in run_jobs(narrow_job, jobs, procs=a.cores, timeout=900):
        if isinstance(r, tuple) and r[1]:
            nar.setdefault(r[0][:2], []).extend(r[1])
    # 3: compose per event (mondrian.replicate pools internally)
    starts = {}
    for ev in evs:
        n, S = ev['n'], ev['S']
        PA = (pieces[(n, 'A')] if a.arm in ('plain', 'both') else []) + nar.get((n, 'A'), [])
        PB = (pieces[(n, 'B')] if a.arm in ('plain', 'both') else []) + nar.get((n, 'B'), [])
        ev['nA'], ev['nB'], ev['narrowedA'], ev['narrowedB'] = len(PA), len(PB), len(nar.get((n, 'A'), [])), len(nar.get((n, 'B'), []))
        comp, raws = [], []
        for Sc in M.comp_sides(S, PA + PB, a.window, a.extra):
            for X, Y in ((PA, PB), (PB, PA)) if ev['A'] != ev['B'] else ((PA, PB),):
                for r in M.replicate(n, Sc, X, Y, procs=a.cores, top=a.top, fq_top=a.fq_top, log=lambda *x: None):
                    comp.append((r['polish'], r['packing']))
                    raws.append((r['polish'], Sc, r['raw']))
        comp.sort(key=lambda c: c[0])
        d = f'{a.out}/n{n}'; os.makedirs(d, exist_ok=True)
        st, seen = [], []
        for sp, pk in comp:
            if all(abs(sp - t) > 1e-9 for t in seen):
                seen.append(sp); f = f'{d}/start{len(st)}.txt'; pk.write(f); st.append(f)
            if len(st) >= 4:
                break
        raws.sort(key=lambda r: r[1])                       # raw composed layouts (before fq), tightest side first
        for k, (_, Sc, p) in enumerate(raws[:8]):
            p.write(f'{d}/raw{k}.txt')
        ev['compose_best'] = comp[0][0] if comp else None
        json.dump(ev, open(f'{d}/event.json', 'w'), indent=1)
        print(f"n={n}: {len(PA)} x {len(PB)} pieces (narrowed {ev['narrowedA']} / {ev['narrowedB']}); composed {len(comp)}, "
              f"best {ev['compose_best']}  posted {S:.10f} ante {ev['ante']:.10f}", flush=True)
        if st:
            starts[n] = st
    # 4: explorers, all at once
    runs = [(ev, starts[ev['n']]) for ev in evs if ev['n'] in starts and a.minutes > 0]
    slots, r = [(ev, st, 0) for ev, st in runs], 1
    while len(slots) < a.cores and runs:
        slots += [(ev, st, r) for ev, st in runs][:a.cores - len(slots)]; r += 1
    P = []
    for ev, st, rep in slots:
        d = f'{a.out}/n{ev["n"]}'
        cmd = [sys.executable, f'{HERE}/explore.py', '--n', str(ev['n']), '--adapt', '--reward', 'value', '--elite-share', '0.3',
               '--polish-gate', '2e-3', '--gate-ref', 'run', '--frontier-s', repr(ev['ante']), '--starts', *st,
               '--minutes', str(a.minutes), '--procs', '1', '--seed', str(ev['n'] * 10 + rep), '--out', f'{d}/explore{rep}']
        P.append(subprocess.Popen(cmd, stdout=open(f'{d}/explore{rep}.log', 'w'), stderr=subprocess.STDOUT, cwd=HERE))
    for p in P:
        p.wait()
    man.finish()
    score(a)


def score(a):
    from collections import Counter
    C = Counter()
    for f in sorted(glob.glob(f'{a.out}/n*/event.json')):
        ev = json.load(open(f)); d = os.path.dirname(f); best = ev.get('compose_best')
        for ad in glob.glob(f'{d}/explore*/archive.jsonl'):
            for l in open(ad):
                e = json.loads(l)
                if not e.get('unpolished') and (best is None or e['s'] < best):
                    best = e['s']
        tier = 'REPLICATED' if best is not None and best <= ev['S'] + 1e-9 else \
            'BEAT-ANTE' if best is not None and best < ev['ante'] - 1e-9 else 'no'
        C[tier] += 1
        print(f"n={ev['n']:3d} composed {ev.get('compose_best') or float('nan'):.10f} final {best or float('nan'):.10f} "
              f"posted {ev['S']:.10f} ante {ev['ante']:.10f} {tier}")
    print(dict(C))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run'); r.add_argument('--set', default='500'); r.add_argument('--out', required=True)
    r.add_argument('--only'); r.add_argument('--arm', default='both', choices=['plain', 'narrow', 'both'])
    r.add_argument('--narrow-deltas', default='0', help="aspect offsets for narrowing, e.g. '0,-0.01,-0.03,0.01,0.03'")
    r.add_argument('--R', type=int, default=3); r.add_argument('--max-pieces', type=int, default=400)
    r.add_argument('--window', type=float, default=0.1); r.add_argument('--extra', type=int, default=2)
    r.add_argument('--top', type=int, default=150); r.add_argument('--fq-top', type=int, default=4)
    r.add_argument('--minutes', type=float, default=15); r.add_argument('--cores', type=int, default=16)
    s = sp.add_parser('score'); s.add_argument('out')
    a = ap.parse_args()
    {'run': run, 'score': score}[a.cmd](a)
