#!/usr/bin/env python3
"""Fresh-packing pass (PLAN-2026-10-12 item 7; 10-13 plan item 1): newly posted packings -> pieces -> composer ->
polish -> certified candidates, in one command.

Why: all three 10-13 records (339, 364, 375) were compositions using crops of packings posted hours earlier (Ostashev's
s(304), BuddenPing's s(146), Xu's s(129)), found by a broad sweep.  This pass does only the new part: every arm uses at
least one piece cut from a new packing, so it is cheap and can run after every sync.

  python3 fresh.py run --since 2026-10-13T08:00 --out runs/fresh/f1 [--ids 291615 ...] [--procs 16] [--budget 30]
  python3 fresh.py new --since ISO          list the new packings only

Steps:
  1. new packings: store packings whose first public sighting (register, pending issue, catalogue, outside pool) is at
     or after --since (or --ids), feasible, n >= 10.  Sync first (pk.py sync-register / sync-pending / sync-kingbird).
  2. pieces: corner crops (mondrian.pieces_of, Pareto-pruned) re-optimised as rectangles (rectlib.reopt, fq), + each
     whole packing as a d = 0 piece; all into the store's piece table (pk.pieces).
  3. compose: compose2 run --new-from (arms must use one of the new packings' pieces), n 20..324 below the grid and
     325..NMAX (--grid), target frontier_now (L-step closed), near-miss probe.
  4. hits (composed side below frontier_now) -> seam polish (fq) -> candidate.py (exactsolve, both checkers, live check).
  summary.json: new packings, pieces (new / dup), composer targets and hits, certified candidates, measured CPU.
Nothing is posted (CLAUDE.md): a certified candidate gets a draft in outreach/ by hand.
"""
from __future__ import annotations
import argparse, json, math, os, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor
from types import SimpleNamespace

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PUBLIC = ('register', 'pending:', 'kingbird', 'pool:', 'watch:')   # watch: = pk/watch.py repo polls (10-14)


def new_packings(since=None, ids=None, nmin=10):
    """[(id, n, s, label)] of packings first seen publicly at or after `since` (or the given ids)."""
    from pk.store import Store
    db = Store().db
    if ids:
        q = f"SELECT id, n, s FROM packing WHERE id IN ({','.join('?' * len(ids))})"
        rows = [(r['id'], r['n'], r['s'], 'ids') for r in db.execute(q, list(ids))]
    else:
        pub = ' OR '.join(f"g.source LIKE '{p}%'" for p in PUBLIC)
        rows = db.execute(f"""SELECT p.id, p.n, p.s, min(g.added) AS first, group_concat(DISTINCT g.source) AS src
            FROM packing p JOIN sighting g ON g.packing = p.id
            WHERE p.status != 'infeasible' AND (p.pen IS NULL OR p.pen <= 1e-9) AND ({pub})
            GROUP BY p.id HAVING first >= ?""", (since,)).fetchall()
        rows = [(r['id'], r['n'], r['s'], r['src'].split(',')[0]) for r in rows]
    return [(i, n, s, f'store:{i}:{src}') for i, n, s, src in rows if n >= nmin]


def make_atlas(rows, path, procs):
    """mondrian atlas (W, H, m, o, sid + sources) of the corner crops of `rows`."""
    import mondrian
    from pk.store import Store
    st = Store()
    src = [dict(label=lab, n=n, s=float(s), sq=np.array(st.get(i).sq)) for i, n, s, lab in rows]
    parts = []
    with ProcessPoolExecutor(procs) as ex:
        for k, r in ex.map(mondrian._pieces_job, [(k, p['s'], p['sq']) for k, p in enumerate(src)]):
            if r is not None:
                parts.append(r + (np.full(len(r[0]), k),))
    if not parts:
        return 0
    W, H, m, o, sid = (np.concatenate(z) for z in zip(*parts))
    W, H, m, o, sid = mondrian.pareto_prune(W, H, m, o, sid)
    np.savez_compressed(path, W=W, H=H, m=m, o=o, sid=sid, labels=np.array([p['label'] for p in src]),
                        src_n=np.array([p['n'] for p in src]), src_s=np.array([p['s'] for p in src]),
                        sq=mondrian._objarr([p['sq'] for p in src]))
    return len(W)


def add_whole(rows):
    """Each new packing as a d = 0 piece (s(327) = s(123) + 7 and 364 = s(146) + 7 used whole packings)."""
    from pk.packing import Packing
    from pk.store import Store
    st = Store()
    new = 0
    for i, n, s, lab in rows:
        p = st.get(i)
        _, c = st.add_piece(Packing(float(p.s), np.array(p.sq), dx=0.0), lab, ref=f'packing:{i}', commit=False)
        new += c
    st.db.commit()
    return new


def sh(cmd, log):
    with open(log, 'w') as f:
        return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=HERE).returncode


def run(a):
    from pk import manifest, pieces, eventbench as E
    import rectlib, scratch
    os.makedirs(a.out, exist_ok=True)
    man = manifest.start(a.out, 'fresh', args=vars(a), procs=a.procs)
    scratch.scope('fresh')
    out = lambda *p: os.path.join(os.path.abspath(a.out), *p)
    t0 = time.time()
    rows = new_packings(a.since, a.ids)
    print(f'{len(rows)} new packings: ' + ' '.join(f'{n}:{s:.6f}' for _, n, s, _ in sorted(rows, key=lambda r: r[1])),
          flush=True)
    S = dict(new=[dict(id=i, n=n, s=s, label=l) for i, n, s, l in rows])
    if not rows:
        json.dump(S, open(out('summary.json'), 'w'), indent=1); man.finish(summary=S); return
    json.dump([r[0] for r in rows], open(out('new_ids.json'), 'w'))
    # 2. pieces
    nat = make_atlas(rows, out('atlas.npz'), a.procs)
    S['crops'] = nat
    if nat:
        rectlib.cmd_build(SimpleNamespace(atlas=out('atlas.npz'), out=out('lib'), mmin=8, mmax=a.mmax, hmax=a.hmax,
                                          per_cell=1, procs=a.procs, limit=0, dmax=8, seed=1))
        S['pieces_seen'], S['pieces_new'] = pieces.import_paths([out('lib', 'pieces.jsonl')])
    S['whole_new'] = add_whole(rows)
    S['t_pieces'] = round(time.time() - t0)
    print(f"pieces: {S.get('pieces_seen', 0)} re-optimised crops ({S.get('pieces_new', 0)} new), "
          f"{S['whole_new']} whole packings new ({S['t_pieces']} s)", flush=True)
    # 3. compose (two ranges: below the grid at n <= 324; n > 324 with --grid, the L-step frontier is the bar there)
    base = [sys.executable, 'compose2.py', 'run', '--new-from', out('new_ids.json'), '--probe', str(a.probe),
            '--budget', str(a.budget), '--procs', str(a.procs), '--target', a.target]
    if a.exclude_src:
        base += ['--exclude-src'] + a.exclude_src
    for tag, ns, extra in (('lo', f'20-{min(324, a.nmax)}', []), ('hi', f'325-{a.nmax}', ['--grid'])):
        if a.nmax < 325 and tag == 'hi':
            continue
        sh(base + ['--ns', ns, '--out', out(f'compose_{tag}')] + extra, out(f'compose_{tag}.log'))
    # 4. hits -> candidate.py
    E.set_public(a.target == 'public', a.exclude_src)
    hits = []
    for tag in ('lo', 'hi'):
        p = out(f'compose_{tag}', 'results.jsonl')
        if not os.path.exists(p):
            continue
        for r in map(json.loads, open(p)):
            if 's_pol' in r and r['s_pol'] < E.frontier_now(r['n']) - 1e-9:
                hits.append(dict(n=r['n'], s_pol=r['s_pol'], A=r.get('A'), B=r.get('B'), file=os.path.join(HERE, r['out'])
                                 if not os.path.isabs(r['out']) else r['out']))
    print(f'{len(hits)} hits: ' + ' '.join(f"{h['n']}:{h['s_pol']:.10f}" for h in hits), flush=True)
    def certify(h):
        name = f"{a.name or os.path.basename(a.out.rstrip('/'))}_n{h['n']}"
        h['rc'] = sh([sys.executable, 'candidate.py', h['file'], '--name', name, '--no-regularize', '--import'],
                     out(f"cand-{h['n']}.log"))
        cs = os.path.join(HERE, 'runs', 'candidates', name, 'summary.json')
        if not os.path.exists(cs):                   # exactsolve failed (seen: under-converged polish): re-polish, retry once
            rp = out(f"repolish-{h['n']}.txt")
            subprocess.run([os.path.join(HERE, 'target', 'release', 'fq'), 'quench', '--in', h['file'], '--out', rp,
                            '--loosen', '1.0', '--no-alm', '--flip-top', '64', '--pit', '2000'], capture_output=True)
            if os.path.exists(rp):
                h['repolished'] = rp
                h['rc'] = sh([sys.executable, 'candidate.py', rp, '--name', name, '--no-regularize', '--import'],
                             out(f"cand-{h['n']}-retry.log"))
        if os.path.exists(cs):
            c = json.load(open(cs))
            h.update(cert=all(v['ok'] for v in c.get('checkers', {}).values()), S_cert=c.get('S_cert'),
                     beaten=c.get('live', {}).get('beaten'), dir=os.path.dirname(cs))
        return h
    from concurrent.futures import ThreadPoolExecutor        # candidate.py is one process each: run them side by side
    with ThreadPoolExecutor(max(1, min(len(hits), a.procs))) as ex:
        for h in ex.map(certify, hits):
            print(f"  {h['n']}: cert {h.get('cert')} S' {h.get('S_cert')} beaten {h.get('beaten')}", flush=True)
    S['hits'] = hits
    S['wall_s'] = round(time.time() - t0)
    json.dump(S, open(out('summary.json'), 'w'), indent=1)
    man.finish(summary=dict(new=len(rows), crops=nat, pieces_new=S.get('pieces_new'), hits=len(hits),
                            certified=sum(1 for h in hits if h.get('cert'))))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run')
    r.add_argument('--since'); r.add_argument('--ids', type=int, nargs='*')
    r.add_argument('--out', required=True); r.add_argument('--name')
    r.add_argument('--procs', type=int, default=16); r.add_argument('--budget', type=float, default=30,
                                                                    help='composer seconds per n')
    r.add_argument('--probe', type=float, default=0.005); r.add_argument('--nmax', type=int, default=480)
    r.add_argument('--target', default='now', help='composer target: now | public (see compose2)')
    r.add_argument('--exclude-src', nargs='*', default=[], help='with --target public (control replays)')
    r.add_argument('--mmax', type=int, default=400); r.add_argument('--hmax', type=float, default=21.5)
    q = sp.add_parser('new'); q.add_argument('--since', required=True)
    a = ap.parse_args()
    if a.cmd == 'new':
        for i, n, s, lab in sorted(new_packings(a.since), key=lambda r: r[1]):
            print(i, n, f'{s:.10f}', lab)
    else:
        if not (a.since or a.ids):
            ap.error('--since or --ids')
        run(a)
