#!/usr/bin/env python3
"""TT2 n = 110: exactsolve every new sub-11 minimum of both arms (analysis/tt2_census.py 'new', 18 rows), compare with
the 939 known exact minima (store table `minimum`, identity.min_key), and place each in the as-of region map (leader
balls, pk/cluster, r = 0.42 / sqrt n) with its nearest known minimum by matching distance.

  python3 analysis/tt2_cert110.py [--procs 6]     -> runs/tt2/cert110/ + runs/tt2/cert110.json
  python3 analysis/tt2_cert110.py --desc           the same on runs/tt2/cert110/desc/ (exact/descend.py outputs) -> cert110_desc.json
"""
from __future__ import annotations
import argparse, contextlib, io, json, os, sqlite3, sys
from multiprocessing import Pool
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PK = os.path.dirname(HERE)
sys.path.insert(0, PK)
sys.path.insert(0, os.path.join(PK, '..', 'exact'))
from pk import cluster as CL, identity as ID, parents as PA       # noqa: E402

N = 110
OUT = os.path.join(PK, 'runs', 'tt2', 'cert110')
DB = os.path.join(PK, 'runs', 'store.sqlite')


def solve(job):
    arm, pid, path = job
    import exactsolve as E, census_exact as CE
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            rep = E.run(path, outdir=os.path.join(OUT, 'exact'), quiet=True)[0]
        v = CE.classify(rep)
        v['S_exact'] = rep.get('S_exact')
        v['vv_milp_dS'] = rep.get('vv_milp_dS')
    except Exception as e:                            # noqa: BLE001
        v = dict(cls='unresolved', status=f'error: {e}'[:120])
    return arm, pid, v


def main(procs, desc=False):
    from pk.store import Store
    os.makedirs(os.path.join(OUT, 'exact'), exist_ok=True)
    D = json.load(open(os.path.join(PK, 'runs', 'tt2', 'census.json')))
    r = [x for x in D['rows'] if x['n'] == N][0]
    fT = r['frontier_T']
    st = Store(); jobs = []
    for a in ('S', 'C'):
        for pid, s, gap in r[a]['sample_pool']:
            path = os.path.join(OUT, f'{a}_{pid}.txt')
            if desc:                                  # --desc: the exact/descend.py outputs of the not-jammed ones
                path = os.path.join(OUT, 'desc', f'{a}_{pid}.txt')
                if not os.path.exists(path):
                    continue
            elif not os.path.exists(path):
                st.get(pid).write(path)
            jobs.append((a, pid, path))
    st.close()
    db = sqlite3.connect(f'file:{DB}?mode=ro', uri=True)
    known = {k: (S, cls) for k, S, cls in db.execute('SELECT key, S, cls FROM minimum WHERE n=?', (N,))}
    # region map: as-of reference (same construction as tt2_census), leaders deepest first
    asof = PA.asof_local('2026-10-10T06:00')
    ref = [(i, s, np.frombuffer(c, '<f8').reshape(-1, 3)) for i, s, c in db.execute(
        "SELECT id, s, coords FROM packing WHERE n=? AND s<11 AND status!='infeasible' AND added<? AND (grid IS NULL OR grid=0)",
        (N, asof))]
    import tt2_census as TC
    ref = TC.dedupe(ref)
    L = CL.Leaders(N)
    for pid, s, sq in ref:
        L.assign(pid, s, sq)
    rank = {pid: i for i, pid in enumerate(sorted(L.ids, key=lambda p: L.s[L.ids.index(p)]))}
    coords = {pid: (s, np.frombuffer(c, '<f8').reshape(-1, 3)) for pid, s, c in db.execute(
        f"SELECT id, s, coords FROM packing WHERE id IN ({','.join(str(j[1]) for j in jobs)})")}
    out = []
    with Pool(procs) as P:
        for a, pid, v in P.imap_unordered(solve, jobs):
            s, sq = coords[pid]
            if desc:                                  # place the descended packing, not the stored one
                from pk.packing import read
                q = read(os.path.join(OUT, 'desc', f'{a}_{pid}.txt')); s, sq = float(q.s), np.asarray(q.sq, float)
            j = TC.lookup(L, s, sq)
            # nearest known minimum (reference side keys) by matching distance, among the 40 closest in side
            near = sorted(ref, key=lambda r: abs(r[1] - s))[:40]
            dm, nid, ns = min((CL.md(q, rs, sq, s), rid, rs) for rid, rs, q in near)
            hit = None
            if v.get('key'):
                for k in ID.same_min_keys(v['S_exact']):
                    if k in known:
                        hit = dict(key=k, S=known[k][0][:30], cls=known[k][1])
            row = dict(arm=a, pid=pid, s=s, gap=s - fT, verdict=v['cls'], S_exact=(v.get('S_exact') or '')[:40],
                       status=v.get('status', '')[:100], known_minimum=hit,
                       region=(f'c110:{L.ids[j]} (leader s {L.s[j]:.10f}, depth rank {rank[L.ids[j]]} of {len(L.ids)})'
                               if j >= 0 else 'outside every known region'),
                       nearest_known=dict(id=nid, s=ns, md=dm))
            out.append(row)
            print(f"{a} {pid} s {s:.13f} gap {s - fT:.2e}  {v['cls']:10s} "
                  f"{'KNOWN ' + hit['key'][:22] if hit else ('NEW exact minimum' if v['cls'] in ('certified', 'weak') else '-'):34s} "
                  f"{row['region']}  nearest md {dm:.4f}", flush=True)
    out.sort(key=lambda r: r['s'])
    json.dump(dict(n=N, frontier_T=fT, known_minima=len(known), ref_leaders=len(L.ids), rows=out),
              open(os.path.join(PK, 'runs', 'tt2', 'cert110_desc.json' if desc else 'cert110.json'), 'w'), indent=1)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--procs', type=int, default=6); ap.add_argument('--desc', action='store_true')
    a = ap.parse_args()
    main(a.procs, a.desc)
