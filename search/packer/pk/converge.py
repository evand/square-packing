"""Converge near-record packings as they enter the store (10-24; Evan: "we should converge packings as they enter the
store, at least for any cheap steps.  Certified high-precision local minimum seems like it's probably too strong a
standard").  PACKER RS3 / RQ1: near-record store packings, polished ones included, often go lower under a fresh
side-free quench (53 % of 3,017; the store best improved at 96 n).

Selection (per sub-grid n): distinct (1e-5) feasible line-free sub-grid packings below best(n) + tol, not tagged
`requenched`, whose familymap code is not G / D / S (45-degree forms only converge onto forms the store holds).  Each:
relax.screen1.0 and relax.screen1.02 from the packing (pk.stages; fq SCREEN = ALM + <= 8 SLP iterations), the better
line-free result kept if below the start by > 1e-9: stored as a new packing (source `converge`, parent = the start,
status screened, tag `requenched`).  The start is tagged `requenched` either way, so nothing is converged twice.

  pk.py converge [--ns 2-1000] [--tol 0.01] [--cap 50] [--procs 8] [--since ISO] [--dry-run]

Run by `pk/jobs.py` after every job's sink with --since = the job's start (one at a time, flock on runs/sink.lock);
the backlog (packings that entered before this module existed) is a separate, explicit run.
"""
from __future__ import annotations
import json, math, os, random, tempfile
from concurrent.futures import ProcessPoolExecutor
import numpy as np

TAG = 'requenched'
SCREENS = ('screen1.0', 'screen1.02')
SKIP_CODES = ('G', 'D', 'S')
FEAS = "status NOT IN ('infeasible') AND (pen IS NULL OR pen<=1e-9) AND (grid IS NULL OR grid=0)"


def select(st, lo=2, hi=1000, tol=0.01, cap=50, since=None):
    """[(id, n, s)] to converge, cheapest checks first (tag, side) and the classifier only on survivors"""
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from familymap import classify
    db = st.db
    out = []
    for n in range(max(2, lo), hi + 1):
        k = math.ceil(math.sqrt(n) - 1e-12)
        b = db.execute(f"SELECT min(s) FROM packing WHERE n=? AND {FEAS} AND s<?", (n, k - 1e-9)).fetchone()[0]
        if b is None:
            continue
        R, seen = [], set()
        q = f"SELECT id, s, coords, tags FROM packing WHERE n=? AND {FEAS} AND s<?" + (" AND added>=?" if since else "")
        for i, s, c, tags in db.execute(q + " ORDER BY s", (n, min(b + tol, k - 1e-9), *([since] if since else []))):
            if round(s, 5) in seen:
                continue
            seen.add(round(s, 5))
            if TAG in json.loads(tags):
                continue
            if classify(n, s, np.frombuffer(c, '<f8').reshape(-1, 3))['code'] in SKIP_CODES:
                continue
            R.append((i, n, s))
        out += R[:cap]
    return out


def _one(args):
    pid, n, s0, coords, store_path = args
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import mcmin
    from layout import full_lines
    from pk import stages as S
    from pk.packing import Packing
    sq0 = np.frombuffer(coords, '<f8').reshape(-1, 3).copy()
    k = math.ceil(math.sqrt(n) - 1e-12)
    st = S.State(Packing(s0, sq0))
    best, cpu = None, 0.0
    for sc in SCREENS:
        st1, row = S.run_stage(f'relax.{sc}', {}, st, random.Random(1), dict(n=n)); cpu += row['cpu']
        if row['info'].get('failed'):
            continue
        with tempfile.TemporaryDirectory(prefix='cv_') as td:
            p = os.path.join(td, 'e.txt'); mcmin.write_deg(p, float(st1.p.s), [tuple(map(float, q)) for q in st1.p.sq])
            g = bool(full_lines(p, k)) or st1.p.s >= k - 1e-9
        if not g and (best is None or st1.p.s < best[0]):
            best = (float(st1.p.s), np.asarray(st1.p.sq, float), sc)
    if best and best[0] < s0 - 1e-9:
        return pid, n, s0, best[0], best[1], best[2], cpu
    return pid, n, s0, None, None, None, cpu


def run(st, lo=2, hi=1000, tol=0.01, cap=50, procs=8, dry=False, log=print, since=None):
    from .packing import Packing
    J = select(st, lo, hi, tol, cap, since)
    log(f'converge: {len(J)} packings over {len({n for _, n, _ in J})} n (n {lo}-{hi}, below best + {tol})')
    if dry or not J:
        return dict(selected=len(J))
    args = [(i, n, s, st.db.execute('SELECT coords FROM packing WHERE id=?', (i,)).fetchone()[0], None) for i, n, s in J]
    improved, new, cpu = 0, 0, 0.0
    with ProcessPoolExecutor(procs) as ex:
        for k, (pid, n, s0, s1, sq1, by, c) in enumerate(ex.map(_one, args, chunksize=2)):
            cpu += c
            if s1 is not None:
                improved += 1
                nid, was_new = st.add(Packing(s1, sq1), 'converge', 'ours', ref=f'{by} from {pid}', status='screened',
                                      parent=pid, tags=(TAG,), commit=False)
                new += bool(was_new)
            st.tag(pid, TAG)
            if k % 200 == 199:
                st.db.commit(); log(f'  {k + 1} / {len(args)}')
    st.db.commit()
    r = dict(selected=len(J), improved=improved, new=new, cpu_h=round(cpu / 3600, 3))
    log(f'converge: {r}')
    return r
