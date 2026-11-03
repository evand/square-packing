"""Converge near-record packings as they enter the store (10-24; Evan: "we should converge packings as they enter the
store, at least for any cheap steps.  Certified high-precision local minimum seems like it's probably too strong a
standard").  PACKER RS3 / RQ1: near-record store packings, polished ones included, often go lower under a fresh
side-free quench (53 % of 3,017; the store best improved at 96 n).

Selection (per sub-grid n): distinct (1e-5) feasible line-free sub-grid packings below best(n) + tol, not tagged
`requenched`, whose familymap code is not G / D / S (45-degree forms only converge onto forms the store holds).  Each:
relax.screen1.0 and relax.screen1.02 from the packing (pk.stages; fq SCREEN = ALM + <= 8 SLP iterations), the better
line-free result kept if below the start by > 1e-9: stored as a new packing (source `converge`, parent = the start,
status screened, tag `requenched`).  The start is tagged `requenched` either way, so nothing is converged twice.

  pk.py converge [--ns 2-1000] [--tol 0.01] [--cap 50] [--serial] [--since ISO] [--dry-run]

New rows reach it as pool units (pk/react.py: `--serial --ids`, a few packings per unit); the backlog (packings that entered
before this module existed) is a separate, explicit run.
"""
from __future__ import annotations
import contextlib, json, math, os, random, tempfile
import numpy as np

TAG = 'requenched'
SCREENS = ('screen1.0', 'screen1.02')
SKIP_CODES = ('G', 'D', 'S')
FEAS = "status NOT IN ('infeasible') AND (pen IS NULL OR pen<=1e-9) AND (grid IS NULL OR grid=0)"


def select(st, lo=2, hi=1000, tol=0.01, cap=50, since=None, ids=None):
    """[(id, n, s)] to converge, cheapest checks first (tag, side) and the classifier only on survivors"""
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from familymap import classify
    db = st.db
    out = []
    if ids is not None:                                           # 10-24: just these packings (a goal job's own)
        ids = list(ids)
        if not ids:
            return out
        q = ','.join('?' * len(ids))
        ns = sorted({n for (n,) in db.execute(f"SELECT DISTINCT n FROM packing WHERE id IN ({q})", ids)})
        lo, hi = max(lo, ns[0]), min(hi, ns[-1])
    for n in range(max(2, lo), hi + 1):
        if ids is not None and n not in ns:
            continue
        k = math.ceil(math.sqrt(n) - 1e-12)
        b = db.execute(f"SELECT min(s) FROM packing WHERE n=? AND {FEAS} AND s<?", (n, k - 1e-9)).fetchone()[0]
        if b is None:
            continue
        R, seen = [], set()
        q = f"SELECT id, s, coords, tags FROM packing WHERE n=? AND {FEAS} AND s<?" + (" AND added>=?" if since else "") + \
            (f" AND id IN ({','.join('?' * len(ids))})" if ids is not None else "")
        for i, s, c, tags in db.execute(q + " ORDER BY s", (n, min(b + tol, k - 1e-9), *([since] if since else []),
                                                             *(ids if ids is not None else []))):
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


def run(st, lo=2, hi=1000, tol=0.01, cap=50, parallel=True, dry=False, log=print, since=None, ids=None):
    """parallel: pool units of a few packings (serial inside a unit)"""
    from .packing import Packing
    from . import tasks
    J = select(st, lo, hi, tol, cap, since, ids)
    log(f'converge: {len(J)} packings over {len({n for _, n, _ in J})} n (n {lo}-{hi}, below best + {tol})')
    if dry or not J:
        return dict(selected=len(J))
    args = [(i, n, s, st.db.execute('SELECT coords FROM packing WHERE id=?', (i,)).fetchone()[0], None) for i, n, s in J]
    improved, new, cpu = 0, 0, 0.0
    ex = tasks.Executor('converge') if parallel else None
    with ex or contextlib.nullcontext():           # ~5 CPU-s per packing: 8 per unit
        for k, (pid, n, s0, s1, sq1, by, c) in enumerate(ex.map(_one, args, chunksize=8) if ex else map(_one, args)):
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
