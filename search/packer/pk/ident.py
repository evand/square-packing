"""Strict identity ("is this the same minimum?") and the conservative store dedupe built on it (10-20; Evan: "the store
tries to do strict dedupe, somewhat lazily and imprecisely"; runs that pick starts do their own fast cluster analysis).

Evidence (PACKER.md "Identity at 274 / 110", "110 against exact minima", "Counting vs the 110 census"):
  * f64 side alone is not an identity: copies of one minimum scatter between 1e-9 and 1e-7, two exact minima at 110 have
    f64 sides equal to ~1e-9, and many 'polished' rows are not minima at all (exactsolve: not jammed).
  * rigid core = squares with no first-order motion at fixed side (fixlp.mobility < 1e-7, 0.1 s at 274).  Two rows are
    compared over the squares locked in BOTH (core membership drifts between copies), best of 8 symmetries + Hungarian.
    Max centre deviation: true copies 3e-11 .. 5e-5 (one 8e-4 at 274, one 1e-3 at 110); distinct minima >= 2e-3 (274),
    >= 8e-3 (110).

RULE (Evan 10-20: "conservative enough for now; revisit before applying it again ... shouldn't apply it automatically
without discussing it again in a future session"):  two rows of one n are the same if
  (a) both carry an exact side and their exact minimum keys (identity.min_key, +-1 quantum) are equal, or
  (b) both are converged candidates (status polished / certified / weak / unresolved; not not-min, screened, raw;
      no grid line), |s1 - s2| <= SIDE_TOL (1e-11) and the rigid-core deviation <= CORE_TOL (1e-4, matched core angles
      within 0.01 deg), and neither side is within GUARD (1e-9) of two distinct known exact minima.
  10-20 tightening (Evan: "conservatism ... those pairs should remain distinct"): two certified minima at 110 differ only
  by which side a corner rests on at near-coincident corners (cores 1e-7 / 6.7e-6 apart, exact sides 7e-10 / 3e-9
  apart; analysis/contactswitch.py), so positions cannot separate them; the side can, given f64 sides of converged rows
  sit within 1e-11 of the exact side (96 %; 99 % within 1e-9) and known minima closer than that are guarded.
Rows are compared with group representatives only (no chaining).  Representative: protected first (policy.protected_ids
+ `extra_protect`), then status rank, then lowest side, then lowest id.  A protected row is never deleted; two
protected rows in one group are linked but both kept.  `check_exact`: a merge whose two rows resolve (s_exact or the
census copy table) to different exact minima is refused and counted.

apply(): per copy, in one transaction: sightings, expansion rows and child parent links move to the representative, tags
union, status = the better; the row goes to table `dedup` (old id, hash, side, rep, method, evidence) and is deleted.
Store.add / Store.get and pk.archive.ensure follow `dedup` (an old id or hash resolves to the representative).
Nothing here runs by itself: `pk.py dedupe --n N` is a dry run; `--apply` needs `--i-discussed-this`.
"""
from __future__ import annotations
import collections, datetime, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE)
SIDE_TOL = 1e-11       # 10-20 (Evan: conservatism; distinct exact minima at 110 come within 4.6e-13 .. 1e-7 of each other)
GUARD = 1e-9           # refuse a merge if a row's side is within GUARD of two or more distinct known exact minima
CORE_TOL = 1e-4
CONVERGED = ('polished', 'certified', 'weak', 'unresolved')
SCHEMA = """CREATE TABLE IF NOT EXISTS dedup (
  old INTEGER PRIMARY KEY, hash TEXT NOT NULL, n INTEGER NOT NULL, s REAL NOT NULL, rep INTEGER NOT NULL,
  method TEXT NOT NULL, evidence TEXT, added TEXT NOT NULL);
CREATE INDEX IF NOT EXISTS dedup_hash ON dedup(hash);
CREATE INDEX IF NOT EXISTS dedup_rep ON dedup(rep);
CREATE INDEX IF NOT EXISTS packing_parent ON packing(parent);"""


def ensure_schema(db):
    db.executescript(SCHEMA)


def _paths():
    for p in (PK, os.path.join(PK, 'rational')):
        if p not in sys.path:
            sys.path.insert(0, p)


def core_mask(args):
    """worker: (pid, s, coords bytes) -> (pid, locked mask)"""
    _paths()
    import fixlp
    pid, s, blob = args
    sq = np.frombuffer(blob, '<f8').reshape(-1, 3).copy()
    try:
        return pid, fixlp.mobility(sq, s) < 1e-7
    except Exception:
        return pid, None


def core_dev(a, b):
    """a, b = (s, sq, mask).  Max centre deviation over squares locked in both, best of 8 symmetries (Hungarian over
    all squares); inf if no symmetry overlaps at least half of the smaller core."""
    _paths()
    from scipy.optimize import linear_sum_assignment
    from core.coremap import D4
    sa, qa, ca = a; sb, qb, cb = b
    best = np.inf
    for Bf in D4(qb, sb):
        C = np.sqrt(((qa[:, None, :2] - Bf[None, :, :2]) ** 2).sum(-1))
        r, c = linear_sum_assignment(C)
        core = ca[r] & cb[c]
        if core.sum() < 0.5 * min(ca.sum(), cb.sum()) or not core.any():
            continue
        best = min(best, float(C[r, c][core].max()))
    return best


def _syms(s, sq):
    """8 box images of sq with a KD tree on each (for core_dev_nn)."""
    _paths()
    from scipy.spatial import cKDTree
    from core.coremap import D4
    return [(Bf, cKDTree(Bf[:, :2])) for Bf in D4(sq, s)]


def core_dev_nn(qa, ca, symsB, cb, tol):
    """core_dev for small deviations, by nearest neighbours: squares are >= ~1 apart, so a deviation <= tol (<< 1) makes
    the nearest neighbour the matched square, and nearest-neighbour distance <= matched distance, so 'dev <= tol' is
    decided exactly; values above tol are lower bounds."""
    best = np.inf
    for Bf, tree in symsB:
        d, j = tree.query(qa[:, :2])
        core = ca & cb[j]
        if core.sum() < 0.5 * min(ca.sum(), cb.sum()) or not core.any():
            continue
        if d[core].max() <= tol:
            ang = np.abs(qa[core, 2] - Bf[j[core], 2]) % 90
            if np.minimum(ang, 90 - ang).max() > 1e-2:     # same centres, different angles: not the same square
                continue
        best = min(best, float(d[core].max()))
        if best <= tol:
            break
    return best


def _group_chain(args):
    """worker: one side chain -> (groups {rep: [(old, method, ev)]}, refused Counter, comparisons)"""
    chain, side_tol, core_tol, known = args
    def guarded(x):                                   # two or more distinct known minima within GUARD of side x
        if len(known) < 2:
            return False
        lo, hi = np.searchsorted(known, [x - GUARD, x + GUARD])
        return hi - lo >= 2
    groups, refused, compared, refused_pairs = collections.defaultdict(list), collections.Counter(), 0, []
    rows = {p: (s, np.frombuffer(c, '<f8').reshape(-1, 3), m, k, e) for p, s, c, m, k, e in chain}
    syms = {}
    reps, rep_of = [], {}
    for p, s, _, m, k, e in chain:
        qa = rows[p][1]; hit = None
        for s2, q in reversed(reps):
            if s - s2 > side_tol:
                break
            if q not in syms:
                syms[q] = _syms(rows[q][0], rows[q][1])
            compared += 1
            d = core_dev_nn(qa, m, syms[q], rows[q][2], core_tol)
            if d <= core_tol:
                hit = (q, d); break
        if hit is None:
            reps.append((s, p)); continue
        q, d = hit
        if guarded(s) or guarded(rows[q][0]):
            refused['near two known minima (guard band)'] += 1; reps.append((s, p))
            refused_pairs.append(dict(a=p, b=q, dside=abs(s - rows[q][0]), core=d, guard=True))
            continue
        if e and rows[q][4] and set(e) != set(rows[q][4]) and len(e) == 1 and len(rows[q][4]) == 1:
            refused['different exact minima'] += 1; reps.append((s, p))
            refused_pairs.append(dict(a=p, b=q, dside=abs(s - rows[q][0]), core=d, ea=list(e), eb=list(rows[q][4])))
            continue
        rep, old = sorted((p, q), key=lambda x: rows[x][3])
        if old in groups:
            for x in groups.pop(old):
                groups[rep].append(x)
        groups[rep].append((old, 'side+core', dict(dside=abs(rows[old][0] - rows[rep][0]), core=d)))
        if rep == p:
            reps = [(s_, x) for s_, x in reps if x != old] + [(s, p)]
    return dict(groups), refused, compared, refused_pairs


def _exact_maps(db, n):
    from .identity import side_key, min_key
    C = collections.defaultdict(set)
    for sk, mk in db.execute('SELECT side_key, min_key FROM copy WHERE n=?', (n,)):
        C[sk].add(mk)
    def resolve(s, s_exact):
        if s_exact:
            return {min_key(s_exact)}
        return C.get(side_key(s), set())
    return resolve


def census_protect(n):
    """Rows the published census data uses (110: runs/landscape110/map/minima.json pids)."""
    P = set()
    if n == 110:
        f = os.path.join(PK, 'runs', 'landscape110', 'map', 'minima.json')
        if os.path.exists(f):
            P.update(m['pid'] for m in json.load(open(f)) if m.get('pid'))
    return P


def plan(st, n, procs=16, side_tol=SIDE_TOL, core_tol=CORE_TOL, check_exact=True, log=print):
    """-> dict(groups=[(rep, [(old, method, evidence)])], stats)"""
    from concurrent.futures import ProcessPoolExecutor
    from . import policy
    from .identity import min_key
    from .store import RANK
    db = st.db
    rows = db.execute("SELECT id, s, status, grid, s_exact, coords FROM packing WHERE n=? AND status != 'infeasible'",
                      (n,)).fetchall()
    prot = policy.protected_ids(st, policy.load()) | census_protect(n)
    prot |= {r[0] for r in db.execute('SELECT pid FROM minimum WHERE n=? AND pid IS NOT NULL', (n,))}
    resolve = _exact_maps(db, n)
    R = {r[0]: r for r in rows}
    cand = [r for r in rows if r[2] in CONVERGED and not (r[3] or 0)]
    # (a) exact keys
    ek = collections.defaultdict(list)
    for r in rows:
        if r[4] and r[2] != 'not-min':
            ek[min_key(r[4])].append(r[0])
    # (b) side windows: only rows with another candidate within side_tol need a core
    cand.sort(key=lambda r: r[1])
    S = np.array([r[1] for r in cand])
    need = set()
    for i in range(len(cand)):
        if (i + 1 < len(cand) and S[i + 1] - S[i] <= side_tol) or (i > 0 and S[i] - S[i - 1] <= side_tol):
            need.add(cand[i][0])
    log(f'n {n}: {len(rows):,} rows, {len(cand):,} converged candidates, {len(need):,} with a neighbour within '
        f'{side_tol:g}; protected {len(prot & set(R)):,}; cores for {len(need):,} rows on {procs} procs')
    cache = os.path.join(PK, 'runs', 'dedup', f'cores_{n}.npz')      # rigid cores keyed by row id (rerun = seconds)
    masks = {}
    if os.path.exists(cache):
        z = np.load(cache)
        masks = {int(p): m.astype(bool) for p, m in zip(z['ids'], z['masks'])}
    jobs = [(p, R[p][1], R[p][5]) for p in need if p not in masks]
    if procs <= 1:
        masks.update(map(core_mask, jobs))
    else:
        with ProcessPoolExecutor(procs) as ex:
            for pid, m in ex.map(core_mask, jobs, chunksize=32):
                masks[pid] = m
    if jobs and procs > 1:
        ok = [(p, m) for p, m in masks.items() if m is not None and len(m) == n]
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        np.savez_compressed(cache, ids=np.array([p for p, _ in ok]), masks=np.array([m for _, m in ok]))
    key = lambda p: (p not in prot, -RANK.get(R[p][2], 0), R[p][1], p)
    rep_of, groups, refused, compared = {}, collections.defaultdict(list), collections.Counter(), 0
    # exact-key groups first
    for k, ids in ek.items():
        ids = sorted(ids, key=key)
        for p in ids[1:]:
            rep_of[p] = ids[0]; groups[ids[0]].append((p, 'exact-key', dict(key=k)))
    # side + core: independent chains (consecutive candidate sides within side_tol), grouped in parallel
    chains, cur, last = [], [], None
    for r in cand:
        p = r[0]
        if p in rep_of or p not in need or masks.get(p) is None:
            continue
        if last is not None and r[1] - last > side_tol and cur:
            chains.append(cur); cur = []
        ex_ = resolve(r[1], r[4]) if check_exact else set()
        cur.append((p, r[1], r[5], masks[p], key(p), tuple(ex_))); last = r[1]
    if cur:
        chains.append(cur)
    chains = [c for c in chains if len(c) > 1]
    log(f'  {len(chains):,} side chains ({sum(map(len, chains)):,} rows; largest {max(map(len, chains), default=0):,})')
    known = np.array(sorted(float(S) for (S,) in db.execute('SELECT S FROM minimum WHERE n=?', (n,))))
    args = [(c, side_tol, core_tol, known) for c in chains]
    done, t0, refused_list = 0, __import__('time').time(), []
    order = sorted(range(len(args)), key=lambda i: -len(args[i][0]))       # big chains first
    it = map(_group_chain, [args[i] for i in order]) if procs <= 1 else None
    ex = None
    if it is None:
        ex = ProcessPoolExecutor(procs); it = ex.map(_group_chain, [args[i] for i in order], chunksize=1)
    try:
        for G, ref, cmp_, rp in it:
            compared += cmp_; refused.update(ref); done += 1; refused_list.extend(rp)
            for rep, L in G.items():
                for old, m, ev in L:
                    rep_of[old] = rep; groups[rep].append((old, m, ev))
                    if old in groups and old != rep:                     # an exact-key representative merged by core
                        for x in groups.pop(old):
                            rep_of[x[0]] = rep; groups[rep].append(x)
            if done % 500 == 0 or done == len(args):
                log(f'  chains {done:,} / {len(args):,}, comparisons {compared:,}, {__import__("time").time() - t0:.0f} s')
    finally:
        if ex is not None:
            ex.shutdown()
    # drop protected rows from deletion (kept, linked only)
    cens = census_protect(n)
    keep_linked = [(rep, old) for rep, L in groups.items() for old, *_ in L if old in prot]
    dele = sum(1 for rep, L in groups.items() for old, *_ in L if old not in prot)
    stats = dict(n=n, rows=len(rows), candidates=len(cand), compared=compared, groups=len(groups), delete=dele,
                 protected_kept=len(keep_linked), refused=dict(refused), cores_failed=sum(m is None for m in masks.values()),
                 by_method=dict(collections.Counter(m for L in groups.values() for _, m, _ in L)),
                 census_hit=sum(1 for L in groups.values() for old, *_ in L if old in cens))
    log(f'  groups {len(groups):,}; rows to delete {dele:,} (by method {stats["by_method"]}); protected rows linked '
        f'but kept {len(keep_linked):,}; refused {dict(refused)}; core comparisons {compared:,}')
    return dict(groups={k: v for k, v in groups.items()}, protected=prot, stats=stats, refused=refused_list)


def apply(st, P, log=print):
    db = st.db
    ensure_schema(db)
    now = datetime.datetime.now().isoformat(timespec='seconds')
    from .store import RANK
    nd = 0
    with db:                                         # one transaction
        for rep, L in P['groups'].items():
            rr = db.execute('SELECT status, tags FROM packing WHERE id=?', (rep,)).fetchone()
            st_rep, tags = rr[0], set(json.loads(rr[1]))
            for old, method, ev in L:
                if old in P['protected']:
                    continue
                o = db.execute('SELECT hash, n, s, status, tags FROM packing WHERE id=?', (old,)).fetchone()
                if o is None:
                    continue
                db.execute('UPDATE OR IGNORE sighting SET packing=? WHERE packing=?', (rep, old))
                db.execute('DELETE FROM sighting WHERE packing=?', (old,))
                db.execute('UPDATE expansion SET packing=? WHERE packing=?', (rep, old))
                db.execute('UPDATE packing SET parent=? WHERE parent=?', (rep, old))
                db.execute('UPDATE dedup SET rep=? WHERE rep=?', (rep, old))
                tags |= set(json.loads(o[4]))
                if RANK.get(o[3], 0) > RANK.get(st_rep, 0):
                    st_rep = o[3]
                db.execute('INSERT OR REPLACE INTO dedup VALUES (?,?,?,?,?,?,?,?)',
                           (old, o[0], o[1], o[2], rep, method, json.dumps(ev), now))
                db.execute('DELETE FROM packing WHERE id=?', (old,))
                nd += 1
            db.execute('UPDATE packing SET status=?, tags=? WHERE id=?', (st_rep, json.dumps(sorted(tags)), rep))
    log(f'  applied: {nd:,} rows deleted into dedup')
    return nd


def resolve_old(db, pid):
    """id of the row that now stands for `pid` (follows dedup links), or None."""
    seen = set()
    while pid is not None and pid not in seen:
        seen.add(pid)
        if db.execute('SELECT 1 FROM packing WHERE id=?', (pid,)).fetchone():
            return pid
        try:
            r = db.execute('SELECT rep FROM dedup WHERE old=?', (pid,)).fetchone()
        except Exception:
            return None
        pid = r[0] if r else None
    return None
