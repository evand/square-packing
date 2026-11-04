"""Family codes cached in the store (10-25, DESIGN.md "Next: the unit loop", step 1: make decisions cheap).

familymap.classify takes 25-80 ms per packing (n 366-980) while the (n, s) index answers a donor query in ~1-10 ms,
so every decision that filters by family (goal donors, the wanted-family best, the reference) was paying a
classification per row it read.  The code is a pure function of the packing, so it is cached in packing.family
(TEXT, NULL = not yet classified) and filled lazily: readers classify only rows that are still NULL and write the
codes back (best effort: a busy store just skips the write; the next reader fills it).  No backfill is needed: 1.7 M
rows lie within 0.06 of their n's best (~19 CPU-h to classify), but decisions read a few hundred per n.

If familymap.classify changes, clear the cache: `UPDATE packing SET family = NULL` (VERSION documents the classifier
the codes came from).
"""
from __future__ import annotations
import os, sqlite3, sys
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STORE = os.path.join(HERE, 'runs', 'store.sqlite')
VERSION = 'familymap.classify 10-25'
FEAS = "status!='infeasible' AND (pen IS NULL OR pen<=1e-9) AND (grid IS NULL OR grid=0)"


def kgrid(n):
    """the grid side for n squares: ceil sqrt n"""
    import math
    return math.ceil(math.sqrt(n) - 1e-12)


def migrate(path=STORE):
    """add packing.family if missing (O(1) schema change; existing readers use rows by name)"""
    db = sqlite3.connect(path, timeout=60)
    try:
        if 'family' not in {r[1] for r in db.execute('PRAGMA table_info(packing)')}:
            db.execute('ALTER TABLE packing ADD COLUMN family TEXT'); db.commit()
            return True
        return False
    finally:
        db.close()


def has_column(db):
    return 'family' in {r[1] for r in db.execute('PRAGMA table_info(packing)')}


def classify(n, s, sq):
    sys.path.insert(0, HERE)
    from familymap import classify as C
    return C(n, s, np.asarray(sq, float))['code']


def fill(rows, n, path=STORE, write=True):
    """rows: [(id, s, coords blob, family or None)] at n -> {id: code}; classifies the NULL ones and writes them back"""
    out, new = {}, []
    for i, s, c, f in rows:
        if f is None:
            f = classify(n, s, np.frombuffer(c, '<f8').reshape(-1, 3)); new.append((f, i))
        out[i] = f
    if new and write:
        try:
            db = sqlite3.connect(path, timeout=2)
            try:
                db.executemany('UPDATE packing SET family=? WHERE id=? AND family IS NULL', new); db.commit()
            finally:
                db.close()
        except sqlite3.OperationalError:                        # busy / read-only: the next reader fills them
            pass
    return out


SHAPE_MOVE = 0.1          # differs(): a centre this far from every centre of the other packing has moved
SHAPE_K = 3               # differs(): more than this many moved squares (best box symmetry) = a different arrangement
SHAPE_SCAN = 64           # near(shape > 0): rows read per pool slot at most


def differs(s, X, Y, move=SHAPE_MOVE, k=SHAPE_K):
    """are X and Y (same n, box side ~s) different arrangements, not one packing up to relax noise and box symmetry?
    Per symmetry of Y, the squares of X whose centre is more than `move` from every centre of Y; different if more than
    k for all 8 symmetries (a rattler or two may wander).  Nearest-centre, not a matching: cheap (a KD tree when scipy
    is there), and enough to tell a noise copy from a rearrangement (11-04: composer rebuilds and lineshift variants
    at one side; relaxed copies of one packing at slightly different sides)."""
    X = np.asarray(X, float)[:, :2]; P0 = np.asarray(Y, float)[:, :2]
    try:
        from scipy.spatial import cKDTree
        T = cKDTree(P0)
        def far(Z):
            return int((T.query(Z, distance_upper_bound=move)[0] > move).sum())
    except ImportError:                                         # brute force, chunked
        def far(Z):
            m = 0
            for a in range(0, len(Z), 256):
                d = np.hypot(Z[a:a + 256, None, 0] - P0[None, :, 0], Z[a:a + 256, None, 1] - P0[None, :, 1]).min(1)
                m += int((d > move).sum())
            return m
    for sw in (False, True):                                    # move X instead of Y: the tree is built once
        for fx in (False, True):
            for fy in (False, True):
                Z = X[:, ::-1].copy() if sw else X.copy()
                if fx: Z[:, 0] = s - Z[:, 0]
                if fy: Z[:, 1] = s - Z[:, 1]
                if far(Z[:64]) <= k and far(Z) <= k:            # a sample first: most wrong symmetries fail there
                    return False
    return True


def near(db, n, s_max, exclude=(), limit=400, path=STORE, dedupe=6, write=True, window=None, before=None, shape=0):
    """feasible packings at n with s < s_max whose family is not in `exclude`, best first, one per side rounded to
    `dedupe` decimals, at most `limit`, and (window) no more than `window` above the first one kept -> [(id, s, sq,
    family)].  before: only packings added before this time (a goal's start).  Classifies lazily, in batches, only as far as it
    reads.  db: any connection to the store (read-only is fine; writes go through their own connection).
    shape (11-04, Evan: equal-side alternatives may behave differently under later moves): dedupe by arrangement
    instead of by side: a packing is skipped if differs() calls it the same arrangement as one already kept (any side),
    and at most `shape` are kept per rounded side; reads at most SHAPE_SCAN * limit rows.  0 = one per rounded side,
    as before (which also keeps relaxed copies of one packing at slightly different sides: 6 of 8 at 745, 11-04)."""
    out, seen, ex = [], set(), set(exclude)
    kept, per = [], {}                                          # shape > 0: coords kept; rounded side -> count

    def first(r):
        key = round(r[1], dedupe)
        if not shape:
            if key in seen:
                return False
            seen.add(key)
            return True
        if per.get(key, 0) >= shape:
            return False
        X = np.frombuffer(r[2], '<f8').reshape(-1, 3)
        if any(differs(r[1], X, Y) is False for Y in kept):
            return False
        kept.append(X); per[key] = per.get(key, 0) + 1
        return True
    q = (f"SELECT id, s, coords, family FROM packing WHERE n=? AND {FEAS} AND s<?"
         + (" AND added<?" if before else "") + " ORDER BY s LIMIT ? OFFSET ?")
    off, B = 0, 64
    while len(out) < limit:
        R = db.execute(q, (n, s_max, *([before] if before else []), B, off)).fetchall()
        if not R:
            break
        off += len(R)
        if shape and window is not None and out:                # compare shapes only inside the window
            if R[0][1] > out[0][1] + window:
                break
            R = [r for r in R if r[1] <= out[0][1] + window]
        R = [r for r in R if first(r)]
        last = bool(shape) and off >= SHAPE_SCAN * limit
        if window is not None and out:                          # classify only what can still be kept
            R = [r for r in R if r[1] <= out[0][1] + window] or R[:1]
        F = fill(R, n, path, write)
        for i, s, c, _ in R:
            if window is not None and out and s > out[0][1] + window:
                return out
            if F[i] not in ex:
                out.append((i, s, np.frombuffer(c, '<f8').reshape(-1, 3).copy(), F[i]))
                if len(out) >= limit:
                    break
        if last:
            break
        B = min(4096, 2 * B)
    return out


def best(db, n, s_max, exclude=(), before=None, path=STORE, write=True):
    """the lowest feasible side at n outside `exclude` (added before `before`) -> (s, id) or (None, None); the cached
    form of the old goal loop's best_match"""
    r = near(db, n, s_max, exclude, limit=1, path=path, write=write, before=before)
    return (r[0][1], r[0][0]) if r else (None, None)
