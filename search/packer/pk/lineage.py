"""Lineage of a store packing, for credit (10-28; after s(172) / s(204) / s(260) for jlevy#845 were traced by hand and
a TODO note had the wrong start for 204).

Parent links live in three places; parent_of() reads them in this order:
  packing.parent    set at intake (explorer archives, <dir>/parents.json, sources.import_paths)
  expansion.child   the move log (pk/store expansion ledger: packing --mechanism--> child)
  file name         a sighting ref whose basename holds `_p<id>_` (drivers that name outputs after their parent)
A walk follows parents up to a root and reports, per step: id, n, side, added, sources / finders, how the link was
found, the move mechanism when the log has it, and flags: `cross-n` (the parent is at another n: insertion, removal,
an L step, a composition donor), `run-dir-missing` (a run: source whose directory is gone, so its log can't be read).

The root is where no link is recorded.  A root sighted only from one of our runs means the lineage is incomplete
(`unrecorded`); a root with a public sighting is a natural end.  For unrecorded roots, donors() lists earlier packings
at other n with the same fractional side (to 1e-9): a block that sets the side (side = k + a fixed fraction) carried
between n by composition leaves exactly this trace.  These are hints, not links.

Credit: the start to credit is the nearest public ancestor (public = store.is_public, or a kaggle: notebook output);
every public ancestor is listed with its finder (SUBMITTING.md: name every start / ancestor found by others).

  pk.py lineage ID|best:N [--json] [--depth 200] [--no-donors]
  pk.py lineage --backfill [--apply]       copy expansion.child / file-name links into packing.parent where unset
"""
from __future__ import annotations
import math, os, re, sqlite3

from .store import is_public

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEAS = "status!='infeasible' AND (pen IS NULL OR pen<=1e-9)"
PREF = re.compile(r'(?:^|_)p(\d+)_')
SHRINK = ('polish', 'converge', 'descend', 'slp2', 'requench', 'exact')     # mechanism words for shrink-only steps


OURS = re.compile(r'^ours$|\(evand\)|^evand$|^Evan Daniel', re.I)
COPY_TOL = 1e-9                  # same n, side equal to this: the same packing re-imported (e.g. an explorer start file)


def public_source(src):
    return is_public(src) or src.startswith('kaggle:')


def others(s):
    """a public sighting found by someone else (our own posted claims are public but not a credit start)"""
    return public_source(s['source']) and not OURS.search(s['finder'] or '') and not s['source'].startswith('pool:ours')


def _copy_of(db, pid, n, s, added):
    """an earlier packing at the same n with the same side (to COPY_TOL), preferring one with a public sighting by others
    or a recorded parent: the walk continues through it ('same-side')"""
    best = None
    for (q,) in db.execute(f'SELECT id FROM packing WHERE n=? AND s BETWEEN ? AND ? AND (added < ? OR (added = ? AND id < ?)) '
                           f'AND {FEAS} ORDER BY added, id', (n, s - COPY_TOL, s + COPY_TOL, added, added, pid)):
        sg = _sightings(db, q)
        if any(others(x) for x in sg):
            return q
        if best is None and (parent_of(db, q, sg)[0] is not None):
            best = q
    return best


def _sightings(db, pid):
    return [dict(source=r[0], finder=r[1], ref=r[2], added=r[3]) for r in
            db.execute('SELECT source, finder, ref, added FROM sighting WHERE packing=? ORDER BY added', (pid,))]


def _ref_parent(db, pid, sights):
    for s in sights:
        m = PREF.search(os.path.basename(str(s['ref'] or '')))
        if m:
            par = int(m.group(1))
            if par != pid and db.execute('SELECT 1 FROM packing WHERE id=?', (par,)).fetchone():
                return par
    return None


def parent_of(db, pid, sights=None):
    """(parent id or None, how, mechanism)"""
    r = db.execute('SELECT parent FROM packing WHERE id=?', (pid,)).fetchone()
    ex = db.execute('SELECT packing, mechanism FROM expansion WHERE child=? ORDER BY added LIMIT 1', (pid,)).fetchone()
    mech = ex[1] if ex else None
    if r and r[0] is not None and r[0] != pid:
        return r[0], 'parent', mech
    if ex and ex[0] != pid:
        return ex[0], 'move-log', mech
    par = _ref_parent(db, pid, sights if sights is not None else _sightings(db, pid))
    return (par, 'file-name', mech) if par else (None, None, mech)


def _run_dir_missing(sights):
    for s in sights:
        if s['source'].startswith('run:'):
            d = s['source'][4:]
            if not os.path.isdir(d if os.path.isabs(d) else os.path.join(HERE, d)):
                return True
    return False


def walk(db, pid, depth=200):
    """steps from pid up to the root (pid first)"""
    out, cur, seen = [], pid, set()
    while cur is not None and cur not in seen and len(out) < depth:
        seen.add(cur)
        row = db.execute('SELECT id, n, s, added FROM packing WHERE id=?', (cur,)).fetchone()
        if not row:
            break
        sights = _sightings(db, cur)
        par, how, mech = parent_of(db, cur, sights)
        pub = [s for s in sights if others(s)]
        if par is None and not pub:
            par = _copy_of(db, cur, row[1], row[2], row[3])
            how = 'same-side' if par is not None else None
        st = dict(id=row[0], n=row[1], s=row[2], added=row[3], sources=sorted({s['source'] for s in sights}),
                  public=[dict(source=s['source'], finder=s['finder'], added=s['added']) for s in pub],
                  parent=par, how=how, mechanism=mech, flags=[])
        if _run_dir_missing(sights):
            st['flags'].append('run-dir-missing')
        if mech and any(w in mech for w in SHRINK):
            st['flags'].append('shrink-only')
        if out and out[-1]['n'] != st['n']:
            out[-1]['flags'].append('cross-n')
        out.append(st)
        cur = par
    if out:
        root = out[-1]
        root['flags'].append('root:public' if root['public'] else
                             ('root:cycle-or-depth' if root['parent'] is not None else 'root:unrecorded'))
    return out


def credit(steps):
    """(nearest public ancestor step or None, [every public ancestor step])"""
    pubs = [s for s in steps if s['public']]
    return (pubs[0] if pubs else None), pubs


def near_above(db, step, gap=1e-4, limit=5):
    """earlier packings at the same n just above step's side: pre-polish / pre-descent ancestors (hints)"""
    rows = db.execute(f'SELECT id, s, added FROM packing WHERE n=? AND s > ? AND s <= ? AND added <= ? AND id != ? AND {FEAS} '
                      'ORDER BY s LIMIT ?', (step['n'], step['s'] + COPY_TOL, step['s'] + gap, step['added'], step['id'], limit))
    return [dict(id=r[0], n=step['n'], s=r[1], added=r[2], sources=sorted({x['source'] for x in _sightings(db, r[0])})[:3])
            for r in rows]


def donors(db, step, tol=1e-9, nmax=1000, limit=8):
    """earlier packings at other n whose side has the same fractional part as step's (composition hints)"""
    f = step['s'] - math.floor(step['s'])
    out = []
    for n in range(2, nmax + 1):
        if n == step['n']:
            continue
        k = math.ceil(math.sqrt(n) - 1e-12)
        lo = k - 1 + f - tol                       # sub-grid sides at n lie in (k - 1, k)
        for r in db.execute(f'SELECT id, s, added FROM packing WHERE n=? AND s BETWEEN ? AND ? AND added < ? AND {FEAS} '
                            'ORDER BY added LIMIT 1', (n, lo, lo + 2 * tol, step['added'])):
            out.append(dict(id=r[0], n=n, s=r[1], added=r[2]))
    out.sort(key=lambda d: d['added'])
    for d in out[:limit]:
        d['sources'] = sorted({s['source'] for s in _sightings(db, d['id'])})[:3]
    return out[:limit]


def backfill(st, apply=False, batch=500, log=print):
    """packing.parent from expansion.child, then from sighting file names, where unset.  Returns (move-log, file-name)."""
    db = st.db
    a = db.execute('SELECT e.child, MIN(e.packing) FROM expansion e JOIN packing p ON p.id=e.child '
                   'WHERE e.child IS NOT NULL AND p.parent IS NULL AND e.packing != e.child GROUP BY e.child').fetchall()
    b = []
    for pid, ref in db.execute("SELECT s.packing, s.ref FROM sighting s JOIN packing p ON p.id=s.packing "
                               "WHERE p.parent IS NULL AND s.source LIKE 'run:%' AND s.ref LIKE '%\\_p%\\_%' ESCAPE '\\'"):
        m = PREF.search(os.path.basename(str(ref)))
        if m and int(m.group(1)) != pid:
            b.append((pid, int(m.group(1))))
    done = {c for c, _ in a}
    b = [(c, p) for c, p in dict(b).items() if c not in done and db.execute('SELECT 1 FROM packing WHERE id=?', (p,)).fetchone()]
    log(f'backfill: {len(a)} links from the move log, {len(b)} from file names' + ('' if apply else ' (dry run)'))
    if apply:
        for i, (c, p) in enumerate(a + b):
            db.execute('UPDATE packing SET parent=? WHERE id=? AND parent IS NULL', (p, c))
            if i % batch == batch - 1:
                db.commit()
        db.commit()
    return len(a), len(b)


def fmt(db, steps, with_donors=True):
    L = []
    for i, s in enumerate(steps):
        who = '; '.join(f"{p['source'][:60]} ({p['finder'] or '?'})" for p in s['public'][:2])
        src = who or ', '.join(x[:50] for x in s['sources'][:2])
        via = f"  <- {s['how']}" + (f" [{s['mechanism']}]" if s['mechanism'] else '') if s['parent'] else ''
        L.append(f"{i:3d}  {s['id']:>8}  n={s['n']:<4} {s['s']:.10f}  {s['added'][:16]}  {src}{via}"
                 + (f"  {{{', '.join(s['flags'])}}}" if s['flags'] else ''))
    near, pubs = credit(steps)
    L.append('')
    if near:
        p = near['public'][0]
        L.append(f"credit start (nearest public ancestor): {near['id']} n={near['n']} {near['s']:.10f}  "
                 f"{p['source']} ({p['finder'] or '?'}), generation {steps.index(near)}")
        others = [s for s in pubs if s is not near]
        if others:
            L.append('other public ancestors: ' + '; '.join(
                f"{s['id']} n={s['n']} ({s['public'][0]['finder'] or s['public'][0]['source'][:40]})" for s in others[:6]))
    else:
        L.append('credit start: none found (no public ancestor on the recorded lineage)')
    gaps = [s for s in steps if 'run-dir-missing' in s['flags']]
    if gaps:
        L.append(f"run dir missing at {', '.join(str(s['id']) for s in gaps)}: those steps can't be checked from logs")
    root = steps[-1] if steps else None
    if root and 'root:unrecorded' in root['flags']:
        L.append(f"lineage INCOMPLETE: root {root['id']} is ours with no recorded parent")
        if with_donors:
            row = lambda d: f"    {d['id']:>8}  n={d['n']:<4} {d['s']:.10f}  {d['added'][:16]}  {', '.join(d['sources'])}"
            na = near_above(db, root)
            if na:
                L.append('  hints, same n just above (pre-polish / pre-descent ancestor?):')
                L += [row(d) for d in na]
            ds = donors(db, root)
            if ds:
                L.append('  hints, same fractional side at other n, added earlier (composition donor?):')
                L += [row(d) for d in ds]
            if na or ds:
                L.append(f"  follow a hint with: pk.py lineage ID")
    return '\n'.join(L)


def connect(path):
    return sqlite3.connect(f'file:{path}?mode=ro', uri=True)
