"""Packing store: every packing we touch (register, pending issues, other people's repos, our runs and candidates) in one
sqlite file, deduplicated up to box symmetry, with provenance (sightings), status and tags.

Tables
  packing(id, hash UNIQUE, n, s, s_str, coords BLOB f64 (n,3), status, s_exact, fp, ntilted, grid, pen, parent, tags, added)
  sighting(packing, source, finder, ref, kind, added)   -- where a packing was seen (many per packing)
  frontier(n, s, s_str, finder, source, ref, pending, updated)  -- best known per n, one row per (n, source)
  basin(n, key, cls, S, status, added)  -- exactsolve verdict per screened f64 side key (f'{s:.9f}'): a lookup cache,
                                           not a minimum (pk/identity.py)
  minimum(n, key, S, cls, status, src, pid, added, updated)  -- one row per exact minimum (10-11): key = identity.min_key(S)
                                           (1e-17), cls certified | weak; census / "new minimum" counts read this table
  copy(n, side_key, min_key, via, added)  -- a screened side that descends / solves to a minimum (copies of one minimum)

status: raw (imported, f64-feasible within tol) < screened < polished < certified; or not-min / infeasible / unresolved.
"""
from __future__ import annotations
import datetime, json, math, os, sqlite3
import numpy as np
from .packing import Packing

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.environ.get('PK_STORE', os.path.join(HERE, '..', 'runs', 'store.sqlite'))
RANK = dict(infeasible=-1, raw=0, screened=1, polished=2, unresolved=2, **{'not-min': 3}, weak=4, certified=5)

SCHEMA = """
CREATE TABLE IF NOT EXISTS packing (
  id INTEGER PRIMARY KEY, hash TEXT UNIQUE NOT NULL, n INTEGER NOT NULL, s REAL NOT NULL, s_str TEXT,
  coords BLOB NOT NULL, status TEXT NOT NULL DEFAULT 'raw', s_exact TEXT, fp TEXT, ntilted INTEGER, grid INTEGER,
  pen REAL, parent INTEGER, tags TEXT NOT NULL DEFAULT '[]', added TEXT NOT NULL);
CREATE INDEX IF NOT EXISTS packing_ns ON packing(n, s);
CREATE TABLE IF NOT EXISTS sighting (
  packing INTEGER NOT NULL REFERENCES packing(id), source TEXT NOT NULL, finder TEXT, ref TEXT, kind TEXT, added TEXT,
  UNIQUE(packing, source, ref));
CREATE INDEX IF NOT EXISTS sighting_src ON sighting(source);
CREATE TABLE IF NOT EXISTS basin (
  n INTEGER NOT NULL, key TEXT NOT NULL, cls TEXT NOT NULL, S TEXT, status TEXT, added TEXT, UNIQUE(n, key));
CREATE TABLE IF NOT EXISTS minimum (
  n INTEGER NOT NULL, key TEXT NOT NULL, S TEXT NOT NULL, cls TEXT NOT NULL, status TEXT, src TEXT, pid INTEGER,
  added TEXT, updated TEXT, UNIQUE(n, key));
CREATE TABLE IF NOT EXISTS copy (
  n INTEGER NOT NULL, side_key TEXT NOT NULL, min_key TEXT NOT NULL, via TEXT, added TEXT, UNIQUE(n, side_key, min_key));
CREATE TABLE IF NOT EXISTS frontier (
  n INTEGER NOT NULL, s REAL NOT NULL, s_str TEXT, finder TEXT, source TEXT NOT NULL, ref TEXT, pending INTEGER NOT NULL,
  packing INTEGER, updated TEXT, UNIQUE(n, source));
"""


def now():
    return datetime.datetime.now().isoformat(timespec='seconds')


class Store:
    def __init__(self, path=DEFAULT):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        self.db = sqlite3.connect(path, timeout=60)
        self.db.row_factory = sqlite3.Row
        self.db.execute('PRAGMA journal_mode=WAL')
        self.db.executescript(SCHEMA)

    def close(self):
        self.db.commit(); self.db.close()

    # ---------- write ----------
    def add(self, p: Packing, source: str, finder=None, ref=None, kind=None, status='raw', parent=None, tags=(),
            check=True, tol=1e-9, fp=None, commit=True, pre=None, s_exact=None):
        """Insert (or find) a packing; record the sighting.  Returns (id, created).  Infeasible (pen > tol) packings are
        stored with status 'infeasible' so the rejection is remembered."""
        pre = pre or {}                                  # precomputed hash / grid / ntilted / pen (parallel importers)
        h = pre.get('hash') or p.canon_hash()
        row = self.db.execute('SELECT id, status, tags FROM packing WHERE hash=?', (h,)).fetchone()
        created = row is None
        if created:
            pen = pre['pen'] if 'pen' in pre else (p.max_pen() if check else None)
            st = 'infeasible' if (pen is not None and pen > tol) else status
            cur = self.db.execute(
                'INSERT INTO packing(hash,n,s,s_str,coords,status,s_exact,fp,ntilted,grid,pen,parent,tags,added) '
                'VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                (h, p.n, float(p.s), p.s_str, p.sq.astype('<f8').tobytes(), st, s_exact, fp,
                 pre.get('ntilted', p.ntilted()) if 'ntilted' in pre else p.ntilted(),
                 pre['grid'] if 'grid' in pre else p.grid_lines(), pen, parent, json.dumps(sorted(set(tags))), now()))
            pid = cur.lastrowid
        else:
            pid = row['id']
            if RANK.get(status, 0) > RANK.get(row['status'], 0) and row['status'] != 'infeasible':
                self.db.execute('UPDATE packing SET status=? WHERE id=?', (status, pid))
            if s_exact:
                self.db.execute('UPDATE packing SET s_exact=? WHERE id=? AND s_exact IS NULL', (s_exact, pid))
            if tags:
                t = sorted(set(json.loads(row['tags'])) | set(tags))
                self.db.execute('UPDATE packing SET tags=? WHERE id=?', (json.dumps(t), pid))
        self.db.execute('INSERT OR IGNORE INTO sighting VALUES (?,?,?,?,?,?)', (pid, source, finder, ref, kind, now()))
        if commit:
            self.db.commit()
        return pid, created

    def set(self, pid, **kw):
        if 'tags' in kw:
            kw['tags'] = json.dumps(sorted(set(kw['tags'])))
        cols = ', '.join(f'{k}=?' for k in kw)
        self.db.execute(f'UPDATE packing SET {cols} WHERE id=?', (*kw.values(), pid)); self.db.commit()

    def tag(self, pid, *tags, remove=()):
        t = set(json.loads(self.db.execute('SELECT tags FROM packing WHERE id=?', (pid,)).fetchone()[0]))
        self.set(pid, tags=(t | set(tags)) - set(remove))

    def set_frontier(self, n, s, s_str, finder, source, ref, pending, packing=None):
        self.db.execute('INSERT OR REPLACE INTO frontier VALUES (?,?,?,?,?,?,?,?,?)',
                        (n, float(s), s_str, finder, source, ref, int(pending), packing, now()))

    def set_basin(self, n, key, cls, S=None, status=None):
        """Certification of a basin, keyed by the polished side to 1e-9 (f'{s:.9f}')."""
        self.db.execute('INSERT OR REPLACE INTO basin VALUES (?,?,?,?,?,?)', (n, key, cls, S, status, now()))
        self.db.commit()

    def basin(self, n, key):
        return self.db.execute('SELECT * FROM basin WHERE n=? AND key=?', (n, key)).fetchone()

    # ---------- minima (pk/identity.py) ----------
    def find_minimum(self, n, S):
        """The minimum row whose key is S's key or a neighbour (1e-17), or None."""
        from .identity import same_min_keys
        ks = same_min_keys(S)
        return self.db.execute(f'SELECT * FROM minimum WHERE n=? AND key IN (?,?,?)', (n, *ks)).fetchone()

    def add_minimum(self, n, S, cls, status=None, src=None, pid=None, commit=True):
        """Record an exact minimum (cls certified | weak).  Merges with an existing row within 1e-17: keeps the longer S
        string, upgrades weak -> certified (never down).  Returns (key, is_new)."""
        from .identity import min_key, RANK as IR, CERTIFIED, WEAK
        assert cls in (CERTIFIED, WEAK), cls
        r = self.find_minimum(n, S)
        t = now()
        if r is None:
            k = min_key(S)
            self.db.execute('INSERT INTO minimum VALUES (?,?,?,?,?,?,?,?,?)', (n, k, str(S), cls, status, src, pid, t, t))
            out = (k, True)
        else:
            S2 = str(S) if len(str(S)) > len(r['S']) else r['S']
            up = IR[cls] > IR[r['cls']]
            self.db.execute('UPDATE minimum SET S=?, cls=?, status=?, src=?, updated=? WHERE n=? AND key=?',
                            (S2, cls if up else r['cls'], status if up else r['status'], src if up else r['src'], t, n, r['key']))
            out = (r['key'], False)
        if commit:
            self.db.commit()
        return out

    def record_verdict(self, n, side_key, r, via, pid=None):
        """One exactsolve verdict r (census_exact.classify dict) for a screened side: basin cache row (side key) and,
        for certified / weak, the minimum row + the copy link side_key -> minimum.  Returns the minimum key or None."""
        from .identity import CERTIFIED, WEAK
        cls = r.get('cls', 'unresolved')
        self.db.execute('INSERT OR REPLACE INTO basin VALUES (?,?,?,?,?,?)', (n, side_key, cls, r.get('S'), r.get('status'), now()))
        mk = None
        if cls in (CERTIFIED, WEAK) and r.get('S'):
            mk, _ = self.add_minimum(n, r['S'], cls, r.get('status'), src=via, pid=pid, commit=False)
            self.add_copy(n, side_key, mk, via, commit=False)
        self.db.commit()
        return mk

    def add_copy(self, n, side_key, mkey, via=None, commit=True):
        self.db.execute('INSERT OR IGNORE INTO copy VALUES (?,?,?,?,?)', (n, side_key, mkey, via, now()))
        if commit:
            self.db.commit()

    def minima(self, n, cls=None):
        q, a = 'SELECT * FROM minimum WHERE n=?', [n]
        if cls:
            cl = [cls] if isinstance(cls, str) else list(cls)
            q += f' AND cls IN ({",".join("?" * len(cl))})'; a += cl
        return self.db.execute(q + ' ORDER BY key', a).fetchall()

    # ---------- read ----------
    def get(self, pid) -> Packing:
        r = self.db.execute('SELECT * FROM packing WHERE id=?', (pid,)).fetchone()
        if r is None:
            raise KeyError(pid)
        return self._pk(r)

    def _pk(self, r) -> Packing:
        sq = np.frombuffer(r['coords'], '<f8').reshape(-1, 3)
        return Packing(r['s'], sq, r['s_str'], meta=dict(id=r['id'], status=r['status'], n=r['n'], fp=r['fp'],
                                                         tags=json.loads(r['tags']), parent=r['parent']))

    def query(self, n=None, below=None, above=None, status=None, source=None, tag=None, nongrid=False, limit=None,
              order='s', feasible=True):
        """Rows (sqlite Row, no coords) matching the filters, ordered by side."""
        w, a = [], []
        if n is not None:
            w.append('p.n=?'); a.append(n)
        if below is not None:
            w.append('p.s<?'); a.append(below)
        if above is not None:
            w.append('p.s>?'); a.append(above)
        if status:
            st = [status] if isinstance(status, str) else list(status)
            w.append(f'p.status IN ({",".join("?" * len(st))})'); a += st
        elif feasible:
            w.append("p.status!='infeasible'")
        if nongrid:
            w.append('p.grid=0')
        if tag:
            w.append('p.tags LIKE ?'); a.append(f'%"{tag}"%')
        if source:
            w.append('p.id IN (SELECT packing FROM sighting WHERE source LIKE ?)'); a.append(source)
        q = ('SELECT p.id,p.n,p.s,p.s_str,p.status,p.s_exact,p.ntilted,p.grid,p.pen,p.parent,p.tags,p.fp,'
             '(SELECT group_concat(DISTINCT source) FROM sighting WHERE packing=p.id) AS sources,'
             '(SELECT group_concat(DISTINCT finder) FROM sighting WHERE packing=p.id) AS finders FROM packing p')
        if w:
            q += ' WHERE ' + ' AND '.join(w)
        q += f' ORDER BY {order}'
        if limit:
            q += f' LIMIT {int(limit)}'
        return self.db.execute(q, a).fetchall()

    def best(self, n, **kw) -> Packing | None:
        r = self.query(n=n, limit=1, **kw)
        return self.get(r[0]['id']) if r else None

    def frontier(self, n=None, pending=True):
        """{n: best frontier row} (register + pending if pending=True)."""
        q = 'SELECT * FROM frontier' + ('' if pending else ' WHERE pending=0')
        out = {}
        for r in self.db.execute(q):
            if (n is None or r['n'] == n) and (r['n'] not in out or r['s'] < out[r['n']]['s']):
                out[r['n']] = r
        return out

    def sightings(self, pid):
        return self.db.execute('SELECT * FROM sighting WHERE packing=? ORDER BY added', (pid,)).fetchall()
