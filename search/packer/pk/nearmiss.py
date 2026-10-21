"""Near-miss queue + finisher (10-14).  A near miss = one of OUR packings (no public sighting) at n <= NMAX, not grid,
feasible, with MIN_GAP (2e-5) < s - frontier_now(n) < THR (1e-3).  The queue is the store table `nearmiss` (one row per packing;
status queued -> running -> done | candidate | failed; a `running` item whose finisher process is gone goes back to queued
at the next `finish`) and the tag 'nearmiss' on the packing (for `pk.py ls --tag`).
A table, not only a tag: the queue needs per-item state (status, gap at enqueue, best result, CPU) that tags cannot hold.

  pk.py nearmiss list [--thr 1e-3] [--per-n 2]      what a scan would enqueue + the queue (read-only)
  pk.py nearmiss tag [--from DIR ...]               scan the store (+ compose2 / harden / fresh output dirs) -> enqueue
  pk.py nearmiss requeue ID ...                     put items back to queued
  pk.py finish [--limit N (default all)] [--stages polish,harden1,harden] [--candidate] [--dry-run]

Sources of near misses: every store packing we made (chain / explore runs once imported (`pk.py import-run`), battery,
corpus, candidates); output dirs (`--from`): compose2 probe arms (`near.jsonl`: layout at the arm's smallest side,
trimmed to n), compose2 / fresh polishes (`results.jsonl` / `polish.jsonl` rows with `s_pol`), harden runs
(`results.jsonl` rows with `ok`).  Those are imported first (source compose2:<dir> / harden:<dir>), then queued.
Grid frontiers (10-20): at an onset n (frontier(n - 1) below ceil(sqrt n), frontier(n) = the grid) non-grid rows within
THR above the grid are near misses too (274, 157, 242, ...): the bar is the grid, a hit is the first sub-grid packing.
Rules: per n at most PER_N queue items (a new one only if it beats every item already there), and none within DISTINCT
(1e-7) in side of an existing item (explorer copies of one minimum).  Finisher outputs (source finish:*) are results,
never new items.

Finisher (`finish`): per item, cheapest first, stopping once below the bar:
  polish   fq quench --loosen 1.0 --no-alm --flip-top 64 --pit 600 (compose2's seam polish)
  harden1  harden.py 1 round, r0 {0.1, 0.2} x noise 0.004, 4 seeds (FH recipe; 8 runs)
  harden   harden.py 4 rounds keep 3, r0 {0.05, 0.1, 0.2} x noise {0.001, 0.002, 0.004}, 6 seeds (NM recipe, PACKER
           "Night 10-14/15"; ~500 runs, ~1 CPU-h); only when the gap is still < --nm-gate (5e-4)
Each improvement goes into the store (source finish:<stage>, parent = the item), the ledger gets one row per item
(Store.record_expansion, mechanism 'finish', measured CPU), the item gets status done / candidate / failed.  Bar =
frontier_now(n) at finisher start (an item that was itself below the frontier at enqueue is compared with the frontier
it was enqueued against).  Below bar - 1e-9 -> candidate.py (`--candidate`; default: the command is printed only).
"""
from __future__ import annotations
import datetime, json, math, os, resource, sqlite3, subprocess, sys, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PK = os.path.dirname(HERE)                                   # search/packer
if PK not in sys.path:
    sys.path.insert(0, PK)                                   # mondrian (near.jsonl trim)
FQ = os.path.join(PK, 'target', 'release', 'fq')
THR = 1e-3            # near-miss threshold: s - frontier_now(n)
TIE = 1e-8            # |gap| <= TIE: a tie with the frontier (eventbench.THR); never an item
MIN_GAP = 2e-5        # scan floor: rows this close above the frontier are mostly the frontier minimum re-found at f64
                      # polish precision.  10-14 drain (runs/nmdrain1014, first 26 items, gaps 1.2e-6..1.6e-5, 6.6 CPU-h):
                      # 15 ended exactly on the frontier minimum (< 4e-11), 11 did not move below their start, 0 gains.
                      # Was 1e-6 (main store 10-14: 1e-8..1e-6 rows at 146, 131, 84, 109, 38, 299, ...)
NMAX = 324            # policy (Evan 10-14): no compute above 324 without a specific reason
NMIN = 10
PER_N = 2             # queue items per n (a new one only if it beats all of them)
DISTINCT = 1e-7       # sides closer than this at one n = the same minimum
BELOW = 1e-9          # candidate: s < bar - BELOW
from .store import PUBLIC, is_public, public_sql   # 10-17: the store's one definition
TERMINAL = ('finish:',)                                       # finisher outputs: results, not items
STAGES = ('polish', 'harden1', 'harden')
NM_GATE = 5e-4
RECIPES = dict(
    polish=['--loosen', '1.0', '--no-alm', '--flip-top', '64', '--pit', '600'],
    harden1=dict(grid='r0=0.1,0.2;noise=0.004', seeds=4, rounds=1, keep=2),
    harden=dict(grid='r0=0.05,0.1,0.2;noise=0.001,0.002,0.004', seeds=6, rounds=4, keep=3),
)


def now():
    return datetime.datetime.now().isoformat(timespec='seconds')


def ro_db(path):
    """read-only connection (no schema writes; safe on a store another process is writing)"""
    db = sqlite3.connect(f'file:{os.path.abspath(path)}?mode=ro', uri=True, timeout=60)
    db.row_factory = sqlite3.Row
    return db


class ReadOnly:
    """a Store stand-in with only a read-only .db (for dry runs on a store another process is writing)"""
    def __init__(self, path):
        self.db = ro_db(path)


def frontier_fn(path=None):
    """eventbench.frontier_now on the store at `path` (default PK_STORE / runs/store.sqlite); caches reset."""
    from . import eventbench as E
    if path:
        E.STORE = os.path.abspath(path)
    E._DB = None; E._LC.clear(); E._BASE.clear()
    return E.frontier_now


def _has_queue(db):
    return db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='nearmiss'").fetchone() is not None


def queue(db, status=None, n=None):
    if not _has_queue(db):
        return []
    w, a = [], []
    if status:
        st = [status] if isinstance(status, str) else list(status)
        w.append(f'status IN ({",".join("?" * len(st))})'); a += st
    if n is not None:
        w.append('n=?'); a.append(n)
    return db.execute('SELECT * FROM nearmiss' + (' WHERE ' + ' AND '.join(w) if w else '') + ' ORDER BY gap, packing',
                      a).fetchall()


def _pub_clause(alias='g'):
    return public_sql(alias)


def _origin(db, pid):
    """first non-public sighting source of a packing (the mechanism that made it)"""
    for r in db.execute('SELECT source FROM sighting WHERE packing=? ORDER BY added', (pid,)):
        if not is_public(r[0]):
            return r[0]
    return None


def select(cands, existing, per_n=PER_N, distinct=DISTINCT):
    """Pure queue rule.  cands: [dict(packing, n, s, gap, ...)]; existing: [dict(n, s, gap)] already queued (any status).
    Keeps, per n in order of gap: not within `distinct` of an existing / kept side, and at most per_n items per n in all,
    except that a candidate better than every item at its n always gets in."""
    byn = {}
    for e in existing:
        byn.setdefault(e['n'], []).append((e['s'], e['gap']))
    out = []
    for c in sorted(cands, key=lambda c: (c['gap'], c['packing'])):
        L = byn.setdefault(c['n'], [])
        if any(abs(c['s'] - s) <= distinct for s, _ in L):
            continue
        if len(L) >= per_n and c['gap'] >= min(g for _, g in L):
            continue
        L.append((c['s'], c['gap']))
        out.append(c)
    return out


def scan(db, frontier, thr=THR, nmin=NMIN, nmax=NMAX, ns=None, ids=None, per_n=PER_N, include_public=False, min_gap=MIN_GAP,
         allow_below=False):
    """Store rows that qualify as near misses and are not queued yet (after `select`).  frontier: n -> side.
    ids: only these packings (fresh imports); allow_below: keep gap < -TIE (a hit: frontier computed before the import)."""
    w = ["p.status!='infeasible'", 'p.grid=0', 'p.n BETWEEN ? AND ?']
    a = [nmin, nmax]
    if not include_public:
        w.append(f'NOT EXISTS (SELECT 1 FROM sighting g WHERE g.packing=p.id AND {_pub_clause()})')
    w.append('NOT EXISTS (SELECT 1 FROM sighting g WHERE g.packing=p.id AND ('
             + ' OR '.join(f"g.source LIKE '{t}%'" for t in TERMINAL) + '))')
    if _has_queue(db):
        w.append('p.id NOT IN (SELECT packing FROM nearmiss)')
    if ids is not None:
        ids = list(ids)
        if not ids:
            return []
        w.append(f'p.id IN ({",".join("?" * len(ids))})'); a += ids
    nset = sorted(set(ns) if ns else {r[0] for r in db.execute(
        'SELECT DISTINCT n FROM packing WHERE n BETWEEN ? AND ?' + (f' AND id IN ({",".join("?" * len(ids))})' if ids else ''),
        [nmin, nmax] + (ids or []))})
    cands = []
    for n in nset:
        if not nmin <= n <= nmax:
            continue
        f = frontier(n)
        if f is None:
            continue
        kg = math.ceil(math.sqrt(n) - 1e-12)
        floor = min_gap
        if f >= kg - 1e-9:
            # the grid is the frontier.  10-20 (274 campaign): at an onset n (n - 1 already below the grid) a non-grid
            # packing just above the grid is a near miss for the first sub-grid packing, so finish towards the grid.
            # No min_gap floor: a non-grid row is never a float copy of the grid minimum.  Elsewhere: nothing to do.
            fp = frontier(n - 1) if n - 1 >= nmin else None
            if fp is None or fp >= kg - 1e-9:
                continue
            floor = 0.0
        lo = -math.inf if allow_below else f + floor
        q = (f'SELECT p.id, p.n, p.s FROM packing p WHERE {" AND ".join(w)} AND p.n=? AND p.s>? AND p.s<? '
             f'ORDER BY p.s LIMIT {max(500, 10 * per_n)}')
        for r in db.execute(q, a + [n, lo, f + thr]):
            g = r['s'] - f
            if abs(g) <= TIE:
                continue
            cands.append(dict(packing=r['id'], n=n, s=r['s'], frontier=f, gap=g))
    existing = [dict(n=r['n'], s=r['s'], gap=r['gap']) for r in queue(db)]
    out = select(cands, existing, per_n)
    for c in out:
        c['origin'] = _origin(db, c['packing'])
    return out


def enqueue(st, items, note=None):
    """items from scan -> nearmiss rows (status queued) + tag 'nearmiss'.  Returns the number added."""
    k = 0
    for c in items:
        cur = st.db.execute('INSERT OR IGNORE INTO nearmiss(packing,n,s,frontier,gap,origin,status,added,updated,note) '
                            "VALUES (?,?,?,?,?,?,'queued',?,?,?)",
                            (c['packing'], c['n'], c['s'], c['frontier'], c['gap'], c.get('origin'), now(), now(), note))
        if cur.rowcount:
            k += 1
            t = set(json.loads(st.db.execute('SELECT tags FROM packing WHERE id=?', (c['packing'],)).fetchone()[0]))
            st.db.execute('UPDATE packing SET tags=? WHERE id=?', (json.dumps(sorted(t | {'nearmiss'})), c['packing']))
    st.db.commit()
    return k


def set_status(st, pid, status, **kw):
    kw.update(status=status, updated=now())
    cols = ', '.join(f'{k}=?' for k in kw)
    st.db.execute(f'UPDATE nearmiss SET {cols} WHERE packing=?', (*kw.values(), pid)); st.db.commit()


# ---------- output dirs (compose2 probe / polish, fresh, harden) -> store rows ----------
def _resolve(path, d):
    for p in (path, os.path.join(d, path), os.path.join(PK, path), os.path.join(d, os.path.basename(path))):
        if os.path.exists(p):
            return os.path.abspath(p)
    return None


def read_outputs(d, thr=THR, nmax=NMAX, frontier=None):
    """[(kind, n, side, Packing)] from one output dir: near.jsonl (compose2 probe arms: layout at S_min, trimmed to n;
    kept when need < thr), results.jsonl / polish.jsonl rows with s_pol + out (compose2, fresh) or ok + out (harden).
    With `frontier`, rows with s - frontier(n) >= thr are dropped before any file is read."""
    from .packing import Packing, read
    out = []
    near = os.path.join(d, 'near.jsonl')
    if os.path.exists(near):
        import mondrian                                        # search/packer on sys.path (pk.py, fresh.py)
        for r in map(json.loads, open(near)):
            if r['n'] > nmax or r.get('need', math.inf) >= thr or not r.get('sq'):
                continue
            if frontier and r['S_min'] - frontier(r['n']) >= thr:
                continue
            X = mondrian.trim(np.array(r['sq'], float), r['n'])
            if len(X) == r['n']:
                out.append(('near', r['n'], float(r['S_min']), Packing(float(r['S_min']), X)))
    for fn in ('results.jsonl', 'polish.jsonl'):
        p = os.path.join(d, fn)
        if not os.path.exists(p):
            continue
        for r in map(json.loads, open(p)):
            if 's_pol' in r and r.get('out'):
                kind, s = 'polish', r['s_pol']
            elif r.get('ok') and r.get('out') and 's' in r:
                kind, s = 'harden', r['s']
            else:
                continue
            if r['n'] > nmax or (frontier and s - frontier(r['n']) >= thr):
                continue
            f = _resolve(r['out'], d)
            q = read(f) if f else None
            if q is not None and q.n == r['n']:
                out.append((kind, r['n'], float(q.s), q))
    return out


def from_outputs(st, dirs, frontier, thr=THR, nmax=NMAX, per_n=PER_N):
    """Import the near rows of output dirs into the store (source compose2:<dir> | harden:<dir>), then scan just them.
    frontier must be evaluated before the import (frontier_now caches on first call; we prefetch here)."""
    rows = []
    for d in dirs:
        R = read_outputs(d, thr, nmax, frontier)
        for _, n, _, _ in R:
            frontier(n)                                         # prefetch (cache) before anything is added
        rows += [(d, *r) for r in R]
    ids = set()
    for d, kind, n, s, p in rows:
        rel = os.path.relpath(os.path.abspath(d), PK)
        src = ('harden:' if kind == 'harden' else 'compose2:') + rel
        pid, _ = st.add(p, src, finder='ours', ref=rel, kind=kind, status='polished' if kind != 'near' else 'raw',
                        tol=1e-7 if kind == 'near' else 1e-9, commit=False)
        ids.add(pid)
    st.db.commit()
    return scan(st.db, frontier, thr, nmax=nmax, ids=ids, per_n=per_n, allow_below=True)


# ---------- finisher ----------
def cpu_now():
    s = resource.getrusage(resource.RUSAGE_SELF); c = resource.getrusage(resource.RUSAGE_CHILDREN)
    return s.ru_utime + s.ru_stime + c.ru_utime + c.ru_stime


def bar(item, fnow):
    """success bar for an item: frontier_now(n), unless the frontier is the item itself (a hit when enqueued)"""
    f = fnow(item['n'])
    if item['gap'] is not None and item['gap'] < 0 and f >= item['s'] - 1e-12:
        return item['frontier']
    return f


def stage_polish(inp, wd, n, target, procs):
    """fq long polish -> [(s, path)] (feasible outputs only), number of proposals"""
    out = os.path.join(wd, 'polish.txt')
    from .engine import QOpt, fq_quench               # the one fq call path (B2)
    q = fq_quench(inp, QOpt.from_args(RECIPES['polish'] + ['--seed', '1']), out=out, timeout=7200)
    try:
        d = json.loads(q.stdout.strip().splitlines()[-1])
    except Exception:
        raise RuntimeError('fq: ' + q.stderr[-300:])
    ok = d.get('min_gap', 0) > -1e-9 and d.get('min_wall', 0) > -1e-9
    return ([(float(d['s']), out)] if ok and os.path.exists(out) else []), 1


def _stage_harden(name):
    def run(inp, wd, n, target, procs):
        r = RECIPES[name]
        hd = os.path.join(wd, name)
        cmd = [sys.executable, os.path.join(PK, 'harden.py'), 'run', '--files', inp, '--out', hd, '--grid', r['grid'],
               '--seeds', str(r['seeds']), '--rounds', str(r['rounds']), '--keep', str(r['keep']), '--procs', str(procs),
               '--target', repr(float(target))]
        with open(os.path.join(wd, f'{name}.log'), 'w') as f:
            rc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=PK).returncode
        res = os.path.join(hd, 'results.jsonl')
        if not os.path.exists(res):
            raise RuntimeError(f'{name}: rc {rc}, no results')
        R = [json.loads(l) for l in open(res)]
        good = sorted((x['s'], _resolve(x['out'], hd)) for x in R if x.get('ok') and x.get('n') == n)
        return [g for g in good if g[1]][:1], len(R)
    return run


RUNNERS = dict(polish=stage_polish, harden1=_stage_harden('harden1'), harden=_stage_harden('harden'))


def candidate_cmd(path, name):
    return [sys.executable, os.path.join(PK, 'candidate.py'), path, '--name', name, '--no-regularize', '--import']


def _pid_alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def reclaim_stale(st, log=print):
    """`running` items whose finisher process is gone -> queued.  The finisher records 'pid <p>' in the note when it starts
    an item; items marked running before that (10-14) carry no pid and are only reported (requeue by hand)."""
    k = 0
    for r in queue(st.db, 'running'):
        m = (r['note'] or '').split()
        if len(m) == 2 and m[0] == 'pid' and m[1].isdigit():
            if not _pid_alive(int(m[1])):
                set_status(st, r['packing'], 'queued', note=f'reclaimed (finisher pid {m[1]} gone)'); k += 1
        else:
            log(f"note: item {r['packing']} is 'running' with no finisher pid (pre-10-15 run?): `pk.py nearmiss requeue "
                f"{r['packing']}` once no finisher runs")
    return k


def finish(st, out, limit=None, stages=STAGES, procs=4, candidate=False, dry=False, fnow=None, runners=None,
           nm_gate=NM_GATE, ids=None, log=print):
    """Drain up to `limit` queued items (None = all; smallest gap first).  Returns [result dict per item]."""
    from .packing import read
    runners = runners or RUNNERS
    fnow = fnow or frontier_fn()
    if not dry:
        k = reclaim_stale(st, log)
        if k:
            log(f'reclaimed {k} stale running item(s)')
    items = queue(st.db, 'queued')
    if ids:
        items = [r for r in items if r['packing'] in set(ids)]
    items = [dict(r) for r in (items if limit is None else items[:limit])]
    for it in items:
        fnow(it['n'])                                       # the bar is fixed before anything is added
    if not dry:
        os.makedirs(out, exist_ok=True)
    results = []
    for it in items:
        pid, n = it['packing'], it['n']
        b = bar(it, fnow)
        log(f"item {pid} n={n} s={it['s']:.12f} gap {it['s'] - b:+.2e} (enqueued {it['gap']:+.2e}, {it['origin']})"
            f" stages {','.join(stages)}")
        if dry:
            results.append(dict(packing=pid, n=n, bar=b, dry=True)); continue
        set_status(st, pid, 'running', note=f'pid {os.getpid()}')
        wd = os.path.join(out, f'nm{pid}')
        os.makedirs(wd, exist_ok=True)
        inp = os.path.join(wd, 'in.txt')
        st.get(pid).write(inp)
        cur_s, cur_path, cur_pid = it['s'], inp, pid
        gridbar = b >= math.ceil(math.sqrt(n) - 1e-12) - 1e-9
        c0, nprop, log_st, err = cpu_now(), 0, [], None
        try:
            for name in stages:
                # 10-17 (Evan): no stop below the bar.  A hit still needs finishing: a deeper / jammed point is lower
                # and often what makes it certifiable (s(229), #674: valid, not first-order jammed)
                if abs(cur_s - b) <= TIE:                  # polished onto the frontier: the record re-found, stop
                    log_st.append(dict(stage=name, skipped='tie with the frontier')); break
                if name == 'harden' and cur_s - b > nm_gate:
                    log_st.append(dict(stage=name, skipped=f'gap {cur_s - b:.1e} > gate {nm_gate:g}')); continue
                swd = os.path.join(wd, name); os.makedirs(swd, exist_ok=True)
                outs, k = runners[name](cur_path, swd, n, min(b, cur_s), procs)   # below the bar: aim below itself
                nprop += k
                if gridbar and outs:
                    # 10-20: at a grid bar, an output with a full line / chain is grid-locked (side >= the grid
                    # whatever it shows): not progress, and never a tie.  Keep the best line-free output.
                    keep = []
                    for s_o, p_o in outs:
                        q = read(p_o)
                        if q is not None and q.n == n and q.grid_lines() == 0:
                            keep.append((s_o, p_o)); break
                    outs = keep
                rec = dict(stage=name, runs=k, best=outs[0][0] if outs else None)
                if outs and outs[0][0] < cur_s - 1e-12:
                    s2, path2 = outs[0]
                    p2 = read(path2)
                    if p2 is not None and p2.n == n:
                        cur_pid, _ = st.add(p2, f'finish:{name}', finder='ours', ref=os.path.relpath(path2, PK),
                                            kind='finish', status='polished', parent=pid, tags=('finish',))
                        cur_s, cur_path = float(p2.s), path2
                        rec['pid'] = cur_pid
                log_st.append(rec)
                log(f"  {name}: {k} runs, best {rec['best']}, now {cur_s:.12f} ({cur_s - b:+.2e})")
        except Exception as ex:
            err = repr(ex)[:300]
        cpu = cpu_now() - c0
        st.record_expansion(pid, n, 'finish', cpu, params=dict(stages=list(stages), log=log_st, nm_gate=nm_gate),
                            run=os.path.relpath(wd, PK), nprop=nprop, best_s=cur_s, frontier_s=b)
        res = dict(packing=pid, n=n, bar=b, s0=it['s'], best_s=cur_s, best_pid=cur_pid, cpu_s=round(cpu, 1),
                   stages=log_st, err=err)
        if err:
            set_status(st, pid, 'failed', best_s=cur_s, best_pid=cur_pid, cpu_s=cpu, note=err)
        elif cur_s < b - BELOW:
            cmd = candidate_cmd(os.path.abspath(cur_path), f'nm{pid}_n{n}')
            res['candidate'] = cmd
            if candidate:
                with open(os.path.join(wd, 'candidate.log'), 'w') as f:
                    res['candidate_rc'] = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=PK).returncode
                note = f"candidate.py rc {res['candidate_rc']}"
            else:
                note = 'candidate (dry run): ' + ' '.join(cmd[1:])
            log(f'  BELOW bar by {b - cur_s:.2e}: {note}')
            set_status(st, pid, 'candidate', best_s=cur_s, best_pid=cur_pid, cpu_s=cpu, note=note)
        else:
            set_status(st, pid, 'done', best_s=cur_s, best_pid=cur_pid, cpu_s=cpu,
                       note=f'best gap {cur_s - b:+.2e}')
        with open(os.path.join(out, 'finish.jsonl'), 'a') as f:
            f.write(json.dumps(res) + '\n')
        results.append(res)
    return results


def print_list(db, frontier, a):
    c = scan(db, frontier, a.thr, nmax=a.nmax, ns=a.ns, per_n=a.per_n, include_public=a.include_public,
             min_gap=a.min_gap)
    print(f'{len(c)} would be queued (gap < {a.thr:g}, n <= {a.nmax}, per n {a.per_n}):')
    for x in c[:a.limit]:
        print(f"  {x['packing']:8d} n={x['n']:3d} s={x['s']:.12f} gap {x['gap']:+.2e}  {x['origin']}")
    Q = queue(db)
    by = {}
    for r in Q:
        by[r['status']] = by.get(r['status'], 0) + 1
    print(f'queue: {len(Q)} items ' + ' '.join(f'{k} {v}' for k, v in sorted(by.items())))
    for r in Q[:a.limit]:
        b = '' if r['best_s'] is None else f" best {r['best_s'] - r['frontier']:+.2e}"
        print(f"  {r['packing']:8d} n={r['n']:3d} gap {r['gap']:+.2e} {r['status']:9s}{b}  {r['origin'] or ''}  {r['note'] or ''}")
