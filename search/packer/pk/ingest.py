"""One ingestion command (10-14): sync + watch + fresh pass, its near misses into the queue.

  pk.py ingest                     dry run: the plan, --since, and what the store already holds since then (read-only)
  pk.py ingest --apply [--no-fresh | --fresh-later] [--since ISO]
  pk.py ingest --fresh-pending              run the fresh pass for every window --fresh-later deferred (10-28)

--apply runs, each as its own process with the same store (PK_STORE), logs in runs/ingest/<stamp>/:
  1. pk.py sync-register   2. pk.py watch --apply (updated issues + watched repos)
  3. fresh.py run --since <last ingest> --nmax 324 --out runs/ingest/<stamp>/fresh   (unless --no-fresh)
  4. near misses of the fresh pass (compose_lo / compose_hi dirs: near.jsonl arms, polished hits) + a store scan of the
     packings added since -> the near-miss queue (pk/nearmiss.py; drained by `pk.py finish`)
State runs/ingest/state.json: `last` = start time of the last ingest whose steps all succeeded (local time, the clock of
sighting.added), so the next --since covers everything first seen after it.  One line per --apply in runs/ingest/log.jsonl.
--fresh-later (10-28): steps 1-2 only, but the window is kept in state `fresh_pending` (the earliest deferred --since),
so a quick refresh need not wait for fresh.py (~20 min with its compositions); `--fresh-pending` later runs one fresh
pass from that time (steps 3-4) and clears it.  --no-fresh drops the window (the state still advances).
Nothing is posted anywhere: fresh.py's hits go to candidate.py (local checks + import), never to an issue.
"""
from __future__ import annotations
import datetime, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PK = os.path.dirname(HERE)
IDIR = os.path.join(PK, 'runs', 'ingest')


def state_path(root=IDIR):
    return os.path.join(root, 'state.json')


def load_state(root=IDIR):
    try:
        return json.load(open(state_path(root)))
    except FileNotFoundError:
        return dict(last=None, runs=[])


def save_state(S, root=IDIR):
    os.makedirs(root, exist_ok=True)
    tmp = state_path(root) + '.tmp'
    json.dump(S, open(tmp, 'w'), indent=1)
    os.replace(tmp, state_path(root))


def plan(since, out, fresh=True, nmax=324, store=None):
    """[(name, argv)] for --apply (cwd search/packer)."""
    py, pk = sys.executable, os.path.join(PK, 'pk.py')
    base = [py, pk] + (['--store', store] if store else [])
    P = [('sync-register', base + ['sync-register']), ('watch', base + ['watch', '--apply'])]
    if fresh:
        P.append(('fresh', [py, os.path.join(PK, 'fresh.py'), 'run', '--since', since, '--nmax', str(nmax),
                            '--out', os.path.join(out, 'fresh'),
                            '--name', 'ingest_' + os.path.basename(out)]))
    return P


def new_since(db, since):
    """rows (id, n, s, first public sighting) first seen publicly at or after since (fresh.new_packings, read-only)"""
    from .store import public_sql
    pub = public_sql('g')
    return db.execute(f"""SELECT p.id, p.n, p.s, min(g.added) AS first FROM packing p JOIN sighting g ON g.packing = p.id
        WHERE p.status != 'infeasible' AND ({pub}) GROUP BY p.id HAVING first >= ? ORDER BY p.n""", (since,)).fetchall()


def run(store=None, apply=False, since=None, fresh=True, nmax=324, root=IDIR, runner=None, log=print,
        later=False):
    from . import nearmiss as NM
    S = load_state(root)
    since = since or S.get('last')
    start = datetime.datetime.now().isoformat(timespec='seconds')
    if later:
        fresh = False
    if (fresh or later) and not since:
        raise SystemExit('ingest: no previous ingest in state; give --since ISO (local time) the first time')
    out = os.path.join(root, start.replace(':', '').replace('-', ''))
    P = plan(since, out, fresh, nmax, store)
    spath = store or os.environ.get('PK_STORE') or os.path.join(PK, 'runs', 'store.sqlite')
    log(f'ingest {"--apply" if apply else "(dry run)"}: since {since} (state last {S.get("last")}); out {out}')
    for name, cmd in P:
        log(f'  {name}: {" ".join(cmd[1:])}')
    if not apply:
        if since and os.path.exists(spath):
            rows = new_since(NM.ro_db(spath), since)
            log(f'store already holds {len(rows)} packings first seen publicly since {since} '
                f'(fresh would use these + whatever the syncs bring): ' +
                ' '.join(f"{r['n']}:{r['s']:.6f}" for r in rows[:30]) + (' ...' if len(rows) > 30 else ''))
        log('dry run: nothing fetched, run or written; --apply to run')
        return dict(since=since, plan=P, apply=False)
    os.makedirs(out, exist_ok=True)
    env = dict(os.environ, PK_STORE=os.path.abspath(spath))
    runner = runner or (lambda cmd, logf: subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT, cwd=PK,
                                                         env=env).returncode)
    rec = dict(start=start, since=since, out=out, steps={})
    ok = True
    for name, cmd in P:
        with open(os.path.join(out, f'{name}.log'), 'w') as f:
            rc = runner(cmd, f)
        rec['steps'][name] = rc
        log(f'  {name}: rc {rc}')
        ok &= rc == 0
    k1, k2 = _nearmiss(spath, out, since, start, nmax, log)
    rec.update(nearmiss_fresh=k1, nearmiss_scan=k2, ok=ok, end=datetime.datetime.now().isoformat(timespec='seconds'))
    if ok:
        S['last'] = start
        if later:
            S['fresh_pending'] = min(x for x in (S.get('fresh_pending'), since) if x)
            rec['fresh_pending'] = S['fresh_pending']
            log(f"fresh pass deferred: pending since {S['fresh_pending']} (pk.py ingest --fresh-pending)")
    else:
        log('a step failed: state not advanced (the next ingest repeats this window)')
    S.setdefault('runs', []).append(dict(start=start, ok=ok, out=out))
    save_state(S, root)
    with open(os.path.join(root, 'log.jsonl'), 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return rec


def _nearmiss(spath, out, since, start, nmax, log):
    """near misses: fresh output dirs + store scan of everything added since `since` -> the queue.  (k_fresh, k_scan)"""
    from . import nearmiss as NM
    from .store import Store
    st = Store(spath)
    fnow = NM.frontier_fn(spath)
    dirs = [d for d in (os.path.join(out, 'fresh', f'compose_{t}') for t in ('lo', 'hi')) if os.path.isdir(d)]
    items = NM.from_outputs(st, dirs, fnow, nmax=nmax) if dirs else []
    k1 = NM.enqueue(st, items, note=f'ingest {start}: fresh outputs')
    k2 = 0
    if since:
        ids = [r[0] for r in st.db.execute('SELECT id FROM packing WHERE added >= ?', (since,))]
        k2 = NM.enqueue(st, NM.scan(st.db, fnow, nmax=nmax, ids=ids), note=f'ingest {start}: store scan')
    log(f'near-miss queue: +{k1} from fresh outputs, +{k2} from the store scan (pk.py finish drains it)')
    st.close()
    return k1, k2


def run_pending(store=None, nmax=324, root=IDIR, runner=None, log=print):
    """the fresh pass for the windows --fresh-later deferred: one fresh.py run from the earliest pending time, then
    its near misses; clears `fresh_pending` on success unless a newer ingest changed it meanwhile."""
    S = load_state(root)
    since = S.get('fresh_pending')
    if not since:
        log('ingest --fresh-pending: nothing pending')
        return None
    start = datetime.datetime.now().isoformat(timespec='seconds')
    out = os.path.join(root, start.replace(':', '').replace('-', '') + '-fresh')
    (name, cmd), = [x for x in plan(since, out, True, nmax, store) if x[0] == 'fresh']
    spath = store or os.environ.get('PK_STORE') or os.path.join(PK, 'runs', 'store.sqlite')
    os.makedirs(out, exist_ok=True)
    env = dict(os.environ, PK_STORE=os.path.abspath(spath))
    runner = runner or (lambda cmd, logf: subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT, cwd=PK,
                                                         env=env).returncode)
    log(f'ingest --fresh-pending: since {since}; out {out}')
    with open(os.path.join(out, 'fresh.log'), 'w') as f:
        rc = runner(cmd, f)
    log(f'  fresh: rc {rc}')
    k1, k2 = _nearmiss(spath, out, since, start, nmax, log)
    S = load_state(root)                                  # re-read: an ingest may have run meanwhile
    if rc == 0 and S.get('fresh_pending') == since:
        S.pop('fresh_pending')
    elif rc:
        log('fresh failed: still pending')
    S.setdefault('runs', []).append(dict(start=start, ok=rc == 0, out=out, fresh_pending=since))
    save_state(S, root)
    with open(os.path.join(root, 'log.jsonl'), 'a') as f:
        f.write(json.dumps(dict(start=start, since=since, out=out, steps={'fresh': rc}, nearmiss_fresh=k1,
                                nearmiss_scan=k2, ok=rc == 0, pending=True)) + '\n')
    return rc

