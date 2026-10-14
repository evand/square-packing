"""One ingestion command (10-14): sync + watch + fresh pass, its near misses into the queue.

  pk.py ingest                     dry run: the plan, --since, and what the store already holds since then (read-only)
  pk.py ingest --apply [--procs 16] [--no-fresh] [--since ISO]

--apply runs, each as its own process with the same store (PK_STORE), logs in runs/ingest/<stamp>/:
  1. pk.py sync-register   2. pk.py sync-pending   3. pk.py watch --apply
  4. fresh.py run --since <last ingest> --nmax 324 --out runs/ingest/<stamp>/fresh   (unless --no-fresh)
  5. near misses of the fresh pass (compose_lo / compose_hi dirs: near.jsonl arms, polished hits) + a store scan of the
     packings added since -> the near-miss queue (pk/nearmiss.py; drained by `pk.py finish`)
State runs/ingest/state.json: `last` = start time of the last ingest whose steps all succeeded (local time, the clock of
sighting.added), so the next --since covers everything first seen after it.  One line per --apply in runs/ingest/log.jsonl.
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


def plan(since, out, procs=16, fresh=True, nmax=324, store=None):
    """[(name, argv)] for --apply (cwd search/packer)."""
    py, pk = sys.executable, os.path.join(PK, 'pk.py')
    base = [py, pk] + (['--store', store] if store else [])
    P = [('sync-register', base + ['sync-register']), ('sync-pending', base + ['sync-pending']),
         ('watch', base + ['watch', '--apply'])]
    if fresh:
        P.append(('fresh', [py, os.path.join(PK, 'fresh.py'), 'run', '--since', since, '--nmax', str(nmax),
                            '--procs', str(procs), '--out', os.path.join(out, 'fresh'),
                            '--name', 'ingest_' + os.path.basename(out)]))
    return P


def new_since(db, since):
    """rows (id, n, s, first public sighting) first seen publicly at or after since (fresh.new_packings, read-only)"""
    from .nearmiss import PUBLIC
    pub = ' OR '.join(f"g.source LIKE '{p}%'" for p in PUBLIC)
    return db.execute(f"""SELECT p.id, p.n, p.s, min(g.added) AS first FROM packing p JOIN sighting g ON g.packing = p.id
        WHERE p.status != 'infeasible' AND ({pub}) GROUP BY p.id HAVING first >= ? ORDER BY p.n""", (since,)).fetchall()


def run(store=None, apply=False, since=None, procs=16, fresh=True, nmax=324, root=IDIR, runner=None, log=print):
    from . import nearmiss as NM
    S = load_state(root)
    since = since or S.get('last')
    start = datetime.datetime.now().isoformat(timespec='seconds')
    if fresh and not since:
        raise SystemExit('ingest: no previous ingest in state; give --since ISO (local time) the first time')
    out = os.path.join(root, start.replace(':', '').replace('-', ''))
    P = plan(since, out, procs, fresh, nmax, store)
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
    # near misses: fresh output dirs + store scan of everything added since `since`
    from .store import Store
    st = Store(spath)
    fnow = NM.frontier_fn(spath)
    dirs = [d for d in (os.path.join(out, 'fresh', f'compose_{t}') for t in ('lo', 'hi')) if os.path.isdir(d)]
    items = NM.from_outputs(st, dirs, fnow, nmax=nmax) if dirs else []
    k1 = NM.enqueue(st, items, note=f'ingest {start}: fresh outputs')
    if since:
        ids = [r[0] for r in st.db.execute('SELECT id FROM packing WHERE added >= ?', (since,))]
        k2 = NM.enqueue(st, NM.scan(st.db, fnow, nmax=nmax, ids=ids), note=f'ingest {start}: store scan')
    else:
        k2 = 0
    rec.update(nearmiss_fresh=k1, nearmiss_scan=k2, ok=ok, end=datetime.datetime.now().isoformat(timespec='seconds'))
    log(f'near-miss queue: +{k1} from fresh outputs, +{k2} from the store scan (pk.py finish drains it)')
    st.close()
    if ok:
        S['last'] = start
    else:
        log('a step failed: state not advanced (the next ingest repeats this window)')
    S.setdefault('runs', []).append(dict(start=start, ok=ok, out=out))
    save_state(S, root)
    with open(os.path.join(root, 'log.jsonl'), 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return rec
