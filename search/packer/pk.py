#!/usr/bin/env python3
"""pk: packing store CLI.  Every packing goes through the store; refer to packings by id or "best of n", never by
copying coordinates.  Store: runs/store.sqlite (env PK_STORE overrides).

  pk.py sync-register                      fetch jlevy register (sparse clone in runs/jlevy), import witnesses + frontier
  pk.py sync-kingbird [--date D] [--all]   crawl Ellsworth's catalogue (all pages, incl. revisions), dated snapshot in
                                           runs/kingbird/<date>/, parse, import with provenance (pk/kingbird.py)
  pk.py sync-kingbird --snapshot DIR       parse + import an existing snapshot dir (no network)
  pk.py watch [--apply]                    jlevy/squares issues updated since the last watch, synced from their
                                           current text (+ pk/overrides.yaml); new commits in watched repos -> every
                                           packing file imported (pk/watch.py; dry run without --apply)
  pk.py sync-issues [N ...] [--all] [--dry]  re-sync the given issues (--all: every issue the store has seen +
                                           overrides.yaml) from their current text
  pk.py import PATH... --source S [--finder F] [--ref R] [--status raw|polished|certified]
  pk.py import-run RUNDIR... [--below K]   explorer archives (parent links kept)
  pk.py frontier [N...] [--lo A --hi B] [--changed]   best known (register / pending / ours) per n
  pk.py ls N [--below S] [--source PAT] [--tag T] [--nongrid] [--limit 20]
  pk.py get ID|best:N [--fmt ours|couzo|ellsworth|json] [-o FILE]    print / write a packing
  pk.py info ID                            row + sightings + parent chain
  pk.py tag ID TAG... [--rm]
  pk.py export N... --dir D [--source PAT]   best packing per n -> n-<n>.txt (inputs for legacy tools)
  pk.py nearmiss list [--thr 1e-3] [--per-n 2]   near misses a scan would queue + the queue (read-only; pk/nearmiss.py)
  pk.py nearmiss tag [--from DIR ...]      scan the store (+ compose2 / harden / fresh output dirs) -> near-miss queue
  pk.py nearmiss requeue ID...             queue items back to queued
  pk.py finish [--limit N] [--stages polish,harden1,harden] [--procs 4] [--candidate] [--dry-run]
                                           drain the queue: fq long polish, harden rounds; below frontier -> candidate.py
                                           (printed only without --candidate)
  pk.py prune --dry-run [--policy F]       what pk/policy.py would keep (tiers, threshold sweep, record lineage)
  pk.py archive-dedupe [ROOT...] [--apply]  drop explorer archive files the store regenerates byte for byte
  pk.py archive-restore RUN...             write them back (reads do this per file automatically)
  pk.py ingest [--apply] [--since ISO] [--no-fresh]   sync-register + watch --apply + fresh.py since
                                           the last ingest -> near-miss queue (pk/ingest.py; dry run without --apply)
"""
import argparse, glob, json, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pk.store import Store
from pk import sources

HERE = os.path.dirname(os.path.abspath(__file__))


def resolve(st, ref, **kw):
    if str(ref).startswith('best:'):
        p = st.best(int(ref[5:]), **kw)
        if p is None:
            sys.exit(f'nothing stored for n={ref[5:]}')
        return p
    return st.get(int(ref))


def fmt_row(r, ref=None):
    gap = '' if ref is None else f' {r["s"] - ref:+.2e}'
    return (f'{r["id"]:7d} n={r["n"]:3d} s={r["s"]:.12f}{gap} {r["status"]:10s} tilt={r["ntilted"] or 0:3d} '
            f'grid={r["grid"]} {r["sources"] or ""} [{r["finders"] or ""}] {",".join(json.loads(r["tags"]))}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--store', default=None)
    sp = ap.add_subparsers(dest='cmd', required=True)
    sp.add_parser('sync-register')
    p = sp.add_parser('sync-kingbird'); p.add_argument('--date'); p.add_argument('--all', action='store_true')
    p.add_argument('--snapshot'); p.add_argument('--procs', type=int, default=4)
    p = sp.add_parser('sync-issues'); p.add_argument('issues', type=int, nargs='*'); p.add_argument('--all', action='store_true')
    p.add_argument('--dry', action='store_true')
    p = sp.add_parser('watch'); p.add_argument('--apply', action='store_true')
    p.add_argument('--max-commits', type=int, default=20)
    p = sp.add_parser('import'); p.add_argument('paths', nargs='+'); p.add_argument('--source', required=True)
    p.add_argument('--finder'); p.add_argument('--ref'); p.add_argument('--status', default='raw')
    p.add_argument('--kind'); p.add_argument('--tol', type=float, default=1e-9)
    p = sp.add_parser('import-run'); p.add_argument('runs', nargs='+'); p.add_argument('--below', type=float)
    p.add_argument('--below-k', action='store_true', help='only entries below ceil(sqrt n)')
    p.add_argument('--procs', type=int, default=12)
    p = sp.add_parser('import-finals', help='battery / corpus finals with a certified / weak verdict -> packing table'); p.add_argument('--n', type=int, default=110)
    p = sp.add_parser('frontier'); p.add_argument('ns', nargs='*', type=int); p.add_argument('--lo', type=int, default=1)
    p.add_argument('--hi', type=int, default=10000); p.add_argument('--pending-only', action='store_true')
    p = sp.add_parser('ls'); p.add_argument('n', type=int); p.add_argument('--below', type=float)
    p.add_argument('--source'); p.add_argument('--tag'); p.add_argument('--status'); p.add_argument('--nongrid', action='store_true')
    p.add_argument('--limit', type=int, default=20)
    p = sp.add_parser('get'); p.add_argument('ref'); p.add_argument('--fmt', default='ours'); p.add_argument('-o')
    p.add_argument('--source')
    p = sp.add_parser('info'); p.add_argument('ref')
    p = sp.add_parser('tag'); p.add_argument('id', type=int); p.add_argument('tags', nargs='+'); p.add_argument('--rm', action='store_true')
    p = sp.add_parser('export'); p.add_argument('ns', nargs='+', type=int); p.add_argument('--dir', required=True)
    p.add_argument('--source'); p.add_argument('--fmt', default='ours')
    p = sp.add_parser('nearmiss'); p.add_argument('action', choices=['list', 'tag', 'requeue'])
    p.add_argument('ids', nargs='*', type=int); p.add_argument('--thr', type=float, default=1e-3)
    p.add_argument('--nmax', type=int, default=324); p.add_argument('--ns', type=int, nargs='*')
    p.add_argument('--per-n', type=int, default=2); p.add_argument('--include-public', action='store_true')
    p.add_argument('--min-gap', type=float, default=2e-5, help='scan floor (closer = the frontier minimum re-found)')
    p.add_argument('--from', dest='dirs', nargs='*', default=[], help='compose2 / fresh / harden output dirs')
    p.add_argument('--no-scan', action='store_true', help='tag: only the --from dirs')
    p.add_argument('--limit', type=int, default=40)
    p = sp.add_parser('finish'); p.add_argument('--limit', type=int, default=None, help='items (default: all queued)'); p.add_argument('--ids', type=int, nargs='*')
    p.add_argument('--stages', default='polish,harden1,harden'); p.add_argument('--procs', type=int, default=4)
    p.add_argument('--nm-gate', type=float, default=5e-4, help='NM harden only while the gap is below this')
    p.add_argument('--candidate', action='store_true', help='run candidate.py on items below the frontier (default: print)')
    p.add_argument('--dry-run', action='store_true', help='print the items and their bars (read-only)')
    p.add_argument('--out')
    p = sp.add_parser('prune', help='store thinning per pk/policy.py (report only for now)')
    p.add_argument('--dry-run', action='store_true', required=True); p.add_argument('--policy')
    p = sp.add_parser('archive-dedupe', help='remove explorer archive files the store regenerates byte for byte (pk/archive.py)')
    p.add_argument('roots', nargs='*'); p.add_argument('--apply', action='store_true')
    p.add_argument('--min-age-h', type=float, default=12.0); p.add_argument('--procs', type=int, default=16)
    p = sp.add_parser('archive-restore', help='write back the archive files dedupe removed'); p.add_argument('runs', nargs='+')
    p = sp.add_parser('ingest'); p.add_argument('--apply', action='store_true'); p.add_argument('--since')
    p.add_argument('--procs', type=int, default=16); p.add_argument('--no-fresh', action='store_true')
    p.add_argument('--nmax', type=int, default=324)
    a = ap.parse_args()
    spath = a.store or os.environ.get('PK_STORE') or os.path.join(HERE, 'runs', 'store.sqlite')
    if a.cmd in ('archive-dedupe', 'archive-restore'):          # store read-only here
        from pk import archive
        if a.cmd == 'archive-dedupe':
            archive.dedupe(a.roots or [os.path.join(HERE, 'runs')], apply=a.apply, min_age_h=a.min_age_h, procs=a.procs)
        else:
            for r in a.runs:
                print(r, archive.restore(r), 'files restored')
        return
    if a.cmd == 'nearmiss' and a.action == 'list':            # read-only: safe while a run writes the store
        from pk import nearmiss
        nearmiss.print_list(nearmiss.ro_db(spath), nearmiss.frontier_fn(spath), a)
        return
    if a.cmd == 'ingest':
        from pk import ingest
        ingest.run(store=a.store, apply=a.apply, since=a.since, procs=a.procs, fresh=not a.no_fresh, nmax=a.nmax)
        return
    if a.cmd == 'finish' and a.dry_run:
        from pk import nearmiss
        nearmiss.finish(nearmiss.ReadOnly(spath), None, a.limit, a.stages.split(','), dry=True,
                        fnow=nearmiss.frontier_fn(spath), ids=a.ids, nm_gate=a.nm_gate)
        return
    st = Store(a.store) if a.store else Store()

    if a.cmd == 'prune':
        from pk import policy, prune
        prune.report(st, policy.load(a.policy))
        return
    if a.cmd == 'sync-register':
        d = os.path.join(HERE, 'runs', 'jlevy')
        sha = sources.fetch_register(d)
        c = sources.import_register(st, f'{d}/packing', sha[:7])
        st.set_sync('register.sha', sha)
        print(f'register {sha[:7]}: {c["rows"]} frontier rows, {c["witnesses"]} witnesses ({c["new"]} new packings); '
              f'failed {c["failed"][:5]}')
    elif a.cmd == 'sync-kingbird':
        from pk import kingbird
        if a.snapshot:
            recs = kingbird.parse_snapshot(a.snapshot, a.procs)
            kingbird.import_snapshot(st, a.snapshot, recs)
        else:
            kingbird.sync(st, a.date, a.all, a.procs)
    elif a.cmd == 'sync-issues':
        from pk import watch
        if not a.issues and not a.all:
            raise SystemExit('sync-issues: give issue numbers or --all')
        watch.run(st, apply=not a.dry, only=a.issues, everything=a.all, repos=False)
    elif a.cmd == 'watch':
        from pk import watch
        watch.run(st, apply=a.apply, max_commits=a.max_commits)
    elif a.cmd == 'import':
        print(sources.import_paths(st, a.paths, a.source, a.finder, a.ref, a.kind, a.status, tol=a.tol))
    elif a.cmd == 'import-run':
        runs = sorted({os.path.dirname(os.path.abspath(f)) for r in a.runs for f in
                       ([f'{r}/archive.jsonl'] if os.path.exists(f'{r}/archive.jsonl') else
                        glob.glob(f'{r}/**/archive.jsonl', recursive=True))})
        print(f'{len(runs)} archives')
        print(sources.import_runs(st, runs, below=a.below, procs=a.procs, below_k=a.below_k))
    elif a.cmd == 'import-finals':
        print(sources.import_finals(st, a.n))
    elif a.cmd == 'frontier':
        F = st.frontier(); R = st.frontier(pending=False)
        ns = a.ns or sorted(n for n in F if a.lo <= n <= a.hi)
        print(f'{"n":>4} {"best known":>16} {"by":24} {"register":>16} {"ours (store)":>16} {"gap":>9}  sub-k nongrid')
        for n in ns:
            f, r = F.get(n), R.get(n)
            if f is None:
                continue
            if a.pending_only and f['pending'] == 0:
                continue
            ours = st.query(n=n, source='run:%', limit=1)
            k = math.ceil(math.sqrt(n) - 1e-12)
            nsub = st.db.execute("SELECT count(*) FROM packing WHERE n=? AND s<? AND grid=0 AND status!='infeasible'",
                                 (n, k)).fetchone()[0]
            o = ours[0]['s'] if ours else None
            print(f'{n:4d} {f["s"]:16.12f} {(f["finder"] or "")[:16] + ("*" if f["pending"] else ""):24s} '
                  f'{r["s"] if r else float("nan"):16.12f} {o if o else float("nan"):16.12f} '
                  f'{(o - f["s"]) if o else float("nan"):+9.2e}  {nsub}')
        print('* = pending (open issue)')
    elif a.cmd == 'ls':
        F = st.frontier().get(a.n)
        for r in st.query(n=a.n, below=a.below, source=a.source, tag=a.tag, status=a.status, nongrid=a.nongrid,
                          limit=a.limit):
            print(fmt_row(r, F['s'] if F else None))
    elif a.cmd == 'get':
        p = resolve(st, a.ref, source=a.source)
        if a.o:
            p.write(a.o, a.fmt); print(f'wrote {a.o}: n={p.n} s={p.s!r} (id {p.meta["id"]})', file=sys.stderr)
        else:
            sys.stdout.write(p.text(a.fmt))
    elif a.cmd == 'info':
        p = resolve(st, a.ref)
        pid = p.meta['id']
        print(fmt_row(st.query(n=p.n, order='id')[0] if False else
                      next(r for r in st.query(n=p.n, limit=None, feasible=False) if r['id'] == pid)))
        for s in st.sightings(pid):
            print(f'   seen: {s["source"]} [{s["finder"]}] ref={s["ref"]} kind={s["kind"]} {s["added"]}')
        chain, q = [], p.meta['parent']
        while q is not None and len(chain) < 50:
            r = st.db.execute('SELECT id, s, parent FROM packing WHERE id=?', (q,)).fetchone()
            chain.append(f'{r["id"]}({r["s"]:.9f})'); q = r['parent']
        if chain:
            print('   parents:', ' <- '.join(chain))
    elif a.cmd == 'tag':
        st.tag(a.id, *([] if a.rm else a.tags), remove=a.tags if a.rm else ())
    elif a.cmd == 'export':
        os.makedirs(a.dir, exist_ok=True)
        for n in a.ns:
            p = st.best(n, source=a.source)
            if p is not None:
                p.write(f'{a.dir}/n-{n}.txt', a.fmt); print(f'{n} {p.s:.12f} id {p.meta["id"]}')
    elif a.cmd == 'nearmiss':
        from pk import nearmiss
        if a.action == 'requeue':
            for i in a.ids:
                nearmiss.set_status(st, i, 'queued')
        else:
            fnow = nearmiss.frontier_fn(spath)
            k1 = nearmiss.enqueue(st, nearmiss.from_outputs(st, a.dirs, fnow, a.thr, a.nmax, a.per_n),
                                  note='tag --from') if a.dirs else 0
            k2 = 0 if a.no_scan else nearmiss.enqueue(
                st, nearmiss.scan(st.db, fnow, a.thr, nmax=a.nmax, ns=a.ns, per_n=a.per_n,
                                  include_public=a.include_public, min_gap=a.min_gap), note='tag scan')
            print(f'queued {k1} from output dirs, {k2} from the store scan')
    elif a.cmd == 'finish':
        from pk import nearmiss
        import datetime, scratch
        scratch.scope('finish')
        out = a.out or os.path.join(HERE, 'runs', 'finish', datetime.datetime.now().strftime('%Y%m%dT%H%M%S'))
        R = nearmiss.finish(st, out, a.limit, a.stages.split(','), a.procs, a.candidate,
                            fnow=nearmiss.frontier_fn(spath), ids=a.ids, nm_gate=a.nm_gate)
        print(f"{len(R)} items: {sum(1 for r in R if r.get('candidate'))} below the frontier, "
              f"{sum(1 for r in R if r.get('err'))} failed; log {out}/finish.jsonl")
    st.close()


if __name__ == '__main__':
    main()
