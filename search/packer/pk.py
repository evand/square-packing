#!/usr/bin/env python3
"""pk: packing store CLI.  Every packing goes through the store; refer to packings by id or "best of n", never by
copying coordinates.  Store: runs/store.sqlite (env PK_STORE overrides).

  pk.py sync-register                      fetch jlevy register (sparse clone in runs/jlevy), import witnesses + frontier
  pk.py sync-pending [--issue 481 ...]     pending issues + repos from pk/pending.yaml
  pk.py import PATH... --source S [--finder F] [--ref R] [--status raw|polished|certified]
  pk.py import-run RUNDIR... [--below K]   explorer archives (parent links kept)
  pk.py frontier [N...] [--lo A --hi B] [--changed]   best known (register / pending / ours) per n
  pk.py ls N [--below S] [--source PAT] [--tag T] [--nongrid] [--limit 20]
  pk.py get ID|best:N [--fmt ours|couzo|ellsworth|json] [-o FILE]    print / write a packing
  pk.py info ID                            row + sightings + parent chain
  pk.py tag ID TAG... [--rm]
  pk.py export N... --dir D [--source PAT] [--best]   n-<n>.txt per n (inputs for legacy tools)
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
    p = sp.add_parser('sync-pending'); p.add_argument('--issue', type=int, nargs='*')
    p = sp.add_parser('import'); p.add_argument('paths', nargs='+'); p.add_argument('--source', required=True)
    p.add_argument('--finder'); p.add_argument('--ref'); p.add_argument('--status', default='raw')
    p.add_argument('--kind'); p.add_argument('--tol', type=float, default=1e-9)
    p = sp.add_parser('import-run'); p.add_argument('runs', nargs='+'); p.add_argument('--below', type=float)
    p.add_argument('--below-k', action='store_true', help='only entries below ceil(sqrt n)')
    p.add_argument('--procs', type=int, default=12)
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
    a = ap.parse_args()
    st = Store(a.store) if a.store else Store()

    if a.cmd == 'sync-register':
        d = os.path.join(HERE, 'runs', 'jlevy')
        sha = sources.fetch_register(d)
        c = sources.import_register(st, f'{d}/packing', sha[:7])
        print(f'register {sha[:7]}: {c["rows"]} frontier rows, {c["witnesses"]} witnesses ({c["new"]} new packings); '
              f'failed {c["failed"][:5]}')
    elif a.cmd == 'sync-pending':
        sources.import_pending(st, only=a.issue)
    elif a.cmd == 'import':
        print(sources.import_paths(st, a.paths, a.source, a.finder, a.ref, a.kind, a.status, tol=a.tol))
    elif a.cmd == 'import-run':
        runs = sorted({os.path.dirname(os.path.abspath(f)) for r in a.runs for f in
                       ([f'{r}/archive.jsonl'] if os.path.exists(f'{r}/archive.jsonl') else
                        glob.glob(f'{r}/**/archive.jsonl', recursive=True))})
        print(f'{len(runs)} archives')
        print(sources.import_runs(st, runs, below=a.below, procs=a.procs, below_k=a.below_k))
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
    st.close()


if __name__ == '__main__':
    main()
