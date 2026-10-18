#!/usr/bin/env python3
"""Drawings for every record in the leaderboard's histories: data/credits.json -> www/data/hist.json.

hist.json is {n: [file or null, one per upper-bound holder in credits.json, same order]}, naming a file in
www/data/p/.  Catalogue holders already carry their drawing (`file`); the rest (records posted since the catalogue,
mostly superseded within days) are found by side length among the drawings the site already has (best-, rec-, alt-,
step-), and failing that in the packing store (search/packer/runs/store.sqlite: register witnesses and pending-issue
packings), which are written out as www/data/p/hist-<n>-<event>.json.  Holders posted without a side (s null) are
matched by the issue that posted them.  A drawing our regularized lists reproduce (index.json dup_of) is replaced by
our copy (best-N).  Grid holders (kind 'grid') get null: the page draws the grid itself.

Not in build.sh `site` (needs the local store); its outputs are committed.  Re-run after `credits` when the
histories change:  ./build.sh history
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.dirname(HERE)
ROOT = os.path.dirname(SITE)
P = os.path.join(SITE, 'www', 'data', 'p')
TOL = 2e-6     # sides as posted are rounded (12 digits, sometimes fewer); distinct records here differ by > 1e-6


def disk(f):
    """credits/timeline file name (.svg, possibly with %27) -> name on disk under www/data/p."""
    return f.replace('.svg', '.json')


def revisions(store, n, ups, row, keep):
    """A catalogue entry improved in place ("Improved by X ... Refined by Y") has one drawing, its latest state, so
    every holder in the run would show it.  The store keeps the catalogue's earlier revisions (square-N_rK.svg,
    kind catalogue-revision, finder = the revision's first named finder, sometimes truncated: "David W", "M"):
    give each earlier holder in the run the next revision, by falling side, whose finder matches its name."""
    import math
    people = json.load(open(os.path.join(SITE, 'data', 'credits.json')))['people']
    i = 0
    while i < len(row):
        j = i
        while j + 1 < len(row) and row[j + 1] and row[j + 1] == row[i]:
            j += 1
        if j == i or not row[i]:
            i = j + 1; continue
        last_s = float(json.load(open(os.path.join(P, row[j])))['s'])
        top = math.ceil(math.sqrt(n) - 1e-12) - 1e-9          # below the holder before the run, if it has a drawing
        if i and row[i - 1]:
            top = min(top, float(json.load(open(os.path.join(P, row[i - 1])))['s']) - 1e-9)
        q = store.db.execute(
            "SELECT DISTINCT p.id, p.s, s.finder FROM packing p JOIN sighting s ON s.packing = p.id "
            "WHERE p.n = ? AND s.kind = 'catalogue-revision' AND p.s > ? AND p.s < ? ORDER BY p.s DESC",
            (n, last_s + 1e-9, top)).fetchall()
        revs, k = [], 0
        seen = set()
        for r in q:
            if r['id'] not in seen:
                seen.add(r['id']); revs.append(r)
        for t in range(i, j):           # the last holder of the run keeps the current drawing
            names = [people.get(w, {}).get('name', w) for w in ups[t].get('who') or []]
            st = ups[t].get('s')
            if st is not None and abs(st - last_s) < TOL:
                continue                # same side as the drawing: it is this holder's packing (or equal to it)
            for m in range(k, len(revs)):
                fz = (revs[m]['finder'] or '').strip()
                if st is not None and abs(revs[m]['s'] - st) > TOL:
                    continue
                if fz and any(nm.startswith(fz) or fz.startswith(nm) for nm in names):
                    pk = store.get(revs[m]['id'])
                    f = f"hist-{n}-r{revs[m]['id']}.json"
                    json.dump({'name': f, 's': repr(float(pk.s)), 'n': n, 'store_id': revs[m]['id'], 'revision': True,
                               'squares': [[float(a), float(b), float(c)] for a, b, c in pk.sq]},
                              open(os.path.join(P, f), 'w'), separators=(',', ':'))
                    keep.add(f); row[t] = f; k = m + 1
                    break
        # a holder with no revision of its own shows the next drawing after it (the page says so), never an older one
        for t in range(j - 1, i - 1, -1):
            if row[t] == row[j] and row[t + 1] != row[j]:
                row[t] = row[t + 1]
        i = j + 1


def main():
    C = json.load(open(os.path.join(SITE, 'data', 'credits.json')))['n']
    have = {}
    for f in os.listdir(P):
        if f.startswith('hist-'):
            continue
        try:
            d = json.load(open(os.path.join(P, f)))
        except Exception:
            continue
        have.setdefault(int(d['n']), []).append((float(d['s']), f))
    store = None
    sys.path.insert(0, os.path.join(ROOT, 'search', 'packer'))
    try:
        from pk.store import Store
        if os.path.exists(os.path.join(ROOT, 'search', 'packer', 'runs', 'store.sqlite')):
            store = Store()
    except Exception as e:
        print('no store:', e)
    if store is None:
        print('history: the packing store is missing; hist.json not rewritten'); return 1
    # our regularized copy stands in for the drawing it was made from (import_lists.py: dup_of), as on Explore/Compare
    idx = json.load(open(os.path.join(SITE, 'www', 'data', 'index.json')))
    reg = {k.replace('.svg', '.json'): v['dup_of'] for k, v in idx['files'].items() if v.get('dup_of')}
    out, made, miss = {}, [], []
    keep = set()
    for n, d in sorted(C.items(), key=lambda x: int(x[0])):
        n_ = int(n); row = []
        for st in d['upper']:
            f = None
            if st['kind'] == 'grid':
                row.append(None); continue
            if st.get('file') and os.path.exists(os.path.join(P, disk(st['file']))):
                row.append(disk(st['file'])); continue
            s = st.get('s')
            if s is not None:
                c = sorted((abs(x - s), g) for x, g in have.get(n_, []) if abs(x - s) < TOL)
                if c:
                    row.append(c[0][1]); continue
            # the store: by side, preferring a sighting from the posting issue; or, side unknown, by the issue alone
            ev = st.get('event') or ''
            rows = store.query(n=n_, above=(s - TOL) if s else None, below=(s + TOL) if s else None)
            if ev:
                tagged = [r for r in rows if f'pending:{ev}' in (r['sources'] or '').split(',')]
                rows = tagged or (rows if s else [])
            if rows:
                r = rows[0]; pk = store.get(r['id'])
                f = f"hist-{n_}-{(ev or str(r['id'])).lstrip('#')}.json"
                if f in keep:
                    f = f"hist-{n_}-{r['id']}.json"
                json.dump({'name': f, 's': repr(float(pk.s)), 'n': n_, 'store_id': r['id'],
                           'squares': [[float(a), float(b), float(c)] for a, b, c in pk.sq]},
                          open(os.path.join(P, f), 'w'), separators=(',', ':'))
                keep.add(f); made.append((n_, ev, f)); row.append(f)
            else:
                miss.append((n_, st.get('who'), st.get('date'), s, ev)); row.append(None)
        revisions(store, n_, d['upper'], row, keep)
        out[n] = [reg.get(f, f) for f in row]
    for f in os.listdir(P):
        if f.startswith('hist-') and f not in keep:
            os.remove(os.path.join(P, f))
    json.dump(out, open(os.path.join(SITE, 'www', 'data', 'hist.json'), 'w'), separators=(',', ':'))
    print(f'history: {sum(len(r) for r in out.values())} holders, {sum(1 for r in out.values() for x in r if x)} drawn, '
          f'{len(made)} exported from the store, {len(miss)} without a drawing')
    for m in miss:
        print('  no drawing:', m)
    return 0


if __name__ == '__main__':
    sys.exit(main())
