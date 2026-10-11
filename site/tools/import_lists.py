#!/usr/bin/env python3
"""Bring the site's record packings up to date from our lists (build.sh step `lists`, after export and timeline).

Sources:
  search/regularize/lists/best.json        the best-choice packing for every n <= 324: a regularized drawing of the best
                                           known packing, certified at its side (REGULARIZE.md); certs/n-NNN.json.gz
  search/regularize/lists/alternates.json  packings at the record side that are truly different: ours (certs/alt/)
                                           and Ellsworth's catalogue drawings (already on the site, annotated here)
  site/data/records2026.json               who posted the records newer than the catalogue, and when (provenance.py)
  search/packer/runs/store.sqlite          coordinates of the history steps that are no longer the record (n <= 100)

Writes, alongside what export.py and timeline.py made from the catalogue:
  www/data/p/best-N.json, alt-N.json, step-N-K.json   squares + analysis (tools/analysis.py), as for the SVGs
  www/data/index.json     files[...] for each (origin 'best' / 'alt' / 'step'); records[n].best for every n <= 324;
                          records[n].credit for records newer than the catalogue; catalogue files at the record side
                          annotated with their relation to the record (alt_relation)
  www/data/timeline.json  one dated entry per record posted since the catalogue (all of them for n <= 100, where
                          Bounds draws the history; the current record elsewhere), the catalogue's entry demoted

Idempotent: entries it added before are removed first.  Analysis is skipped when the output exists and its input
digest is unchanged (--force redoes it).

  import_lists.py [--procs 8] [--force] [--only N,...]
"""
import argparse, gzip, hashlib, json, math, os, sqlite3, sys, time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.abspath(os.path.join(HERE, '..'))
ROOT = os.path.abspath(os.path.join(SITE, '..'))
LISTS = os.path.join(ROOT, 'search', 'regularize', 'lists')
STORE = os.path.join(ROOT, 'search', 'packer', 'runs', 'store.sqlite')
W = os.path.join(SITE, 'www', 'data')
NMAX = 324
MONTHS = ['', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October',
          'November', 'December']
sys.path.insert(0, HERE)


def gz_json(rel):
    return json.load(gzip.open(os.path.join(ROOT, rel.replace('.cert.gz', '.json.gz'))))


def snap(sq):
    """Axis squares carry angles like -1e-70 in the 50-digit configurations; the analysis groups angles mod 90, which
    files those at 89.999...: an extra 'rotation'.  Snap |angle| < 1e-30 (mod 90) to exactly 0 (a move far below the
    certificate's own rounding)."""
    out = []
    for x, y, t in sq:
        f = float(t) % 90.0
        out.append([x, y, '0' if min(f, 90.0 - f) < 1e-30 else t])
    return out


def digest(s, sq):
    return hashlib.sha256(json.dumps([s, sq]).encode()).hexdigest()[:16]


def analyse(job):
    key, s, sq = job
    from export import compact
    from analysis import analyze
    out = os.path.join(W, 'p', key)
    t = time.time()
    try:
        a = analyze(float(s), sq)
    except Exception as e:
        return key, f'ERROR {type(e).__name__}: {e}'
    rec = {'name': key, 's': s, 'n': len(sq), 'squares': [[float(x) for x in q] for q in sq], 'analysis': compact(a),
           'src': digest(s, sq)}
    json.dump(rec, open(out, 'w'), separators=(',', ':'))
    return key, f'ok {time.time() - t:.1f}s'


def summary(a):
    return {'category': a['category'], 'n_angles': a['n_angles'], 'n_rotated': a['n_rotated'],
            'symmetry': a['symmetry']['class'], 'rigid': a['rigid'], 'free': len(a['free']),
            'wedged': len(a['wedged']), 'mobile': sum(a['mobile']), 'waste': a['waste'],
            'angles': [[round(g['theta'], 6), g['count']] for g in a['angle_groups']],
            'min_gap': a['near_misses'][0]['gap'] if a['near_misses'] else None}


SAME_TOL = 1e-6


def same_drawing(A, B):
    """True when two exported packings are the same drawing: equal side (to SAME_TOL) and, under one of the 8
    symmetries of the container, every square of A has a square of B within SAME_TOL in centre and in angle mod 90."""
    s = float(A['s'])
    if abs(s - float(B['s'])) > SAME_TOL or len(A['squares']) != len(B['squares']):
        return False
    from scipy.spatial import cKDTree
    PA = [(float(x), float(y), float(t) % 90) for x, y, t in A['squares']]
    c = s / 2
    for g in range(8):
        PB = []
        for x, y, t in B['squares']:
            x, y, t = float(x) - c, float(y) - c, float(t)
            for _ in range(g % 4):
                x, y, t = -y, x, t + 90
            if g >= 4:
                x, t = -x, -t
            PB.append((x + c, y + c, t % 90))
        tree = cKDTree([(x, y) for x, y, _ in PB])
        ok = True
        for x, y, t in PA:
            d, j = tree.query((x, y))
            da = abs(t - PB[j][2]); da = min(da, 90 - da)
            if d > SAME_TOL or da > SAME_TOL:
                ok = False
                break
        if ok:
            return True
    return False


def ev_date(e):
    y, m, d = map(int, e['date'].split('-'))
    return {'year': y, 'month': m, 'day': d, 'approx': False, 'text': f'{MONTHS[m]} {d}, {y}'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--only', default=None, help='comma list of n to (re)analyse; the index is always rebuilt in full')
    a = ap.parse_args()
    only = {int(x) for x in a.only.split(',')} if a.only else None

    best = {b['n']: b for b in json.load(open(os.path.join(LISTS, 'best.json')))}
    alts = json.load(open(os.path.join(LISTS, 'alternates.json')))
    prov = json.load(open(os.path.join(SITE, 'data', 'records2026.json')))
    EV = prov['events']
    idx = json.load(open(os.path.join(W, 'index.json')))
    TL = json.load(open(os.path.join(W, 'timeline.json')))

    # ---- strip what a previous run added ----
    idx['files'] = {k: v for k, v in idx['files'].items() if not v.get('origin')}
    for v in idx['files'].values():
        v.pop('alt_relation', None)
        v.pop('dup_of', None)
    for n, r in list(idx['records'].items()):
        if r.get('origin'):
            del idx['records'][n]
            continue
        for k in ('best', 'best_s', 'credit', 'newer'):
            r.pop(k, None)
    for n in list(TL):
        TL[n] = [d for d in TL[n] if not d.get('origin')]
        for d in TL[n]:
            if 'was_record' in d:
                d['is_record'] = d.pop('was_record')
            if 'catalogue_file' in d:
                d['file'] = d.pop('catalogue_file')

    # ---- the packings ----
    jobs, files = [], {}
    for n in range(1, NMAX + 1):
        b = best[n]
        d = gz_json(b['cert'])
        P = d['packing']
        key = f'best-{n}.json'
        files[key] = dict(n=n, s=P['S'], sq=snap(P['sq']), origin='best', variant=b['variant'], level=b['level'],
                          level_text=b['level_text'], groups=b['groups'], orientation=b['orientation'],
                          shape=b['shape'], cert=b['cert'], s_exact=d['S_exact'])
    for r in alts:
        if r['source'].startswith('ours'):
            d = gz_json(r['cert'])
            P = d['packing']
            key = f"alt-{r['n']}.json"
            files[key] = dict(n=r['n'], s=P['S'], sq=snap(P['sq']), origin='alt', relation=r['relation'], groups=r['groups'],
                              parent=r.get('parent'), cert=r['cert'], s_exact=d['S_exact'])
    # history steps that are no longer the record: coordinates from the store (the poster's frontier row)
    db = sqlite3.connect(STORE)
    sys.path.insert(0, os.path.join(ROOT, 'search', 'packer'))
    from pk.store import Store
    st = Store(STORE)
    def stored(n, ev, side):
        src = 'register' if not ev.startswith('#') else f'pending:{ev}'
        row = db.execute('select packing, s_str from frontier where n=? and source=? and abs(s-?)<1e-8',
                         (n, src, side)).fetchone()
        if not row or row[0] is None:       # registered since: the register row carries the same packing
            row = db.execute('select packing, s_str from frontier where n=? and source=? and abs(s-?)<1e-8',
                             (n, 'register', side)).fetchone()
        return row

    def dec(s_str):
        if '/' in s_str:
            from fractions import Fraction
            from decimal import Decimal, getcontext
            getcontext().prec = 40
            fr = Fraction(s_str)
            return str(Decimal(fr.numerator) / Decimal(fr.denominator))
        return s_str

    # records the lists could not certify (provenance 'raw'): drawn from the posting, unregularized
    raw_key = {}
    for ns, e in prov['n'].items():
        if not e.get('raw'):
            continue
        n = int(ns)
        drawn = e.get('drawn') or e['first']
        row = stored(n, drawn, float(e['s']))
        if not row or row[0] is None:
            print(f'n={n}: no stored packing for the raw record {drawn}', file=sys.stderr)
            continue
        p = st.get(row[0])
        key = raw_key[n] = f'rec-{n}.json'
        cert = os.path.join(LISTS, 'certs', 'raw', f'n-{n:03d}.cert.gz')     # search/regularize/register_cert.py
        files[key] = dict(n=n, s=dec(row[1]), sq=[[repr(float(x)), repr(float(y)), repr(float(t))] for x, y, t in p.sq],
                          origin='raw', event=drawn, store_id=row[0],
                          cert=os.path.relpath(cert, ROOT) if os.path.exists(cert) else None,
                          credit={k: EV[drawn][k] for k in ('who', 'date', 'label', 'url')})
    for ns, e in prov['n'].items():
        n = int(ns)
        for k, step in enumerate(e.get('steps', [])):
            if step['current']:
                continue
            ev = step['event']
            row = stored(n, ev, step['s'])
            if not row or row[0] is None:
                print(f'n={n}: no stored packing for step {ev} s={step["s"]}', file=sys.stderr)
                continue
            p = st.get(row[0])
            s_str = dec(row[1])
            key = f'step-{n}-{k}.json'
            step['file'] = key
            files[key] = dict(n=n, s=s_str, sq=[[repr(float(x)), repr(float(y)), repr(float(t))] for x, y, t in p.sq],
                              origin='step', event=ev, store_id=row[0],
                              credit={k: EV[ev][k] for k in ('who', 'date', 'label', 'url')})

    # ---- analysis (skipped when unchanged) ----
    todo = []
    for key, f in files.items():
        out = os.path.join(W, 'p', key)
        if only is not None and f['n'] not in only and os.path.exists(out):
            continue
        if not a.force and os.path.exists(out):
            try:
                if json.load(open(out)).get('src') == digest(f['s'], f['sq']):
                    continue
            except ValueError:
                pass
        todo.append((key, f['s'], f['sq']))
    todo.sort(key=lambda j: -len(j[2]))           # big ones first
    print(f'analysing {len(todo)} of {len(files)} packings on {a.procs} processes', flush=True)
    if todo:
        with Pool(a.procs) as pool:
            for key, msg in pool.imap_unordered(analyse, todo):
                print(key, msg, flush=True)

    # ---- index ----
    cat_alt = {}
    for r in alts:
        if r['source'].startswith('Ellsworth'):
            cat_alt[r['source'].split()[-1]] = r['relation']
    for f, rel in cat_alt.items():
        if f in idx['files']:
            idx['files'][f]['alt_relation'] = rel
    for key, f in files.items():
        out = os.path.join(W, 'p', key)
        if not os.path.exists(out):
            print('missing analysis', key, file=sys.stderr)
            continue
        an = json.load(open(out))['analysis']
        e = {k: v for k, v in f.items() if k not in ('sq', 's')}
        e.update(svg=key, pages=[], n_parsed=f['n'], s=f['s'] if f['origin'] in ('step', 'raw') else f['s_exact'][:40],
                 errors=[], summary=summary(an))
        idx['files'][key] = e

    for n in range(1, NMAX + 1):
        key = f'best-{n}.json'
        r = idx['records'].get(str(n))
        if r is None:
            r = idx['records'][str(n)] = {'n': n, 'svg': None, 'svg_link': None, 'pictured_n': n, 's_tex': None,
                                          's_dec': None, 'polys': [], 'prose': None, 'rigid': None, 'origin': 'lists'}
        e = prov['n'].get(str(n))
        r['best'] = raw_key.get(n, key)
        r['best_s'] = e['s'] if n in raw_key else best[n]['side']
        if e:
            ev = EV[e.get('drawn') or e['first']] if n in raw_key else EV[e['first']]
            r['newer'] = True
            r['credit'] = {'who': ev['who'], 'date': ev['date'], 'label': ev['label'], 'url': ev['url'],
                           # lineage (provenance.py LINEAGE): the construction it started from, and the improvements
                           'found_by': ({k: EV[e['found_by']][k] for k in ('who', 'date', 'label', 'url')}
                                        if e.get('found_by') else None),
                           'improved': [{k: EV[x][k] for k in ('who', 'date', 'label', 'url')}
                                        for x in reversed(e.get('improved_in') or [])],
                           'also': [dict(who=EV[x]['who'], date=EV[x]['date'], label=EV[x]['label'], url=EV[x]['url'])
                                    for x in e['also']],
                           'registered': e['registered'], 'refined_by': e['refined_by'],
                           'register_url': f'https://jlevy.github.io/squares/cases/{n}.html'}

    # ---- timeline: the records posted since the catalogue ----
    for ns, e in prov['n'].items():
        n = int(ns)
        L = TL.setdefault(ns, [])
        for d in L:
            if d.get('is_record'):
                d['was_record'] = True
                d['is_record'] = False
        posted = (e.get('drawn') or e['first']) if e.get('raw') else e['first']
        steps = e.get('steps') or [{'s': float(e['s']), 'event': posted, 'current': True}]
        for k, step in enumerate(steps):
            ev = EV[posted if step['current'] else step['event']]
            key = (raw_key.get(n) or f'best-{n}.json') if step['current'] else step.get('file')
            if not key or key not in idx['files']:
                continue
            dt = ev_date(ev)
            who = ev['who']
            verb = 'Found'
            text = f"{verb} by {who}, posted {dt['text']} ({ev['label']})."
            if step['current'] and e.get('found_by'):     # credit the construction it started from (LINEAGE)
                root = EV[e['found_by']]
                who = root['who']
                text = (f"Found by {who} ({root['label']}, {ev_date(root)['text']}); improved by "
                        + '; '.join(f"{EV[x]['who']} ({EV[x]['label']}, {ev_date(EV[x])['text']})"
                                    for x in reversed(e['improved_in'])) + '.')
            L.append({'file': key, 'n': n, 's': float(idx['files'][key]['s']) if not step['current'] else float(e['s']),
                      'events': [{'verb': verb, 'who': who, 'date': dt, 'text': text}], 'prose': text, 's_tex': None,
                      'rigid': False, 'date': [dt['year'], dt['month']], 'is_record': bool(step['current']),
                      'origin': 'lists', 'url': ev['url'], 'label': ev['label']})

    # ---- one packing, one entry: where our record drawing is the catalogue's record drawing ----
    # (all squares within SAME_TOL in position and angle mod 90, under the best of the 8 container symmetries), the
    # catalogue's dated timeline entry points at our certified file and the catalogue SVG is marked dup_of (Explore
    # and Compare skip it; a ?p= link still opens it).  The regularization level is relative to the lists' input
    # (the exact batch), not to the drawing a reader sees, so it can say "rattlers only" for an identical drawing.
    ndup = 0
    for n in range(1, NMAX + 1):
        r = idx['records'].get(str(n)) or {}
        key = r.get('best')
        if not key or not key.startswith('best-'):
            continue
        # the catalogue's record entry first, then any other catalogue drawing of n at the record side (the plain
        # grids, e.g. square-12.svg, are linked from Ellsworth's pages without an entry)
        side = float(idx['files'][key]['s'])
        cands = [r['svg']] if r.get('svg') in idx['files'] else []
        cands += sorted(k for k, f in idx['files'].items() if not f.get('origin') and f.get('n_parsed') == n
                        and k not in cands and f.get('s') and abs(float(f['s']) - side) < SAME_TOL)
        cat = None
        for c in cands:
            if idx['files'][c].get('n_parsed') != n:
                continue
            try:
                A = json.load(open(os.path.join(W, 'p', c.replace('.svg', '.json'))))
                B = json.load(open(os.path.join(W, 'p', key)))
            except (OSError, ValueError):
                continue
            if same_drawing(A, B):
                cat = c
                break
        if not cat:
            continue
        idx['files'][cat]['dup_of'] = key
        idx['files'][key]['same_as'] = cat
        for d in TL.get(str(n), []):
            if d['file'] == cat:
                d['catalogue_file'] = cat
                d['file'] = key
        ndup += 1
    print(f'{ndup} record drawings identical to the catalogue drawing: listed once')

    json.dump(idx, open(os.path.join(W, 'index.json'), 'w'), separators=(',', ':'), ensure_ascii=False)
    json.dump(TL, open(os.path.join(W, 'timeline.json'), 'w'), separators=(',', ':'), ensure_ascii=False)
    print(f"index: {len(idx['records'])} records, {len(idx['files'])} files "
          f"({sum(1 for v in idx['files'].values() if v.get('origin'))} from the lists); "
          f"timeline: {sum(len(v) for v in TL.values())} entries")


if __name__ == '__main__':
    main()
