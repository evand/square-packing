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
        row = stored(n, e['first'], float(e['s']))
        if not row or row[0] is None:
            print(f'n={n}: no stored packing for the raw record {e["first"]}', file=sys.stderr)
            continue
        p = st.get(row[0])
        key = raw_key[n] = f'rec-{n}.json'
        cert = os.path.join(LISTS, 'certs', 'raw', f'n-{n:03d}.cert.gz')     # search/regularize/register_cert.py
        files[key] = dict(n=n, s=dec(row[1]), sq=[[repr(float(x)), repr(float(y)), repr(float(t))] for x, y, t in p.sq],
                          origin='raw', event=e['first'], store_id=row[0],
                          cert=os.path.relpath(cert, ROOT) if os.path.exists(cert) else None,
                          credit={k: EV[e['first']][k] for k in ('who', 'date', 'label', 'url')})
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
            ev = EV[e['first']]
            r['newer'] = True
            r['credit'] = {'who': ev['who'], 'date': ev['date'], 'label': ev['label'], 'url': ev['url'],
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
        steps = e.get('steps') or [{'s': float(e['s']), 'event': e['first'], 'current': True}]
        for k, step in enumerate(steps):
            ev = EV[step['event']]
            key = (raw_key.get(n) or f'best-{n}.json') if step['current'] else step.get('file')
            if not key or key not in idx['files']:
                continue
            dt = ev_date(ev)
            who = ev['who']
            verb = 'Found'
            text = f"{verb} by {who}, posted {dt['text']} ({ev['label']})."
            L.append({'file': key, 'n': n, 's': float(idx['files'][key]['s']) if not step['current'] else float(e['s']),
                      'events': [{'verb': verb, 'who': who, 'date': dt, 'text': text}], 'prose': text, 's_tex': None,
                      'rigid': False, 'date': [dt['year'], dt['month']], 'is_record': bool(step['current']),
                      'origin': 'lists', 'url': ev['url'], 'label': ev['label']})

    json.dump(idx, open(os.path.join(W, 'index.json'), 'w'), separators=(',', ':'), ensure_ascii=False)
    json.dump(TL, open(os.path.join(W, 'timeline.json'), 'w'), separators=(',', ':'), ensure_ascii=False)
    print(f"index: {len(idx['records'])} records, {len(idx['files'])} files "
          f"({sum(1 for v in idx['files'].values() if v.get('origin'))} from the lists); "
          f"timeline: {sum(len(v) for v in TL.values())} entries")


if __name__ == '__main__':
    main()
