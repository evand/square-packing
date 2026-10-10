#!/usr/bin/env python3
"""Who found the records that are newer than Ellsworth's catalogue, and when: site/data/records2026.json.

The site's record packings come from search/regularize/lists/best.json (regularized drawings of the best known packing
for every n <= 324).  For most n that packing is the catalogue's record and the catalogue's own wording credits it.  For
the rest (found since mid-September 2026) this script reads who posted it first:

  - the packing store (search/packer/runs/store.sqlite): jlevy/squares register rows and pending registration requests
    (open issues, pk/pending.yaml), each with its side;
  - the register's case files (sparse clone search/packer/runs/jlevy, refreshed by `pk.py sync-register`): finder,
    refinements, source keys;
  - EVENTS below: the dated public postings (issue creation dates from the GitHub API, repo dates from the register's
    source keys).

A packing counts as "the same" as a posting if the sides agree to SAME_TOL (exact re-solves of one arrangement move the
side by < 1e-9).  The credited finder is the earliest posting of the record packing; later postings of the same
packing are listed as `also`.  For n <= 100 (the Bounds page's range) every posting that improved on the best before
it is kept as a history step.

  provenance.py            write site/data/records2026.json
  provenance.py --show     print a table, write nothing
"""
import argparse, glob, json, math, os, re, sqlite3, sys
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.abspath(os.path.join(HERE, '..'))
ROOT = os.path.abspath(os.path.join(SITE, '..'))
STORE = os.path.join(ROOT, 'search', 'packer', 'runs', 'store.sqlite')
REG = os.path.join(ROOT, 'search', 'packer', 'runs', 'jlevy')
BEST = os.path.join(ROOT, 'search', 'regularize', 'lists', 'best.json')
SAME_TOL = 1e-8
HISTORY_MAX = 100

ISSUE = 'https://github.com/jlevy/squares/issues/'
# Dated public postings.  Issue dates: GitHub createdAt (UTC), fetched 2026-10-10.
EVENTS = {
    'dewinter': dict(who='Joost de Winter', date='2026-09-16', label='JoostdeWinter/square-packing-211',
                     url='https://github.com/JoostdeWinter/square-packing-211'),
    'casson': dict(who='Griff Casson (griffcass)', date='2026-09-23', label='griffcass/square-packing',
                   url='https://github.com/griffcass/square-packing'),
    'couzo0927': dict(who='Francisco Couzo', date='2026-09-27', label='franciscouzo/square-packing',
                      url='https://github.com/franciscouzo/square-packing'),
    'couzo1003': dict(who='Francisco Couzo', date='2026-10-03', label='franciscouzo/square-packing (10-03 update)',
                      url='https://github.com/franciscouzo/square-packing'),
    '#399': dict(who='Evan Daniel (this project)', date='2026-10-07', label='jlevy/squares#399', url=ISSUE + '399'),
    '#401': dict(who='Nate Chaoweeraprasit (itsnaka, SQUISH)', date='2026-10-07', label='jlevy/squares#401', url=ISSUE + '401'),
    '#432': dict(who='ry-xu', date='2026-10-08', label='jlevy/squares#432', url=ISSUE + '432'),
    '#451': dict(who='Francisco Couzo', date='2026-10-08', label='jlevy/squares#451', url=ISSUE + '451'),
    '#465': dict(who='Evan Daniel (this project)', date='2026-10-09', label='jlevy/squares#465', url=ISSUE + '465'),
    '#470': dict(who='Mishapolk', date='2026-10-09', label='jlevy/squares#470', url=ISSUE + '470'),
    '#476': dict(who='Francisco Couzo', date='2026-10-09', label='jlevy/squares#476', url=ISSUE + '476'),
    '#481': dict(who='Nate Chaoweeraprasit (itsnaka, SQUISH)', date='2026-10-09', label='jlevy/squares#481', url=ISSUE + '481'),
    '#484': dict(who='Kevin Fang (TheSnakeFang)', date='2026-10-10', label='jlevy/squares#484', url=ISSUE + '484'),
    '#488': dict(who='Francisco Couzo', date='2026-10-10', label='jlevy/squares#488', url=ISSUE + '488'),
    '#489': dict(who='Evan Daniel (this project)', date='2026-10-10', label='jlevy/squares#489', url=ISSUE + '489'),
}
# Same-day order (issue creation times), so that "earliest" is well defined.
ORDER = ['dewinter', 'casson', 'couzo0927', 'couzo1003', '#399', '#401', '#432', '#451', '#465', '#470', '#476', '#481',
         '#484', '#488', '#489']
# Our own posts are not in the store's pending table (pk/pending.yaml lists other people's): sides from the issue titles.
OURS = {'#465': {132: 11.98709933224506, 155: 12.95249894401401},
        '#489': {132: 11.985680198846, 308: 17.998269879526}}
# Register finders -> the posting their registered packing came from.  Couzo's 10-03 update changed these n.
COUZO_1003 = {208, 209, 228, 263, 272, 303, 306}


def register_case(n):
    f = os.path.join(REG, 'packing', 'frontier', f'n-{n:03d}.md')
    return yaml.safe_load(open(f).read().split('---')[1])['packing']


def register_event(n, finders):
    who = ' '.join(finders)
    if 'Couzo' in who:
        return 'couzo1003' if n in COUZO_1003 else 'couzo0927'
    if 'Chaoweeraprasit' in who or 'itsnaka' in who:
        return '#401'
    if 'ry-xu' in who:
        return '#432'
    if 'de Winter' in who:
        return 'dewinter'
    if 'Evan Daniel' in who:
        return '#399'
    return None


def side(v):
    v = str(v)
    if '/' in v:
        a, b = v.split('/')
        return int(a) / int(b)
    return float(v)


def build():
    best = {b['n']: b for b in json.load(open(BEST))}
    idx = json.load(open(os.path.join(SITE, 'www', 'data', 'index.json')))
    db = sqlite3.connect(STORE)
    reg_ref = os.popen(f'git -C {REG} log -1 --format="%h %cs"').read().strip()
    out = {'note': 'generated by site/tools/provenance.py; see its docstring', 'register': f'jlevy/squares@{reg_ref}',
           'events': EVENTS, 'n': {}}
    for n in range(1, 325):
        s = list_s = float(best[n]['side'])
        rec = idx['records'].get(str(n))
        cat = float(rec['s_dec']) if rec and rec.get('s_dec') else None
        grid = math.ceil(math.sqrt(n) - 1e-12)
        ref = cat if cat is not None else grid
        case = register_case(n)
        ru = case.get('reported_upper_bound') or {}
        reg_s = side(ru['value']) if ru.get('value') else None
        reg_verified = case.get('verified_upper_bound') or {}
        # every posting below the catalogue: (side, event)
        posts = []
        if reg_s is not None and reg_s < ref - 1e-9:
            ev = register_event(n, ru.get('found_by') or [])
            if ev:
                posts.append((reg_s, ev))
        for src, finder, ps in db.execute('select source, finder, s from frontier where n=? and pending=1', (n,)):
            ev = src.split(':', 1)[1]
            if ev in EVENTS and ps < ref - 1e-9:
                posts.append((ps, ev))
        for ev, d in OURS.items():
            if n in d:
                posts.append((d[n], ev))
        # this project's exact optimum of someone's packing is theirs; dedupe by (event, side)
        posts = sorted(set(posts), key=lambda p: ORDER.index(p[1]))
        # A posting below the list's packing: the list could not certify it (exactsolve fails; REGULARIZE.md, n = 261).
        # The site then draws that posting as posted ('raw') instead of a regularized copy.
        raw = None
        if posts and min(p[0] for p in posts) < list_s - SAME_TOL:
            raw = min(posts, key=lambda p: (round(p[0], 8), ORDER.index(p[1])))
            s = raw[0]
            print(f'n={n}: posted {s:.12f} ({raw[1]}) is below the list packing {list_s:.12f}: drawn as posted',
                  file=sys.stderr)
        if s > ref - 1e-9:
            continue                       # the catalogue's record (or the grid) is still the best
        same = [p for p in posts if abs(p[0] - s) < SAME_TOL]
        if not same:
            print(f'n={n}: no posting matches the list side {s} (posts {posts})', file=sys.stderr)
        first = same[0][1] if same else None
        e = {'s': best[n]['side'] if raw is None else repr(raw[0]), 'catalogue_s': cat, 'first': first, 'raw': raw is not None,
             'also': [p[1] for p in same[1:] if EVENTS[p[1]]['who'] != EVENTS[first]['who']] if first else [],
             'registered': reg_s is not None and abs(reg_s - s) < SAME_TOL,
             'register_s': ru.get('value'), 'register_verified': bool(reg_verified.get('value')) and
             reg_s is not None and abs(side(reg_verified['value']) - reg_s) < 1e-12,
             'refined_by': ru.get('improved_by') or []}
        if n <= HISTORY_MAX:
            steps, cur = [], ref
            for ps, ev in posts:               # in posting order: keep each that beats everything before it
                if ps < cur - 1e-9 and not (abs(ps - s) < SAME_TOL and ev != first):
                    steps.append({'s': ps, 'event': ev, 'current': abs(ps - s) < SAME_TOL})
                    cur = ps
            e['steps'] = steps
        out['n'][str(n)] = e
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--show', action='store_true')
    a = ap.parse_args()
    out = build()
    for n, e in out['n'].items():
        ev = EVENTS.get(e['first'], {})
        print(f"{n:>4} {float(e['s']):.10f} {ev.get('who', '?'):<40} {ev.get('date', ''):<11} {e['first'] or '?':<10}"
              f" also={e['also']} reg={'yes' if e['registered'] else 'no'}{'(verified)' if e['register_verified'] else ''}"
              + (f" steps={[(round(x['s'], 9), x['event']) for x in e['steps']]}" if 'steps' in e else ''))
    if not a.show:
        p = os.path.join(SITE, 'data', 'records2026.json')
        json.dump(out, open(p, 'w'), indent=1, ensure_ascii=False)
        print('wrote', p, len(out['n']), 'n')


if __name__ == '__main__':
    main()
