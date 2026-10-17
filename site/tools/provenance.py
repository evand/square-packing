#!/usr/bin/env python3
"""Who found the records that are newer than Ellsworth's catalogue, and when: site/data/records2026.json.

The site's record packings come from search/regularize/lists/best.json (regularized drawings of the best known packing
for every n <= 324).  For most n that packing is the catalogue's record and the catalogue's own wording credits it.  For
the rest (found since mid-September 2026) this script reads who posted it first:

  - the packing store (search/packer/runs/store.sqlite): jlevy/squares register rows and pending registration requests
    (open issues, synced by pk/watch.py), each with its side;
  - the register's case files (sparse clone search/packer/runs/jlevy, refreshed by `pk.py sync-register`): finder,
    refinements, source keys;
  - the dated public postings: jlevy/squares issues from the store's `issue` table (creation time, credit name; 10-17:
    pk/watch.py syncs every issue), and data/lineage.yaml `events` for postings that are not issues (repos the
    register's source keys point at).

A packing counts as "the same" as a posting if the sides agree to SAME_TOL and search/packer/analysis/same_minimum.py
has not found them to be distinct local minima (it exact-solves every two postings at one n within SAME_TOL).  The credited finder is the earliest posting of the record packing; later postings of the same
packing are listed as `also`.  Where a posting says it started from someone else's packing (LINEAGE, from its own
lineage section; data/lineage.yaml `lineage`), the root of that chain is credited as the finder and the later
postings as improvements only when same_minimum.py certifies the two are the same local minimum (a path between them
with the side never above the higher one's); otherwise the later poster is the finder, `based_on` the root.  The
register's found_by, when its own posting is a different minimum from the register's value, is `based_on` too.
Undecided cases default in two directions on purpose: an undecided lineage link counts as distinct (the poster is
credited), while two undecided postings within SAME_TOL count as one packing (the earlier is credited; n = 261).  A record
the lists could not certify, or posted after them, is drawn as posted (`raw`): those are the n that still need
regularizing (search/regularize/REGULARIZE.md).  For every n <= 324 (10-17; was n <= 100) every posting that improved on the best before
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
SAME_TOL = 1e-7     # 10-17: exactsolve's KKT point of a float-posted packing can sit ~2e-8 below it (#668's 227, 260)
HISTORY_MAX = 324     # 10-17 (Evan): every n's 2026 postings are Compare cards, not just n <= 100

ISSUE = 'https://github.com/jlevy/squares/issues/'
LINEAGE_FILE = os.path.join(SITE, 'data', 'lineage.yaml')
# 10-17 (Evan): a lineage link means "found from" (how the search started).  Its root is credited as the FINDER only
# when the two postings are the same local minimum (a refinement: descent or more exact computation), as decided by
# search/packer/analysis/same_minimum.py; otherwise the poster is the finder and the chain is `based_on`.
MINIMA_FILE = os.path.join(SITE, 'data', 'lineage_minima.json')


def same_minimum_links():
    """(SAME, DISTINCT): SAME = lineage links "posting|n" decided the same minimum; DISTINCT = {(n, frozenset({a, b}))}
    for every pair (lineage link or two close postings) decided distinct minima"""
    if not os.path.exists(MINIMA_FILE):
        return set(), set()
    P = json.load(open(MINIMA_FILE))['pairs']
    same, dist = set(), set()
    for k, v in P.items():
        f = k.split('|')
        a, b, n = (f[0], v['root'], int(f[1])) if len(f) == 2 else (f[0], f[1], int(f[2]))
        if v['verdict'] == 'same' and len(f) == 2:
            same.add(k)
        if v['verdict'] == 'distinct':
            dist.add((n, frozenset((a, b))))
    return same, dist


def load_events(db):
    """(EVENTS, ORDER, LINEAGE).  10-17: postings that are jlevy/squares issues come from the packing store's `issue`
    table (pk/watch.py syncs every issue: creation time, the author's credit name via pk/people.py); postings that are
    not issues and the lineage come from data/lineage.yaml.  ORDER = posting time (issues to the second)."""
    Y = yaml.safe_load(open(LINEAGE_FILE))
    ev, when = {}, {}
    for k, v in (Y.get('events') or {}).items():
        ev[k] = dict(v); when[k] = str(v['date']) + 'T00:00:00Z'
    for num, who, created, title in db.execute('SELECT number, credit, created, title FROM issue WHERE created IS NOT NULL'):
        k = f'#{num}'
        ev[k] = dict(who=who, date=created[:10], label=f'jlevy/squares#{num}', url=ISSUE + str(num))
        when[k] = created
    order = sorted(ev, key=lambda k: (when[k], k))
    lineage = {str(k): {int(n): str(e) for n, e in v.items()} for k, v in (Y.get('lineage') or {}).items()}
    return ev, order, lineage


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
    if 'ry-xu' in who or 'Ryan Xu' in who:
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
    global EVENTS, ORDER, LINEAGE
    EVENTS, ORDER, LINEAGE = load_events(db)
    SAME, DISTINCT = same_minimum_links()
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
        # every issue's claim at n, open or closed (10-17: closed = registered; before, only open issues were read
        # and a registered posting counted only if register_event knew its finder)
        for src, finder, ps in db.execute("select source, finder, s from frontier where n=? and source like 'pending:%'", (n,)):
            ev = src.split(':', 1)[1]
            if ev in EVENTS and ps < ref - 1e-9:
                posts.append((ps, ev))
        # the register's value, filed under its found_by posting (register_event), unless that posting's own claim at n
        # is a different minimum from the register's value (10-17: s(146)'s register value is #606's side, found_by
        # Ryan Xu: #432's claim is a distinct minimum, so the value goes to the posting that carries it)
        reg_ev = reg_from = None
        if reg_s is not None and reg_s < ref - 1e-9:
            ev = register_event(n, ru.get('found_by') or [])
            if ev is None:
                # 10-17: a finder register_event doesn't know (Ostashev, PingYou, Mishapolk, Fang ...): their earliest
                # posting at n, matched by name
                names = [f.split('(')[0].strip() for f in ru.get('found_by') or [] if f.strip()]
                mine = [p for p in posts if any(nm and nm in EVENTS[p[1]]['who'] for nm in names)]
                if mine:
                    ev = min(mine, key=lambda p: ORDER.index(p[1]))[1]
            own = [p for p in posts if p[1] == ev]
            carry = sorted((p for p in posts if abs(p[0] - reg_s) < 1e-10), key=lambda p: ORDER.index(p[1]))
            # found_by's own claim is not the register's packing: farther than SAME_TOL, or decided a distinct minimum
            away = own and all(abs(p[0] - reg_s) >= SAME_TOL or any((n, frozenset((p[1], c[1]))) in DISTINCT for c in carry)
                               for p in own)
            if carry and (ev is None or away):
                if away:
                    reg_from = ev                        # the register's found_by: what the carrier was based on
                ev = carry[0][1]
            elif away:
                reg_from, ev = ev, None                  # nobody we know posted the register's value (211): own event
            if not ev:
                # 10-17: a register record no posting we know carries (e.g. 211: Couzo's T-188 was a note inside #607):
                # its own dated event, credited as the register credits it
                ev = f'register-{n}'
                # found_by's own posting is another minimum (211): the value is the improvers' (Couzo's T-188)
                makers = (ru.get('improved_by') or []) if reg_from else (ru.get('found_by') or [])
                EVENTS[ev] = dict(who=', '.join(makers) or 'jlevy/squares register',
                                  date=str(ru.get('source_date') or case.get('source_reviewed') or ''),
                                  label=f'jlevy/squares register, s({n})',
                                  url=f'https://jlevy.github.io/squares/cases/{n}.html')
                ORDER.append(ev)
            if not any(p[1] == ev and abs(p[0] - reg_s) < 1e-10 for p in posts):
                posts.append((reg_s, ev))
            reg_ev = ev
        # this project's exact optimum of someone's packing is theirs; dedupe by (event, side)
        posts = sorted(set(posts), key=lambda p: ORDER.index(p[1]))
        # A posting below the list's packing: the list could not certify it (exactsolve fails; REGULARIZE.md, n = 261).
        # The site then draws that posting as posted ('raw') instead of a regularized copy.
        # The smallest posting is drawn (n = 261: #505's refinement certifies, #432's posting does not); the credit
        # still goes to the earliest posting of that arrangement (`first`, within SAME_TOL) and its lineage.
        raw = None
        if posts and min(p[0] for p in posts) < list_s - SAME_TOL:
            raw = min(posts, key=lambda p: (p[0], ORDER.index(p[1])))
            s = raw[0]
            print(f'n={n}: posted {s:.12f} ({raw[1]}) is below the list packing {list_s:.12f}: drawn as posted',
                  file=sys.stderr)
        if s > ref - 1e-9:
            continue                       # the catalogue's record (or the grid) is still the best
        # one packing = sides within SAME_TOL and not decided distinct minima (10-17: 146's #505 / #606 are 1.4e-8
        # apart); anchored at the posting nearest the drawn side
        def one(p, q):
            return abs(p[0] - q[0]) < SAME_TOL and (n, frozenset((p[1], q[1]))) not in DISTINCT
        anchor = raw or (min(posts, key=lambda p: (abs(p[0] - s), p[1].startswith('register-'), ORDER.index(p[1])))
                         if posts else None)
        same = [p for p in posts if anchor and abs(anchor[0] - s) < SAME_TOL and one(p, anchor)]
        if not same:
            print(f'n={n}: no posting matches the list side {s} (posts {posts})', file=sys.stderr)
        first = same[0][1] if same else None
        drawn = raw[1] if raw else first          # the posting whose packing the site draws
        # lineage: from the drawn posting back along links to the SAME minimum (refinements): the root is the finder
        chain, ev = [], drawn
        while ev and n in LINEAGE.get(ev, {}) and ev not in chain and f'{ev}|{n}' in SAME:
            chain.append(ev)
            ev = LINEAGE[ev][n]
        found = ev if chain else None             # the root posting, credited as the finder
        based, b = [], ev                          # then the links to other minima: what it was found from
        while b and n in LINEAGE.get(b, {}) and LINEAGE[b][n] not in based and LINEAGE[b][n] != drawn:
            b = LINEAGE[b][n]
            based.append(b)
        registered = reg_ev is not None and any(p[1] == reg_ev for p in same)
        if registered and reg_from and reg_from not in based and reg_from not in (first, found) and reg_from not in chain:
            based.append(reg_from)
        # "based on" someone's own earlier packing says nothing (179: itsnaka's #422 from itsnaka's #401)
        credited = EVENTS[found or first]['who'] if (found or first) in EVENTS else None
        based = [b for b in based if b in EVENTS and EVENTS[b]['who'] != credited]
        based = [b for i, b in enumerate(based) if all(EVENTS[c]['who'] != EVENTS[b]['who'] for c in based[:i])]   # one per person (180: #422, #401)
        shown = [EVENTS[x]['who'] for x in [found or first] + (chain if found else []) if x in EVENTS]
        e = {'s': best[n]['side'] if raw is None else repr(raw[0]), 'catalogue_s': cat, 'first': first, 'raw': raw is not None,
             'drawn': drawn, 'found_by': found, 'improved_in': chain if found else [], 'based_on': based,
             'also': [p[1] for i, p in enumerate(same) if i and EVENTS[p[1]]['who'] != EVENTS[first]['who'] and p[1] not in chain
                      and all(EVENTS[q[1]]['who'] != EVENTS[p[1]]['who'] for q in same[1:i])]   # one per person (108)
                     if first else [],
             'registered': registered,
             'register_s': ru.get('value'), 'register_verified': bool(reg_verified.get('value')) and
             reg_s is not None and abs(side(reg_verified['value']) - reg_s) < 1e-12,
             # the register's improved_by describes the register's packing: only when that is the one drawn
             # less anyone already credited on this line (the finder, the refinements in `improved_in`)
             'refined_by': [r for r in ru.get('improved_by') or []
                            if not any(r.split('(')[0].strip() in w for w in shown)] if registered else []}
        steps, cur = [], ref
        for ps, ev in posts:                   # in posting order: keep each that beats everything before it
            cur_pk = (ps, ev) in same
            if ps < cur - 1e-9 and not (cur_pk and ev != first):
                steps.append({'s': ps, 'event': ev, 'current': cur_pk})
                cur = ps
        for st in steps:                       # each step's own same-minimum root (credits.py credits the root)
            r, seen = st['event'], []
            while n in LINEAGE.get(r, {}) and r not in seen and f'{r}|{n}' in SAME:
                seen.append(r)
                r = LINEAGE[r][n]
            if seen:
                st['root'] = r
            bo, b = [], r                      # and what it was found from (other minima)
            while n in LINEAGE.get(b, {}) and LINEAGE[b][n] not in bo and LINEAGE[b][n] != st['event']:
                b = LINEAGE[b][n]
                bo.append(b)
            if bo:
                st['based_on'] = bo
        if n <= HISTORY_MAX:
            e['steps'] = steps                 # the Bounds page draws these (import_lists.py fetches their packings)
        elif len(steps) > 1:
            e['history'] = steps               # every n: the leaderboard's credit chain (site/tools/credits.py)
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
              + (f" drawn={e['drawn']}(raw)" if e['raw'] else '') + (f" found_by={e['found_by']} via {e['improved_in']}" if e['found_by'] else '')
              + f" also={e['also']} reg={'yes' if e['registered'] else 'no'}{'(verified)' if e['register_verified'] else ''}"
              + (f" steps={[(round(x['s'], 9), x['event']) for x in e['steps']]}" if 'steps' in e else ''))
    if not a.show:
        p = os.path.join(SITE, 'data', 'records2026.json')
        json.dump(out, open(p, 'w'), indent=1, ensure_ascii=False)
        print('wrote', p, len(out['n']), 'n')


if __name__ == '__main__':
    main()
