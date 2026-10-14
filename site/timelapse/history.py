#!/usr/bin/env python3
"""Dated record history for n = 1..324 -> out/history.json (input of reg.py / render.py in this directory).

Per n, a list of states (time, side, who, drawing), strictly decreasing in side:
  1. grid (s = ceil(sqrt n)), from the start;
  2. catalogue era (to 2026-10-05): site/data/credits.json upper-bound holders (Ellsworth's catalogue prose and
     site/data/records2026.json postings), dated (year, month) -> the 15th, year only -> July 1, `day` when given;
     drawing = site/www/data/hist.json (catalogue SVG parses, or superseded 2026 records from the store);
  3. the 10-05 register (search/exact/batch/inputs/n-N.txt) where it is below the catalogue chain (a register record
     with no dated posting in credits): dated by the register case file's source_date (in the sparse register clone),
     drawing = the store packing at that side;
  4. record events 2026-10-07 .. 10-13 00:42 UTC: search/packer/bench/events.json (pk/eventbench.py; issue post times),
     drawing = the store packing at the event's side (eventbench.store_match);
  5. later public postings in the store (register / pending:#NNN sightings) below everything above, dated by the
     pending issue's createdAt (issue_dates.json, fetched read-only with `gh issue view`), else the sighting time.
Everything else (precision-only claims, n > 324, our unposted results) is left out.

  history.py [--fetch-issue-dates]
"""
import argparse, json, math, os, re, subprocess, sys
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SITE = os.path.join(ROOT, 'site')
PKDIR = os.path.join(ROOT, 'search', 'packer')
sys.path.insert(0, PKDIR)
from pk import eventbench as eb  # noqa: E402

OUT = os.path.join(HERE, 'out')
NMAX = 324
REG_T = '2026-10-05T00:00Z'
ISSUE_DATES = os.path.join(HERE, 'issue_dates.json')
ALIAS = {'dostashev': 'Daniil Ostashev', 'ry-xu': 'Ryan Xu'}
DROP = 1e-9          # a new state must be below the current one by more than this


def iso(y, m=None, d=None):
    if m is None:
        return f'{y:04d}-07-01T00:00Z'
    return f'{y:04d}-{m:02d}-{d or 15:02d}T00:00Z'


def drawing_side(path):
    return float(json.load(open(path))['s'])


def catalogue(n, C, H):
    """credits.json chain (non-grid holders) as states"""
    out = []
    for i, u in enumerate(C[str(n)]['upper']):
        if u['kind'] == 'grid':
            continue
        f = H[str(n)][i]
        path = os.path.join(SITE, 'www', 'data', 'p', f) if f else None
        if u.get('day'):
            t = u['day'] + 'T00:00Z'
        else:
            y, m = u['date']
            t = iso(y, m)
        s = u.get('s')
        draw = dict(kind='file', path=os.path.relpath(path, ROOT)) if path else None
        if s is None and path:
            later = [H[str(n)][j] for j in range(i + 1, len(H[str(n)]))]
            if f in later:
                # the catalogue redrew this entry in place: the drawing is a later holder's.  Use the largest-side
                # catalogue revision in the store below the previous drawing's side, if any; else skip the holder.
                s, draw = revision(n, i, C, H)
            else:
                s = drawing_side(path)
        who = ', '.join(C_people.get(w, w) for w in u['who'])
        out.append(dict(t=t, s=s, who=who, via=u.get('event') or u.get('file') or u['kind'],
                        text=u.get('text', ''), src='catalogue', draw=draw))
    return out


def revision(n, i, C, H):
    """(side, drawing) for holder i of n whose own drawing was replaced in place: the largest-side kingbird
    catalogue revision (store sighting kind catalogue-revision) strictly between the next distinct holder's side and
    the side before holder i; (None, None) if there is none"""
    up = C[str(n)]['upper']
    hi = math.ceil(math.sqrt(n) - 1e-12)
    for j in range(i - 1, -1, -1):
        if up[j].get('s'):
            hi = up[j]['s']; break
    lo = min([u['s'] for u in up[i + 1:] if u.get('s')] or [0.0])
    r = eb.db().execute("SELECT p.id, p.s FROM packing p JOIN sighting g ON g.packing = p.id WHERE p.n = ? AND "
                        "g.kind = 'catalogue-revision' AND p.s < ? AND p.s > ? ORDER BY p.s DESC LIMIT 1",
                        (n, hi - DROP, lo + DROP)).fetchone()
    return (float(r[1]), dict(kind='store', id=int(r[0]))) if r else (None, None)


def register_meta(n):
    """register case file (sparse clone): finder list and the dated claims of its priority notes"""
    p = os.path.join(PKDIR, 'runs', 'jlevy', 'packing', 'frontier', f'n-{n:03d}.md')
    if not os.path.exists(p):
        return {}
    txt = open(p).read()
    fb = re.search(r'found_by:\n((?:\s+- .*\n)+)', txt)
    claims = []
    for m in re.finditer(r"- claim: (.*)\n(?: {6}.*\n)*\s+claimed_by:\n((?:\s+- .*\n)+)\s+published: (.*)\n", txt):
        v = re.search(r's\(\d+\) <= ([\d.]+)', m[1])
        d = re.search(r'(\d{4}-\d\d-\d\d)(?:T(\d\d:\d\d))?', m[3])
        if v and d and 'exact optim' not in m[1]:        # skip "Prior reported ... from [evand exact optima]" refinements
            claims.append(dict(s=float(v[1]), t=f'{d[1]}T{d[2] or "00:00"}Z',
                               who=', '.join(x.strip()[2:] for x in m[2].splitlines())))
    return dict(found_by=[x.strip()[2:] for x in fb[1].splitlines()] if fb else [], claims=claims)


def store_packing(n, s, tol=1e-7):
    p = eb.store_match(n, s, tol)
    return None if p is None else int(p.meta['id'])


def issue_dates(fetch=False, issues=()):
    D = json.load(open(ISSUE_DATES)) if os.path.exists(ISSUE_DATES) else {}
    if fetch:
        for i in issues:
            if str(i) in D:
                continue
            r = subprocess.run(['gh', 'issue', 'view', str(i), '-R', 'jlevy/squares', '--json', 'createdAt', '-q', '.createdAt'],
                               capture_output=True, text=True)
            if r.returncode == 0 and r.stdout.strip():
                D[str(i)] = r.stdout.strip().replace(':00Z', 'Z')[:16] + 'Z'
        json.dump(D, open(ISSUE_DATES, 'w'), indent=1, sort_keys=True)
    return D


def later_public(n, below, fetch):
    """public store packings at n below `below`: [(t, s, id, who, via)] earliest posting per packing"""
    rows = eb.db().execute(
        "SELECT p.id, p.s, g.source, g.finder, g.added FROM packing p JOIN sighting g ON g.packing = p.id WHERE p.n = ? "
        "AND p.s < ? AND p.status != 'infeasible' AND (p.pen IS NULL OR p.pen <= 1e-9) AND (g.source LIKE 'register%' "
        "OR g.source LIKE 'pending:%') ORDER BY p.s", (n, below - DROP)).fetchall()
    iss = [int(r['source'].split('#')[1]) for r in rows if r['source'].startswith('pending:#')]
    D = issue_dates(fetch, iss)
    best = {}
    for r in rows:
        if r['source'].startswith('pending:#'):
            i = r['source'].split('#')[1]
            t, via = D.get(i), '#' + i
            approx = t is None
            if t is None:      # sighting times are local; good to a day
                t = r['added'][:16] + 'Z'
        else:
            t, via, approx = r['added'][:16] + 'Z', 'register', True
        k = r['id']
        if k not in best or t < best[k][0]:
            best[k] = (t, r['s'], k, r['finder'], via, approx)
    return sorted(best.values())


def envelope(hist, notes):
    """s(n) <= s(n+1): drop a square.  Where n+1's best known packing (own or itself derived) is below n's own record,
    n's best known packing is n+1's minus one square: add 'derived' states (drawing = the base state's drawing with
    squares removed, chosen in reg.py).  Processed per change time, from the top of each row down."""
    times = sorted({st['t'] for v in hist.values() for st in v if st['t']})
    own = {n: list(v) for n, v in hist.items()}
    eff = {}                                   # n -> (s, base n, base state index) currently shown
    out = {n: [v[0]] for n, v in own.items()}
    cnt = 0
    for T in [None] + times:
        for n in range(NMAX, 0, -1):
            i = max(j for j, st in enumerate(own[n]) if st['t'] is None or (T is not None and st['t'] <= T))
            cand = (own[n][i]['s'], n, i)
            if n < NMAX and eff[n + 1][0] < cand[0] - DROP:
                cand = eff[n + 1]
            if eff.get(n) != cand:
                eff[n] = cand
                if T is None:
                    continue
                s, m, j = cand
                if m == n:
                    out[n].append(own[n][j])
                else:
                    b = own[m][j]
                    out[n].append(dict(t=T, s=s, who=b['who'], via=f'from s({m})', text=f'{m - n} square(s) removed from '
                                       f'the s({m}) record ({b["via"]})', src='derived', base_n=m, base_t=b['t'],
                                       draw=dict(kind='derived', n=m, remove=m - n, base=b['draw'])))
                    cnt += 1
    # own states are shown at their own time or never (eff(n+1) only decreases), so out[n] is in time order
    hist.update(out)
    notes.append(f'envelope: {cnt} derived states (n+1 minus a square below n\'s own record) in '
                 f'{sum(any(st["src"] == "derived" for st in v) for v in hist.values())} n')


def build(fetch=False):
    global C_people, BATCH_BY
    BATCH_BY = {r['n']: r.get('by') for r in json.load(open(os.path.join(eb.REG, '..', 'results.json')))}
    CR = json.load(open(os.path.join(SITE, 'data', 'credits.json')))
    C, C_people = CR['n'], {k: v.get('name', k) for k, v in CR['people'].items()}
    H = json.load(open(os.path.join(SITE, 'www', 'data', 'hist.json')))
    EV = [e for e in eb.events() if e['n'] <= NMAX]
    hist, notes = {}, []
    for n in range(1, NMAX + 1):
        k = math.ceil(math.sqrt(n) - 1e-12)
        S = [dict(t=None, s=float(k), who='', via='grid', text='grid packing', src='grid', draw=dict(kind='grid'))]
        cur = float(k)

        def push(st):
            nonlocal cur
            if st['s'] is not None and st['s'] < cur - DROP:
                S.append(st); cur = st['s']
                return True
            return False
        for st in sorted(catalogue(n, C, H), key=lambda x: x['t']):
            if st['t'] >= REG_T:
                continue
            if st['s'] is None:
                notes.append(f'n={n}: catalogue holder {st["who"]} {st["t"][:10]} has no side and no drawing of its own '
                             '(entry redrawn in place, no revision in the store): skipped')
                continue
            if not push(st):
                notes.append(f'n={n}: catalogue holder {st["who"]} {st["t"][:10]} s={st["s"]:.10f} not below {cur:.10f}: skipped')
        reg = eb.base(n)[0] if n > 1 else 1.0
        if reg < cur - 1e-8:
            # a register record with no dated posting in credits: date it by the register's priority-note claim at
            # that side, else by the finders' known posting (Couzo's repo report 2026-09-27); drawing = the frozen
            # 10-05 batch input (search/exact/batch/inputs, read by pk.packing.read)
            m = register_meta(n)
            by = BATCH_BY.get(n) or []          # finders as of the 10-05 batch (the case file's are today's)
            hit = [c for c in m.get('claims', []) if abs(c['s'] - reg) < 1e-8
                   and any(b.split()[-1] in c['who'] for b in by)]      # the finder's own claim, not a refinement
            if hit:
                t, how = min(c['t'] for c in hit), 'register priority note'
            elif any('Couzo' in b for b in by):
                t, how = '2026-09-27T00:00Z', 'assumed Couzo 09-27 report'
            else:
                t, how = '2026-10-04T00:00Z', 'undated: placed just before the 10-05 register'
            notes.append(f'n={n}: 10-05 register {reg:.10f} ({", ".join(by)}) below the catalogue chain {cur:.10f}; '
                         f'dated {t[:16]} ({how})')
            push(dict(t=t, s=reg, who=', '.join(by), via='register 10-05', text=how, src='register',
                      approx_date=not hit, draw=dict(kind='file', path=os.path.relpath(f'{eb.REG}/n-{n}.txt', ROOT))))
        elif reg > cur + 1e-8:
            notes.append(f'n={n}: catalogue chain {cur:.10f} below the 10-05 register {reg:.10f} (chain kept)')
        for e in sorted((e for e in EV if e['n'] == n), key=lambda e: e['date']):
            pid = store_packing(n, e['new'])
            ok = push(dict(t=e['date'], s=e['new'], who=e['author'], via='#' + e['issue'], text=e['family'] + ': ' + e['parents'],
                           src='event', draw=dict(kind='store', id=pid) if pid else None))
            if not ok:
                notes.append(f'n={n}: event {e["issue"]} {e["new"]} not below {cur}: skipped')
        for t, s, pid, who, via, approx in later_public(n, cur, fetch):
            if push(dict(t=t, s=s, who=ALIAS.get(who, who) or '', via=via, text='', src='store', approx_date=approx,
                         draw=dict(kind='store', id=pid))):
                notes.append(f'n={n}: post-events record {s:.10f} {via} {t}')
        hist[n] = S
    envelope(hist, notes)
    os.makedirs(OUT, exist_ok=True)
    json.dump(dict(built=datetime.now(timezone.utc).isoformat(timespec='seconds'), notes=notes,
                   hist={str(n): v for n, v in hist.items()}), open(os.path.join(OUT, 'history.json'), 'w'), indent=0)
    return hist, notes


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fetch-issue-dates', action='store_true')
    a = ap.parse_args()
    hist, notes = build(a.fetch_issue_dates)
    print('\n'.join(notes))
    print('states', sum(len(v) for v in hist.values()), 'non-grid', sum(len(v) - 1 for v in hist.values()))
