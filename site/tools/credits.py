#!/usr/bin/env python3
"""Who holds, and who held, each bound: site/data/credits.json (the leaderboard's data; same sources as the Bounds page).

Upper bounds (packings), per n <= 324:
  - Ellsworth's catalogue prose (data/records_raw.json): the record's own text (main page) and, for the 33 n on the
    "compared" page, the earlier records with their side lengths.  Each sentence is parsed into a credit:
      found / improved          a new record: the holder changes at that date
      refined / optimized / fixed / refound   credit for the refiner, the holder stays
      proved                    the bound is optimal (a lower-bound credit, below)
      construction              a sentence like 'Extends the s(65) found by Frits Göbel in early 1979' with no finder
                                of its own: credited as a construction on that person's packing, dated by it
      alternative / nonrecord   not a record (same side or worse): kept out of the holder chain
      improvement               a 2026 posting that improved the root of its own lineage chain (records2026 found_by):
                                credited, but the root keeps the record, as on the rest of the site
    plus `based_on`: the packings a credit says it started from ('based on / extending / adding an "L" to the s(N)
    found by ...').
  - data/records2026.json (tools/provenance.py): postings newer than the catalogue, with their lineage.
  - data/credit_overrides.yaml: hand-written credits for sentences the parser gets wrong (keyed by file).
  - n with no catalogue entry and side ceil(sqrt n): the grid ("trivial").
Lower bounds: data/lower_bounds.json histories (n <= 324) with its `sources` table, names mapped through people.yaml.

Holders are the chain of found/improved credits in date order; a holder holds from its date to the next one's.
Dates are (year, month); a missing month counts as mid-year (early = Feb, late = Nov, as in timeline.py).

  credits.py            write data/credits.json and print a summary
  credits.py --show N   print the parsed chain for n = N (several allowed)
  credits.py --review   write notes/credits-review.md: every n's chain, for reading through
"""
import argparse, datetime, json, math, os, re, sys
from collections import Counter, defaultdict
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.abspath(os.path.join(HERE, '..'))
D = os.path.join(SITE, 'data')
NMAX = 324
TODAY = datetime.date.today()  # record-years run to this date; written to credits.json as `generated`
MONTHS = {m: i + 1 for i, m in enumerate(['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August',
                                          'September', 'October', 'November', 'December'])}
MON = '|'.join(MONTHS)
DATE_RE = r'(?:(?:early|late|mid)[- ])?(?:(?:%s)(?:/(?:%s))? )?\d{4}' % (MON, MON)

# ---------------------------------------------------------------- people
PEOPLE = yaml.safe_load(open(os.path.join(D, 'people.yaml')))
ALIAS = {}
for pid, p in PEOPLE.items():
    for a in p.get('aliases', []) + [p['name']]:
        ALIAS[a] = pid
SOURCE_PEOPLE = defaultdict(list)  # lower_bounds.json source key -> [ids]
for pid, p in PEOPLE.items():
    for k in p.get('sources', []):
        SOURCE_PEOPLE[k].append(pid)
NAME_RE = '|'.join(re.escape(a) for a in sorted(ALIAS, key=len, reverse=True))
NAMES_RE = r'(?:%s)(?:(?:, | and |, and )(?:%s))*' % (NAME_RE, NAME_RE)


def people_of(text):
    """'A, B, and C' -> [ids]."""
    out = []
    for m in re.finditer(NAME_RE, text):
        pid = ALIAS[m.group(0)]
        if pid not in out:
            out.append(pid)
    return out


def parse_date(t):
    if not t:
        return None
    m = re.search(DATE_RE, t)
    if not m:
        return None
    t = m.group(0)
    y = int(re.search(r'\d{4}', t).group(0))
    mon = None
    for name, i in MONTHS.items():
        if name in t:
            mon = i  # a range 'November/December' keeps the last
    if mon is None:
        mon = 2 if t.startswith('early') else 11 if t.startswith('late') else 6 if t.startswith('mid') else None
    return [y, mon]


def dkey(d):
    return (d[0], d[1] or 6.5) if d else (9999, 0)


def years(a, b):
    return (dkey(b)[0] + (dkey(b)[1] - 0.5) / 12) - (dkey(a)[0] + (dkey(a)[1] - 0.5) / 12)


# ---------------------------------------------------------------- catalogue prose
VERBS = [  # (regex at sentence start, kind)
    (r'(?:Rigid\. )?Found first', 'found'), (r'Originally found', 'found'), (r'Found and improved', 'found'),
    (r'Apparently found', 'found'), (r'Found', 'found'),
    (r'Improved and optimized', 'improved'), (r'Improved again', 'improved'), (r'Improved', 'improved'),
    (r'Refound and refined', 'refined'), (r'Refound', 'refined'), (r'Refined', 'refined'),
    (r'Optimized', 'refined'), (r'Fixed', 'refined'), (r'Drafted', 'refined'), (r'Polynomial root solution found', 'refined'),
    (r'Side length found', 'found'), (r'Proved', 'proved'),
    (r'(?:Pictured |Semi-rigid |Rigid |Symmetric |Mirror-symmetric |Rotationally symmetric )?(?:[Aa]lternative|version|[Rr]earrangement)[\w ,-]*? (?:found|arranged)', 'alternative'),
    (r'Rotationally symmetric rearrangement[\w ,-]*?', 'alternative'),
    (r'Experiment', 'alternative'),
    (r'Improved independently', 'improved'),
    (r'Invalid; incorrectly calculated', 'invalid'),
    (r'This alternative, converting [^,]*?, found', 'alternative'),
    (r'Alternative constructed', 'alternative'),
]
VERB_RE = [(re.compile(r'^(?:Rigid\. )?(%s)(?: \(.*?\))? by (%s)(?: \(.*?\))?(?: (?:in|between) (%s))?' % (v, NAMES_RE, DATE_RE)), k)
           for v, k in VERBS]
# a reference to another packing: '... the $s(37)$ [that ...] found [and improved] by X [in D]'
REF_RE = re.compile(r'\$s\((\d+)\)\$(?:[^$.]{0,60}?)(?:found|improved|found/improved|found and improved|augmented|converted)'
                    r'(?: and improved)? by (%s)(?: in (%s))?' % (NAMES_RE, DATE_RE))
SELF_REF_RE = re.compile(r'\$s\((\d+)\)\$[^$.]{0,40}? he (?:found|improved|found/improved)')
PACKING_RE = re.compile(r'(?:the |)packing found by (%s)(?: in (%s))?' % (NAMES_RE, DATE_RE))
PATTERN_RE = re.compile(r'pattern found by (%s)(?: in (%s))?' % (NAMES_RE, DATE_RE))
LINEAGE_START = re.compile(r'^(?:Extends|Extended|Adds|Now adds|Continues|Based on|Combines|Unextends|Improves upon|'
                           r'Improved upon|Fits|Similar|Resembles|Uses|Is an alternative|Can also|Originally added|'
                           r'With these improvements|Unimproved form|Converted)')
NONRECORD_RE = re.compile(r"Didn't set|Did not set|didn't set a record|Not optimal|Invalid")
SPLIT_RE = re.compile(r'(?<=[.!?])(?<! [A-Z]\.)(?<!\b[A-Z]\.[A-Z]\.)\s+(?=[A-Z])')


def sentences(prose):
    t = (prose or '').replace('\xa0', ' ').replace('\n', ' ')
    t = re.sub(r'\s+', ' ', t).strip()
    return [s for s in SPLIT_RE.split(t) if s]


def refs_in(text, n, main_who):
    out = []
    for m in REF_RE.finditer(text):
        out.append({'n': int(m.group(1)), 'who': people_of(m.group(2)), 'date': parse_date(m.group(3))})
    for m in SELF_REF_RE.finditer(text):
        out.append({'n': int(m.group(1)), 'who': list(main_who), 'date': None})
    for m in PACKING_RE.finditer(text):
        out.append({'n': n, 'who': people_of(m.group(1)), 'date': parse_date(m.group(2)), 'packing': True})
    for m in PATTERN_RE.finditer(text):
        out.append({'n': None, 'who': people_of(m.group(1)), 'date': parse_date(m.group(2)), 'pattern': True})
    return out


def parse_prose(prose, n):
    """-> (credits, flags).  credits: [{kind, who, date, based_on, text}]; flags: {'nonrecord', 'invalid', ...}"""
    credits, flags, unparsed = [], set(), []
    for s in sentences(prose):
        if s.startswith('Alternative constructed'):
            flags.add('alternative')
            continue
        if NONRECORD_RE.search(s):
            flags.add('invalid' if s.startswith('Invalid') else 'nonrecord')
        m = None
        for rx, kind in VERB_RE:
            m = rx.match(s)
            if m:
                break
        if m:
            who = people_of(m.group(2))
            rest = s[m.end():]
            date = parse_date(m.group(3))
            if date is None:  # 'Found by X (by adding ...) in 2002 or slightly earlier'
                mm = re.search(r' in (%s)' % DATE_RE, rest)
                date = parse_date(mm.group(1)) if mm else None
            based = refs_in(rest, n, who)
            # 'Improves upon / Improved upon the $s(n)$ found by Y': the previous holder, same n
            c = {'kind': kind, 'who': who, 'date': date, 'based_on': based, 'text': s}
            # 'Found by X in 2024. Improved by Y in 2025, ...' : a second clause 'and improved by Z in D'
            for m2 in re.finditer(r'(?:,| and) (improved|fixed|optimized|refined)(?: again)? by (%s) in (%s)' % (NAMES_RE, DATE_RE), rest):
                k2 = 'improved' if m2.group(1) == 'improved' else 'refined'
                credits.append(c)
                c = {'kind': k2, 'who': people_of(m2.group(2)), 'date': parse_date(m2.group(3)), 'based_on': [], 'text': s}
            credits.append(c)
            if kind in ('alternative', 'invalid'):
                flags.add(kind)
            continue
        if LINEAGE_START.match(s):
            refs = refs_in(s, n, [])
            if s.startswith(('Improves upon', 'Improved upon')):
                for r in refs:
                    if r['n'] == n:
                        credits.append({'kind': 'previous', 'who': r['who'], 'date': r['date'], 'based_on': [], 'text': s})
                continue
            if refs:
                own = re.match(r'Converted in (%s)' % DATE_RE, s)
                credits.append({'kind': 'construction', 'who': [], 'date': parse_date(own.group(1)) if own else None,
                                'own_date': bool(own), 'based_on': refs, 'text': s})
                continue
        if re.search(r'\bby (?:%s)' % NAME_RE, s) and not re.search(r'pending\.$|^Didn\'t set an overall record, but proved', s) \
                and not s.startswith(('Note that', 'It was foreshadowed', 'It is the first',
                                                                         'This is', "It's", 'The ', 'See also')):
            unparsed.append(s)
    # a lineage sentence in a file that names its own finder is that finder's lineage, not a separate construction
    main = [c for c in credits if c['kind'] in ('found', 'improved', 'proved')]
    if main:
        for c in [c for c in credits if c['kind'] == 'construction']:
            main[0]['based_on'] += c['based_on']
            credits.remove(c)
    return credits, flags, unparsed


# ---------------------------------------------------------------- assemble
def load():
    R = json.load(open(os.path.join(D, 'records_raw.json')))
    REC = json.load(open(os.path.join(D, 'records2026.json')))
    LB = json.load(open(os.path.join(D, 'lower_bounds.json')))
    OV = {}
    p = os.path.join(D, 'credit_overrides.yaml')
    if os.path.exists(p):
        OV = yaml.safe_load(open(p)) or {}
    return R, REC, LB, OV


def catalogue_files(R):
    """n -> [{file, s, prose, page, record}] (main-page record first, then the compared page, oldest first by s)."""
    out = defaultdict(list)
    seen = set()
    for e in R['squares_in_squares.html']:
        if e.get('hidden') or not e.get('n') or e['svg'] in seen:
            continue
        seen.add(e['svg'])
        for n in e.get('ns') or [e['n']]:  # one picture can be the record for several n (s(147): s(148) less a square)
            out[n].append({'file': e['svg'], 's': float(e['s_dec']) if e.get('s_dec') else None, 'prose': e.get('prose'),
                           'page': 'main', 'record': True, 'pictured_n': e['n']})
    for e in R['squares_in_squares__compared.html']:
        if not e.get('n') or e['svg'] in seen:
            continue
        seen.add(e['svg'])
        out[e['n']].append({'file': e['svg'], 's': float(e['s_dec']) if e.get('s_dec') else None, 'prose': e.get('prose'),
                            'page': 'compared', 'record': False})
    return out


def upper_chain(n, files, rec2026, events, OV, problems):
    """-> {holders: [...], credits: [...], files: [...]}"""
    all_credits = []
    for f in files:
        if f['file'] in OV:
            cr, flags, unp = OV[f['file']].get('credits', []), set(OV[f['file']].get('flags', [])), []
            for c in cr:
                c['who'] = [ALIAS.get(w, w) for w in c.get('who', [])]
                c.setdefault('based_on', [])
                c.setdefault('text', '(override)')
                c['date'] = parse_date(c['date']) if isinstance(c.get('date'), str) else c.get('date')
        else:
            cr, flags, unp = parse_prose(f['prose'], n)
        for s in unp:
            problems.append((n, f['file'], 'unparsed', s))
        last = max((i for i, c in enumerate(cr) if c['kind'] in ('found', 'improved', 'construction')), default=None)
        for i, c in enumerate(cr):
            c['file'] = f['file']
            c['s'] = f['s'] if i == last else None
            c['page'] = f['page']
            if flags & {'alternative', 'nonrecord', 'invalid'} and c['kind'] in ('found', 'improved'):
                c['kind'] = 'alternative' if 'alternative' in flags else 'nonrecord'
            for w in c['who']:
                if w not in PEOPLE:
                    problems.append((n, f['file'], 'unknown person', w))
            all_credits.append(c)
    # 'construction' credits take their date and origin from what they build on
    for c in all_credits:
        if c['kind'] == 'construction' and not c.get('own_date'):
            dated = [r for r in c['based_on'] if r.get('date')]
            c['date'] = max((r['date'] for r in dated), key=dkey) if dated else None
            c['origin'] = sorted({w for r in c['based_on'] for w in r['who']})
    # 2026 postings
    r = rec2026.get(str(n))
    if r:
        first = r['first']
        root = r.get('found_by') or first
        chain = [root] + [e for e in r.get('improved_in', []) if e != root]
        for i, ev in enumerate(chain):
            E = events[ev]
            # the root of a same-minimum chain holds the record (10-17: provenance.py credits a root only when the
            # postings are the same local minimum); later postings of that minimum are refinements.  A record found
            # from another minimum carries it as based_on ('based on a packing by ...')
            bo = [{'who': people_of(events[b]['who']), 'date': [int(x) for x in events[b]['date'].split('-')[:2]]}
                  for b in r.get('based_on', [])] if i == 0 else []
            all_credits.append({'kind': 'found' if i == 0 else 'refined', 'who': people_of(E['who']),
                                'date': [int(x) for x in E['date'].split('-')[:2]], 'day': E['date'],
                                'based_on': bo, 'text': E['label'], 'url': E['url'], 'event': ev,
                                # the record's side belongs to the last posting of the chain
                                's': float(r['s']) if i == len(chain) - 1 else None, 'page': '2026'})
        for name in r.get('refined_by', []):  # register `improved_by`: names, undated
            all_credits.append({'kind': 'refined', 'who': people_of(name), 'date': None, 'based_on': [],
                                'text': 'refinement in the jlevy/squares register', 'page': '2026'})
        for step in r.get('steps') or r.get('history') or []:  # earlier 2026 records superseded since
            ev = step['event']
            if ev in chain:
                continue
            E = events[ev]
            c = {'kind': 'found', 'who': people_of(E['who']), 'date': [int(x) for x in E['date'].split('-')[:2]],
                 'day': E['date'], 'based_on': [], 'text': E['label'], 'url': E['url'], 'event': ev, 's': step['s'], 'page': '2026'}
            if step.get('root'):  # the same minimum as the root's: the root holds it, the poster refined it
                R0 = events[step['root']]
                all_credits.append({'kind': 'refined', 'who': c['who'], 'date': c['date'], 'day': E['date'], 'based_on': [],
                                    'text': E['label'], 'url': E['url'], 'event': ev, 'page': '2026'})
                c.update(who=people_of(R0['who']), text=R0['label'], url=R0['url'])
            c['based_on'] = [{'who': people_of(events[b]['who']), 'date': [int(x) for x in events[b]['date'].split('-')[:2]]}
                             for b in step.get('based_on', [])]
            all_credits.append(c)
    # holders: found / improved / construction / previous, in date order, deduplicated on (who, date)
    hs = []
    for c in all_credits:
        if c['kind'] not in ('found', 'improved', 'construction', 'previous'):
            continue
        who = c['who'] if c['kind'] != 'construction' else []
        key = (tuple(who), tuple(c['date'] or ()), c.get('day'))
        if c['kind'] == 'construction' and any(h['kind'] != 'construction' for h in hs if h['date'] == c['date']):
            continue
        dup = [h for h in hs if (tuple(h['who']), tuple(h['date'] or ()), h.get('day')) == key]
        if dup:
            if c.get('s') and not dup[0].get('s'):
                dup[0].update({k: c[k] for k in ('s', 'file') if k in c})
            continue
        hs.append(c)
    hs.sort(key=lambda c: (dkey(c['date']), c.get('day') or '', -(c.get('s') or 0)))
    # a construction that is not the latest credit for its file does not hold on its own
    holders = []
    for c in hs:
        if holders and c['kind'] == 'construction' and holders[-1]['kind'] != 'construction':
            continue
        holders.append({k: c[k] for k in ('kind', 'who', 'date', 'based_on', 'text', 'file', 's', 'url', 'event', 'day', 'origin')
                        if k in c})
    merged = []  # one entry per run of 2026 records by the same people (someone improving their own record)
    for h in holders:
        if merged and h['who'] and h['who'] == merged[-1]['who'] and h.get('day') and merged[-1].get('day'):
            m = merged[-1]
            m.setdefault('labels', [m['text']])
            if h['text'] not in m['labels']:
                m['labels'].append(h['text'])
            m.update(s=h.get('s') or m.get('s'), last_day=h['day'], **({'url': h['url']} if h.get('url') else {}))
            m['text'] = ' → '.join([m['labels'][0]] + [t.replace('jlevy/squares', '') for t in m['labels'][1:]])
            continue
        merged.append(h)
    holders = merged
    for a, b in zip(holders, holders[1:]):
        a['until'] = b['date']
    known = [h for h in holders if h.get('s')]
    for a, b in zip(known, known[1:]):
        if b['s'] > a['s'] + 1e-9:
            problems.append((n, b.get('file') or b.get('event'), 'side goes up', '%s %.9f -> %s %.9f' % (fmt_date(a['date']), a['s'], fmt_date(b['date']), b['s'])))
    if not holders:
        rec = [f for f in files if f['record']]
        k = math.isqrt(n - 1) + 1
        proved = [c for c in all_credits if c['kind'] == 'proved']
        if rec and rec[0]['s'] is not None and abs(rec[0]['s'] - k) < 1e-12 or not rec and not r:
            holders = [{'kind': 'grid', 'who': [], 'date': None, 's': float(k), 'text': 'grid packing'}]
        elif proved:  # 'Proved by X': the packing and its proof are one credit
            c = proved[0]
            holders = [{'kind': 'found', 'who': c['who'], 'date': c['date'], 'based_on': [], 'text': c['text'],
                        'file': c['file'], 's': rec[0]['s'] if rec else None, 'proved': True}]
        else:
            problems.append((n, None, 'no holder', ''))
    return holders, [c for c in all_credits if c['kind'] not in ('found', 'improved', 'construction', 'previous')]


def lower_chain(n, LB, problems):
    """The floors the Bounds page counts (status proved / preprint, not `reported`), each one that beats the last; the
    rest (Nagamochi 2005's gap, Green's claims, reported-only figures) as `uncounted` credits.  Entries marked
    `awaiting` are held to the stricter bar Evan set for our own s(k^2-4) family (10-12): listed, not counted, until
    jlevy/squares verifies them."""
    e = LB.get(str(n))
    if not e or 'history' not in e:
        return [], []
    S = LB.get('sources', {})
    def item(h):
        src = S.get(h['source'], {})
        return {'source': h['source'], 'who': SOURCE_PEOPLE.get(h['source'], []), 'authors': src.get('authors'),
                'year': h.get('year'), 'value': h['value'], 'status': 'awaiting' if h.get('awaiting') else h.get('status'),
                'url': h.get('url') or src.get('url')}
    hist = [h for h in e['history'] if h['source'] != 'area-bound']
    hist.sort(key=lambda h: (h.get('year') or 0, h.get('date') or '', h['value']))
    out, other, best = [], [], -1
    for h in hist:
        if h.get('status') not in ('proved', 'preprint') or h.get('reported') or h.get('awaiting'):
            other.append(item(h))
            continue
        if h['value'] <= best + 1e-12:
            continue
        best = h['value']
        out.append(item(h))
    for a, b in zip(out, out[1:]):
        a['until'] = b['year']
    if out and abs(out[-1]['value'] - e['best']['value']) > 1e-12 and not e['best'].get('awaiting'):
        problems.append((n, None, 'lower chain', 'ends at %r, Bounds page floor %r' % (out[-1]['value'], e['best']['value'])))
    return out, other


def build():
    R, REC, LB, OV = load()
    files = catalogue_files(R)
    problems = []
    N = {}
    for n in range(1, NMAX + 1):
        fs = files.get(n, [])
        if not fs and str(n) not in REC['n']:
            k = math.isqrt(n - 1) + 1
            N[n] = {'upper': [{'kind': 'grid', 'who': [], 'date': None, 's': float(k), 'text': 'grid packing'}], 'other': []}
        else:
            holders, other = upper_chain(n, fs, REC['n'], REC['events'], OV, problems)
            N[n] = {'upper': holders, 'other': other}
        lc, lo = lower_chain(n, LB, problems)
        if lc:
            N[n]['lower'] = lc
        if lo:
            N[n]['lower_other'] = lo
    return N, problems


def held(N, today=None):
    """person -> tallies.  Upper bounds: `now` / `ever` (n where they found or improved the record), `years` (record-years:
    each holder's years, split evenly between co-finders), `built` (n whose current record is a construction on their
    packing, e.g. Göbel's L-additions: credited, but not counted as finds).  Lower bounds (n <= 324): `lb_now`, `lb_ever`."""
    today = today or [TODAY.year, TODAY.month]
    T = defaultdict(lambda: {'now': [], 'ever': set(), 'years': 0.0, 'built': [], 'lb_now': [], 'lb_ever': set()})
    for n, e in N.items():
        hs = e['upper']
        for i, h in enumerate(hs):
            last = i == len(hs) - 1
            if h['kind'] == 'construction':
                if last:
                    for w in h.get('origin', []):
                        T[w]['built'].append(n)
                continue
            if not h['who'] or not h.get('date'):
                continue
            y = max(0.0, years(h['date'], h.get('until') or list(today))) / len(h['who'])
            for w in h['who']:
                T[w]['ever'].add(n)
                T[w]['years'] += y
                if last:
                    T[w]['now'].append(n)
        lc = e.get('lower', [])
        for i, l in enumerate(lc):
            for w in l['who']:
                T[w]['lb_ever'].add(n)
                if i == len(lc) - 1:
                    T[w]['lb_now'].append(n)
    return T


def fmt_date(d):
    if not d:
        return '?'
    return '%d' % d[0] if not d[1] else '%d-%02d' % (d[0], d[1])


def show(N, n):
    e = N[n]
    print('n = %d' % n)
    for h in e['upper']:
        who = ', '.join(PEOPLE[w]['name'] for w in h['who']) or ('(construction on %s)' % ', '.join(PEOPLE[w]['name'] for w in h.get('origin', [])) if h.get('origin') else '-')
        print('  %-12s %-8s → %-8s %-12s %s  | %s' % (h['kind'], fmt_date(h.get('date')), fmt_date(h.get('until')),
                                                    ('%.9f' % h['s']) if h.get('s') else '', who, h.get('text', '')[:90]))
    for c in e.get('other', []):
        print('    (%s) %s %s' % (c['kind'], fmt_date(c.get('date')), ', '.join(PEOPLE[w]['name'] for w in c['who'] if w in PEOPLE)))
    for l in e.get('lower', []):
        print('  lower %.6f %s %s %s' % (l['value'], l['year'], l['source'], l['status']))


def review(N, T):
    L = ['# Credits review (generated by `site/tools/credits.py --review`; do not edit)', '',
         'Every n: the record holders in order (kind, from → until, side, who, source sentence), then other credits',
         '(refinements, alternatives, 2026 improvements) and the lower-bound chain (n <= 324).  Read for wrong names,',
         'wrong dates, or a holder that should be a refinement; fix in `data/credit_overrides.yaml` (keyed by catalogue',
         'file) or the parser.', '']
    def nm(ids):
        return ', '.join(PEOPLE[w]['name'] for w in ids)
    for n in sorted(N):
        e = N[n]
        if [h['kind'] for h in e['upper']] == ['grid'] and not e.get('lower'):
            continue
        L.append('## %d' % n)
        for h in e['upper']:
            who = nm(h['who']) or ('construction on %s' % nm(h.get('origin', [])) if h.get('origin') else h['kind'])
            L.append('- **%s** %s → %s%s: %s. _%s_' % (h['kind'], fmt_date(h.get('date')), fmt_date(h.get('until')) if h.get('until') else 'now',
                                                   (', s = %.9f' % h['s']) if h.get('s') else '', who, h.get('text', '').replace('$', '')))
        for c in e.get('other', []):
            L.append('  - (%s) %s %s' % (c['kind'], fmt_date(c.get('date')), nm([w for w in c['who'] if w in PEOPLE])))
        for l in e.get('lower', []):
            L.append('  - lower %.6f (%s, %s): %s' % (l['value'], l['year'], l['status'], nm(l['who']) or l['source']))
        L.append('')
    open(os.path.join(SITE, 'notes', 'credits-review.md'), 'w').write('\n'.join(L))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--show', type=int, nargs='*')
    ap.add_argument('--review', action='store_true')
    a = ap.parse_args()
    N, problems = build()
    if a.show:
        for n in a.show:
            show(N, n)
        return
    for p in problems:
        print('PROBLEM n=%s %s: %s %s' % p, file=sys.stderr)
    T = held(N)
    if a.review:
        review(N, T)
    out = {'note': 'generated by site/tools/credits.py; see its docstring', 'generated': TODAY.isoformat(),
           'people': {k: {kk: v[kk] for kk in ('name', 'handle', 'note') if kk in v} for k, v in PEOPLE.items()},
           'n': N}
    json.dump(out, open(os.path.join(D, 'credits.json'), 'w'), ensure_ascii=False, separators=(',', ':'))
    print('credits: %d n, %d problems' % (len(N), len(problems)))
    print('  %-24s %5s %5s %8s %6s   %6s %6s' % ('', 'now', 'ever', 'rec-yrs', 'built', 'lb now', 'lb ever'))
    for w, t in sorted(T.items(), key=lambda x: (-len(x[1]['now']), -x[1]['years'])):
        print('  %-24s %5d %5d %8.1f %6d   %6d %6d' % (PEOPLE[w]['name'], len(t['now']), len(t['ever']), t['years'],
                                                   len(t['built']), len(t['lb_now']), len(t['lb_ever'])))


if __name__ == '__main__':
    main()
