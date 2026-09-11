#!/usr/bin/env python3
"""Extract dated events from the scraped prose and build data/timeline.json:
   { n: { 'record': {file, s, date}, 'events': [ {file, s, date, verb, who, text} ... sorted ] } }
Dates are (year, month) with month guessed from 'early'/'late' when absent.  Every event keeps the
sentence it came from, so the site can show the source wording."""
import json, re, os
HERE = os.path.dirname(os.path.abspath(__file__)); D = HERE + '/../data/'
MONTHS = {m: i + 1 for i, m in enumerate(['January','February','March','April','May','June','July','August','September','October','November','December'])}
DATE = re.compile(r'\b(?:in|In|of)\s+(?:(early|late|mid)-?\s*)?(?:(' + '|'.join(MONTHS) + r')(?:[-–](' + '|'.join(MONTHS) + r'))?\s*)?(\d{4})(?:[-–](\d{4}))?\b')
VERB = re.compile(r'\b(Found first|Found|Improved|Proved optimal|Proved|Optimized|Refined|Refound|Fixed|Drafted|Converted|Bounded|Extended|Side length found|Improvement|Rigid alternative found|Discovered)\b(?:\s+(?:by|and \w+ by)\s+(.+?))?(?=\s+(?:in|In)\s+(?:early|late|mid|January|February|March|April|May|June|July|August|September|October|November|December|\d{4})|[.,;]|$)')

def parse_date(m):
    qual, m1, m2, y1, y2 = m
    year = int(y2 or y1)
    month = MONTHS.get(m2 or m1) if (m1 or m2) else (2 if qual == 'early' else 11 if qual == 'late' else 6 if qual == 'mid' else None)
    return {'year': year, 'month': month, 'approx': month is None and not qual, 'text': ' '.join(x for x in m if x)}

def events_of(prose):
    prose = (prose or '').replace('\xa0', ' ').replace('\n', ' ')
    out = []
    for sent in re.split(r'(?<=[.!?])\s+(?=[A-Z])', prose):
        dm = DATE.search(sent)
        vm = VERB.search(sent)
        if not dm: continue
        d = parse_date(dm.groups())
        who = (vm.group(2) or '').strip() if vm else ''
        who = re.sub(r'\s*\(.*?\)\s*', ' ', who).strip(' ,')
        out.append({'verb': vm.group(1) if vm else '', 'who': who, 'date': d, 'text': sent.strip()})
    return out

def key(d): return (d['year'], d['month'] or 6)

def main():
    R = json.load(open(D + 'records_raw.json'))
    P = json.load(open(D + 'packings.json'))
    files = {}
    for page, es in R.items():
        for e in es:
            f = e.get('svg')
            if not f or f in files: continue
            ev = events_of(e.get('prose'))
            s = P.get(f, {}).get('s') or e.get('s_dec')
            files[f] = {'file': f, 'n': e.get('n'), 's': float(s) if s else None, 'events': ev, 'prose': e.get('prose'),
                        's_tex': e.get('s_tex'), 'rigid': e.get('rigid', False)}
    # per n
    T = {}
    rec_files = {r['svg'] for r in R['squares_in_squares.html'] if not r['hidden']}
    for f, d in files.items():
        n = d['n'] or (P.get(f, {}).get('n'))
        if not n: continue
        dates = [key(e['date']) for e in d['events'] if e['verb'] in ('Found', 'Found first', 'Improved', 'Optimized', 'Refined', 'Refound', 'Fixed', 'Side length found', 'Improvement', 'Discovered', '')]
        d['date'] = max(dates) if dates else None
        d['is_record'] = f in rec_files
        T.setdefault(str(n), []).append(d)
    for n, L in T.items():
        L.sort(key=lambda d: (d['s'] is None, -(d['s'] or 0), d['date'] or (0, 0)))
    json.dump(T, open(D + 'timeline.json', 'w'), separators=(',', ':'), ensure_ascii=False)
    nd = sum(1 for L in T.values() for d in L if d['date'] is None and d['s'])
    print('timeline: %d n, %d files, %d undated' % (len(T), sum(len(L) for L in T.values()), nd))

if __name__ == '__main__': main()
