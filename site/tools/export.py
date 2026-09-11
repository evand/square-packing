#!/usr/bin/env python3
"""Run the analysis on every parsed packing and write the site's data files:
   www/data/p/<name>.json   one per SVG: squares (double precision), analysis
   www/data/index.json      catalogue: per n the record + alternatives, with scraped attribution
Usage: export.py [names...]   (default: all; skips ones already exported unless --force)
"""
import os, sys, json, re, time, traceback
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis import analyze
HERE = os.path.dirname(os.path.abspath(__file__))
D = HERE + '/../data/'; W = HERE + '/../www/data/'
os.makedirs(W + 'p', exist_ok=True)

def compact(a):
    """Trim the analysis for the browser."""
    out = dict(a)
    out['contacts'] = [{k: v for k, v in c.items() if k in ('i','j','wall','type','points','normal','exact','gap_mp','overlap')} for c in a['contacts']]
    for c in out['contacts']:
        c['points'] = [[round(x, 6) for x in p] for p in c['points']]     # drawing precision; exactness lives in gap_mp / exact
        c['normal'] = [round(x, 6) for x in c['normal']]
        if 'gap_mp' in c: c['gap_mp'] = float('%.3g' % c['gap_mp'])
    out['near_misses'] = [{k: (float('%.6g' % v) if isinstance(v, float) else v) for k, v in m.items()} for m in a['near_misses'][:25]]
    out['regions'] = {}
    for k, r in a['regions'].items():
        r = dict(r)
        if 'mask' in r:
            # run-length encode: per row, alternating run lengths starting with a run of '0'
            rows = []
            for row in r['mask']:
                runs, cur, cnt = [], '0', 0
                for ch in row:
                    if ch == cur: cnt += 1
                    else: runs.append(cnt); cur = ch; cnt = 1
                runs.append(cnt); rows.append(','.join(map(str, runs)))
            r['rle'] = ';'.join(rows); del r['mask']
        for kk in ('area', 'cell'):
            if kk in r: r[kk] = float('%.5g' % r[kk])
        for kk in ('dx', 'dy', 'rot'):
            if kk in r: r[kk] = [float('%.5g' % v) for v in r[kk]]
        if 'slide' in r:
            # A channel slide has zero clearance sideways, so a coarsely rounded direction cannot be
            # replayed without overlapping something: keep enough digits that it can be.
            r['slide'] = {'dir': [float('%.12g' % v) for v in r['slide']['dir']], 'len': float('%.12g' % r['slide']['len'])}
        out['regions'][k] = r
    out.pop('motion_dirs', None)
    return out

def main(argv):
    force = '--force' in argv; names = [a for a in argv if not a.startswith('--')]
    P = json.load(open(D + 'packings.json'))
    todo = names or sorted(P.keys(), key=lambda k: (int(re.match(r'square-(\d+)', k).group(1)) if re.match(r'square-(\d+)', k) else 0, k))
    log = open(D + 'export.log', 'a')
    for name in todo:
        outp = W + 'p/' + name.replace('.svg', '.json')
        if os.path.exists(outp) and not force: continue
        p = P[name]
        if [e for e in p['errors'] if not e.startswith('warn')]:
            log.write('%s SKIP %s\n' % (name, p['errors'][:2])); log.flush(); continue
        t = time.time()
        try:
            a = analyze(float(p['s']), p['squares'])
        except Exception as e:
            log.write('%s ERROR %s\n%s\n' % (name, e, traceback.format_exc())); log.flush(); continue
        rec = {'name': name, 's': p['s'], 'n': p['n'],
               'squares': [[float(x) for x in q] for q in p['squares']],
               'analysis': compact(a)}
        json.dump(rec, open(outp, 'w'), separators=(',', ':'))
        log.write('%s ok n=%d %.1fs free=%d wedged=%d mobile=%d rigid=%s\n' % (name, p['n'], time.time() - t, len(a['free']), len(a['wedged']), sum(a['mobile']), a['rigid'])); log.flush()
    build_index(P)

def usable(P, f):
    """A file we can actually show: parsed, and not rejected as an invalid packing."""
    p = P.get(f)
    return bool(p) and not [e for e in p['errors'] if not e.startswith('warn')]

def is_start(f):
    """An optimiser's starting configuration, not a record attempt (its s is above the record)."""
    return '_start' in f

def build_index(P):
    R = json.load(open(D + 'records_raw.json'))
    main_page = R['squares_in_squares.html']
    idx = {'records': {}, 'files': {}}
    for r in main_page:
        if r['hidden']: continue
        # the box's <embed> is the picture and the <a href> the link; for the "2, 3" box the embed
        # points at square-3.svg, which does not exist -- fall back to the link.
        svg = r['svg']
        if not usable(P, svg) and usable(P, r.get('svg_link')): svg = r['svg_link']
        # what the picture actually contains, which need not be the largest n on the label
        pictured = P[svg]['n'] if usable(P, svg) else r['n']
        for n in r['ns']:
            idx['records'][n] = {'n': n, 'svg': svg, 'svg_link': r['svg_link'], 'pictured_n': pictured,
                                 's_tex': r.get('s_tex'), 's_dec': r.get('s_dec'), 'polys': r.get('polys', []),
                                 'prose': r.get('prose'), 'rigid': r.get('rigid')}
    # every file: metadata from any page that shows it
    for page, entries in R.items():
        for e in entries:
            f = e.get('svg')
            if not f: continue
            d = idx['files'].setdefault(f, {'svg': f, 'pages': [], 'n': e.get('n')})
            d['pages'].append(page.replace('squares_in_squares', 'sis').replace('.html', ''))
            for k in ('s_tex', 's_dec', 'polys', 'prose', 'rigid'):
                if e.get(k) and not d.get(k): d[k] = e[k]
            # Ellsworth's rigid page presents each rigid packing next to the better non-rigid one it
            # loses to, so membership alone says nothing; his wording does.
            if page == 'squares_in_squares__rigid.html':
                pr = (e.get('prose') or '').lower()
                d['rigid_claim'] = 'no' if ('not rigid' in pr or 'semi-rigid' in pr) else 'yes'
    # analysis summaries
    for f, p in P.items():
        outp = W + 'p/' + f.replace('.svg', '.json')
        d = idx['files'].setdefault(f, {'svg': f, 'pages': [], 'n': p['n']})
        d['n_parsed'] = p['n']; d['s'] = p['s']; d['errors'] = p['errors']
        if is_start(f): d['start'] = True
        if os.path.exists(outp):
            a = json.load(open(outp))['analysis']
            d['summary'] = {'category': a['category'], 'n_angles': a['n_angles'], 'n_rotated': a['n_rotated'],
                            'symmetry': a['symmetry']['class'], 'rigid': a['rigid'], 'free': len(a['free']),
                            'wedged': len(a['wedged']), 'mobile': sum(a['mobile']), 'waste': a['waste'],
                            'angles': [[round(g['theta'], 6), g['count']] for g in a['angle_groups']],
                            'min_gap': a['near_misses'][0]['gap'] if a['near_misses'] else None}
    # Cross-check: our rigidity analysis is written independently of Ellsworth's page, so how his
    # rigid page reads and how our analysis comes out is the one external check the site has.
    # Both directions matter -- agreeing on what is *not* rigid is the harder half.
    claim = {f: d['rigid_claim'] for f, d in idx['files'].items() if d.get('rigid_claim') and d.get('summary')}
    ours = lambda f: idx['files'][f]['summary']['rigid']
    yes = sorted(f for f, c in claim.items() if c == 'yes')
    no = sorted(f for f, c in claim.items() if c == 'no')
    idx['checks'] = {'rigid_yes': len(yes), 'rigid_yes_disagree': [f for f in yes if not ours(f)],
                     'rigid_no': len(no), 'rigid_no_disagree': [f for f in no if ours(f)]}
    json.dump(idx, open(W + 'index.json', 'w'), separators=(',', ':'), ensure_ascii=False)
    c = idx['checks']
    print('index: %d records, %d files; rigid check: %d claimed rigid (%d disagree), %d claimed not (%d disagree)'
          % (len(idx['records']), len(idx['files']), c['rigid_yes'], len(c['rigid_yes_disagree']),
             c['rigid_no'], len(c['rigid_no_disagree'])))

if __name__ == '__main__':
    main(sys.argv[1:])
