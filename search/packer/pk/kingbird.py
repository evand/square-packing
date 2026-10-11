"""David Ellsworth's squares-in-squares catalogue (kingbird.myphotos.cc/packing) -> dated snapshot -> store (10-10).

Why: the catalogue's former records (`square-N_rK.svg`) are seeds and provenance we did not have; the 09-30 site
mirror held 549 files, the catalogue links 862 (476 revisions), most of them only from the `compared2` / `compared3`
pages.  Every packing is kept with its source file, fetch date, finder line and full header comment.

  crawl(): follow every squares_in_squares*.html page linked from the main page (same directory only), collect
           every linked .svg name.  Polite: one request at a time, UA with a contact URL, DELAY between requests.
  fetch(): snapshot dir runs/kingbird/<date>/ (pages/, svg/, manifest.json).  Revision files are reused from the
           newest earlier snapshot that has them (they are history); current-record files and new names are fetched.
  parse(): site/tools/parse_svg.parse (the site's validated parser, exact arithmetic), f64 overlap check, header
           comment -> finder / date / full text.
  import_snapshot(): store rows (source 'kingbird', ref '<file>@<date>', kind 'catalogue-current' / -revision /
           -variant), tags 'kingbird', 'kingbird:current' | 'kingbird:revision' | 'kingbird:variant'.

Third-party data (attribute "David Ellsworth, kingbird.myphotos.cc"); raw files stay in gitignored runs/.
"""
from __future__ import annotations
import hashlib, json, os, re, shutil, subprocess, sys, time, urllib.parse
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
ROOT = os.path.abspath(os.path.join(HERE, '..', 'runs', 'kingbird'))
BASE = 'https://kingbird.myphotos.cc/packing/'
UA = 'square-packing-explorer (+https://github.com/evand/square-packing)'
DELAY = 0.4
NAME = re.compile(r'^square-(\d+)(.*)\.svg$')


def _get(url, dest):
    r = subprocess.run(['curl', '-sS', '-f', '-A', UA, '-o', dest, url], capture_output=True, text=True)
    time.sleep(DELAY)
    return r.returncode == 0


def crawl(pages_dir):
    os.makedirs(pages_dir, exist_ok=True)
    todo, seen, svgs, failed = ['squares_in_squares.html'], [], set(), []
    while todo:
        p = todo.pop(0)
        if p in seen:
            continue
        seen.append(p)
        dest = os.path.join(pages_dir, p)
        if not _get(BASE + urllib.parse.quote(p), dest):
            failed.append(p); continue
        s = open(dest, encoding='utf-8', errors='replace').read()
        for h in re.findall(r'href="([^"#?]+)"', s):
            h = urllib.parse.unquote(h)
            if h.startswith(BASE):
                h = h[len(BASE):]
            if '/' in h or h.startswith('http'):
                continue
            if h.endswith('.svg') and NAME.match(h):
                svgs.add(h)
            elif h.startswith('squares_in_squares') and h.endswith('.html') and h not in seen:
                todo.append(h)
    return seen, sorted(svgs), failed


def kind_of(name):
    m = NAME.match(name)
    suf = m.group(2) if m else ''
    if suf == '':
        return 'current'
    if re.fullmatch(r'_r\d+[a-z]?', suf):
        return 'revision'
    return 'variant'               # 104s, 1102b, _r4_start, _no_min_rot, __invalid, ...


def snapshots():
    return sorted(d for d in os.listdir(ROOT) if re.fullmatch(r'\d{4}-\d{2}-\d{2}', d)) if os.path.isdir(ROOT) else []


def fetch(date=None, refetch_all=False, log=print):
    date = date or time.strftime('%Y-%m-%d')
    snap = os.path.join(ROOT, date)
    pages, names, failed = crawl(os.path.join(snap, 'pages'))
    log(f'crawl: {len(pages)} pages ({len(failed)} failed), {len(names)} svg names, '
        f'{sum(kind_of(n) == "revision" for n in names)} revisions')
    svgdir = os.path.join(snap, 'svg')
    os.makedirs(svgdir, exist_ok=True)
    older = [os.path.join(ROOT, d, 'svg') for d in snapshots() if d < date][::-1]
    got, reused, fail = 0, 0, []
    for nm in names:
        dest = os.path.join(svgdir, nm)
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            continue
        if not refetch_all and kind_of(nm) == 'revision':
            src = next((os.path.join(d, nm) for d in older if os.path.exists(os.path.join(d, nm))), None)
            if src:
                shutil.copy2(src, dest); reused += 1; continue
        if _get(BASE + nm, dest):
            got += 1
        else:
            fail.append(nm)
            if os.path.exists(dest):
                os.remove(dest)
    log(f'fetch: {got} fetched, {reused} reused from older snapshots, {len(fail)} failed')
    json.dump(dict(date=date, pages=pages, failed_pages=failed, names=names, failed=fail),
              open(os.path.join(snap, 'crawl.json'), 'w'), indent=1)
    return snap


def header(path):
    """Full header comment, finder (latest 'Improved by', else first 'Found by'), and that line's date text."""
    s = open(path, encoding='utf-8', errors='replace').read(20000)
    m = re.search(r'<!--(.*?)-->', s, re.S)
    text = re.sub(r'[ \t]+', ' ', m.group(1)).strip() if m else ''
    fin, when = None, None
    DATE = r'(?:on |in )?((?:early |late |mid-)?(?:[A-Z][a-z]+ )?(?:\d{1,2}, )?\d{4})'
    imp = list(re.finditer(r'Improved by ([^,.(\n]+?)(?= on | in |,|\.|\n|$)', text))
    fnd = list(re.finditer(r'[Ff]ound by ([^,.(\n]+?)(?= on | in |,|\.|\n|$)', text))
    m = imp[-1] if imp else (fnd[0] if fnd else None)
    if m:
        fin = m.group(1).strip()
        d = re.match(r'\s*' + DATE, text[m.end():m.end() + 60])
        when = d.group(1) if d else None
    return text, fin, (when or None)


def _parse_job(path):
    sys.path.insert(0, os.path.join(REPO, 'site', 'tools'))
    from parse_svg import parse
    sys.path.insert(0, os.path.join(HERE, '..'))
    from pk.packing import Packing
    nm = os.path.basename(path)
    rec = dict(name=nm, kind=kind_of(nm), sha256=hashlib.sha256(open(path, 'rb').read()).hexdigest())
    m = NAME.match(nm)
    rec['n_name'] = int(m.group(1)) if m else None
    try:
        s, sq = parse(path)
        rec['s_str'] = str(s)
        sq = np.array([[float(a) for a in r] for r in sq])
        rec['n'] = len(sq)
        p = Packing(float(s), sq, s_str=str(s))
        rec['pen'] = p.max_pen() if len(sq) <= 4000 else None
        rec['coords'] = sq.tolist()
        rec['s'] = float(s)
    except Exception as e:                      # keep the failure in the manifest
        rec['error'] = repr(e)[:300]
    rec['comment'], rec['finder'], rec['when'] = header(path)
    return rec


def parse_snapshot(snap, procs=4, log=print):
    paths = sorted(os.path.join(snap, 'svg', f) for f in os.listdir(os.path.join(snap, 'svg')) if f.endswith('.svg'))
    with ProcessPoolExecutor(procs) as ex:
        recs = list(ex.map(_parse_job, paths, chunksize=4))
    bad = [r for r in recs if 'error' in r]
    over = [r for r in recs if r.get('pen') and r['pen'] > 1e-9]
    mism = [r for r in recs if 'n' in r and r['n_name'] and r['n'] != r['n_name']]
    log(f'parse: {len(recs)} files, {len(bad)} parse errors, {len(over)} overlapping (pen > 1e-9), '
        f'{len(mism)} with square count != file-name n')
    man = [{k: v for k, v in r.items() if k != 'coords'} for r in recs]
    json.dump(man, open(os.path.join(snap, 'manifest.json'), 'w'), indent=1)
    return recs


def import_snapshot(st, snap, recs, log=print):
    from pk.packing import Packing
    date = os.path.basename(snap)
    made = 0
    for r in recs:
        if 'coords' not in r or not r['coords']:
            continue
        p = Packing(r['s'], np.array(r['coords']), s_str=r.get('s_str'))
        status = 'infeasible' if (r.get('pen') or 0) > 1e-9 or 'invalid' in r['name'] else 'raw'
        pid, created = st.add(p, 'kingbird', finder=r.get('finder'), ref=f"{r['name']}@{date}",
                              kind=f"catalogue-{r['kind']}", status=status,
                              tags=('kingbird', f"kingbird:{r['kind']}"), commit=False)
        r['store_id'] = pid
        made += created
    st.db.commit()
    man = json.load(open(os.path.join(snap, 'manifest.json')))
    ids = {r['name']: r.get('store_id') for r in recs}
    for m in man:
        m['store_id'] = ids.get(m['name'])
    json.dump(man, open(os.path.join(snap, 'manifest.json'), 'w'), indent=1)
    log(f'store: {sum(1 for r in recs if r.get("store_id"))} catalogue packings, {made} new rows')


def latest(name):
    """(sq, s) of a catalogue file in the newest snapshot that has it (name with or without .svg), or None."""
    nm = name if name.endswith('.svg') else name + '.svg'
    for d in snapshots()[::-1]:
        man = os.path.join(ROOT, d, 'manifest.json')
        if not os.path.exists(man):
            continue
        for m in json.load(open(man)):
            if m['name'] == nm and m.get('store_id'):
                from pk.store import Store
                p = Store().get(m['store_id'])
                return np.array(p.sq), float(p.s)
    return None


def sync(st, date=None, refetch_all=False, procs=4, log=print):
    snap = fetch(date, refetch_all, log)
    recs = parse_snapshot(snap, procs, log)
    import_snapshot(st, snap, recs, log)
    return snap
