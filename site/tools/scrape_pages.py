#!/usr/bin/env python3
"""Scrape David Ellsworth's 'Squares in Squares' HTML pages into structured JSON.

Main page: one <div class="box"> per record (label like "14, 15" means the picture shows the
largest n; smaller ones come from removing squares).  Compared/rigid/group pages: <table
class="svg-table"> with a row of boxes and a row of descriptions.  Commented-out boxes are kept
with hidden=true (they are superseded/merged entries).
"""
import re, json, html, sys, os
D = os.path.dirname(os.path.abspath(__file__)) + '/../data/ellsworth/'
PAGES = ['squares_in_squares.html', 'squares_in_squares__compared.html', 'squares_in_squares__rigid.html',
         'squares_in_squares__n^2-n-1.html', 'squares_in_squares__Göbel_squares.html',
         'squares_in_squares__Göbel_strips.html', 'squares_in_squares__triangular_table.html']

def clean(s):
    s = re.sub(r'<br\s*/?>', '\n', s)
    s = re.sub(r'<[^>]+>', '', s)
    s = html.unescape(s)
    s = re.sub(r'[ \t]+', ' ', s)
    return '\n'.join(l.strip() for l in s.split('\n') if l.strip())

def parse_desc(d):
    """Split a description cell into s-decimal, s-latex, polynomials, prose."""
    out = {}
    fr1 = re.search(r'<span class="frame1">(.*?)</span>', d, re.S)
    polys = re.findall(r'<span class="frames">(.*?)</span>', d, re.S)
    if fr1:
        s_tex = clean(fr1.group(1)); out['polys'] = [clean(p).strip('$ ') for p in polys]
        d = re.sub(r'<span class="toggle">.*?</span>\s*</span>', '', d, flags=re.S)
    else:
        m = re.search(r'\$(.*?)\$', d, re.S)
        s_tex = clean(m.group(0)) if m else ''
        if m: d = d.replace(m.group(0), '', 1)
    out['s_tex'] = s_tex.strip('$ ')
    m = re.search(r'\\Nn\{([0-9.]+)\}', s_tex)
    if m: out['s_dec'] = m.group(1)
    else:
        m = re.search(r's\s*=\s*([0-9]+(?:\.[0-9]+)?)\s*$', out['s_tex'])
        if m: out['s_dec'] = m.group(1)
    out['prose'] = clean(d)
    out['rigid'] = 'Rigid' in out['prose']
    return out

def box_svgs(b):
    href = re.search(r'<a href="([^"]+\.svg)"', b); emb = re.search(r'<embed src="([^"]+\.svg)"', b)
    return (href.group(1) if href else None), (emb.group(1) if emb else None)

def scrape_main(t):
    recs = []
    # mark commented boxes
    for hidden, chunk in [(True, c) for c in re.findall(r'<!--(.*?)-->', t, re.S)] + [(False, re.sub(r'<!--.*?-->', '', t, flags=re.S))]:
        for b in re.split(r'<div class="box">', chunk)[1:]:
            lab = re.search(r'<font size="\+3">([^<]*)<br>', b)
            if not lab: continue
            ns = [int(x) for x in re.findall(r'\d+', lab.group(1))]
            if not ns: continue
            href, emb = box_svgs(b)
            desc = re.search(r'<div align="center">(.*?)</div></div>', b, re.S)
            r = {'ns': ns, 'n': max(ns), 'svg': emb or href, 'svg_link': href, 'hidden': hidden}
            if desc: r.update(parse_desc(desc.group(1)))
            recs.append(r)
    return recs

def scrape_tables(t, page):
    out = []
    t = re.sub(r'<!--.*?-->', '', t, flags=re.S)
    for tb in re.findall(r'<table class="svg-table".*?</table>', t, re.S):
        rows = re.findall(r'<tr[^>]*>(.*?)</tr>', tb, re.S)
        if len(rows) < 2: continue
        boxes = re.findall(r'<td class="box">(.*?)</td>', rows[0], re.S)
        descs = re.findall(r'<td class="desc">(.*?)</td>', rows[1], re.S)
        n = None
        for i, b in enumerate(boxes):
            lab = re.search(r'<font size="\+3">\s*([0-9, ]*)\.?<br>', b)
            if lab and re.findall(r'\d+', lab.group(1)): n = max(int(x) for x in re.findall(r'\d+', lab.group(1)))
            href, emb = box_svgs(b)
            r = {'n': n, 'svg': emb or href, 'svg_link': href, 'page': page, 'col': i, 'ncols': len(boxes)}
            if i < len(descs): r.update(parse_desc(descs[i]))
            out.append(r)
    return out

if __name__ == '__main__':
    result = {}
    for p in PAGES:
        t = open(D + p, encoding='utf-8').read()
        result[p] = scrape_main(t) if p == 'squares_in_squares.html' else scrape_tables(t, p)
        sys.stderr.write('%s: %d entries\n' % (p, len(result[p])))
    json.dump(result, open(D + '../records_raw.json', 'w'), indent=1, ensure_ascii=False)
