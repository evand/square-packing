"""Rediscovery report for a corpus whose parents exclude a withheld neighbourhood of the record.

  rediscovery.py CORPUS VARIANT --dist runs/dist110.json [--withhold 0.02] [--by arm|band]

Per group: proposals, CPU s, distinct certified basins, new (not in the corpus' frozen known set), record hits (proposals
landing in the record basin), withheld basins rediscovered (certified basins within --withhold of the record), per CPU-h.
"""
import argparse, collections, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pk.store import Store

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs', 'corpus')
key = lambda s: f'{s:.9f}'

ap = argparse.ArgumentParser(); ap.add_argument('corpus'); ap.add_argument('variant')
ap.add_argument('--dist', required=True); ap.add_argument('--withhold', type=float, default=0.02); ap.add_argument('--by', default='arm')
a = ap.parse_args()
d = os.path.join(ROOT, a.corpus)
meta = json.load(open(f'{d}/meta.json')); n, k = meta['n'], meta['k']
known = {key(s) for s in meta['known']}
D = json.load(open(a.dist))
rec_key = key(D['rec_s'])
withheld = {key(x['s']) for x in D['basins'] if x['d'] < a.withhold}
pdist = {x['id']: x['d'] for x in D['basins']}
st = Store()
B = {r['key']: r['cls'] for r in st.db.execute('SELECT key, cls FROM basin WHERE n=?', (n,))}
items = {it['i']: it for it in map(json.loads, open(f'{d}/items.jsonl'))}
R = [json.loads(l) for l in open(f'{d}/out_{a.variant}.jsonl')]
G = collections.defaultdict(list)
for r in R:
    it = items[r['i']]
    if a.by == 'arm':
        g = it['arm']
    else:
        dd = pdist.get(it['parent'], -1)
        g = '<0.05' if dd < 0.05 else '<0.1' if dd < 0.1 else '>=0.1'
    G[g].append(r)
print(f'{a.corpus}/{a.variant}: {len(R)} proposals; withheld basins (d < {a.withhold}): {len(withheld)}; record {rec_key}')
print(f'{"group":10s} {"N":>5s} {"CPU s":>6s} {"cert":>5s} {"new":>4s} {"rec hits":>8s} {"withheld":>8s} {"new/CPUh":>8s} {"wh/CPUh":>7s}')
for g, rs in sorted(G.items()):
    cpu = sum(r['sec'] for r in rs)
    cert, new, wh, rh = set(), set(), set(), 0
    for r in rs:
        if r['st'] != 'new?' or r.get('unpolished') or r['s'] >= k:
            continue
        kk = key(r['s'])
        if B.get(kk) != 'certified':
            continue
        cert.add(kk)
        if kk not in known: new.add(kk)
        if kk in withheld: wh.add(kk)
        if kk == rec_key: rh += 1
    print(f'{g:10s} {len(rs):5d} {cpu:6.0f} {len(cert):5d} {len(new):4d} {rh:8d} {len(wh):8d} {len(new) / cpu * 3600:8.1f} {len(wh) / cpu * 3600:7.1f}')
