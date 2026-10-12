#!/usr/bin/env python3
"""Pack the explainer data into one JSON (v2, 10-11).

  python3 landscape/export.py [--out ../../site/www/data/landscape110.json]      (run from packer/, after mapdata.py,
                                                                                  layout.py and the tour2 build)
Inputs (runs/landscape110/):
  map/minima.json, map/fpD.npy   certified minima below 11 (store snapshot), display orientation, tilt counts (mapdata.py)
  map/layout.json                map coordinates (layout.py)
  tour2/reps.json, cycle.json    tour v2 (tour2_prep.py, build.py barrier / cycle), display orientation
  settle/*.json                  settling runs (settle.py)
Census rows: S (16 digits), tilted (= 110 - axis), counts [larger band, smaller band, axis, junction, other], island (fp
cluster at 0.02, ranked by best side), flat (exact flat motions), xy (map), sq.  Forces: tour stops with a v1 exactsolve
report (forces.py on the v1 file, mapped to the display orientation by the stop's v1 symmetry).
Coordinates are integers: x, y in 1e-4 units, theta in 1e-2 degrees.  Sides as strings."""
import argparse, glob, json, math, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, PK)
import numpy as np
from scipy.sparse.csgraph import connected_components

RUN = os.path.join(PK, 'runs/landscape110')
SITE = os.path.join(PK, '..', '..', 'site', 'www', 'data', 'landscape110.json')


def enc_rad(sq):
    out = []
    for x, y, t in sq:
        out += [round(x * 1e4), round(y * 1e4), round(math.degrees(t) * 100)]
    return out


def nflat(st):
    m = re.search(r'modulo (\d+)', st or '') or re.search(r'PSD with (\d+) zero', st or '')
    return int(m.group(1)) if m else 0


def tsym_pt(x, y, g, s):
    """box symmetry (transpose, flip y, flip x) on a point, as mapdata.sym."""
    t, fy, fx = g
    if t: x, y = y, x
    if fy: y = s - y
    if fx: x = s - x
    return x, y


def tsym_ang(th, g):
    t, fy, fx = g
    if t: th = math.pi / 2 - th
    if fy: th = -th
    if fx: th = -th
    return th


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default=SITE); ap.add_argument('--no-forces', action='store_true')
    a = ap.parse_args()
    M = json.load(open(f'{RUN}/map/minima.json')); D = np.load(f'{RUN}/map/fpD.npy'); L = json.load(open(f'{RUN}/map/layout.json'))
    S = np.array([float(m['S']) for m in M]); assert len(L['xy']) == len(M)
    k, lab = connected_components(D < .02, directed=False)
    rank = {c: r for r, c in enumerate(sorted(set(lab), key=lambda c: S[lab == c].min()))}
    from pk.store import Store
    st = Store(); status = {r[0]: r[1] for r in st.db.execute("SELECT key, status FROM minimum WHERE n = 110")}
    census = [dict(S=m['S'], tilted=110 - int(m['fam'][2]), counts=[int(v) for v in m['fam']], island=rank[lab[i]],
                   flat=nflat(status.get(m['key'])), xy=[round(v, 4) for v in L['xy'][i]], sq=m['sq']) for i, m in enumerate(M)]
    # ---- tour v2
    R = json.load(open(f'{RUN}/tour2/reps.json')); cyc = json.load(open(f'{RUN}/tour2/cycle.json'))
    stops = [dict(S=r['S'], census_index=r['census_index'], tilted=r['tilted'], island=rank[lab[r['census_index']]]) for r in R]
    # one copy per minimum: tour stops are drawn everywhere (census, map, tour) from the tour's own file; copies of a
    # minimum can differ by flat motions / slack (up to ~0.5 at 110), which would show as a jump between sections
    for r in R:
        Lr = [l.split() for l in open(os.path.join(PK, r['path']))]
        n = int(Lr[0][0]); census[r['census_index']]['sq'] = [round(float(v) * f) for l in Lr[1:n + 1] for v, f in zip(l[:3], (1e4, 1e4, 100))]
    legs = [dict(i=l['i'], j=l['j'], ok=l.get('ok'), s=l.get('s'), frames=[[round(s_, 6), enc_rad(q)] for s_, q in l.get('frames', [])])
            for l in cyc['legs']]
    # ---- forces (tour stops with a v1 exactsolve report)
    forces_out = []
    if not a.no_forces:
        import forces
        v1 = {r['S'][:13]: r['path'] for r in json.load(open(f'{RUN}/reps.json'))}
        for idx in cyc['tour']:
            r = R[idx]; p = v1.get(r['S'][:13]); g = r.get('v1_sym')
            if not p or g is None:
                continue
            f = forces.forces(os.path.join(PK, p) if not os.path.isabs(p) else p); s = f['S']
            sq = [(*tsym_pt(x, y, g, s), tsym_ang(t, g)) for x, y, t in f['sq']]
            cts = [[c[0], c[1], c[2], *[round(v * 1e4) for v in tsym_pt(c[3], c[4], g, s)], round(c[5], 5)] for c in f['contacts']]
            forces_out.append(dict(S=r['S'], stop=idx, sq=enc_rad(sq), free=f['free'], flat=f['flat'], contacts=cts))
        forces_out.sort(key=lambda d: float(d['S']))
    # ---- settle (unchanged runs)
    settle = []
    for fn in sorted(glob.glob(f'{RUN}/settle/*.json')):
        d = json.load(open(fn)); tr = d['traj']
        if len(tr) > 160:
            idx = sorted(set(np.linspace(0, len(tr) - 1, 160).round().astype(int).tolist())); tr = [tr[i] for i in idx]
        settle.append(dict(kind=d['kind'], seed=d['seed'], sigma=d.get('sigma'), s_final=d['s_final'], jammed=d['jammed'],
                           frames=[[s_, enc_rad(q)] for s_, q in tr]))
    out = dict(n=110, snapshot=dict(minima=len(census), store_sha=__import__('hashlib').sha256(open(os.path.join(PK, 'runs', 'store.sqlite'), 'rb').read()).hexdigest()[:12], layout=dict((k_, L[k_]) for k_ in ('p', 'k', 'min_dist', 'spread', 'seed'))),
               census=census, stops=stops, tour=cyc['tour'], legs=legs, forces=forces_out, settle=settle)
    json.dump(out, open(a.out, 'w'), separators=(',', ':'))
    print(a.out, os.path.getsize(a.out) // 1024, 'KB;', len(census), 'census,', len(stops), 'stops,', len(legs), 'legs,',
          len(forces_out), 'force figures,', len(settle), 'settle')


if __name__ == '__main__':
    main()
