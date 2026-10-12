#!/usr/bin/env python3
"""Tour v2 set-up (10-11): reps in display orientation (mapdata.py tree propagation), island coverage, reuse of the
tour-v1 expansions whose relative orientation is unchanged.

  python3 landscape/tour2_prep.py [--min-island 5] [--pocket2 6]      (run from packer/; seconds)

Reps = the v1 reps still certified + the best-side minimum of every fp island (cut 0.02) with >= min-island members not
yet covered + farthest-point picks in the second deep pocket up to --pocket2 reps.  Files: tour2/reps/*.txt (f64, deg,
display orientation: v1 files mapped by the box symmetry matching their display image, new reps = store packing + g);
tour2/reps.json (census_index = row of map/minima.json, for `build.py barrier --near K --dist map/fpD.npy`);
tour2/barriers.json seeded with v1 pairs whose two files got the same symmetry, identity or transpose: the v1 walk is
the display-frame walk (frame 0) or the transposed-frame walk (frame 1, build.py --frame 1)."""
import argparse, collections, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, PK)
import numpy as np
from scipy.sparse.csgraph import connected_components
from mapdata import sym, SYMS
from pk import fingerprint as fp
from pk.store import Store
RUN = os.path.join(PK, 'runs/landscape110'); OUT = f'{RUN}/tour2'


def loaddeg(p):
    L = [l.split() for l in open(p) if l.strip() and not l.startswith('#')]
    n, s = int(L[0][0]), float(L[0][1]); return s, np.array([[float(v) for v in l[:3]] for l in L[1:n + 1]])


def write(p, s, sq):
    with open(p, 'w') as f:
        f.write(f'{len(sq)} {float(s)!r}\n' + ''.join(f'{float(x)!r} {float(y)!r} {float(t)!r}\n' for x, y, t in sq))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--min-island', type=int, default=5); ap.add_argument('--pocket2', type=int, default=6)
    a = ap.parse_args(); os.makedirs(f'{OUT}/reps', exist_ok=True)
    M = json.load(open(f'{RUN}/map/minima.json')); D = np.load(f'{RUN}/map/fpD.npy'); S = np.array([float(m['S']) for m in M])
    k, lab = connected_components(D < .02, directed=False)
    rank = {c: r for r, c in enumerate(sorted(set(lab), key=lambda c: S[lab == c].min()))}; isl = np.array([rank[c] for c in lab])
    byS = {m['S'][:13]: i for i, m in enumerate(M)}
    disp = lambda i: np.array(M[i]['sq']).reshape(-1, 3) / [1e4, 1e4, 100]
    st = Store(); reps = []                    # (census index, source, g applied to the source file or None)
    R1 = json.load(open(f'{RUN}/reps.json')); drop = []
    for r in R1:
        i = byS.get(r['S'][:13])
        if i is None: drop.append(r['S'][:13]); continue
        s, sq = loaddeg(os.path.join(PK, r['path']) if not os.path.isabs(r['path']) else r['path'])
        T = fp.image(S[i], disp(i)); g = SYMS[int(np.argmin([np.abs(fp.image(s, sym(sq, s, *e)) - T).sum() for e in SYMS]))]
        reps.append(dict(ci=i, src='v1', v1S=r['S'], g=g, s=s, sq=sym(sq, s, *g)))
    have = {isl[r['ci']] for r in reps}
    big = [r_ for r_, n in sorted(collections.Counter(isl).items()) if n >= a.min_island]
    def add(i, why):
        p = st.get(M[i]['pid']); s = float(p.s); sq = sym(np.array([tuple(map(float, q)) for q in p.sq]), s, *M[i]['g'])
        reps.append(dict(ci=i, src=why, g=None, s=s, sq=sq))
    for r_ in big:
        if r_ not in have:
            m = np.where(isl == r_)[0]; add(int(m[np.argmin(S[m])]), f'island {r_}')
    m = np.where(isl == 1)[0]                  # second deep pocket: farthest-point picks
    while sum(isl[r['ci']] == 1 for r in reps) < a.pocket2:
        cur = [r['ci'] for r in reps if isl[r['ci']] == 1]
        add(int(m[np.argmax(D[np.ix_(m, cur)].min(1))]), 'pocket 2')
    reps.sort(key=lambda r: S[r['ci']])
    out = []
    for r in reps:
        p = f'{OUT}/reps/m{M[r["ci"]]["S"][:13]}.txt'; write(p, r['s'], r['sq'])
        out.append(dict(S=M[r['ci']]['S'], census_index=int(r['ci']), path=os.path.relpath(p, PK), src=r['src'], island=int(isl[r['ci']]),
                        tilted=110 - int(M[r['ci']]['fam'][2]), v1_sym=r['g']))
    json.dump(out, open(f'{OUT}/reps.json', 'w'), indent=0)
    v1 = {r['v1S']: r for r in reps if r['src'] == 'v1'}; seed = []
    for b in json.load(open(f'{RUN}/barriers.json')):
        A, B = v1.get(b['S_i']), v1.get(b['S_j'])
        if A and B and A['g'] == B['g']:
            if A['g'] not in ((0, 0, 0), (1, 0, 0)):
                continue                       # other symmetries move the spread origin: not a pure frame change
            seed.append(dict(b, S_i=M[A['ci']]['S'], S_j=M[B['ci']]['S'], frame=int(A['g'] == (1, 0, 0))))
    json.dump(seed, open(f'{OUT}/barriers.json', 'w'), indent=0)
    cov = collections.Counter(r['island'] for r in out)
    print(f'{len(out)} reps ({sum(r["src"] == "v1" for r in out)} from v1, dropped {drop}); islands >= {a.min_island}: {len(big)}, covered '
          f'{sum(1 for r_ in big if cov[r_])}; pocket 1 / 2: {cov[0]} / {cov[1]}; tilted counts {dict(sorted(collections.Counter(r["tilted"] for r in out).items()))}')
    print(f'seeded {len(seed)} v1 pairs (same symmetry on both files) of {len(json.load(open(f"{RUN}/barriers.json")))}')


if __name__ == '__main__':
    main()
