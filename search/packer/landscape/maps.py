#!/usr/bin/env python3
"""Connectivity-map variants for the s(110) explainer (10-11; exploratory, figures to compare and iterate on).

  ./landscape/mapenv.sh python3 landscape/maps.py [--dir runs/landscape110/map]        (run from packer/)

Inputs: mapdata.py output (minima.json, fpD.npy); runs/landscape110/{reps,barriers}.json for the disconnectivity graph.
Outputs (PNG in --dir):
  layouts.png   UMAP (2 settings), t-SNE (2 perplexities), classical MDS; colour = family, grey = MST of fp distance
  layouts_S.png same layouts coloured by side
  families.png  family graph: one node per family (>= 3 members), MDS of medoid distances, edges = MST of the
                single-link (closest pair) distance between families
  disconn.png   disconnectivity graph of the 48 tour minima from the stored walk_carrot expansions (upper bounds)
"""
import argparse, collections, json, os, sys
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.sparse.csgraph import minimum_spanning_tree
HERE = os.path.dirname(os.path.abspath(__file__)); PK = os.path.dirname(HERE)
RUN = os.path.join(PK, 'runs/landscape110')
REC, TOP = 10.99678339663, 11.0          # colour scale: record to 11, absolute side, all digits
from matplotlib.ticker import FormatStrFormatter


def sidebar(sc, ax, label='side'):
    cb = plt.colorbar(sc, ax=ax, shrink=.8, label=label)
    cb.set_ticks([REC, 10.9975, 10.998, 10.9985, 10.999, 10.9995, TOP])
    cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.7f'))
    return cb


def famcolours(fam, k=10):
    top = [f for f, _ in collections.Counter(fam).most_common(k)]
    cm = plt.get_cmap('tab10')
    return {f: cm(i) for i, f in enumerate(top)}


def flab(f):
    return f'{f[0]}/{f[1]} tilted, {f[2]} axis, {f[3]} junction' + (f', {f[4]} other' if f[4] else '')


def mds(D, k=2):
    n = len(D); J = np.eye(n) - 1 / n; B = -.5 * J @ (D ** 2) @ J
    w, v = np.linalg.eigh(B); return v[:, -k:][:, ::-1] * np.sqrt(np.maximum(w[-k:][::-1], 0))


def scatter(a, E, fam, col, S, T, by, title):
    for i, j in zip(T.row, T.col):
        a.plot(E[[i, j], 0], E[[i, j], 1], color='0.82', lw=.5, zorder=0)
    if by == 'S':
        sc = a.scatter(E[:, 0], E[:, 1], c=S, cmap='viridis_r', s=9, vmin=REC, vmax=TOP); sidebar(sc, a)
    else:
        c = [col.get(f, (.7, .7, .7, 1)) for f in fam]
        a.scatter(E[:, 0], E[:, 1], c=c, s=9)
    b = int(np.argmin(S)); a.scatter(*E[b], marker='*', s=220, c='red', edgecolor='k', zorder=5)
    a.set_title(title, fontsize=10); a.set_xticks([]); a.set_yticks([])


def layouts(D, fam, S, out):
    import umap
    from sklearn.manifold import TSNE
    V = [('UMAP k=15, min_dist 0.1', lambda: umap.UMAP(n_neighbors=15, min_dist=.1, metric='precomputed', random_state=0).fit_transform(D)),
         ('UMAP k=60, min_dist 0.5', lambda: umap.UMAP(n_neighbors=60, min_dist=.5, metric='precomputed', random_state=0).fit_transform(D)),
         ('t-SNE perplexity 30', lambda: TSNE(perplexity=30, metric='precomputed', init='random', random_state=0).fit_transform(D)),
         ('t-SNE perplexity 150', lambda: TSNE(perplexity=150, metric='precomputed', init='random', random_state=0).fit_transform(D)),
         ('classical MDS', lambda: mds(D))]
    E = {name: f() for name, f in V}
    np.savez(f'{out}/layouts.npz', **{k.replace(' ', '_'): v for k, v in E.items()})
    T = minimum_spanning_tree(D).tocoo(); col = famcolours(fam)
    for by, fn in (('fam', 'layouts.png'), ('S', 'layouts_S.png')):
        fig, ax = plt.subplots(2, 3, figsize=(18, 12))
        for a, (name, Y) in zip(ax.flat, E.items()):
            scatter(a, Y, fam, col, S, T, by, name)
        a = ax.flat[-1]; a.axis('off')
        if by == 'fam':
            cnt = collections.Counter(fam)
            for f, c in col.items():
                a.scatter([], [], color=c, label=f'{flab(f)} ({cnt[f]})')
            a.scatter([], [], color=(.7, .7, .7), label=f'other families ({sum(v for f, v in cnt.items() if f not in col)})')
            a.legend(loc='center', fontsize=9, title=f'{len(S)} certified minima, {len(cnt)} families')
        fig.suptitle('s(110) certified minima: fingerprint distance; grey = minimum spanning tree; star = record')
        plt.tight_layout(); plt.savefig(f'{out}/{fn}', dpi=72); plt.close()


def families(D, fam, S, out, kmin=3):
    cnt = collections.Counter(fam); F = [f for f in cnt if cnt[f] >= kmin]
    idx = {f: np.array([i for i, g in enumerate(fam) if g == f]) for f in F}
    med = {f: idx[f][np.argmin(D[np.ix_(idx[f], idx[f])].sum(1))] for f in F}
    Dm = np.array([[D[med[f], med[g]] for g in F] for f in F])
    Ds = np.array([[0 if f == g else D[np.ix_(idx[f], idx[g])].min() for g in F] for f in F])
    E = mds(Dm); T = minimum_spanning_tree(Ds).tocoo()
    best = np.array([S[idx[f]].min() for f in F])
    fig, a = plt.subplots(figsize=(13, 11))
    for i, j, w in zip(T.row, T.col, T.data):
        a.plot(E[[i, j], 0], E[[i, j], 1], color='0.6', lw=1, zorder=0)
        a.text(*(E[i] + E[j]) / 2, f'{w:.3f}', fontsize=6, color='0.4', ha='center')
    sc = a.scatter(E[:, 0], E[:, 1], s=[12 * cnt[f] ** .8 for f in F], c=best, cmap='viridis_r', edgecolor='k', lw=.5, vmin=REC, vmax=TOP)
    for k, f in enumerate(F):
        a.annotate(f'{f[0]}/{f[1]}/{f[2]}/{f[3]}' + (f'/{f[4]}' if f[4] else '') + f' ({cnt[f]})', E[k], fontsize=7,
                   xytext=(4, 4), textcoords='offset points')
    sidebar(sc, a, 'best side in family'); a.set_xticks([]); a.set_yticks([])
    a.set_title(f'Families (L/B/axis/J, count) with >= {kmin} certified minima: position = MDS of medoid distances; '
                'edges = MST of closest-pair distance (labels)', fontsize=9)
    plt.tight_layout(); plt.savefig(f'{out}/families.png', dpi=80); plt.close()


def disconn(out, col):
    R = json.load(open(f'{RUN}/reps.json')); B = json.load(open(f'{RUN}/barriers.json'))
    S = [float(r['S']) for r in R]; n = len(R)
    def fam_of(path):
        L = open(path).read().split('\n'); m = int(L[0].split()[0])
        sq = [tuple(float(v) for v in l.split()[:3]) for l in L[1:m + 1]]
        sys.path.insert(0, PK); from mapdata import family; return family(sq)
    fam = [fam_of(r['path'] if os.path.isabs(r['path']) else os.path.join(PK, r['path'])) for r in R]
    E = sorted((b['s'], b['i'], b['j']) for b in B if b.get('s'))
    par = list(range(n)); node = {i: ('leaf', i) for i in range(n)}
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for s, i, j in E:
        a, b = find(i), find(j)
        if a != b:
            par[b] = a; node[a] = ('join', s, node[a], node[b])
    roots = {find(i) for i in range(n)}
    order = []
    def leaves(t):
        if t[0] == 'leaf': order.append(t[1]); return
        l, r = t[2], t[3]
        # lower subtree (by best side) first
        lo = lambda u: min(S[k] for k in leafset(u))
        for u in sorted((l, r), key=lo): leaves(u)
    def leafset(t): return [t[1]] if t[0] == 'leaf' else leafset(t[2]) + leafset(t[3])
    for r in sorted(roots, key=lambda r: min(S[k] for k in leafset(node[r]))): leaves(node[r])
    xpos = {k: x for x, k in enumerate(order)}
    fig, a = plt.subplots(figsize=(15, 8))
    def draw(t):
        if t[0] == 'leaf': return xpos[t[1]], S[t[1]]
        (xl, yl), (xr, yr) = draw(t[2]), draw(t[3]); h = t[1]
        a.plot([xl, xl], [yl, h], color='0.3', lw=.8); a.plot([xr, xr], [yr, h], color='0.3', lw=.8)
        a.plot([xl, xr], [h, h], color='0.3', lw=.8)
        return (xl + xr) / 2, h
    for r in roots: draw(node[r])
    for k in range(n):
        a.scatter(xpos[k], S[k], color=col.get(fam[k], (.6, .6, .6, 1)), s=30, zorder=3, edgecolor='k', lw=.3)
    a.axhline(11, ls='--', color='0.5', lw=.8)
    a.set_yscale('function', functions=(lambda y: np.sign(y - REC) * np.abs(y - REC) ** .5,
                                        lambda u: REC + np.sign(u) * u ** 2))
    a.set_ylim(REC - 2e-5, None)
    a.set_yticks([REC, 10.998, 10.999, 11.0, 11.01, 11.02, 11.05, 11.1, 11.15, 11.2])
    a.yaxis.set_major_formatter(FormatStrFormatter('%.7f'))
    a.set_ylabel('side (square-root scale above the record)'); a.set_xticks([])
    for f, c in col.items():
        if f in fam: a.scatter([], [], color=c, label=flab(f))
    a.scatter([], [], color=(.6, .6, .6), label='other families')
    a.legend(fontsize=7, loc='upper left')
    a.set_title(f'Disconnectivity graph of the {n} tour minima: two groups join at the smallest expansion (walk_carrot upper '
                f'bound) linking them; {len(E)} of {n * (n - 1) // 2} pairs computed', fontsize=9)
    plt.tight_layout(); plt.savefig(f'{out}/disconn.png', dpi=80); plt.close()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--dir', default=f'{RUN}/map'); a = ap.parse_args()
    M = json.load(open(f'{a.dir}/minima.json')); D = np.load(f'{a.dir}/fpD.npy')
    fam = [tuple(m['fam']) for m in M]; S = np.array([float(m['S']) for m in M])
    nn = np.argsort(D, 1)[:, 1:6]
    print('5-NN same family:', round(np.mean([[fam[j] == fam[i] for j in nn[i]] for i in range(len(M))]), 3))
    disconn(a.dir, famcolours(fam)); families(D, fam, S, a.dir); layouts(D, fam, S, a.dir)
    print('->', a.dir)


if __name__ == '__main__':
    main()
