"""`pk.py prune --dry-run`: what pk/policy.py would keep, per tier / reason / n, with a threshold sweep and the lineage
check (how far above the record were the ancestors of records we found?).  Read-only.  Applying (stubbing rows) is a
separate step, not built yet: pick thresholds from this report first."""
from __future__ import annotations
import collections, math
import numpy as np
from . import policy


def _rows(st):
    cur = st.db.execute('SELECT id, n, s, status, grid, parent, length(coords) FROM packing')
    I, N, S, ST, G, P, B = [], [], [], [], [], [], []
    for r in cur:
        I.append(r[0]); N.append(r[1]); S.append(r[2]); ST.append(r[3]); G.append(r[4] or 0)
        P.append(r[5] if r[5] is not None else -1); B.append(r[6] or 0)
    return (np.array(I), np.array(N), np.array(S), np.array(ST, dtype=object), np.array(G), np.array(P), np.array(B))


def report(st, cfg=None, sweep=(0.05, 0.1, 0.2, 0.5, 1, 2, 4), log=print):
    cfg = cfg or policy.load()
    I, N, S, ST, G, P, B = _rows(st)
    F = policy.references(st)
    prot = policy.protected_ids(st, cfg)
    isprot = np.fromiter((i in prot for i in I), bool, len(I))
    Fn = np.array([F.get(int(n), np.nan) for n in N])
    A = S * S - Fn * Fn                                      # area excess (unit squares)
    screened = np.isin(ST, ['screened', 'raw'])
    rng = np.array([(n <= cfg['nmax']) or (n in cfg['extra_n']) or (n in cfg['targets']) for n in N])
    inf = ST == 'infeasible'
    grid = G > 0
    gb = lambda m: B[m].sum() / 1e9
    log(f'store: {len(I):,} rows, coordinates {B.sum() / 1e9:.2f} GB; protected {isprot.sum():,} ({gb(isprot):.2f} GB)')
    log(f'cfg: nmax {cfg["nmax"]}, full {cfg["full"]}, screened_full {cfg["screened_full"]}, drop (intake) {cfg["drop"]}, '
        f'targets {sorted(cfg["targets"])}')

    def tiers(full, sfull):
        tier = np.full(len(I), 'stub', dtype=object); why = np.full(len(I), 'far', dtype=object)
        lim = np.where(screened, sfull, full)
        near = ~np.isnan(A) & (A <= lim)
        for m, t, w in ((near, 'full', 'near'), (np.isnan(A), 'full', 'no reference'), (grid, 'stub', 'grid line'),
                        (~rng, 'stub', 'n range'), (inf, 'stub', 'infeasible'), (isprot, 'full', 'protected')):
            tier[m] = t; why[m] = w                          # later masks win: same order as policy.classify, reversed
        # per-n target thresholds
        for n, t in cfg['targets'].items():
            m = (N == n) & ~isprot & ~inf & ~grid
            l2 = np.where(screened[m], t.get('screened_full', sfull), t.get('full', full))
            ok = A[m] <= l2
            tier[np.where(m)[0][ok]] = 'full'; why[np.where(m)[0][ok]] = 'near (target)'
        return tier, why

    tier, why = tiers(cfg['full'], cfg['screened_full'])
    log('\nby reason (rows, coordinate GB):')
    for w in ('protected', 'near', 'near (target)', 'no reference', 'far', 'grid line', 'n range', 'infeasible'):
        m = why == w
        if m.any():
            log(f'  {w:14s} -> {tier[m][0]:4s}  {m.sum():>10,}  {gb(m):6.2f} GB')
    full = tier == 'full'
    log(f'  total full {full.sum():,} ({gb(full):.2f} GB), stub {(~full).sum():,} ({gb(~full):.2f} GB freed)')
    log('\nby status (full / stub rows):')
    for s_ in collections.Counter(ST).most_common():
        m = ST == s_[0]
        log(f'  {s_[0]:11s} {(m & full).sum():>10,} / {(m & ~full).sum():>10,}')

    log('\nsweep (polished limit `full`; screened kept at min(full, screened_full)): full rows / GB')
    for f in sweep:
        t, _ = tiers(f, min(f, cfg['screened_full']))
        m = t == 'full'
        log(f'  full <= {f:<5g} {m.sum():>10,} rows  {gb(m):6.2f} GB')
    log('  (screened sweep at full = cfg)')
    for f in (0.01, 0.05, 0.1, 0.5, 1):
        t, _ = tiers(cfg['full'], f)
        m = t == 'full'
        log(f'  screened <= {f:<5g} {m.sum():>10,} rows  {gb(m):6.2f} GB')

    log('\nn with the most full rows (cfg):')
    c = collections.Counter(N[full].tolist())
    log('  ' + ', '.join(f'{n}: {k:,}' for n, k in c.most_common(12)))

    # lineage: ancestors of our records (rows at their n's reference, found by us = a sighting not from a public source)
    pos = {int(i): k for k, i in enumerate(I)}
    ours = {r[0] for r in st.db.execute("SELECT DISTINCT packing FROM sighting WHERE finder='ours'")}
    rec = [k for k in range(len(I)) if not np.isnan(A[k]) and abs(A[k]) <= 1e-9 * max(1.0, S[k]) and int(I[k]) in ours
           and P[k] >= 0]
    anc = []
    for k in rec:
        seen, j, d = set(), pos.get(int(P[k])), 1
        while j is not None and j not in seen:
            seen.add(j); anc.append((int(N[j]), d, float(A[j]), str(ST[j]), str(tier[j]), bool(isprot[j])))
            j = pos.get(int(P[j])) if P[j] >= 0 else None; d += 1
    log(f'\nlineage: {len(rec)} of our record-side rows have a stored parent chain; {len(anc)} ancestors')
    if anc:
        a = np.array([x[2] for x in anc])
        log('  ancestor area excess quantiles (50/90/99/max): ' + ' / '.join(f'{v:.3g}' for v in np.quantile(a, [.5, .9, .99, 1])))
        for f in sweep:
            log(f'  full <= {f:<5g} keeps coordinates of {np.mean(a <= f) * 100:5.1f} % of ancestors')
        log(f'  under cfg: {sum(x[4] == "stub" for x in anc)} ancestors would be stubs (lineage kept, coordinates dropped)')
    from . import archive
    back = archive.backing_ids()
    bk = np.fromiter((int(i) in back for i in I), bool, len(I))
    log(f'\narchive files removed by dedupe are backed by {bk.sum():,} rows; {(bk & ~full).sum():,} of them would be stubs '
        f'under cfg (their run files then lose their coordinates: pruning must restore or accept that per run)')
    return dict(tier=tier, why=why)
