#!/usr/bin/env python3
"""Stratified extrapolation of a Lean-tree census sample to the whole D4 region, priced by costmodel.MODELS.

  python3 search/golf/estimate.py TAG [TAG ...]     (reads search/golf/data/lc_TAG_{cells,roots}.jsonl)

Light cells (strata S1-S5 of search/golf/data/strata.json) are sampled as whole 1/10 cells (8 u-bins); the 8 heavy
cells (S6) as zm_mixed roots (1/20 x 1/20 x 1/32 in u) in strata R1-R5 of strata_roots.json.  Expansion estimator
N/n * sum per stratum.  A root that timed out is counted as its partial tree (generator running stats): a LOWER
bound for that root, flagged in the output.
"""
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from costmodel import MODELS, cpu, feats, RSS_BASE, RSS_PER_DIGIT  # noqa: E402

KEYS = ['Z', 'E', 'adm', 'ch', 'L', 'T', 'S', 'digits', 'boxes', 'cpu', 'ptclaims', 'pairs', 'segclaims', 'UNCERT']
D = 'search/golf/data'


def partial_feats(r):
    p = r['partial']
    z = p.get('ADM', 0) + p.get('CHAIN1', 0) + p.get('CHAIN2', 0) + p.get('PIECE', 0)
    # leafL counts Z leaves with an L-block; the completed heavy roots average ~1.5 lines per such leaf
    return dict(Z=z, E=p.get('E', 0), adm=p.get('ptclaims', 0), ch=0, L=1.5 * p.get('leafL', 0), T=p.get('leafT', 0),
                S=0, digits=None, boxes=p.get('boxes', 0), cpu=r.get('cpu', 0), UNCERT=p.get('UNCERT', 0),
                ptclaims=p.get('ptclaims', 0), pairs=0, segclaims=0)


def load(path):
    rows = {}
    if not os.path.exists(path):
        return rows
    for l in open(path):
        r = json.loads(l)
        k = r['cell'] if isinstance(r['cell'], str) else ','.join(map(str, r['cell']))
        if k in rows and not (rows[k].get('timeout') or rows[k].get('error')):
            continue
        rows[k] = r
    return rows


def strat(units, rows, tot, var, flags, per_unit):
    for name, N, samp in units:
        fs = []
        for k in samp:
            r = rows.get(k)
            if r is None:
                flags.append((name, k, 'absent'))
                continue
            if r.get('timeout'):
                flags.append((name, k, f"timeout at {r['cpu']:.0f} s: partial tree counted (lower bound)"))
                f = partial_feats(r)
            elif r.get('error'):
                flags.append((name, k, r['error']))
                continue
            else:
                f = feats(r)
            fs.append(f)
            per_unit[k] = f
        if not fs:
            continue
        n = len(fs)
        for key in KEYS:
            tot[key] += N * sum((f[key] or 0) for f in fs) / n
        for mk, m in MODELS.items():
            cs = [cpu(m, f) for f in fs]
            mu = sum(cs) / n
            if 1 < n < N:
                s2 = sum((c - mu) ** 2 for c in cs) / (n - 1)
                var[mk] += N * N * (1 - n / N) * s2 / n


def estimate(tag, quiet=False):
    st = json.load(open(f'{D}/strata.json'))
    sr = json.load(open(f'{D}/strata_roots.json'))
    tot, var, flags, per = defaultdict(float), {k: 0.0 for k in MODELS}, [], {}
    strat([(n, v['N'], v['sample']) for n, v in st.items() if n != 'S6'], load(f'{D}/lc_{tag}_cells.jsonl'),
          tot, var, flags, per)
    light = dict(tot)
    strat([(n, v['N'], [x['label'] for x in v['sample']]) for n, v in sr.items()], load(f'{D}/lc_{tag}_roots.jsonl'),
          tot, var, flags, per)
    heavy = {k: tot[k] - light.get(k, 0) for k in tot}
    out = dict(tag=tag, tot=dict(tot), light=light, heavy=heavy, flags=flags, per=per)
    for mk, m in MODELS.items():
        out[mk] = dict(cpu_h=cpu(m, tot) / 3600, se_h=var[mk] ** 0.5 / 3600, heavy_h=cpu(m, heavy) / 3600)
    if not quiet:
        t = tot
        print(f"== {tag}: estimated whole D4 region")
        print(f"  Z leaves {t['Z']:.0f} (heavy cells {heavy.get('Z', 0):.0f}), E {t['E']:.0f}, "
              f"boxes {t['boxes']:.0f}, generator CPU {t['cpu'] / 3600:.1f} h")
        print(f"  point entries {(t['adm'] + t['ch']) / 1e6:.1f} M (ADM {t['adm'] / 1e6:.1f} M, chain "
              f"{t['ch'] / 1e6:.2f} M), {(t['adm'] + t['ch']) / max(t['Z'], 1):.0f} per Z leaf; "
              f"L-lines {t['L']:.0f}, T-groups {t['T']:.0f}")
        if t['digits']:
            print(f"  digits (completed units only) {t['digits'] / 1e6:.1f} M; at 15 GB/process "
                  f"~{t['digits'] / ((15 - RSS_BASE) / RSS_PER_DIGIT):.0f} part files")
        for mk in MODELS:
            print(f"  kernel CPU [{mk}]: {out[mk]['cpu_h']:.0f} h (heavy cells {out[mk]['heavy_h']:.0f} h; "
                  f"sampling s.e. {out[mk]['se_h']:.0f} h)")
        for f in flags:
            print('   ', f)
    return out


def paired(t0, t1):
    """per-unit ratios of the central-model cost, candidate t1 vs t0, on units completed in both."""
    a, b = estimate(t0, True)['per'], estimate(t1, True)['per']
    m = MODELS['central']
    print(f"== paired {t1} / {t0} (central model, per sampled unit)")
    for k in sorted(set(a) & set(b)):
        ca, cb = cpu(m, a[k]), cpu(m, b[k])
        if ca > 1:
            print(f"  {k:8s} {ca:9.0f} s -> {cb:9.0f} s  x{cb / ca:.2f}   Z {a[k]['Z']:.0f} -> {b[k]['Z']:.0f}, "
                  f"pts/leaf {a[k]['adm'] / max(a[k]['Z'], 1):.0f} -> {b[k]['adm'] / max(b[k]['Z'], 1):.0f}")


if __name__ == '__main__':
    if sys.argv[1] == 'paired':
        paired(sys.argv[2], sys.argv[3])
    else:
        for tag in sys.argv[1:]:
            estimate(tag)
