#!/usr/bin/env python3
"""Per-cell census of the Lean mixed tree (the search of `lean/scripts/gen_zmmtree.py`, unchanged, used as a
library), for cost estimation.  No Lean is emitted.

For each requested root cell (pitch 1/10 over the D4 region, all 8 u-bins, as gen_zmmtree) it runs the
generator's own search (exact mirror included) and walks the resulting tree, writing one JSON line with the
leaf census and the per-leaf features the kernel cost depends on: Z leaves by point kind (none / ADM / CHAIN1 /
CHAIN2), claimed point entries (ADM witnesses, chain entries, pivots), S-blocks, T-groups and coupled pairs,
L-blocks and L-lines, segment claims, E leaves, C nodes, splits, UNCERT boxes, max depth, CPU.

Usage (from the repo root):
  taskset -c 6-11 python3 search/golf/lean_census.py COVER --n 21 --cells "13,13;14,12" \
      --nproc 6 --timeout 1800 --out search/golf/data/lc_current.jsonl
  --cells all        every cell;  --cells-file F  one "i,j" per line
Cells already in --out are skipped (resumable).
"""
import argparse
import json
import os
import signal
import sys
import time
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'lean', 'scripts'))
sys.path.insert(0, os.path.join(HERE, '..'))
import gen_zmmtree as GM  # noqa: E402
import gen_zmtree as G  # noqa: E402
import zm_mixed as ZX  # noqa: E402
import mixed_cover as MC  # noqa: E402

_S = None


class Timeout(Exception):
    pass


def _alarm(sig, frm):
    raise Timeout()


def walk(t, acc, depth=0):
    k = t[0]
    if k == 'E':
        acc['E'] += 1
        acc['digits'] += 1
    elif k == 'UNC':
        acc['UNCERT'] += 1
    elif k == 'C':
        acc['C'] += 1
        acc['digits'] += 3
        walk(t[2], acc, depth)
    elif k == 'Z':
        claims, sb, tg, lb, val = t[1]
        acc['Z'] += 1
        # base-2^20 digits exactly as gen_zmmtree.replay/leaf_digits_pc emit them
        dg = 1 + len(GM.leaf_digits_pc(claims, sb, tg, lb))
        if len(t) > 2:
            _, chA, chB, emp, ents = t[2]
            dg += 1 + sum(1 if tg_[0] == 4 else 3 for _, tg_ in ents) + 2 + 2 * (len(chA) + len(chB)) \
                + 1 + len(emp)
        else:
            dg += 1
        acc['digits'] += dg
        acc['segclaims'] += len(claims)
        acc['Sblk'] += len(sb)
        acc['Tgrp'] += len(tg)
        acc['pairs'] += sum(len(g[8]) for g in tg)
        acc['Lblk'] += len(lb)
        acc['Llines'] += sum(len(B[0]) for B in lb)
        if len(t) > 2:
            kind, chA, chB, emp, ents = t[2]
            acc['pt_' + kind] += 1
            nch = sum(1 for _, tg_ in ents if tg_[0] != 4)
            acc['adm_ents'] += len(ents) - nch
            acc['chain_ents'] += nch
            acc['pivots'] += len(chA) + len(chB)
            acc['ptclaims_hist'].append(len(ents))
        else:
            acc['pt_none'] += 1
            acc['ptclaims_hist'].append(0)
    else:
        acc['split'] += 1
        acc['digits'] += 1
        walk(t[1], acc, depth + 1)
        walk(t[2], acc, depth + 1)


def work(cell_root):
    cell, root, timeout, nub = cell_root
    t0 = time.process_time()
    for k in _S.stat:
        _S.stat[k] = 0
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(timeout)
    try:
        def rec(box, n):
            if n == 1:
                return _S._build(box, 0)
            l, r = G.split(box, 2)
            return (2, rec(l, n // 2), rec(r, n // 2))
        tree = rec(root, nub)
        signal.alarm(0)
    except Timeout:
        return dict(cell=cell, timeout=True, cpu=time.process_time() - t0, partial=dict(_S.stat))
    except RecursionError:
        signal.alarm(0)
        return dict(cell=cell, error='recursion', cpu=time.process_time() - t0)
    from collections import Counter
    acc = Counter()
    acc['ptclaims_hist'] = []
    walk(tree, acc)
    hist = acc.pop('ptclaims_hist')
    out = dict(cell=cell, cpu=time.process_time() - t0, maxdepth=_S.stat['maxdepth'],
               boxes=_S.stat['boxes'], inherit=_S.stat['inherit'], **acc)
    out['ptclaims_max'] = max(hist) if hist else 0
    out['ptclaims_sum'] = sum(hist)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cover')
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--cells', default=None)
    ap.add_argument('--cells-file', default=None)
    ap.add_argument('--roots-file', default=None,
                    help='zm_mixed-style roots, one per line: LABEL x0 x1 y0 y1 u0 u1 (fractions); no u-bin split')
    ap.add_argument('--nproc', type=int, default=1)
    ap.add_argument('--timeout', type=int, default=1800)
    ap.add_argument('--depth', type=int, default=30)
    ap.add_argument('--out', required=True)
    ap.add_argument('--pitch', default='1/10')
    ap.add_argument('--ubins', type=int, default=8)
    args = ap.parse_args()
    cv = MC.load(args.cover)
    MC.validate(cv)
    s = F(cv['s_num'], cv['s_den'])
    D, W = cv['D'], cv['W']
    Pts = sorted(tuple(p) for p in cv['points'] if p[2] > 0)
    Sg = []
    for (X0, Y0, X1, Y1, w) in cv['segments']:
        if w == 0:
            continue
        a, b = (X0, Y0), (X1, Y1)
        if b < a:
            a, b = b, a
        Sg.append((a[0], a[1], b[0], b[1], w))
    Sg.sort()
    tot = sum(p[2] for p in Pts) + sum(e[4] for e in Sg)
    assert tot < args.n * W
    S = 2 ** 12
    Q = D * S
    Mq = int(s * D)
    M = Mq * S
    R = 2 ** 32
    ctx = G.Ctx(S, Q, R, W, 3 * Q // 4)
    cells = int((s / 2) / F(args.pitch))
    cov = ZX.Cover(cv)
    assert cov.symmetric_d4()
    mchk = ZX.MixedChecker(cov, max_depth=args.depth, cert_mode=True)
    P = GM.Pieces(ctx, cov, Sg)
    global _S
    _S = GM.MSearch(ctx, P, mchk, max_depth=args.depth, use_thr=True, use_lin=True, pts=Pts, D=D, m=s)
    _S.ubins = args.ubins
    roots = G.roots_of(ctx, M, cells, args.ubins)
    wd = (M // 2) // cells
    byc = {(r[0] // wd, r[2] // wd): r for r in roots}
    if args.roots_file:
        want, byc = [], {}
        for l in open(args.roots_file):
            if not l.strip():
                continue
            lab, *fr = l.split()
            x0, x1, y0, y1, u0, u1 = (F(v) for v in fr)
            box = tuple(int(v * Q) for v in (x0, x1, y0, y1)) + (int(u0 * R), int(u1 * R))
            assert all(F(b, Q) == v for b, v in zip(box[:4], (x0, x1, y0, y1)))
            assert F(box[4], R) == u0 and F(box[5], R) == u1
            want.append(lab)
            byc[lab] = box
        nub = 1
    elif args.cells == 'all':
        want = sorted(byc)
    elif args.cells_file:
        want = [tuple(map(int, l.split(','))) for l in open(args.cells_file) if l.strip()]
    else:
        want = [tuple(map(int, c.split(','))) for c in args.cells.split(';')]
    if not args.roots_file:
        nub = args.ubins
    done = set()
    if os.path.exists(args.out):
        for l in open(args.out):
            r = json.loads(l)
            if not r.get('timeout') and not r.get('error'):
                done.add(r['cell'] if isinstance(r['cell'], str) else tuple(r['cell']))
    todo = [(c, byc[c], args.timeout, nub) for c in want if c not in done]
    print(f"{len(todo)} cells to do ({len(done)} done), {len(Pts)} points, {len(Sg)} segments",
          file=sys.stderr, flush=True)
    if args.nproc > 1:
        import multiprocessing as mp
        pool = mp.get_context('fork').Pool(args.nproc)
        it = pool.imap_unordered(work, todo, chunksize=1)
    else:
        it = map(work, todo)
    with open(args.out, 'a') as fo:
        for r in it:
            fo.write(json.dumps(r) + '\n')
            fo.flush()
            print(json.dumps(r), file=sys.stderr, flush=True)


if __name__ == '__main__':
    main()
