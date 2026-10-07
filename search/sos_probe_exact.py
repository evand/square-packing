#!/usr/bin/env python3
"""sos_probe_exact.py -- round a float infeasibility certificate to an EXACT rational identity (2026-10-03).

Same arguments as sos_probe_sdp.py (infeasibility mode, fixed T), plus --bits K.

    1. solve for lam* ; re-solve with lam = lam*/2 maximising the common interior margin mu;
    2. round every Gram matrix / scalar to the grid 2^-K;
    3. Peyrl-Parrilo projection: correct x by A_R^T y with (A_R A_R^T) y = residual_R solved EXACTLY
       (python-flint), R a maximal independent set of monomial rows;
    4. re-check the identity  -lam = sum_blocks <Gram, basis-products * g> + sum_scalars tau g_a g_b
       in exact rational arithmetic on ALL monomial rows (normal form mod s_i^2 + c_i^2 - 1), and
       every Gram matrix PSD by exact LDL^T over Q; scalars >= 0;
    5. independent float check: evaluate the certificate at random points of the variety using the
       UNREDUCED generators.

A pass is a proof (modulo the encoding in sos_probe_enc.py) that the separation type has no
configuration with container side T.
"""
import json
import math
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np
import scipy.sparse as sps

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sos_probe_sdp as S                                               # noqa: E402

import flint                                                            # noqa: E402


def ldl_psd(Q):
    """exact PSD test of a symmetric rational matrix (list of lists of flint.fmpq).
    returns (ok, min positive pivot, rank)."""
    k = len(Q)
    M = [row[:] for row in Q]
    minpiv, rank = None, 0
    for i in range(k):
        d = M[i][i]
        if d < 0:
            return False, None, rank
        if d == 0:
            if any(M[i][j] != 0 for j in range(i + 1, k)):
                return False, None, rank
            continue
        rank += 1
        minpiv = d if minpiv is None or d < minpiv else minpiv
        inv = 1 / d
        for j in range(i + 1, k):
            if M[i][j] == 0:
                continue
            f = M[j][i] * inv
            for l in range(i + 1, k):
                if M[i][l] != 0:
                    M[j][l] -= f * M[i][l]
    return True, minpiv, rank


def to_fmpq(x, K):
    return flint.fmpq(int(round(x * 2 ** K)), 2 ** K)


def exact_round(P, vals, lam, bits, bound, out):
    """round vals, project exactly onto the identity, verify; fills `out`, returns rounded blocks."""
    t1 = time.time()
    # exact data
    emats = P.build(exact=True)
    Mrows = len(P.mono_index)
    zero = (0,) * P.nv
    r0 = P.mono_index[zero]
    # rounded variables, flattened in block order
    xs, blocks = [], []
    for (kind, k, rows, cols, v, name), val in zip(emats, vals):
        if kind == 'psd':
            Q = np.asarray(val)
            Q = (Q + Q.T) / 2
            R = [[to_fmpq(Q[i, j], bits) for j in range(k)] for i in range(k)]
            blocks.append((kind, k, name, R))
        else:
            R = [max(flint.fmpq(0), to_fmpq(x, bits)) for x in np.asarray(val).ravel()] if k else []
            blocks.append((kind, k, name, R))

    def residual():
        res = [flint.fmpq(0)] * Mrows
        res[r0] = -flint.fmpq(lam.numerator, lam.denominator)
        if bound:
            res[P.mono_index[tuple(1 if i == 0 else 0 for i in range(P.nv))]] += 1
        for (kind, k, rows, cols, v, name), (_k2, _kk, _nm, R) in zip(emats, blocks):
            for rr, cc, vv in zip(rows, cols, v):
                x = R[cc % k][cc // k] if kind == 'psd' else R[cc]
                if x != 0:
                    res[rr] -= flint.fmpq(vv.numerator, vv.denominator) * x
        return res          # res = (-lam) - sum(...): must become 0

    res = residual()
    out['resid_rounded_max'] = float(max(abs(float(x)) for x in res))

    # independent rows (float), then exact projection
    fm = P.build(exact=False)
    Afull = sps.hstack([m[2] for m in fm]).tocsr()
    import scipy.linalg as sla
    AAT = (Afull @ Afull.T).toarray()
    _q, rr_, piv = sla.qr(AAT, mode='economic', pivoting=True)
    dg = np.abs(np.diag(rr_))
    rank = int((dg > 1e-12 * dg[0]).sum())
    Rrows = sorted(piv[:rank].tolist())
    out['rows'] = Mrows
    out['rank'] = rank
    # exact A restricted to rows R, as sparse dict per column (global column index)
    colmap = []          # global col -> (block index, local col)
    offs = []
    off = 0
    for bi, (kind, k, rows, cols, v, name) in enumerate(emats):
        offs.append(off)
        off += k * k if kind == 'psd' else k
    ncol = off
    Rpos = {r: i for i, r in enumerate(Rrows)}
    AR = {}              # (i_row_in_R, global col) -> Fr
    for bi, (kind, k, rows, cols, v, name) in enumerate(emats):
        for rr, cc, vv in zip(rows, cols, v):
            if rr in Rpos:
                key = (Rpos[rr], offs[bi] + cc)
                AR[key] = AR.get(key, 0) + vv
    # AAT exact
    bycol = {}
    for (i, c), vv in AR.items():
        bycol.setdefault(c, []).append((i, vv))
    nR = len(Rrows)
    AATx = [[flint.fmpq(0)] * nR for _ in range(nR)]
    for c, lst in bycol.items():
        for i, vi in lst:
            fi = flint.fmpq(vi.numerator, vi.denominator)
            for j, vj in lst:
                AATx[i][j] += fi * flint.fmpq(vj.numerator, vj.denominator)
    t2 = time.time()
    Mq = flint.fmpq_mat(nR, nR, [x for row in AATx for x in row])
    rhs = flint.fmpq_mat(nR, 1, [res[r] for r in Rrows])
    y = Mq.solve(rhs)
    t3 = time.time()
    yv = [y[i, 0] for i in range(nR)]
    # x += A_R^T y   (res = b - A x, so A dx = res -> dx = A^T y with AAT y = res)
    for c, lst in bycol.items():
        dx = flint.fmpq(0)
        for i, vi in lst:
            dx += flint.fmpq(vi.numerator, vi.denominator) * yv[i]
        if dx == 0:
            continue
        bi = max(b for b in range(len(offs)) if offs[b] <= c)
        kind, k, name, R = blocks[bi]
        loc = c - offs[bi]
        if kind == 'psd':
            R[loc % k][loc // k] += dx
        else:
            R[loc] += dx
    res2 = residual()
    nonzero = sum(1 for x in res2 if x != 0)
    out['identity_exact'] = (nonzero == 0)
    out['identity_nonzero_rows'] = nonzero
    # PSD checks
    bad, minpiv, sizes, maxden = 0, None, [], 1
    for kind, k, name, R in blocks:
        if kind == 'psd':
            sym = all(R[i][j] == R[j][i] for i in range(k) for j in range(i))
            ok, mp, rk = ldl_psd(R)
            if not (ok and sym):
                bad += 1
            elif mp is not None:
                minpiv = mp if minpiv is None or mp < minpiv else minpiv
            for row in R:
                for x in row:
                    maxden = max(maxden, int(x.q))
        else:
            if any(x < 0 for x in R):
                bad += 1
            for x in R:
                maxden = max(maxden, int(x.q))
    out['blocks_not_psd'] = bad
    out['min_pivot'] = float(minpiv) if minpiv is not None else None
    out['max_denominator_bits'] = maxden.bit_length()
    out['lam_exact'] = str(lam)
    out['t_exact'] = round(time.time() - t1, 1)
    out['t_solve_AAT'] = round(t3 - t2, 1)

    return blocks


def main():
    ap = S.parser()
    ap.add_argument('--bits', type=int, default=30)
    ap.add_argument('--lamfrac', type=float, default=0.5)
    a = ap.parse_args()
    t0 = time.time()
    n, typ, G, names, idx = S.setup(a)
    P = S.make_program(n, G, names, a.order, a.mult, a.prod, a.pdeg, a.csp, a.csp0)
    info1, _ = S.solve_infeas(P, a.solver)
    lam_star = info1['lam']
    out = dict(cfg=a.cfg, sub=idx, classes=a.classes, orthant=a.orthant, order=a.order, T=a.T,
               lam_star=lam_star, status1=info1['status'])
    if not lam_star or lam_star <= 0:
        out['verdict'] = 'no certificate (lam* <= 0)'
        print(json.dumps(out)); return
    lam = Fr(lam_star * a.lamfrac).limit_denominator(10 ** 6)
    info2, vals = S.solve_infeas(P, a.solver, lamfix=float(lam))
    out.update(mu=info2['mu'], status2=info2['status'])
    t1 = time.time()

    out['t_sdp'] = round(time.time() - t0, 1)
    blocks = exact_round(P, vals, lam, a.bits, False, out)
    # independent float check with UNREDUCED generators at random variety points
    rng = np.random.default_rng(7)
    worst = 0.0
    gens = {str(nm): g for g, nm in zip(G, names)}
    for _ in range(20):
        z = np.zeros(P.nv); z[0] = a.T
        for i in range(n):
            th = rng.uniform(0, math.pi / 4) * (1 if not a.orthant or a.orthant[i] != '-' else -1)
            z[1 + 4 * i], z[2 + 4 * i] = rng.uniform(0, 3, 2)
            z[3 + 4 * i], z[4 + 4 * i] = math.cos(th), math.sin(th)
        if a.classes:            # respect the class substitution
            lab = a.classes.split(',')
            rep = {}
            for i, l in enumerate(lab):
                if l == '0':
                    z[3 + 4 * i], z[4 + 4 * i] = 1.0, 0.0
                elif l in rep:
                    z[3 + 4 * i], z[4 + 4 * i] = z[3 + 4 * rep[l]], z[4 + 4 * rep[l]]
                else:
                    rep[l] = i
        tot = float(lam)
        for (name_b, B, g), (kind, k, nm, R) in zip(P.blocks, [b for b in blocks if b[0] == 'psd']):
            v = np.array([np.prod([z[i] ** e for i, e in enumerate(m) if e]) for m in B])
            Qf = np.array([[float(x) for x in row] for row in R])
            tot += float(v @ Qf @ v) * g.ev(z)
        sc = [b for b in blocks if b[0] == 'scal'][0][3]
        for (nm, g), x in zip(P.scal, sc):
            tot += float(x) * g.ev(z)
        worst = max(worst, abs(tot))
    out['float_check_max'] = worst
    out['verdict'] = ('EXACT CERTIFICATE' if out['identity_exact'] and bad == 0 else 'rounding failed')
    print(json.dumps(out), flush=True)


if __name__ == '__main__':
    main()
