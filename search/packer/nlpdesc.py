#!/usr/bin/env python3
"""Second-order descent for slow-valley points of min S (squares of side 1 in a box [0, S]^2).

Fixed-line NLP: for each near pair, one edge line of one square (the best SAT separator at the current point) must keep
all 4 corners of the other square outside: smooth, and a sufficient condition for disjointness, so every NLP-feasible
point is a packing (to IPOPT's tolerance; the repair at the end scales about the origin to clear residual overlap).
Walls: every corner within [0, S].  Objective S.  Exact Hessians.  Trust box |v - v_k| <= R.
Outer loop: re-pick lines at the new point; when the NLP stops moving, try the line switches the exact corner-corner MILP
(exactsolve.vv_branch_test, via descend.milp_dir) asks for, if it reports a descent.

  ./ipenv.sh python3 nlpdesc.py   (ipenv.sh: nix python with cyipopt) in.txt out.txt [--rounds 30] [--R 0.02]
"""
import argparse, math, sys, time
import numpy as np
import cyipopt

H = 0.5
U = np.array([(H, H), (-H, H), (-H, -H), (H, -H)])
E4 = np.array([(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)])


def load(path):
    L = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    n, s = int(L[0][0]), float(L[0][1])
    return s, np.array([(float(a), float(b), math.radians(float(c))) for a, b, c, *_ in L[1:n + 1]])


def save(path, s, sq):
    with open(path, 'w') as f:
        f.write(f'{len(sq)} {float(s)!r}\n')
        for x, y, t in sq:
            f.write(f'{float(x)!r} {float(y)!r} {math.degrees(t) % 90!r}\n')


def rot(t):
    c, s = math.cos(t), math.sin(t)
    return np.array([[c, -s], [s, c]])


def line_gap(sq, s_, o, k):
    """min over corners of o of the signed distance outside edge k of square s_ (minus 1/2)."""
    xs, ys, ts = sq[s_]; xo, yo, to = sq[o]
    n = rot(ts) @ E4[k]
    P = np.array([xo, yo]) + U @ rot(to).T
    return float(((P - [xs, ys]) @ n).min() - H)


def best_line(sq, i, j):
    best = None
    for (a, b) in ((i, j), (j, i)):
        for k in range(4):
            g = line_gap(sq, a, b, k)
            if best is None or g > best[0]:
                best = (g, a, b, k)
    return best


def sat_gap(sq, i, j):
    return best_line(sq, i, j)[0]


def worst(s, sq, near):
    w = -1.0
    for i, (x, y, t) in enumerate(sq):
        e = H * (abs(math.cos(t)) + abs(math.sin(t)))
        w = max(w, e - x, x + e - s, e - y, y + e - s)
    for (i, j) in near:
        w = max(w, -sat_gap(sq, i, j))
    return w


def near_pairs(sq, cut=1.5):
    n = len(sq)
    d = np.hypot(sq[:, 0][:, None] - sq[:, 0][None], sq[:, 1][:, None] - sq[:, 1][None])
    return [(i, j) for i in range(n) for j in range(i + 1, n) if d[i, j] < cut]


class Prob:
    """Constraints: ('P', s_, o, k, a): corner a of o outside edge k of s_; ('W', i, a, axis, side)."""

    def __init__(self, n, cons):
        self.n, self.cons = n, cons
        self.nv = 3 * n + 1
        rows, cols = [], []
        hr, hc = [], []
        self.hidx = {}
        for r, c in enumerate(cons):
            vs = self.vars(c)
            for v in vs:
                rows.append(r); cols.append(v)
            for a in vs:
                for b in vs:
                    if a >= b and (a, b) not in self.hidx:
                        self.hidx[(a, b)] = len(hr); hr.append(a); hc.append(b)
        self.jr, self.jc = np.array(rows), np.array(cols)
        self.hr, self.hc = np.array(hr), np.array(hc)

    def vars(self, c):
        if c[0] == 'P':
            _, s_, o, k, a = c
            return [3 * s_, 3 * s_ + 1, 3 * s_ + 2, 3 * o, 3 * o + 1, 3 * o + 2]
        _, i, a, ax, side = c
        return [3 * i + ax, 3 * i + 2] + ([3 * self.n] if side == 1 else [])

    def objective(self, v): return v[-1]

    def gradient(self, v):
        g = np.zeros(self.nv); g[-1] = 1.0; return g

    def _eval(self, v, c, order):
        if c[0] == 'P':
            _, s_, o, k, a = c
            xs, ys, ts, xo, yo, to = v[3 * s_], v[3 * s_ + 1], v[3 * s_ + 2], v[3 * o], v[3 * o + 1], v[3 * o + 2]
            e = E4[k]; u = U[a]
            d = np.array([xo - xs, yo - ys])
            A = rot(-ts); Ap = np.array([[-math.sin(ts), math.cos(ts)], [-math.cos(ts), -math.sin(ts)]])
            ph = to - ts; Rp = rot(ph); Rq = rot(ph + math.pi / 2)
            val = e @ A @ d + e @ Rp @ u - H
            if order == 0:
                return val
            n = rot(ts) @ e                      # = A^T e
            dts = e @ Ap @ d - e @ Rq @ u
            dto = e @ Rq @ u
            grad = [-n[0], -n[1], dts, n[0], n[1], dto]
            if order == 1:
                return val, grad
            r2 = -(e @ Rp @ u)                   # e R''(ph) u
            eAp = e @ Ap
            Hm = np.zeros((6, 6))
            Hm[2, 2] = -(e @ A @ d) + r2         # e A'' d + e R'' u
            Hm[5, 5] = r2
            Hm[2, 5] = Hm[5, 2] = -r2
            Hm[2, 0] = Hm[0, 2] = -eAp[0]; Hm[2, 1] = Hm[1, 2] = -eAp[1]
            Hm[2, 3] = Hm[3, 2] = eAp[0]; Hm[2, 4] = Hm[4, 2] = eAp[1]
            return val, grad, Hm
        _, i, a, ax, side = c
        t = v[3 * i + 2]; ru = rot(t) @ U[a]; rq = rot(t + math.pi / 2) @ U[a]
        p = v[3 * i + ax] + ru[ax]
        sg = 1.0 if side == 0 else -1.0
        val = sg * p + (v[-1] if side == 1 else 0.0)
        if order == 0:
            return val
        grad = [sg, sg * rq[ax]] + ([1.0] if side == 1 else [])
        if order == 1:
            return val, grad
        Hm = np.zeros((len(grad), len(grad))); Hm[1, 1] = -sg * ru[ax]
        return val, grad, Hm

    def constraints(self, v): return np.array([self._eval(v, c, 0) for c in self.cons])

    def jacobian(self, v): return np.concatenate([self._eval(v, c, 1)[1] for c in self.cons])

    def jacobianstructure(self): return self.jr, self.jc

    def hessianstructure(self): return self.hr, self.hc

    def hessian(self, v, lam, objf):
        out = np.zeros(len(self.hr))
        for r, c in enumerate(self.cons):
            if lam[r] == 0.0:
                continue
            _, _, Hm = self._eval(v, c, 2)
            vs = self.vars(c)
            for p, a in enumerate(vs):
                for q, b in enumerate(vs):
                    if a >= b and Hm[p, q] != 0.0:
                        out[self.hidx[(a, b)]] += lam[r] * Hm[p, q]
        return out


def build(s, sq, R, lines_override=None, wall_cut=0.2):
    n = len(sq)
    cons = []
    near = near_pairs(sq)
    for (i, j) in near:
        g, a, b, k = best_line(sq, i, j)
        if lines_override and (i, j) in lines_override:
            a, b, k = lines_override[(i, j)]
        if g > 0.3:
            continue
        cons += [('P', a, b, k, c) for c in range(4)]
    for i, (x, y, t) in enumerate(sq):
        for a in range(4):
            p = np.array([x, y]) + rot(t) @ U[a]
            for ax in range(2):
                if p[ax] < wall_cut: cons.append(('W', i, a, ax, 0))
                if s - p[ax] < wall_cut: cons.append(('W', i, a, ax, 1))
    return cons, near


def solve(s, sq, R, lines_override=None, tol=1e-13, iters=300, verbose=0):
    n = len(sq)
    cons, near = build(s, sq, R, lines_override)
    P = Prob(n, cons)
    v0 = np.concatenate([sq.ravel(), [s]])
    lb, ub = v0 - R, v0 + R
    nlp = cyipopt.Problem(n=P.nv, m=len(cons), problem_obj=P, lb=lb, ub=ub, cl=np.zeros(len(cons)),
                          cu=np.full(len(cons), 1e20))
    for k, val in (('tol', tol), ('constr_viol_tol', 1e-14), ('acceptable_tol', 1e-11), ('max_iter', iters),
                   ('print_level', verbose), ('mu_strategy', 'adaptive'), ('bound_push', 1e-12),
                   ('bound_frac', 1e-12), ('honor_original_bounds', 'yes'), ('sb', 'yes'),
                   ('nlp_scaling_method', 'none')):
        nlp.add_option(k, val)
    v, info = nlp.solve(v0)
    return v[-1], v[:-1].reshape(n, 3), info, near


def repair(s, sq, near):
    p = worst(s, sq, near)
    if p <= 0:
        return s, sq, p
    f = 1.0 + 1.01 * p
    for _ in range(40):
        s2, sq2 = s * f, sq.copy(); sq2[:, :2] *= f
        if worst(s2, sq2, near) <= 0:
            return s2, sq2, p
        f = 1.0 + 2.0 * (f - 1.0)
    return float('inf'), sq, p


def milp_switch(s, sq):
    """Exact corner-corner MILP at the point (descend.milp_dir); returns dS* and the direction."""
    sys.path.insert(0, '/home/evand/math/square-packing/public/search/exact')
    sys.path.insert(0, '/home/evand/math/square-packing/public/search/packer')
    import descend
    dS, z = descend.milp_dir(s, [tuple(r) for r in sq], 1e-7)
    return dS, z


def run(s, sq, rounds=30, R=0.02, log=print, use_milp=True):
    s0 = s
    near = near_pairs(sq)
    s, sq, _ = repair(s, sq, near)
    stall = 0
    for rd in range(rounds):
        t0 = time.time()
        s1, sq1, info, near = solve(s, sq, R)
        s2, sq2, p = repair(s1, sq1, near)
        st = info['status']
        log(f'round {rd}: S {s:.13f} -> nlp {s1:.13f} (status {st}, {info["status_msg"][:40]!r}) overlap {p:.1e} '
            f'-> {s2:.13f}  ({time.time() - t0:.1f}s)')
        if s2 < s - 1e-14:
            moved = s - s2
            s, sq = s2, sq2
            if moved < 1e-12:
                stall += 1
            else:
                stall = 0
        else:
            stall += 1
        if stall >= 1 and use_milp:
            dS, z = milp_switch(s, sq)
            log(f'  MILP dS* = {dS:.3e}')
            if dS is None or dS > -1e-11:
                break
            # step along the MILP direction just enough to change the best lines, then re-solve from there
            moved = False
            for t in (1e-4, 1e-5, 1e-6, 1e-7):
                sqt = sq + t * z[:-1].reshape(-1, 3)
                st_, sqt, _ = repair(s + t * z[-1], sqt, near_pairs(sqt))
                s3, sq3, info3, near3 = solve(st_, sqt, R)
                s3, sq3, _ = repair(s3, sq3, near3)
                log(f'  kick t={t:.0e}: repaired {st_:.13f} -> nlp {s3:.13f}')
                if s3 < s - 1e-14:
                    s, sq, moved = s3, sq3, True
                    break
            if not moved:
                break
            stall = 0
        elif stall >= 2:
            break
    log(f'done: {s0:.13f} -> {s:.13f} (gain {s0 - s:.3e})')
    return s, sq


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('input'); ap.add_argument('output')
    ap.add_argument('--rounds', type=int, default=30); ap.add_argument('--R', type=float, default=0.02)
    ap.add_argument('--no-milp', action='store_true')
    a = ap.parse_args()
    s, sq = load(a.input)
    s2, sq2 = run(s, sq, a.rounds, a.R, use_milp=not a.no_milp)
    save(a.output, s2, sq2)
