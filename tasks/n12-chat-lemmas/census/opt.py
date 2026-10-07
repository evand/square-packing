"""Epigraph margin maximisation for n unit squares in [0, L]^2 with per-square centre boxes.

    max delta  s.t.  for every pair (i, j) a line n(phi_ij).p = c_ij with
                       n.v <= c - delta/2  for the 4 vertices v of square i,
                       n.v >= c + delta/2  for the 4 vertices v of square j,
                     each square delta away from every wall,
                     centre of square i in its box (bounds).

Variables: x(n), y(n), theta(n) in [0, pi/2], phi(P), c(P), delta.  SLSQP with analytic Jacobians.
For a fixed pose set the optimum over (phi, c) of the pair gap is the Euclidean distance (disjoint) or
minus the penetration depth (overlapping); the reported margin is recomputed independently by `sat_margin`
(edge-normal separating axes, plus wall distances), which is what RESULTS.md quotes.
"""
import numpy as np
from scipy.optimize import minimize

L = 4.0
U = 0.5 * np.array([[1, 1], [1, -1], [-1, -1], [-1, 1]], float)


def verts(P):
    x, y, t = P[:, 0], P[:, 1], P[:, 2]
    c, s = np.cos(t), np.sin(t)
    vx = x[:, None] + U[None, :, 0] * c[:, None] - U[None, :, 1] * s[:, None]
    vy = y[:, None] + U[None, :, 0] * s[:, None] + U[None, :, 1] * c[:, None]
    return np.stack([vx, vy], -1)          # n x 4 x 2


def pair_sat(P):
    """SAT gap matrix (edge normals of both squares); negative = penetration depth along best axis."""
    V = verts(P)
    n = len(P)
    G = np.full((n, n), np.inf)
    for i in range(n):
        for j in range(i + 1, n):
            best = -np.inf
            for t in (P[i, 2], P[i, 2] + np.pi / 2, P[j, 2], P[j, 2] + np.pi / 2):
                d = np.array([np.cos(t), np.sin(t)])
                a, b = V[i] @ d, V[j] @ d
                best = max(best, b.min() - a.max(), a.min() - b.max())
            G[i, j] = G[j, i] = best
    return G


def wall_gap(P):
    V = verts(P)
    return np.minimum(np.minimum(V[..., 0].min(1), L - V[..., 0].max(1)),
                      np.minimum(V[..., 1].min(1), L - V[..., 1].max(1)))


def sat_margin(P):
    G = pair_sat(P)
    return float(min(G[np.isfinite(G)].min(), wall_gap(P).min()))


class Problem:
    def __init__(self, boxes):
        self.n = n = len(boxes)
        self.boxes = np.array(boxes, float)
        self.I, self.J = np.triu_indices(n, 1)
        self.P = len(self.I)
        self.nv = 3 * n + 2 * self.P + 1

    def split(self, z):
        n, P = self.n, self.P
        return z[:n], z[n:2 * n], z[2 * n:3 * n], z[3 * n:3 * n + P], z[3 * n + P:3 * n + 2 * P], z[-1]

    def cons(self, z):
        x, y, t, phi, c, d = self.split(z)
        I, J = self.I, self.J
        cp, sp = np.cos(phi), np.sin(phi)
        ai, aj = phi - t[I], phi - t[J]
        # projections n.v for the 4 vertices: x cos phi + y sin phi + ux cos(a) + uy sin(a)
        pi_ = (cp * x[I] + sp * y[I])[:, None] + U[None, :, 0] * np.cos(ai)[:, None] + U[None, :, 1] * np.sin(ai)[:, None]
        pj_ = (cp * x[J] + sp * y[J])[:, None] + U[None, :, 0] * np.cos(aj)[:, None] + U[None, :, 1] * np.sin(aj)[:, None]
        gi = (c - d / 2)[:, None] - pi_
        gj = pj_ - (c + d / 2)[:, None]
        w = 0.5 * (np.cos(t) + np.sin(t))
        gw = np.concatenate([x - w - d, L - x - w - d, y - w - d, L - y - w - d])
        return np.concatenate([gi.ravel(), gj.ravel(), gw])

    def jac(self, z):
        x, y, t, phi, c, d = self.split(z)
        n, Pn, I, J = self.n, self.P, self.I, self.J
        cp, sp = np.cos(phi), np.sin(phi)
        m = 8 * Pn + 4 * n
        Jm = np.zeros((m, self.nv))
        r = np.arange(Pn)
        for side, K, sign, off in ((0, I, -1.0, 0), (1, J, 1.0, 4 * Pn)):
            a = phi - t[K]
            ca, sa = np.cos(a), np.sin(a)
            for k in range(4):
                ux, uy = U[k]
                rows = off + 4 * r + k
                # d proj / d var
                dx, dy = cp, sp
                dth = ux * sa - uy * ca                      # d/dtheta of ux cos(phi-t)+uy sin(phi-t)
                dph = -sp * x[K] + cp * y[K] - ux * sa + uy * ca
                Jm[rows, K] += sign * dx
                Jm[rows, n + K] += sign * dy
                Jm[rows, 2 * n + K] += sign * dth
                Jm[rows, 3 * n + r] += sign * dph
                Jm[rows, 3 * n + Pn + r] = -sign            # gi: +c ; gj: -c
                Jm[rows, -1] = -0.5
        o = 8 * Pn
        ii = np.arange(n)
        dw = 0.5 * (-np.sin(t) + np.cos(t))
        Jm[o + ii, ii] = 1; Jm[o + ii, 2 * n + ii] = -dw
        Jm[o + n + ii, ii] = -1; Jm[o + n + ii, 2 * n + ii] = -dw
        Jm[o + 2 * n + ii, n + ii] = 1; Jm[o + 2 * n + ii, 2 * n + ii] = -dw
        Jm[o + 3 * n + ii, n + ii] = -1; Jm[o + 3 * n + ii, 2 * n + ii] = -dw
        Jm[o:, -1] = -1
        return Jm

    def bounds(self):
        b = self.boxes
        bx = [(lo, hi) for lo, hi, _, _ in b]
        by = [(lo, hi) for _, _, lo, hi in b]
        bt = [(0.0, np.pi / 2)] * self.n
        return bx + by + bt + [(None, None)] * (2 * self.P) + [(-3.0, 1.0)]

    def lines_from(self, P):
        """Best separating line per pair among edge normals and the centre direction."""
        V = verts(P)
        phi = np.zeros(self.P); c = np.zeros(self.P); g = np.zeros(self.P)
        for p, (i, j) in enumerate(zip(self.I, self.J)):
            cand = [P[i, 2], P[j, 2]]
            cand = [a + q * np.pi / 2 for a in cand for q in range(4)]
            dd = P[j, :2] - P[i, :2]
            cand.append(np.arctan2(dd[1], dd[0]))
            best = (-np.inf, 0, 0)
            for a in cand:
                dv = np.array([np.cos(a), np.sin(a)])
                pa, pb = V[i] @ dv, V[j] @ dv
                gap = pb.min() - pa.max()
                if gap > best[0]:
                    best = (gap, a, 0.5 * (pb.min() + pa.max()))
            g[p], phi[p], c[p] = best
        return phi, c, g

    def pack(self, P):
        phi, c, g = self.lines_from(P)
        d = min(g.min(), wall_gap(P).min())
        return np.concatenate([P[:, 0], P[:, 1], P[:, 2], phi, c, [max(d, -3.0)]])

    def poses(self, z):
        x, y, t = z[:self.n], z[self.n:2 * self.n], z[2 * self.n:3 * self.n]
        return np.column_stack([x, y, t])

    def clip(self, P):
        P = P.copy()
        b = self.boxes
        P[:, 0] = np.clip(P[:, 0], b[:, 0], b[:, 1])
        P[:, 1] = np.clip(P[:, 1], b[:, 2], b[:, 3])
        P[:, 2] = np.mod(P[:, 2], np.pi / 2)
        return P

    def solve(self, P0, rounds=4, maxiter=400, tol=1e-12):
        P = self.clip(P0)
        best = (sat_margin(P), P)
        for _ in range(rounds):
            z0 = self.pack(P)
            res = minimize(lambda z: -z[-1], z0, jac=lambda z: np.eye(1, self.nv, self.nv - 1).ravel() * -1,
                           method='SLSQP', bounds=self.bounds(),
                           constraints=[{'type': 'ineq', 'fun': self.cons, 'jac': self.jac}],
                           options={'maxiter': maxiter, 'ftol': tol})
            P = self.clip(self.poses(res.x))
            m = sat_margin(P)
            improved = m > best[0] + 1e-10
            if m > best[0]:
                best = (m, P)
            if not improved:
                break
        return best


def soft_relax(P0, boxes, iters=300, rng=None):
    """Penalty pre-relaxation (seed generator only; never used as a final answer)."""
    prob = Problem(boxes)
    from scipy.optimize import minimize as mn
    n = len(P0)

    def energy(v):
        P = v.reshape(n, 3)
        G = pair_sat(P)
        e = np.sum(np.minimum(G[np.triu_indices(n, 1)], 0) ** 2)
        e += np.sum(np.minimum(wall_gap(P), 0) ** 2)
        return e
    b = prob.boxes
    bnds = []
    for i in range(n):
        bnds += [(b[i, 0], b[i, 1]), (b[i, 2], b[i, 3]), (0, np.pi / 2)]
    r = mn(energy, prob.clip(P0).ravel(), method='L-BFGS-B', bounds=bnds, options={'maxiter': iters})
    return r.x.reshape(n, 3)
