"""The staircase end of lean/Sqpack/RyuStairCount.lean (our variant of Ryu's Section 6), in floats.

Construction for one end region (strip [0, w], ceiling g(X) = y + (X - sigma) mu, mu = tan(theta)):
  m = floor(g(w)) + 1;  column tilt tan(alpha/2) = (sqrt 2 + 1/v)/v (RyuStairEnd.end_bound);
  columns i < I = floor((w - m tau)/sec) with J_i = floor((g(i sec) - H0)/(1 - cos alpha)), m - J_i squares;
  layers j < J1 = max(J_{I-1}, floor(g(xi_I - g(w) tau))) of floor(w - a_j) squares, a_j = xi_{p(j)} - j tau;
  wedge rows r < floor(g0) of floor((m - r - 1) tau) squares.
Prints, per v and b in {v^3/2, v^3} and lifts y over the window [3v^2/4, 3v^2/4 + 2]: the worst end cost
E = A0(y) - N, the closed-form bound A0(y) - lb of RyuStairCount.EndP.count, and Eub(v) = 4.33 v^2 + 11 v + 30.

  python3 search/ryu/staircase.py [v ...]
"""
import bisect, math, sys


def end_cost(b, v, y, k=10 ** 9):
    tb = 1 / b + 1 / (8 * k * b ** 3)                       # the band angle of RyuStairBound.packs_stair
    gam, sig = (1 - tb * tb) / (1 + tb * tb), 2 * tb / (1 + tb * tb)
    w, mu = b * gam + sig, sig / gam
    g = lambda X: y + (X - sig) * mu
    g0, gb = g(0), g(w)
    m = math.floor(gb) + 1
    t = (math.sqrt(2) + 1 / v) / v
    ca, sa = (1 - t * t) / (1 + t * t), 2 * t / (1 + t * t)
    ta, se = sa / ca, 1 / ca
    H0 = m * ca + sa
    assert H0 <= g0, 'tilt condition'
    I = math.floor((w - m * ta) / se)
    J = [math.floor((g(i * se) - H0) / (1 - ca)) for i in range(I)]
    xi = lambda i: m * ta + i * se
    J1 = max(max(J) if J else 0, math.floor(g(xi(I) - gb * ta)))
    a = lambda j: xi(bisect.bisect_right(J, j)) - j * ta
    R0 = math.floor(g0)
    N = (sum(m - Ji for Ji in J) + sum(math.floor(w - a(j)) for j in range(J1))
         + sum(math.floor((m - r - 1) * ta) for r in range(R0)))
    lb = (I * m + J1 * (w - m * ta - I * se - 1) + ta * J1 * (J1 - 1) / 2
          + (se - 1) * ((I * (g0 - H0) + se * mu * (I * (I - 1) / 2)) / (1 - ca) - I)
          + R0 * (m - 1) * ta - ta * R0 * (R0 - 1) / 2 - R0)
    A0 = y * w + gam * sig * (1 + b * b) / 2
    assert lb <= N + 1e-6
    return A0 - N, A0 - lb


if __name__ == '__main__':
    vs = [float(x) for x in sys.argv[1:]] or [10, 12, 15, 20, 30, 50, 100]
    for v in vs:
        for b in sorted({math.ceil(v ** 3 / 2), math.floor(v ** 3)}):
            ys = [0.75 * v * v + 2 * i / 24 for i in range(25)]
            E = max(end_cost(b, v, y)[0] for y in ys)
            Elb = max(end_cost(b, v, y)[1] for y in ys)
            ub = 4.33 * v * v + 11 * v + 30
            print(f'v={v:g} b={b}: E <= {E:.1f}, A0 - lb <= {Elb:.1f}, Eub = {ub:.1f}  ({E / v ** 2:.3f}, '
                  f'{Elb / v ** 2:.3f}, {ub / v ** 2:.3f} x v^2)')
