"""Conjectured strip frontier min(n, F_n(H)), 1 <= H < 2 (Friedman, Math Magic June 2010, n <= 10; even n >= 12
by the extrapolated Morandi family, checked by explicit construction 10-22)."""
import math
import numpy as np
from scipy.optimize import brentq

R2 = math.sqrt(2)

def even_ab(m, t):
    c = math.sqrt(1 - t * t); D = m - (2 * m - 1) * t * t
    a = 2 * (1 - t * t - (m - 1) * t ** 3 + (m - 1) * c ** 3) / D
    b = 2 * (m - (m - 1) * (m - 2) * t - (2 * m - 1) * t * t - (m - 1) * t ** 3 + (m - 1) * (t * t - 2 * t + m) * c) / D
    return a, b

def F(n, H):
    if n % 2:
        return 4 + (2 * n - 4) * R2 - (n - 1) * H
    m = n // 2
    if m == 1: return math.inf
    # branch t in (1/sqrt2, t_pole) with a decreasing from 2; find a(t) = H
    tp = math.sqrt(m / (2 * m - 1))
    lo, hi = 1 / R2 + 1e-12, tp - 1e-9
    f = lambda t: even_ab(m, t)[0] - H
    if f(lo) * f(hi) > 0:
        return math.inf
    return even_ab(m, brentq(f, lo, hi, xtol=1e-15))[1]

def target(n, H):
    return min(n, F(n, H))
