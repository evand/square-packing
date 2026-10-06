"""Occupancy classes for 12 unit squares in side L = 4 with strip width h.

Regions (closed boxes in centre coordinates):
  corners  K0=BL K1=BR K2=TR K3=TL    (within h of two walls)
  edges    E0=B  E1=R  E2=T  E3=L     (within h of exactly that wall)
  centre   Z                          (farther than h from every wall)
Wall w's strip holds E_w + K_w + K_{w+1}  (wall B: K0, K1; R: K1, K2; T: K2, K3; L: K3, K0).
A class is (k0..k3, e0..e3, z) up to D4.
"""
import itertools

L = 4.0


def strips(c):
    k, e = c[:4], c[4:8]
    return tuple(e[w] + k[w] + k[(w + 1) % 4] for w in range(4))


def _rot(c):          # rotate the picture by 90 deg: corner i -> i+1, wall w -> w+1
    k, e, z = c[:4], c[4:8], c[8]
    return tuple(k[(i - 1) % 4] for i in range(4)) + tuple(e[(i - 1) % 4] for i in range(4)) + (z,)


def _refl(c):         # reflect x -> L - x: BL<->BR, TL<->TR, walls B,T fixed, L<->R
    k, e, z = c[:4], c[4:8], c[8]
    return (k[1], k[0], k[3], k[2], e[0], e[3], e[2], e[1], z)


def orbit(c):
    out, x = [], c
    for _ in range(4):
        out += [x, _refl(x)]
        x = _rot(x)
    return out


def canon(c):
    return max(orbit(c))


def enumerate_classes(n=12, maxstrip=3):
    seen = set()
    for k in itertools.product((0, 1), repeat=4):
        for e in itertools.product(range(n + 1), repeat=4):
            z = n - sum(k) - sum(e)
            if z < 0:
                continue
            c = k + e + (z,)
            if max(strips(c)) > maxstrip:
                continue
            seen.add(canon(c))
    return sorted(seen, key=lambda c: (-c[8], -sum(c[:4]), c))


def regions(h):
    """Centre boxes (xlo, xhi, ylo, yhi) for the 9 regions, index 0..3 corners, 4..7 edges, 8 centre."""
    lo, hi = 0.5, L - 0.5
    a, b = h, L - h
    return [
        (lo, a, lo, a), (b, hi, lo, a), (b, hi, b, hi), (lo, a, b, hi),     # K0..K3
        (a, b, lo, a), (b, hi, a, b), (a, b, b, hi), (lo, a, a, b),         # E_B, E_R, E_T, E_L
        (a, b, a, b),                                                        # Z
    ]


def assignment(c):
    """List of 12 region indices realising class c."""
    out = []
    for r in range(9):
        out += [r] * c[r]
    return out


def classify(P, h):
    """Class vector of poses P (n x 3) at strip width h (closed: within h means dist <= h)."""
    c = [0] * 9
    for x, y, _ in P:
        nb, nt, nl, nr = int(y <= h), int(L - y <= h), int(x <= h), int(L - x <= h)
        m = nb + nt + nl + nr
        if m >= 2:
            idx = {(1, 0, 1, 0): 0, (1, 0, 0, 1): 1, (0, 1, 0, 1): 2, (0, 1, 1, 0): 3}.get((nb, nt, nl, nr))
            c[idx] += 1
        elif m == 1:
            c[4 + [nb, nr, nt, nl].index(1)] += 1
        else:
            c[8] += 1
    return tuple(c)


def label(c):
    k, e, z = sum(c[:4]), sum(c[4:8]), c[8]
    return f"k{k} e{e} z{z} K{''.join(map(str, c[:4]))} E{''.join(map(str, c[4:8]))} S{''.join(map(str, strips(c)))}"


if __name__ == "__main__":
    cs = enumerate_classes()
    print(len(cs), "classes")
    from collections import Counter
    print(Counter((sum(c[:4]), c[8]) for c in cs))
