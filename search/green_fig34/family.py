"""Parametric families read off Figure 34 of Friedman's DS7 (pixel coordinates in sets.py)."""
import math

def green17(t, a, b, s, v):
    """16 points, 4 rows of 4; rows alternate type A / type B, set invariant under 180-deg rotation.
    type A row: x = a, a+e, a+e+s, t-a ; type B row: x = a, a+s, a+2s, t-a ; e = t-2a-2s.
    rows y = b (A), b+v (B), t-b-v (A), t-b (B)."""
    e = t - 2 * a - 2 * s
    A = [a, a + e, a + e + s, t - a]
    B = [a, a + s, a + 2 * s, t - a]
    P = []
    for y, row in ((b, A), (b + v, B), (t - b - v, A), (t - b, B)):
        P += [(x, y) for x in row]
    return P

# pixel data (GIF 171x171, frame pixels 0 and 170; dot centres)
L17_PX = [(37, 35), (56, 35), (94, 35), (132, 35), (37, 68), (75, 68), (113, 68), (132, 68),
          (37, 101), (56, 101), (94, 101), (132, 101), (37, 135), (75, 135), (113, 135), (132, 135)]
L19_PX = [(38, 38), (69, 38), (100, 38), (132, 38), (38, 69), (59, 69), (90, 69), (121, 69), (132, 69),
          (38, 100), (48, 100), (80, 100), (111, 100), (132, 100), (38, 132), (69, 132), (100, 132), (132, 132)]

T17 = (40 * math.sqrt(2) + 19) / 17
T19 = 6 * math.sqrt(2) - 4

def from_px(px, t, frame=170):
    k = t / frame
    return [(x * k, (frame - y) * k) for (x, y) in px]
