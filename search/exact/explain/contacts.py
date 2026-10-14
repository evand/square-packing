"""Contacts of a regularized record: corner-on-wall and corner-on-side (side labels u+/u-/v+/v- in each square's own
raw-angle frame; a corner at a corner yields both side labels)."""
import math, render, jstyle
def contacts(S, sq, tol=1e-9):
    out = []
    C = [jstyle.corners(*q) for q in sq]
    def inside_edge(p, q, k):
        x, y, t = q; c, s = math.cos(t), math.sin(t)
        u = (p[0]-x)*c + (p[1]-y)*s; v = -(p[0]-x)*s + (p[1]-y)*c
        out = []
        if abs(abs(u)-.5) < tol and abs(v) <= .5+tol: out.append('u+' if u > 0 else 'u-')
        if abs(abs(v)-.5) < tol and abs(u) <= .5+tol: out.append('v+' if v > 0 else 'v-')
        return out
    for i, q in enumerate(sq):
        for k, p in enumerate(C[i]):
            for nm, val in (('L', p[0]), ('R', S-p[0]), ('B', p[1]), ('T', S-p[1])):
                if abs(val) < tol: out.append(('wall', i, k, nm))
            for j, r in enumerate(sq):
                if j != i:
                    for e in inside_edge(p, r, k): out.append(('corner-on-edge', i, k, j, e))
    return out
if __name__ == '__main__':
    import sys; n = int(sys.argv[1]); S, sq = render.load(n)
    for i, q in enumerate(sq): print(i, [round(v, 6) for v in q[:2]], round(math.degrees(q[2]), 6))
    for c in contacts(S, sq): print(c)
