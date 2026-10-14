"""Each s(11) contact in the record and in the real phantom root S = -1.853: (along, across) of the corner in the
other square's frame.  One coordinate is +-1/2 (the equation); the other is where along the side line it lands."""
import os
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'runs', 'explain', 'phantom11.json')
import json, math, render, contacts as C, jstyle
Sr, sq = render.load(11)
cr, sr = math.cos(sq[6][2]), math.sin(sq[6][2])          # record frame for tilted squares
P = json.load(open(OUT))
def frames(c, s):
    return {i: ((1, 0), (0, 1)) if i < 6 else ((c, s), (-s, c)) for i in range(11)}
def config(o):
    S = o['S'][0]; g = lambda k: o[k][0]
    cen = [(.5, .5), (.5, S - .5), (S - .5, g('y2')), (S - 1.5, .5), (S - .5, .5), (S - .5, 1.5)]
    cen += [(g(f'x{i}'), g(f'y{i}')) for i in range(6, 11)]
    return S, cen, g('c0'), g('s0')
# symsys tilted frame: fold angle a in (0,90): c0 = cos(a); record raw angle -40.18 => a = 49.82
def local_signs(i, k):
    x, y, t = sq[i]; p = jstyle.corners(*sq[i])[k]
    a = t % (math.pi/2) if i >= 6 else 0.0
    du = (p[0]-x)*math.cos(a) + (p[1]-y)*math.sin(a); dv = -(p[0]-x)*math.sin(a) + (p[1]-y)*math.cos(a)
    return (.5 if du > 0 else -.5), (.5 if dv > 0 else -.5)
def run(o, label):
    S, cen, c, s = config(o); F = frames(c, s)
    print(label, 'S =', round(S, 4))
    for ct in C.contacts(Sr, sq):
        if ct[0] != 'corner-on-edge': continue
        _, i, k, j, ed = ct
        A, B = local_signs(i, k); (u, v) = F[i]
        P_ = (cen[i][0] + A*u[0] + B*v[0], cen[i][1] + A*u[1] + B*v[1])
        (uj, vj) = F[j]; d = (P_[0]-cen[j][0], P_[1]-cen[j][1])
        du = d[0]*uj[0] + d[1]*uj[1]; dv = d[0]*vj[0] + d[1]*vj[1]
        print(f'  corner of {i+1} on square {j+1}: (along, across) = ({du:+.3f}, {dv:+.3f})')
run(P[1], 'record'); run(P[0], 'phantom')
