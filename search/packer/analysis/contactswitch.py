"""Where do two nearby minima differ in contacts?  (10-20; Evan: "places where a corner contacts a different side, or a
transition would require a corner-crossing motion?").  Inputs: exactsolve's <base>.exact.txt + <base>.contacts.json for
two packings.  Each contact is located geometrically (corner point P, outward side normal / wall), B is mapped into A's
frame (pk.distinct.best_match symmetry + square matching), and per square pair the contacts are compared:
  same        a contact in both (same corner owner, P within 0.05, normal within 5 deg)
  open/close  the pair (or square-wall) touches in only one packing: a gap opens; continuous, no corner crossing
  switch-side same owner and corner, but on a different side of the other square (normal turned ~90 deg): the corner
              must pass the other square's corner, i.e. a corner-crossing (through a corner-corner touch)
  switch-role the other square's corner now does the touching (roles swapped): also through a corner-corner touch
  other       anything else on a pair touching in both (e.g. a different corner of the same owner on the same side)
Wall contacts: corner a of i on wall w; 'switch-corner' = a different corner of i on the same wall (the square turned
through having a side flat on the wall), 'switch-wall' = the same corner on a different wall.
  python3 analysis/contactswitch.py BASE_A BASE_B"""
import sys, os, json, math, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
from pk.distinct import best_match

CU = ((1, 1), (-1, 1), (-1, -1), (1, -1)); EK = ((1, 0), (0, 1), (-1, 0), (0, -1)); H = 0.5


def load(base):
    L = [l.split() for l in open(base + '.exact.txt') if l.strip()]
    n, s = int(L[0][0]), float(L[0][1])
    sq = np.array([[float(v) for v in l[:3]] for l in L[1:n + 1]])
    return s, sq, json.load(open(base + '.contacts.json'))


def geom(s, sq, K):
    """contacts -> list of dicts in this packing's own frame"""
    out = []
    for e in K['equations']:
        if e[0] == 'W':
            _, i, a, w = e; c, si = math.cos(math.radians(sq[i, 2])), math.sin(math.radians(sq[i, 2]))
            su, sv = CU[a]; P = (sq[i, 0] + H * (c * su - si * sv), sq[i, 1] + H * (si * su + c * sv))
            out.append(dict(kind='W', i=i, P=np.array(P), w=w))
        else:
            _, j, a, i, k = e
            cj, sj = math.cos(math.radians(sq[j, 2])), math.sin(math.radians(sq[j, 2]))
            ci, si = math.cos(math.radians(sq[i, 2])), math.sin(math.radians(sq[i, 2]))
            su, sv = CU[a]; P = (sq[j, 0] + H * (cj * su - sj * sv), sq[j, 1] + H * (sj * su + cj * sv))
            ex, ey = EK[k]; nrm = (ci * ex - si * ey, si * ex + ci * ey)
            out.append(dict(kind='C', owner=j, other=i, P=np.array(P), n=np.array(nrm)))
    return out


def symp(P, s, k):
    """box symmetry k (pk.cluster._sym convention) on a point"""
    x, y = P
    if k >= 4:
        x = s - x
    for _ in range(k % 4):
        x, y = s - y, x
    return np.array([x, y])


def symv(v, k):
    x, y = v
    if k >= 4:
        x = -x
    for _ in range(k % 4):
        x, y = -y, x
    return np.array([x, y])


def symw(w, k):
    m = dict(L='L', R='R', B='B', T='T')
    if k >= 4:
        m = {kk: dict(L='R', R='L').get(v, v) for kk, v in m.items()}
    for _ in range(k % 4):                    # (x, y) -> (s - y, x): L(x=0) -> B, B(y=0) -> R, R -> T, T -> L
        m = {kk: dict(L='B', B='R', R='T', T='L')[v] for kk, v in m.items()}
    return m[w]


def compare(A_base, B_base):
    sA, qA, KA = load(A_base); sB, qB, KB = load(B_base)
    _, k, perm, _ = best_match(qA, sA, qB, sB)                # A[i] ~ Bk[perm[i]]
    inv = {int(perm[i]): i for i in range(len(perm))}
    GA = geom(sA, qA, KA)
    GB = []
    for g in geom(sB, qB, KB):
        h = dict(g); h['P'] = symp(g['P'], sB, k) * sA / sB
        if g['kind'] == 'W':
            h['i'] = inv[g['i']]; h['w'] = symw(g['w'], k)
        else:
            h['owner'] = inv[g['owner']]; h['other'] = inv[g['other']]; h['n'] = symv(g['n'], k)
        GB.append(h)
    def key(g):
        return ('W', g['i']) if g['kind'] == 'W' else ('C',) + tuple(sorted((g['owner'], g['other'])))
    byA, byB = collections.defaultdict(list), collections.defaultdict(list)
    for g in GA: byA[key(g)].append(g)
    for g in GB: byB[key(g)].append(g)
    def same(g, h):
        if g['kind'] == 'W':
            return g['w'] == h['w'] and np.linalg.norm(g['P'] - h['P']) < 0.05
        return g['owner'] == h['owner'] and np.linalg.norm(g['P'] - h['P']) < 0.05 and g['n'] @ h['n'] > math.cos(math.radians(5))
    C = collections.Counter(); ex = []; SW = []
    vvA = [tuple(v) for v in (KA.get('vv') or [])]
    vvB = [tuple(inv[x] for x in v) for v in (KB.get('vv') or [])]
    for kk in set(byA) | set(byB):
        a, b = byA.get(kk, []), byB.get(kk, [])
        ua = [g for g in a if not any(same(g, h) for h in b)]; ub = [h for h in b if not any(same(g, h) for g in a)]
        if not ua and not ub:
            continue
        if not a or not b:
            C['open/close'] += 1; continue
        for g in ua:                                           # A-only contact on a pair touching in both
            if kk[0] == 'W':
                kind = 'switch-wall' if any(np.linalg.norm(g['P'] - h['P']) < 0.05 for h in ub) else 'switch-corner'
            else:
                hs = [h for h in ub]
                if any(h['owner'] == g['owner'] and np.linalg.norm(g['P'] - h['P']) < 0.05 for h in hs):
                    kind = 'switch-side'
                elif any(h['owner'] != g['owner'] for h in hs):
                    kind = 'switch-role'
                elif hs:
                    kind = 'other'
                else:
                    kind = 'drop (pair still touches elsewhere)'
            C[kind] += 1; ex.append((kind, kk))
            if kind.startswith('switch') and kk[0] == 'C':
                i1, j1 = kk[1], kk[2]
                isvv = any(set(v) == {i1, j1} for v in vvA) or any(set(v) == {i1, j1} for v in vvB)
                # distance from the touching corner to the nearest corner of the other square (in A)
                oth = g['other']; c_, s_ = math.cos(math.radians(qA[oth, 2])), math.sin(math.radians(qA[oth, 2]))
                dmin = min(np.linalg.norm(g['P'] - (qA[oth, :2] + H * np.array([c_ * su - s_ * sv, s_ * su + c_ * sv])))
                           for su, sv in CU)
                SW.append((kind, isvv, dmin))
        for h in ub:
            if not ua:
                C['drop (pair still touches elsewhere)'] += 1
    vv = (len(KA.get('vv') or []), len(KB.get('vv') or []))
    sw = collections.Counter((kd, 'at a corner-corner touch' if v else 'not vv') for kd, v, _ in SW)
    far = [round(d, 4) for kd, v, d in SW if not v]
    return dict(symmetry=k, kinds=dict(C), switches=dict((f'{a} / {b}', c) for (a, b), c in sw.items()),
                corner_dist_not_vv=sorted(far)[:10], corner_dist_max=round(max((d for *_, d in SW), default=0), 5), vv=vv)


if __name__ == '__main__':
    print(json.dumps(compare(sys.argv[1], sys.argv[2]), default=str, indent=1))
