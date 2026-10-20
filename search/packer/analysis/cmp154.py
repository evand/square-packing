"""154: contact-graph diff + expansion barrier between certified packings (Evan 10-19: no magic move thresholds)."""
import json, sys, os, time, math
from multiprocessing import Pool
PK = '/home/evand/math/square-packing/public/search/packer'
sys.path.insert(0, PK); sys.path.insert(0, PK + '/landscape')
import morph
W = PK + '/runs/candidates/'
P = {'552': W + 'chk154_552', '583': W + 'chk154_583', '586': W + 's154_sw/prior', 'hit': W + 's154_sw'}
PAIRS = [('hit', '586'), ('552', '583'), ('583', '586'), ('552', '586')]
H = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]

def ld(k):
    from fractions import Fraction
    L = [l.split() for l in open(P[k] + '/cand.txt') if l.strip()]
    s = float(Fraction(L[0][1]))
    return s, [(float(r[0]), float(r[1]), math.radians(float(r[2]))) for r in L[1:]]

SWAP = [dict(L='B', B='L', R='T', T='R'), dict(L='R', R='L'), dict(B='T', T='B')]   # transpose, x-flip, y-flip

def sym(X, k):
    s, out = X[0], []
    for x, y, t in X[1]:
        if k & 1: x, y, t = y, x, math.pi / 2 - t
        if k & 2: x, t = s - x, -t
        if k & 4: y, t = s - y, -t
        out.append((x, y, t))
    return s, out

def side(lab, k):
    for b, sw in zip((1, 2, 4), SWAP):
        if k & b: lab = sw.get(lab, lab)
    return lab

def frame_of(A, B):
    return min(range(8), key=lambda k: morph.match(A, sym(B, k))[1])

def graph(k):
    d = json.load(open(P[k] + '/cand.contacts.json'))
    E = set()
    for e in d['equations']:
        E.add(('W', e[1], e[3]) if e[0] == 'W' else ('C', min(e[1], e[3]), max(e[1], e[3])))
    return E, set(d['free'])

def cdiff(a, b):
    A, B0 = ld(a), ld(b)
    k = frame_of(A, B0); B = sym(B0, k)
    perm = morph.match2(A, B)
    Ea, Fa = graph(a); Eb0, Fb = graph(b)
    Eb = {('W', e[1], side(e[2], k)) if e[0] == 'W' else e for e in Eb0}
    m = lambda e: ('W', perm[e[1]], e[2]) if e[0] == 'W' else ('C', *sorted((perm[e[1]], perm[e[2]])))
    Ma = {m(e) for e in Ea}
    return dict(pair=f'{a}-{b}', sym=k, cost=round(morph.match(A, B)[1], 3), edges_a=len(Ea), edges_b=len(Eb), shared=len(Ma & Eb), a_only=len(Ma - Eb),
                b_only=len(Eb - Ma), rattlers_a=len(Fa), rattlers_b=len(Fb),
                rattler_shared=len({perm[i] for i in Fa} & Fb))

def walk(job):
    a, b, frame = job
    A, B0 = ld(a), ld(b)
    B = sym(B0, frame_of(A, B0))
    if frame:
        T = lambda X: (X[0], [(y, x, morph.Q - t) for x, y, t in X[1]])
        A, B = T(A), T(B)
    perm = morph.match2(A, B)
    base = max(A[0], B[0]); out = []
    for h in H:
        t0 = time.time(); ok, d, _ = morph.walk_carrot(base + h, A, B, perm)
        out.append((h, bool(ok), float(d), round(time.time() - t0)))
        print(a, '->', b, 'frame', frame, h, ok, f'{d:.2e}', flush=True)
        if ok:
            break
    return dict(a=a, b=b, frame=frame, base=base, tries=out, conn=next((h for h, ok, _, _ in out if ok), None))

if __name__ == '__main__':
    for a, b in PAIRS:
        print(json.dumps(cdiff(a, b)), flush=True)
    jobs = [(a, b, f) for a, b in PAIRS for f in (0, 1)] + [(b, a, f) for a, b in PAIRS for f in (0, 1)]
    with Pool(8) as pool:
        R = pool.map(walk, jobs, chunksize=1)
    json.dump(R, open(PK + '/runs/cmp154.json', 'w'), indent=1)
    for r in R:
        print('RESULT', r['a'], '->', r['b'], 'frame', r['frame'], 'connects at base +', r['conn'])
