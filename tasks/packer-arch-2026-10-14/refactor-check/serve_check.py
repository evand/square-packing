"""fq_quench mode='subprocess' vs mode='serve' (one long-lived worker, requests in two different orders)."""
import os, sys, json, random, tempfile
TREE = sys.argv[1]
sys.path.insert(0, TREE); os.chdir(TREE)
from pk.engine import QOpt, SCREEN, POLISH, POLISH_DEEP, fq_quench, FQ, read_cfg
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
OPTS = {'screen': SCREEN.but(loosen=1.02), 'polish': POLISH, 'deep': POLISH_DEEP, 'default': QOpt(),
        'kick': POLISH.but(kick=1e-4, seed=12345), 'kickalm': QOpt(loosen=1.0, kick=3e-3, seed=7),
        'soft05': QOpt(loosen=1.0, round_r0=0.05, hold_f=1.0, round_steps=40, round_kick=0.03, seed=11),
        'soft15c': QOpt(loosen=1.0, round_r0=0.15, round_compress=True, round_steps=40, round_kick=0.02, seed=11),
        'shrink': QOpt(loosen=1.0, kick=0.02, shrink=True, seed=3), 'rect': QOpt(loosen=1.0, rect_dx=0.25, flip_top=32),
        'softsq': QOpt(loosen=1.0, soft_a0=0.05, soft_compress=True, hold_f=1.0, soft_steps=20, soft_kick=0.01, seed=5),
        'stag': POLISH.but(stag_tol=1e-7)}
INP = []
for n in (40, 110, 237):
    s, sq = read_cfg(f'{DATA}/n{n}.txt')
    rng = random.Random(n)
    INP.append((f'n{n}file', f'{DATA}/n{n}.txt'))
    INP.append((f'n{n}kick', (s, [(x + rng.gauss(0, 0.01), y + rng.gauss(0, 0.01), a + rng.gauss(0, 0.2) + 90 * rng.choice((-1, 0, 1)))
                                  for x, y, a in sq])))
JOBS = [(i, o) for i in range(len(INP)) for o in OPTS if not (INP[i][0].startswith('n237') and o in ('soft05', 'soft15c', 'softsq', 'shrink', 'deep', 'default'))]
T = tempfile.mkdtemp()
def run(mode, order, tag):
    res = {}
    for i, o in order:
        name, inp = INP[i]
        d = os.path.join(T, f'{tag}_{name}_{o}'); os.makedirs(d)
        r = fq_quench(inp, OPTS[o], tmp=d, mode=mode, timeout=600, mod90=True)
        st = json.loads(r.stdout) if r.stdout else {'err': r.stderr}
        for k in ('sec', 't_alm'): st.pop(k, None)
        st['trace'] = [[q[0], *q[2:]] for q in st.get('trace', [])]
        res[(i, o)] = (open(r.out, 'rb').read() if r.out else None, st)
    return res
A = run('subprocess', JOBS, 'sub')
B = run('serve', JOBS, 'srv1')
rev = list(reversed(JOBS)); random.Random(1).shuffle(rev)
C = run('serve', rev, 'srv2')               # same worker process, different request order (state across requests?)
bad = [(INP[k[0]][0], k[1], w) for k in JOBS for w, X in (('serve', B), ('serve-shuffled', C))
       if X[k][0] != A[k][0] or X[k][1] != A[k][1]]
print(f'{len(JOBS)} (input, option) pairs x 3 runs; output files + JSON (timing masked) differing from subprocess: {len(bad)}', bad[:20])
print('worker calls', FQ and __import__('pk.engine', fromlist=['fq']).fq().calls)
print('subprocess runs with no output file / error:', [(INP[k[0]][0], k[1]) for k in JOBS if A[k][0] is None or 'err' in A[k][1]])
print('feasible outputs:', sum(A[k][1].get('min_gap', -1) >= 0 and A[k][1].get('min_wall', -1) >= 0 for k in JOBS), 'of', len(JOBS))
