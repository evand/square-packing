"""Old-vs-new call-site driver.  python3 drive.py TREE OUT  (TREE = .../search/packer).  Runs every refactored call
site on fixed inputs + seeds; writes fq argv log (FQLOG), output files and results.pkl into OUT."""
import os, sys, pickle, random, shutil, subprocess
TREE, OUT = sys.argv[1], os.path.abspath(sys.argv[2])
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
os.makedirs(OUT, exist_ok=True)
os.environ['FQLOG'] = os.path.join(OUT, 'argv.log')
open(os.environ['FQLOG'], 'w').close()
sys.path.insert(0, TREE); os.chdir(TREE)
import numpy as np
import mcmin, hop, explore, compose2, harden, relax, rectlib, chain

R = {}
def note(key, *msg):
    with open(os.environ['FQLOG'], 'a') as f:
        f.write(f'### {key} {" ".join(map(str, msg))}\n')

def d(name):
    p = os.path.join(OUT, name); os.makedirs(p, exist_ok=True); return p

def strip(x):
    if isinstance(x, dict):
        return {k: strip(v) for k, v in x.items() if k not in ('sec', 't_alm', 'q_sec', 'h_sec')}
    if isinstance(x, (list, tuple)):
        return type(x)(strip(v) for v in x)
    return x

L = {n: mcmin.load_deg(f'{DATA}/n{n}.txt') for n in (40, 41, 110, 111, 237, 238)}

def perturb(sq, seed, sig=0.02):
    rng = random.Random(seed)
    return [(x + rng.gauss(0, sig), y + rng.gauss(0, sig), a + rng.gauss(0, 20 * sig)) for x, y, a in sq]

# A. hop.quench (legacy EXTRA + extra)
CASES = [(['--loosen', '1.02'], ()), (['--loosen', '1.0'], ('--no-alm',)), (['--loosen', '1.02'], ('--pit', '8', '--flip-top', '0')),
         (['--loosen', '1.02'], ('--no-alm', '--loosen', '1.0'))]
for n in (40, 110, 237):
    s, sq = L[n]
    for ci, (ex, extra) in enumerate(CASES):
        if n == 237 and ci not in (1, 2):
            continue
        key = f'hop_n{n}_c{ci}'; note(key)
        hop.EXTRA[:] = ex
        t = d(key)
        r = hop.quench(s, perturb(sq, ci) if ci != 1 else sq, t, extra=extra)
        R[key] = r

# B. explore._screen_polish (+ kick_retry) with hop.propose proposals
for n, seed in ((40, 1), (40, 2), (110, 3), (110, 4), (237, 5)):
    s, sq = L[n]
    rng = random.Random(seed)
    prop, kind = hop.propose(sq, rng)
    k = int(np.ceil(np.sqrt(n) - 1e-12))
    for pe, kr in (((), (0, 1e-4, True, False)), (('--stag-tol', '1e-7'), (2, 1.234567e-4, True, False)),
                   ((), (2, 1e-4, False, False))):
        if n == 237 and kr[0]:
            continue
        key = f'sp_n{n}_s{seed}_{len(pe)}_{kr[0]}{kr[2]}'; note(key)
        t = d(key)
        r = explore._screen_polish(s, prop, kind, 'desc', rng, k, [], 0.0, k + 0.05, float('inf'), pe, t, None, kr=kr)
        R[key] = (strip(r), rng.random())

# C. compose2.polish
for n in (40, 110):
    s, sq = L[n]
    key = f'c2_n{n}'; note(key)
    R[key] = compose2.polish(n, s + 0.001, perturb(sq, 9, 0.0005), d(key), tag='x')

# D. harden.make_input + harden.job
mk = {}
for n in (40, 110, 237):
    s, sq = L[n]
    for xf in ('none', 'insert:1', 'insert:3', 'remove:all', 'remove:5'):
        for seed in range(1, 5):
            for k in ((0, 7, n - 1) if xf == 'remove:all' else (0,)):
                mk[(n, xf, seed, k)] = harden.make_input(sq, s, xf, seed, k)
R['harden_make_input'] = mk
for n in (40, 110):
    s, sq = L[n]
    key = f'hj_n{n}'; note(key); t = d(key)
    n2, s2, sq2 = harden.make_input(sq, s, 'insert:1', 3, 0)
    inp = os.path.join(t, 'in.txt'); harden.write_cfg(inp, n2, s2, sq2)
    R[key] = strip(harden.job((inp, os.path.join(t, 'raw.txt'), os.path.join(t, 'pol.txt'),
                               dict(pre='200', hold='500', release='400', settle='300'), 3)))

# E. relax.job
for n, ops in ((40, list(relax.OPS)), (110, ['polish', 'loose'])):
    for op in ops:
        key = f'rx_n{n}_{op}'; note(key); t = d(key)
        R[key] = strip(relax.job((f'{DATA}/n{n}.txt', os.path.join(t, 'out.txt'), op, 11)))

# F. rectlib._fq
for n in (40, 110):
    for extra in ([], ['--no-alm']):
        key = f'rl_n{n}_{len(extra)}'; note(key); t = d(key)
        r, err = rectlib._fq(f'{DATA}/n{n}.txt', os.path.join(t, 'out.txt'), np.float64(0.25), list(extra))
        R[key] = (None if r is None else (r[0], r[2]), err)

# G/H. fresh re-polish, rectlib compose polish (old: inline argv; new: helper / POLISH_DEEP)
for n in (110,):
    key = f'fresh_n{n}'; note(key); t = d(key); rp = os.path.join(t, 'rp.txt')
    import fresh
    if hasattr(fresh, 'repolish'):
        fresh.repolish(f'{DATA}/n{n}.txt', rp)
    else:
        subprocess.run([os.path.join(TREE, 'target', 'release', 'fq'), 'quench', '--in', f'{DATA}/n{n}.txt', '--out', rp,
                        '--loosen', '1.0', '--no-alm', '--flip-top', '64', '--pit', '2000'], capture_output=True)
    key = f'rlc_n{n}'; note(key); t = d(key); outp = os.path.join(t, 'n.txt')
    from pk import engine
    if hasattr(engine, 'POLISH_DEEP'):
        q = engine.fq_quench(f'{DATA}/n{n}.txt', engine.POLISH_DEEP, out=outp, binary=rectlib.FQ)
    else:
        q = subprocess.run([rectlib.FQ, 'quench', '--in', f'{DATA}/n{n}.txt', '--out', outp, '--loosen', '1.0', '--no-alm', '--flip-top', '64',
                            '--pit', '600'], capture_output=True, text=True)

# I. chain.seed_job (TRIES = 4)
chain.TRIES = 4
for src, k, n, seed in ((41, 1, 40, 1), (41, 1, 40, 2), (238, 1, 237, 3), (40, -1, 41, 4), (110, -2, 112, 5), (110, -1, 111, 6)):
    key = f'seed_{src}_{k}_{seed}'; note(key); t = d(key)
    R[key] = chain.seed_job((f'{DATA}/n{src}.txt', k, n, os.path.join(t, 'seed.txt'), seed))

pickle.dump(R, open(os.path.join(OUT, 'results.pkl'), 'wb'))
print('done', len(R))
