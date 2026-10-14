"""compare.py OLD NEW: byte-compare every file (except argv.log / results.pkl), results.pkl (exact), and the fq argv
logs per call: raw-identical, or equal after fq's own parse (last occurrence wins; values parsed as fq parses them)."""
import os, sys, pickle, filecmp
import numpy as np
A, B = sys.argv[1], sys.argv[2]

files = sorted(os.path.relpath(os.path.join(r, f), A) for r, _, fs in os.walk(A) for f in fs)
filesB = sorted(os.path.relpath(os.path.join(r, f), B) for r, _, fs in os.walk(B) for f in fs)
assert files == filesB, set(files) ^ set(filesB)
data = [f for f in files if f not in ('argv.log', 'results.pkl')]
bad = [f for f in data if not filecmp.cmp(os.path.join(A, f), os.path.join(B, f), shallow=False)]
print(f'files: {len(data)} compared, {len(bad)} differ', bad[:10])

def eq(x, y):
    if isinstance(x, np.ndarray) or isinstance(y, np.ndarray):
        return isinstance(x, np.ndarray) and isinstance(y, np.ndarray) and x.dtype == y.dtype and x.shape == y.shape \
            and x.tobytes() == y.tobytes()
    if type(x) != type(y):
        return False
    if isinstance(x, dict):
        return x.keys() == y.keys() and all(eq(x[k], y[k]) for k in x)
    if isinstance(x, (list, tuple)):
        return len(x) == len(y) and all(eq(a, b) for a, b in zip(x, y))
    if isinstance(x, float):
        return x.hex() == y.hex() or (x != x and y != y)
    return x == y

TIMING = ('sec', 't_alm', 'q_sec', 'h_sec')        # wall-clock fields of fq / harden JSON; trace[i][1] = step time
def norm(x, root):
    if isinstance(x, dict):
        return {k: (norm([[*r[:1], None, *r[2:]] for r in v], root) if k == 'trace' else norm(v, root))
                for k, v in x.items() if k not in TIMING}
    if isinstance(x, (list, tuple)):
        return type(x)(norm(v, root) for v in x)
    if isinstance(x, str):
        return x.replace(os.path.abspath(root), 'RUN')
    return x
RA = norm(pickle.load(open(f'{A}/results.pkl', 'rb')), A)
RB = norm(pickle.load(open(f'{B}/results.pkl', 'rb')), B)
assert RA.keys() == RB.keys()
rb = [k for k in RA if not eq(RA[k], RB[k])]
print(f'results (timing fields and run dir masked; all floats compared by bits): {len(RA)} keys, {len(rb)} differ', rb)
mk = RA['harden_make_input']
print(f'  harden make_input cases: {len(mk)}; sample', list(mk)[:2])
none = [k for k in RA if RA[k] is None]
print('  results that are None (both):', none)

BOOL = {'--no-alm', '--no-polish', '--shrink', '--round-compress', '--soft-compress', '--sqp', '--finish', '--ident', '--v'}
INT = {'--outer', '--pit', '--flip-top', '--stag-w', '--seed', '--round-steps', '--soft-steps'}
def parse(argv):
    """fq's view: (cmd, {flag: value}) with last occurrence winning; --in / --out paths dropped (they are per-tree)."""
    t = argv.split()
    cmd, o, i = t[0], {}, 1
    while i < len(t):
        f = t[i]
        if f in BOOL:
            o[f] = True; i += 1; continue
        v = t[i + 1]
        if f not in ('--in', '--out'):
            o[f] = int(v) if f in INT else float(v).hex()
        i += 2
    return cmd, o

def calls(p):
    out, key = {}, None
    for l in open(p):
        l = l.rstrip('\n')
        if l.startswith('### '):
            key = l.split()[1]; out[key] = []
        else:
            out[key].append(l)
    return out
CA, CB = calls(f'{A}/argv.log'), calls(f'{B}/argv.log')
assert CA.keys() == CB.keys()
n = raw = sem = 0
for k in CA:
    assert len(CA[k]) == len(CB[k]), (k, len(CA[k]), len(CB[k]))
    for x, y in zip(CA[k], CB[k]):
        n += 1
        strip = lambda s: ' '.join(w for w in s.split() if not w.startswith('/'))
        if strip(x) == strip(y):
            raw += 1
        elif parse(x) == parse(y):
            sem += 1
            if '-v' in sys.argv: print(f'  {k}: argv order/format differs, same fq parse:\n     old {strip(x)}\n     new {strip(y)}')
        else:
            print(f'  {k}: ARGV DIFFER\n     old {x}\n     new {y}')
print(f'argv: {n} fq calls; {raw} identical (paths aside), {sem} same after fq parse, {n - raw - sem} differ')
