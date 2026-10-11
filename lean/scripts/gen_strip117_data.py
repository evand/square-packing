#!/usr/bin/env python3
"""Generate lean/Sqpack/Strip117Data.lean from the wall-strip certificate at h = 117/100.

    python3 lean/scripts/gen_strip117_data.py [cert.json] [out.lean]

Input: tasks/n12-chat-lemmas/strip-demo/certs/h1.170_fine.json (PROOF.md section 7).  Output: the
certificate in the integer format of `Sqpack/StripCheck.lean` (`Cert`), plus the coverage witness and
the min-plus chain tables, which are computed here.  Before writing, every check of
`StripCheck.Cert` (paramsOK, ruleOK, coverOK, flagsOK, leavesOK, chainOK) is replayed in Python with
the same list arithmetic, so the `decide +kernel` proofs in `Sqpack/Strip117.lean` are expected to
succeed; Lean re-checks everything and trusts nothing produced here.
"""
import json, sys, os
from fractions import Fraction as F
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CERT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    ROOT, 'tasks/n12-chat-lemmas/strip-demo/certs/h1.170_fine.json')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, 'lean/Sqpack/Strip117Data.lean')

K = 10 ** 6
TR = 80000
DEPTH = 6

# ---- BernZ list arithmetic, exactly as in Sqpack/BernsteinZ.lean -------------------------------
def add(p, q):
    if not p:
        return list(q)
    if not q:
        return list(p)
    return [p[0] + q[0]] + add(p[1:], q[1:])

def smul(c, p):
    return [c * a for a in p]

def mul(p, q):
    if not p:
        return []
    return add(smul(p[0], q), [0] + mul(p[1:], q))

def hom(U0, H, R, n, p):
    if not p:
        return []
    if len(p) == 1:
        return [p[0] * R ** n]
    return add([p[0] * R ** n], mul([U0, H], hom(U0, H, R, n - 1, p[1:])))

def check(p, n, U0, U1, R):
    if len(p) > n + 1:
        return False
    b = hom(U0, U1 - U0, R, n, p)
    if len(b) > n + 1:
        return False
    bb = lambda k: b[k] if k < len(b) else 0
    return all(sum(comb(n - k, i - k) * bb(k) for k in range(i + 1)) >= 0 for i in range(n + 1))

def chkSub(p, n, U0, U1, R, d):
    if check(p, n, U0, U1, R):
        return True
    if d == 0:
        return False
    return chkSub(p, n, 2 * U0, U0 + U1, 2 * R, d - 1) and chkSub(p, n, U0 + U1, 2 * U1, 2 * R, d - 1)

def chk(p, T0, T1):
    return chkSub(p, len(p) - 1, T0, T1, TR, DEPTH)

pN, pM, pW, pMW = [1, 0, 1], [1, 0, -1], [0, 2], [1, 2, -1]

def comboOK(A, B, T0, T1, C0, C1):
    if not B:
        return chk(A, T0, T1)
    return chk(add(A, smul(C1, B)), T0, T1) and (
        chk(add(A, smul(C0, B)), T0, T1) or chk(add(mul([2, 0, 2], A), mul(smul(K, pMW), B)), T0, T1))

A11 = lambda Zl, Zr, Q: add(smul(K, pN), add(smul(Zl - Zr, pW), smul(-Q, pM)))
A22 = lambda Zl, Zr, Q: add(smul(K, pN), add(smul(-(Zl - Zr), pM), smul(-Q, pW)))
A12 = lambda Zl, Zr, Q: add(mul(pW, add(smul(K, pN), smul(2 * Zl, pW))),
                            add(mul(pM, add(smul(2 * Zr, pM), smul(K, pN))), smul(-2 * Q, mul(pM, pW))))
B12 = smul(-2, add(mul(pW, pW), mul(pM, pM)))
A21 = lambda Zl, Zr, Q: add(mul(pM, add(smul(K, pN), smul(-2 * Zl, pM))),
                            add(mul(pW, add(smul(K, pN), smul(-2 * Zr, pW))), smul(-2 * Q, mul(pM, pW))))
B21 = smul(2, add(mul(pM, pM), mul(pW, pW)))
AE1 = lambda Zr, Q: add(smul(K, mul(pMW, pM)), add(mul(pN, add(smul(K, pN), smul(-2 * Zr, pW))),
                                                   smul(-2 * Q, mul(pN, pM))))
BE1 = smul(2, mul(pN, pW))
AE2 = lambda Zr, Q: add(smul(K, mul(pMW, pW)), add(mul(pN, add(smul(2 * Zr, pM), smul(K, pN))),
                                                   smul(-2 * Q, mul(pN, pW))))
BE2 = smul(-2, mul(pN, pM))
A1E = lambda Zl, Q: add(mul(pN, add(smul(K, pN), smul(2 * Zl, pW))), add(smul(K, mul(pMW, pM)),
                                                                        smul(-2 * Q, mul(pN, pM))))
B1E = smul(-2, mul(pN, pW))
A2E = lambda Zl, Q: add(mul(pN, add(smul(K, pN), smul(-2 * Zl, pM))), add(smul(K, mul(pMW, pW)),
                                                                         smul(-2 * Q, mul(pN, pW))))
B2E = smul(2, mul(pN, pM))

def leafOK(menu, bx, l, r, Q):
    T0, T1, C0, C1 = bx
    # Lean: `j - 1` is natural subtraction, so side 0 reads menu[0]; only used when l, r != 0
    Zl = menu[max(l - 1, 0)] if max(l - 1, 0) < len(menu) else 0
    Zr = menu[max(r - 1, 0)] if max(r - 1, 0) < len(menu) else 0
    if l == 0:
        if r == 0:
            return False
        return comboOK(AE1(Zr, Q), BE1, *bx) and comboOK(AE2(Zr, Q), BE2, *bx)
    if r == 0:
        return comboOK(A1E(Zl, Q), B1E, *bx) and comboOK(A2E(Zl, Q), B2E, *bx)
    return (comboOK(A11(Zl, Zr, Q), [], *bx) and comboOK(A22(Zl, Zr, Q), [], *bx) and
            comboOK(A12(Zl, Zr, Q), B12, *bx) and comboOK(A21(Zl, Zr, Q), B21, *bx))

def flagPoly(aN, C1):
    return add(smul(2 * (aN - C1), mul(pN, pN)), add(smul(K, mul(pN, pMW)), smul(-2 * K, mul(pM, pW))))

def pLow(C0, S0, S1):
    return 2 * C0 <= K or chk(add(smul(K, pMW), smul(-2 * C0, pN)), S0, S1)

# ---- load -----------------------------------------------------------------------------------------
d = json.load(open(CERT))
h, a = F(d['h']), F(d['a'])
hN, aN = int(h * K), int(a * K)
assert F(hN, K) == h and F(aN, K) == a
menu = []
for z in d['menu']:
    Z = F(z) * K
    assert Z.denominator == 1
    menu.append(int(Z))
boxes = []
for b in d['boxes']:
    t0, t1, c0, c1 = (F(x) for x in b)
    T0, T1, C0, C1 = t0 * TR, t1 * TR, c0 * K, c1 * K
    assert all(x.denominator == 1 for x in (T0, T1, C0, C1))
    boxes.append((int(T0), int(T1), int(C0), int(C1)))
m = len(boxes)
rule = [[x + 1 for x in row] for row in d['rule']]
claims = []
for L in d['leaves']:
    row = []
    for l, r, q in L:
        Q = F(q) * K
        assert Q.denominator == 1 and Q >= 0
        row.append((l + 1, r + 1, int(Q)))
    assert len({(l, r) for l, r, _ in row}) == len(row)
    claims.append(row)

# ---- params, rule ---------------------------------------------------------------------------------
assert (2 * aN + K) ** 2 <= 8 * K ** 2 and 2 * hN < 2 * aN + K and aN < K
for b in range(m):
    for b2 in range(m):
        R = rule[b][b2]
        assert R != 0 and 2 * hN < 2 * menu[R - 1] + K and menu[R - 1] < K

# ---- flags: unflagged boxes must pass the quartic test ----------------------------------------------
flags = []
for (T0, T1, C0, C1) in boxes:
    flags.append(not chk(flagPoly(aN, C1), T0, T1))
print('flagged boxes:', [i for i, f in enumerate(flags) if f])

# ---- leaves ---------------------------------------------------------------------------------------
bad = [(b, e) for b in range(m) for e in claims[b] if not leafOK(menu, boxes[b], *e)]
assert not bad, bad[:10]
print('leaves: all', sum(map(len, claims)), 'pass')

# ---- coverage -------------------------------------------------------------------------------------
ts = sorted({0, TR} | {b[0] for b in boxes} | {b[1] for b in boxes})
cover = []
for S0, S1 in zip(ts, ts[1:]):
    cand = [i for i, b in enumerate(boxes) if b[0] <= S0 and S1 <= b[1]]
    first = [i for i in cand if pLow(boxes[i][2], S0, S1)]
    assert first, (S0, S1)
    b0 = max(first, key=lambda i: boxes[i][3])
    chain, cur = [b0], boxes[b0][3]
    while cur < hN:
        nxt = [i for i in cand if boxes[i][2] <= cur]
        bn = max(nxt, key=lambda i: boxes[i][3])
        assert boxes[bn][3] > cur, ('coverage gap', S0, S1, cur)
        chain.append(bn)
        cur = boxes[bn][3]
    cover.append((S0, S1, chain))
print('coverage:', len(cover), 't-strips')

# ---- chain tables ---------------------------------------------------------------------------------
v = [{(l, r): q for l, r, q in L} for L in claims]
phi1 = {}
for b1 in range(m):
    for b2 in range(m):
        z = rule[b1][b2]
        key = (b2, z, flags[b1])
        phi1[key] = min(phi1.get(key, 10 ** 30), v[b1][(0, z)])

def step(phi):
    nxt = {}
    for (b, z, f), s in phi.items():
        for b2 in range(m):
            z2 = rule[b][b2]
            key = (b2, z2, f or flags[b])
            nxt[key] = min(nxt.get(key, 10 ** 30), s + v[b][(z, z2)])
    return nxt

phi2 = step(phi1)
phi3 = step(phi2)
mn = min(s + v[b][(z, 0)] for (b, z, f), s in phi3.items() if f or flags[b])
assert mn >= 4 * K, mn
print('min chain value', F(mn, K), '=', mn / K)

def tab(phi):
    rows = [[] for _ in range(m)]
    for (b, z, f), s in sorted(phi.items()):
        rows[b].append((z, f, s))
    return rows

# ---- emit -----------------------------------------------------------------------------------------
def lst(xs):
    return '[' + ', '.join(xs) + ']'

def tri(e):
    return f'({e[0]}, {e[1]}, {e[2]})'

def trib(e):
    return f'({e[0]}, {"true" if e[1] else "false"}, {e[2]})'

import hashlib
sha = hashlib.sha256(open(CERT, 'rb').read()).hexdigest()
o = []
o.append('import Sqpack.StripCheck\n')
o.append('/-!\n# The wall-strip certificate at `h = 117/100`, as data (generated; do not edit)\n')
o.append(f'Generated by `lean/scripts/gen_strip117_data.py` from\n'
         f'`tasks/n12-chat-lemmas/strip-demo/certs/h1.170_fine.json`, sha256\n`{sha}`.\n'
         f'{m} boxes, {len(menu)} menu heights, {sum(map(len, claims))} leaf claims, '
         f'{sum(flags)} H-flagged boxes, {len(cover)} coverage strips; chain tables with '
         f'{len(phi1)} / {len(phi2)} / {len(phi3)} entries; minimum chain value {mn}/{K}.\n'
         'Heights and claims are over `K = 10^6`, `t = tan(θ/2)` over `TR = 80000`; side `0` is the wall,\n'
         'side `j ≥ 1` the menu height `j`.\n-/\n')
o.append('set_option linter.style.longLine false\nset_option maxRecDepth 20000\n')
o.append('namespace SquarePacking.Strip117Data\n')
o.append('open StripCheck\n')
o.append(f'def menu : List ℕ := {lst(map(str, menu))}\n')
o.append('def boxes : List (ℕ × ℕ × ℕ × ℕ) :=\n  ' + lst(f'({a}, {b}, {c}, {e})' for a, b, c, e in boxes) + '\n')
o.append('def rule : List (List ℕ) :=\n  [' + ',\n   '.join(lst(map(str, row)) for row in rule) + ']\n')
o.append('def flags : List Bool := ' + lst('true' if f else 'false' for f in flags) + '\n')
for b, row in enumerate(claims):
    o.append(f'def claims{b} : List (ℕ × ℕ × ℕ) :=\n  ' + lst(map(tri, row)) + '\n')
o.append('def claims : List (List (ℕ × ℕ × ℕ)) :=\n  ' + lst(f'claims{b}' for b in range(m)) + '\n')
o.append('def cover : List (ℕ × ℕ × List ℕ) :=\n  [' + ',\n   '.join(
    f'({S0}, {S1}, {lst(map(str, ch))})' for S0, S1, ch in cover) + ']\n')
for name, phi in (('phi1', phi1), ('phi2', phi2), ('phi3', phi3)):
    for b, row in enumerate(tab(phi)):
        o.append(f'def {name}_{b} : List (ℕ × Bool × ℕ) := ' + lst(map(trib, row)) + '\n')
    o.append(f'def {name} : List (List (ℕ × Bool × ℕ)) :=\n  ' + lst(f'{name}_{b}' for b in range(m)) + '\n')
o.append(f'''/-- The certificate. -/
def cert : Cert where
  hN := {hN}
  aN := {aN}
  menu := menu
  boxes := boxes
  rule := rule
  flags := flags
  claims := claims
  cover := cover
  phi1 := phi1
  phi2 := phi2
  phi3 := phi3
''')
o.append('end SquarePacking.Strip117Data\n')
open(OUT, 'w').write('\n'.join(o))
print('wrote', OUT)
