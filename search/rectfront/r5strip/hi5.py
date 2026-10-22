#!/usr/bin/env python3
"""Final5Hi generator (2026-10-21): the old Final5 pipeline (../r5: cert5, gen5, tight5) restricted to the high strip
2r - 1 <= H < 2 (r = sqrt 2; low strips are `rect5_low`).

Changes from ../r5:
* H ranges: 'lo' (target 5) is [2r - 1, h0], 'hi' (target 4 + 6r - 4H) is [h0, 2], h0 = (6r - 1)/4.
* Non-tight leaves are linear in H: each is proved from its two H-endpoint polynomials (variables u, v, r only),
  with edge factors (u, 1 - u, 1 + u^2, ...) stripped before the LP certificate (cert5.certify_split).
* Tight leaves (the zigzag): tight5's (1 + u)^m G = R + k Q decomposition, generators adapted.

Usage:  hi5.py collect | certify NPROC | emit
Cache: runs/rectfront/r5hi/plans.pkl.  Output: lean/Sqpack/Rect5Hi{Leaves,Tight,CasesFF,..,TT}.lean, Rect5Final.lean.
"""
import itertools, os, pickle, re, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'r5'))
import sympy as sp
import cert5
from cert5 import parts, certify, certify_split, lean_expr, lean_hint, hyp_name, u, v, H, r, H0, strip_edges, R_LO, R_HI
import gen5
from gen5 import L, tag, positive_den, sub_r, emit_poly_leaf

ROOT = '/home/evand/math/square-packing/public'
CACHE = f'{ROOT}/runs/rectfront/r5hi/plans.pkl'
LEAN = f'{ROOT}/lean/Sqpack'
HL = 2 * r - 1
CASES = [(e1, e2, b1, b2, c) for e1, e2, b1, b2 in itertools.product((False, True), (False, True), 'UP', 'UP')
         for c in ('up0', 'up90', 'a1', 'a2')]
TIGHT_CORES = []
for _e1, _e2, _mp in ((True, False, False), (False, True, True)):
    for _cand in ('a1', 'a2'):
        for _kp in (True, False):
            TIGHT_CORES.append(((_e1, _e2, 'P', 'P', _cand), _mp, _kp))


def gens_hi(region, Hint):
    G = [g for g in cert5.gens_for(region, 'lo') if g[2] != 'H']
    if Hint == 'lo':
        H_G = [('H-hl', H - HL, 'H'), ('h0-H', H0 - H, 'H')]
    else:
        H_G = [('H-h0', H - H0, 'H'), ('2-H', 2 - H, 'H')]
    rG = [g for g in G if g[2] == 'r']
    uvG = [g for g in G if g[2] == 'uv']
    return uvG + H_G + rG


gen5.gens_for = gens_hi           # emit_case uses it for both the region facts and the leaf generators

def _positive_den(D):
    c, fl = sp.factor_list(D)
    sign = 1 if c > 0 else -1
    terms, proofs = [abs(c)], ['(by norm_num)']
    for f, m in fl:
        fp, pr = gen5.pos_factor_proof(f)
        pr = {u: 'hu', v: 'hv'}.get(sp.expand(fp), '(by first | positivity | linarith | nlinarith [hu, hv, mul_pos hu hv])')
        if sp.expand(fp - f) != 0: sign *= (-1) ** m
        terms.append(fp ** m if m > 1 else fp); proofs.append(pr if m == 1 else f'(pow_pos {pr} {m})')
    Dpos = sp.Mul(*terms, evaluate=False)
    pf = proofs[0]
    for p_ in proofs[1:]: pf = f'(mul_pos {pf} {p_})'
    return sign, Dpos, pf
gen5.positive_den = _positive_den

_parts0 = cert5.parts
def _parts_exact(*a):
    """cert5.parts with exact coefficients (its upright candidate (1, 0) gives Python floats in (1 + 0)/2)."""
    E, region = _parts0(*a)
    q = lambda e: sp.nsimplify(e, rational=True) if e.atoms(sp.Float) else e
    return q(E), [q(g) for g in region]
cert5.parts = _parts_exact
gen5.parts = _parts_exact


def hrange(Hint):
    return (HL, H0) if Hint == 'lo' else (H0, sp.Integer(2))


def tight_names():
    out = set()
    for case, mp, kp in TIGHT_CORES:
        t = tag(case)
        sub = f"{'p' if mp else 'n'}{'p' if kp else 'n'}"
        out.add(f'f5h_{t}_{sub}_hi'); out.add(f'f5h_{t}_{sub}_lo')
    return out


def fix_case_text(txt):
    """repairs of gen5.emit_case's text (never compiled before 10-21)."""
    txt = txt.replace('psi_branch_pos hC10 ', 'psi_branch_pos hC₁0 ').replace('psi_branch_pos hC20 ', 'psi_branch_pos hC₂0 ')
    for w in 'uv':
        txt = re.sub(rf"(  have h{w}1' : {w} < 1 := half_lt_one [^\n]*\n)",
                     rf"\1  have h{w}n1 : (1:ℝ) - {w} ≠ 0 := (by linarith : (0:ℝ) < 1 - {w}).ne'\n"
                     rf"  have h{w}n2 : (1:ℝ) - {w} ^ 2 ≠ 0 := (by nlinarith : (0:ℝ) < 1 - {w} ^ 2).ne'\n", txt)
    return txt


def collect():
    """all leaves (lname, P, Gh) of all case theorems, the case texts, and the tight set."""
    leaves, texts = [], {}
    for case in CASES:
        certs = {}
        txt = gen5.emit_case(case, certs, 'f5h')
        txt = txt.replace('(hH1 : 1 ≤ H) (hH2 : H < 2)\n', '(hH1 : 1 ≤ H) (hH2 : H < 2) (hHl : 2 * Real.sqrt 2 - 1 ≤ H)\n')
        txt = re.sub(r'^theorem f5_', 'theorem f5h_', txt, flags=re.M)
        texts[case] = fix_case_text(txt)
        for lname, P, Gh, c, sub, Hint in certs['need']:
            leaves.append((lname, sp.expand(P), Gh, Hint))
    return leaves, texts


def endpoint_polys(P, Hint):
    assert sp.degree(P, H) <= 1, 'not linear in H'
    Ha, Hb = hrange(Hint)
    return sp.expand(P.subs(H, Ha)), sp.expand(P.subs(H, Hb))


def base_gens(Gh):
    return [g for g in Gh if g[2] != 'H']


def ekey(core, G0):
    return (sp.srepr(core), tuple((g[0], sp.srepr(g[1])) for g in G0))


def endpoint_jobs(leaves):
    tn = tight_names(); jobs = {}
    for lname, P, Gh, Hint in leaves:
        if lname in tn: continue
        G0 = base_gens(Gh)
        for Pe in endpoint_polys(P, Hint):
            fac, core = strip_edges(Pe)
            if core.is_number:
                continue
            jobs.setdefault(ekey(core, G0), (core, G0))
    return jobs


def work(item):
    key, (core, G0) = item
    t = time.time()
    plan = certify_split(core, G0)
    return key, plan, time.time() - t


def cmd_certify(nproc):
    from multiprocessing import Pool
    leaves, _ = collect()
    jobs = endpoint_jobs(leaves)
    cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
    todo = [(k, v) for k, v in jobs.items() if k not in cache]
    print('endpoint cores:', len(jobs), ' to do:', len(todo), flush=True)
    with Pool(nproc, maxtasksperchild=1) as pool:
        for i, (key, plan, dt) in enumerate(pool.imap_unordered(work, todo)):
            cache[key] = plan
            pickle.dump(cache, open(CACHE, 'wb'))
            print(f'{i + 1}/{len(todo)} {"OK " + str(len(plan)) + " box" if plan else "FAIL"} {dt:.0f}s', flush=True)


ROBUST = set(x for x in os.environ.get('HI5_ROBUST', 'f5h_FTPU_a1_pp_lo_b_core,f5h_FTPU_a1_pp_hi_b_core').split(',') if x)

def robust_plan(core, plan):
    """merge the supports of several LP solutions per box (degrees 6, 8, 10, random costs) for linarith."""
    import numpy as np
    out = []
    for box, used, Gb in plan:
        seen = {lab for lab, _ in used}; allu = list(used)
        for deg in (6, 8, 10):
            cols, labels, target = cert5.certify_cols(core, Gb, deg)
            monos = sorted({m for c in cols for m in c} | set(target)); idx = {m: i for i, m in enumerate(monos)}
            A = np.zeros((len(monos), len(cols)))
            for j, c in enumerate(cols):
                for m, val in c.items(): A[idx[m], j] = val
            b = np.zeros(len(monos))
            for m, val in target.items(): b[idx[m]] = val
            rng = np.random.default_rng(deg)
            for trial in range(3):
                cost = np.ones(len(cols)) if trial == 0 else rng.uniform(0.5, 2, len(cols))
                res = cert5.linprog(cost, A_eq=A, b_eq=b, bounds=[(0, None)] * len(cols), method='highs')
                if res.status != 0: break
                for j in range(len(cols)):
                    if res.x[j] > 1e-12 and labels[j] not in seen:
                        seen.add(labels[j]); allu.append((labels[j], res.x[j]))
        out.append((box, allu, Gb))
    return out

def nonneg_fac_proof(fac):
    """explicit proof term of 0 <= fac, fac a product of powers of u, v, 1-u, 1-v, 1+u, 1+v, 1+u^2, 1+v^2."""
    c, fl = sp.factor_list(sp.expand(fac))
    terms = []
    pf = {u: 'g_u', v: 'g_v', 1 - u: '(by linarith : (0:ℝ) ≤ 1 - u)', u - 1: None,
          1 - v: '(by linarith : (0:ℝ) ≤ 1 - v)', v - 1: None, u + 1: '(by linarith : (0:ℝ) ≤ u + 1)',
          v + 1: '(by linarith : (0:ℝ) ≤ v + 1)', u**2 + 1: '(by positivity : (0:ℝ) ≤ u ^ 2 + 1)',
          v**2 + 1: '(by positivity : (0:ℝ) ≤ v ^ 2 + 1)'}
    sign = 1
    for f, m in fl:
        f = sp.expand(f)
        if f in (u - 1, v - 1):
            sign *= (-1) ** m; f = -f
        t = pf[f]
        terms.append(t if m == 1 else f'(pow_nonneg {t} {m})')
    assert sign * c > 0, fac
    out = f'(by norm_num : (0:ℝ) ≤ {lean_expr(abs(c))})'
    for t in terms: out = f'(mul_nonneg {out} {t})'
    return f'(by have h := {out}; convert h using 1; ring)'



def emit_endpoint_leaf(lname, P, Gh, Hint, cache, cores):
    """the leaf lemma (statement as gen5 expects) from two endpoint lemmas."""
    G0 = base_gens(Gh)
    Ha, Hb = hrange(Hint)
    out = []
    names = []
    for side, Pe in zip('ab', endpoint_polys(P, Hint)):
        fac, core = strip_edges(Pe)
        nm = f'{lname}_{side}'
        names.append((nm, Pe))
        hyps = ' '.join(f'({hyp_name(g[0])} : 0 ≤ {lean_expr(g[1])})' for g in G0)
        if core.is_number:
            assert core >= 0
            out.append(f'theorem {nm} (u v H r : ℝ) {hyps} :\n    0 ≤ {lean_expr(Pe)} := by\n'
                       f'  have e : {lean_expr(Pe)} = {lean_expr(core)} * ({lean_expr(fac)}) := by ring\n'
                       f'  rw [e]; exact mul_nonneg (by norm_num) ({nonneg_fac_proof(fac)})\n')
            continue
        k = ekey(core, G0)
        regs = [g for g in G0 if g[0].startswith('reg')]
        diag = any(sp.expand(a[1] + b[1]) == 0 and sp.expand(a[1] - (v - u)) == 0 for a in regs for b in regs)
        if diag:
            # degenerate sign region v = u: certify the one-variable polynomial
            core1 = sp.expand(core.subs(v, u))
            G1 = [g for g in G0 if g[0] in ('u', '1-u', 'r-rlo', 'rhi-r')]
            k1 = ('diag', sp.srepr(core1))
            plan1 = cache.get(k1) or certify_split(core1, G1)
            if not plan1: return None
            cache[k1] = plan1
            cn = cores.get(k)
            if cn is None:
                cn = f'{nm}_core'; cores[k] = cn
                out.append(emit_poly_leaf(cn + '_1', core1, G1, plan1))
                hyps0 = ' '.join(f'({hyp_name(g[0])} : 0 ≤ {lean_expr(g[1])})' for g in G0)
                out.append(f'theorem {cn} (u v H r : ℝ) {hyps0} :\n    0 ≤ {lean_expr(core)} := by\n'
                           f'  have huv : v = u := by linarith\n'
                           f'  have h := {cn}_1 u u H r g_u g_1u g_rlo g_rhi\n  rw [huv]; linarith [h]\n')
            args = ' '.join(hyp_name(g[0]) for g in G0)
            if sp.expand(fac - 1) == 0:
                out.append(f'theorem {nm} (u v H r : ℝ) {hyps} :\n    0 ≤ {lean_expr(Pe)} := by\n'
                           f'  have h := {cn} u v H r {args}\n  have e : {lean_expr(Pe)} = {lean_expr(core)} := by ring\n'
                           f'  rw [e]; exact h\n')
            else:
                out.append(f'theorem {nm} (u v H r : ℝ) {hyps} :\n    0 ≤ {lean_expr(Pe)} := by\n'
                           f'  have h := {cn} u v H r {args}\n'
                           f'  have e : {lean_expr(Pe)} = ({lean_expr(fac)}) * ({lean_expr(core)}) := by ring\n'
                           f'  rw [e]; exact mul_nonneg ({nonneg_fac_proof(fac)}) h\n')
            continue
        plan = cache.get(k)
        Gc, extra = G0, ''
        if not plan:
            # degenerate sign region: reg0 >= 0 together with -reg0 = (u + v)(1 - u v) >= 0 on the box
            reg0 = [g for g in G0 if g[0] == 'reg0']
            if not reg0 or sp.expand(-reg0[0][1] - (u + v) * (1 - u * v)) != 0:
                return None
            Gc = G0 + [('nreg0', sp.expand(-reg0[0][1]), 'uv')]
            k = ('nreg0',) + k
            plan = cache.get(k) or certify_split(core, Gc)
            if not plan:
                return None
            cache[k] = plan
            extra = (f'  have g_nreg0 : 0 ≤ {lean_expr(Gc[-1][1])} := by\n'
                     f'    have e : {lean_expr(Gc[-1][1])} = (u + v) * (1 - u * v) := by ring\n'
                     f'    rw [e]; exact mul_nonneg (by linarith) (by nlinarith)\n')
        cn = cores.get(k)
        if cn is None:
            cn = f'{nm}_core'; cores[k] = cn
            if cn in ROBUST: plan = robust_plan(core, plan)
            out.append(emit_poly_leaf(cn, core, Gc, plan))
        args = ' '.join(hyp_name(g[0]) for g in Gc)
        if sp.expand(fac - 1) == 0:
            out.append(f'theorem {nm} (u v H r : ℝ) {hyps} :\n    0 ≤ {lean_expr(Pe)} := by\n' + extra +
                       f'  have h := {cn} u v H r {args}\n  have e : {lean_expr(Pe)} = {lean_expr(core)} := by ring\n'
                       f'  rw [e]; exact h\n')
        else:
            out.append(f'theorem {nm} (u v H r : ℝ) {hyps} :\n    0 ≤ {lean_expr(Pe)} := by\n' + extra +
                       f'  have h := {cn} u v H r {args}\n'
                       f'  have e : {lean_expr(Pe)} = ({lean_expr(fac)}) * ({lean_expr(core)}) := by ring\n'
                       f'  rw [e]; exact mul_nonneg ({nonneg_fac_proof(fac)}) h\n')
    hyps = ' '.join(f'({hyp_name(g[0])} : 0 ≤ {lean_expr(g[1])})' for g in Gh)
    args = ' '.join(hyp_name(g[0]) for g in G0)
    (na, Pa), (nb, Pb) = names
    hlo, hhi = [hyp_name(g[0]) for g in Gh if g[2] == 'H']
    out.append(f'theorem {lname} (u v H r : ℝ) {hyps} :\n    0 ≤ {lean_expr(P)} := by\n'
               f'  have ha := {na} u v H r {args}\n  have hb := {nb} u v H r {args}\n'
               f'  have e : ({lean_expr(Hb - Ha)}) * ({lean_expr(P)}) = ({lean_expr(Hb)} - H) * ({lean_expr(Pa)}) + '
               f'(H - ({lean_expr(Ha)})) * ({lean_expr(Pb)}) := by ring\n'
               f'  have hd : 0 < {lean_expr(Hb - Ha)} := by linarith\n'
               f'  have h1 : 0 ≤ {lean_expr(Hb)} - H := by linarith\n  have h2 : 0 ≤ H - ({lean_expr(Ha)}) := by linarith\n'
               f'  have : 0 ≤ ({lean_expr(Hb - Ha)}) * ({lean_expr(P)}) := by rw [e]; exact add_nonneg (mul_nonneg h1 ha) (mul_nonneg h2 hb)\n'
               f'  exact (mul_nonneg_iff_of_pos_left hd).1 this\n')
    return '\n'.join(out)


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'collect'
    if cmd == 'collect':
        leaves, texts = collect()
        tn = tight_names()
        print('leaves:', len(leaves), ' distinct:', len({(sp.srepr(P), str(Gh)) for _, P, Gh, _ in leaves}),
              ' tight:', sum(1 for l in leaves if l[0] in tn))
        jobs = endpoint_jobs(leaves)
        print('endpoint cores:', len(jobs))
    elif cmd == 'certify':
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        cmd_certify(int(sys.argv[2]) if len(sys.argv) > 2 else 4)
    elif cmd == 'shard':                 # shard i n: certify jobs i mod n sequentially into plans.<i>.pkl
        i, n = int(sys.argv[2]), int(sys.argv[3])
        leaves, _ = collect(); jobs = sorted(endpoint_jobs(leaves).items(), key=lambda kv: kv[0][0])
        base = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
        out = f'{CACHE[:-4]}.{i}.pkl'
        mine = pickle.load(open(out, 'rb')) if os.path.exists(out) else {}
        todo = [kv for j, kv in enumerate(jobs) if j % n == i and kv[0] not in base and kv[0] not in mine]
        print('shard', i, 'todo', len(todo), flush=True)
        for j, (key, (core, G0)) in enumerate(todo):
            t = time.time(); plan = certify_split(core, G0); mine[key] = plan
            pickle.dump(mine, open(out, 'wb'))
            print(f'{j + 1}/{len(todo)} {"OK " + str(len(plan)) + " box" if plan else "FAIL"} {time.time() - t:.0f}s', flush=True)
    elif cmd == 'merge':
        import glob
        base = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
        for fn in glob.glob(f'{CACHE[:-4]}.*.pkl'):
            base.update(pickle.load(open(fn, 'rb')))
        pickle.dump(base, open(CACHE, 'wb'))
        print('cache', len(base), 'ok', sum(1 for v in base.values() if v), 'fail', [1 for v in base.values() if not v].__len__())


def emit_corner_leaf(lname, P, Gh):
    """degenerate sign region reg0 = (u + v)(u v - 1) >= 0: only the corners u = v = 0, u = v = 1, where P = 0."""
    reg0 = [g for g in Gh if g[0] == 'reg0']
    if not reg0 or sp.expand(reg0[0][1] - (u + v) * (u * v - 1)) != 0: return None
    if sp.expand(P.subs({u: 0, v: 0})) != 0 or sp.expand(P.subs({u: 1, v: 1})) != 0: return None
    hyps = ' '.join(f'({hyp_name(g[0])} : 0 ≤ {lean_expr(g[1])})' for g in Gh)
    Pt = lean_expr(P)
    return (f'theorem {lname} (u v H r : ℝ) {hyps} :\n    0 ≤ {Pt} := by\n'
            f'  have hz : (u + v) * (u * v - 1) = 0 := by\n'
            f'    have h1 : 0 ≤ (u + v) * (u * v - 1) := by\n'
            f'      have e : {lean_expr(reg0[0][1])} = (u + v) * (u * v - 1) := by ring\n'
            f'      rw [← e]; exact g_reg0\n'
            f'    have h2 : (u + v) * (u * v - 1) ≤ 0 := mul_nonpos_of_nonneg_of_nonpos (by linarith) (by nlinarith)\n'
            f'    linarith\n'
            f'  rcases mul_eq_zero.1 hz with h | h\n'
            f'  · have hu0 : u = 0 := by linarith\n    have hv0 : v = 0 := by linarith\n'
            f'    have e : {Pt} = 0 := by rw [hu0, hv0]; ring\n    rw [e]\n'
            f'  · have hu0 : u = 1 := by nlinarith\n    have hv0 : v = 1 := by rw [hu0] at h; linarith\n'
            f'    have e : {Pt} = 0 := by rw [hu0, hv0]; ring\n    rw [e]\n')


# ======================================================================== emission
HDR = '/-! Generated by `search/rectfront/r5strip/hi5.py`; do not edit by hand. -/\n\n'
NLEAF = 4


def cmd_emit():
    import tight5
    tight5.gens_for = gens_hi
    leaves, texts = collect()
    cache = pickle.load(open(CACHE, 'rb'))
    tn = tight_names()
    canon, leaf_txt, cores, missing = {}, [], {}, []
    for lname, P, Gh, Hint in leaves:
        key = (sp.srepr(P), tuple((g[0], sp.srepr(g[1])) for g in Gh))
        if key in canon or lname in tn:
            continue
        canon[key] = lname
        t = emit_endpoint_leaf(lname, P, Gh, Hint, cache, cores)
        if t is None: t = emit_corner_leaf(lname, P, Gh)
        if t is None: missing.append(lname)
        else: leaf_txt.append(t)
    pickle.dump(cache, open(CACHE, 'wb'))
    # case texts: deduplicated leaf names, tight leaves get `r ^ 2 = 2`
    case_text = {}
    rename = {}
    for lname, P, Gh, Hint in leaves:
        key = (sp.srepr(P), tuple((g[0], sp.srepr(g[1])) for g in Gh))
        if lname not in tn: rename[lname] = canon[key]
    for case in CASES:
        txt = texts[case]
        txt = re.sub(r'\bf5h_\w+', lambda m: rename.get(m.group(0), m.group(0)), txt)
        for nm in tn:
            txt = re.sub(rf'(have hP := {nm} u v H \(Real\.sqrt 2\)[^\n]*)', r'\1 (Real.sq_sqrt (by norm_num))', txt)
        txt = re.sub(r'^theorem f5h_', 'set_option maxHeartbeats 4000000 in\ntheorem f5h_', txt, flags=re.M)
        grp = ('T' if case[0] else 'F') + ('T' if case[1] else 'F')
        case_text.setdefault(grp, []).append(txt)
    # tight
    tight_txt, tbad = [], []
    for case, mp, kp in TIGHT_CORES:
        t = tag(case); sub = f"{'p' if mp else 'n'}{'p' if kp else 'n'}"
        txt, err = tight5.emit_tight(f'f5h_{t}_{sub}_hi', f'f5h_{t}_{sub}_lo', case, mp, kp)
        if txt is None: tbad.append((t, sub, err))
        else: tight_txt.append(txt)
    print('leaves', len(leaf_txt), 'missing', missing, 'tight bad', tbad, flush=True)
    pre = 'namespace SquarePacking\n\nset_option linter.unusedVariables false\n\n'
    nper = -(-len(leaf_txt) // NLEAF)
    opts = 'set_option maxRecDepth 20000\nset_option maxHeartbeats 1000000\n\n'
    for i in range(NLEAF):
        imp = 'import Mathlib\n' if i == 0 else f'import Sqpack.Rect5HiLeaves{i - 1}\n'
        with open(f'{LEAN}/Rect5HiLeaves{i}.lean', 'w') as f:
            f.write(imp + '\n' + HDR + pre + opts + '\n'.join(leaf_txt[i * nper:(i + 1) * nper]) + '\nend SquarePacking\n')
    with open(f'{LEAN}/Rect5HiTight.lean', 'w') as f:
        f.write('import Mathlib\n\n' + HDR + '/-! Tight leaves (the zigzag): `(1 + u)ᵐ G = R + k Q`, see `../r5/tight5.py`. -/\n\n'
                + pre + '\n'.join(tight_txt) + '\nend SquarePacking\n')
    imports = f'import Sqpack.Rect5HiLeaves{NLEAF - 1}\n'
    for grp, txts in case_text.items():
        with open(f'{LEAN}/Rect5HiCases{grp}.lean', 'w') as f:
            f.write('import Sqpack.Rect5HiGlue\n' + imports + 'import Sqpack.Rect5HiTight\n\n' + HDR + pre +
                    'set_option linter.unusedSimpArgs false\n\n' + '\n'.join(txts) + '\nend SquarePacking\n')
    order = [(e1, e2, b1, b2, c) for e1 in (False, True) for e2 in (False, True) for b1 in 'UP' for b2 in 'UP'
             for c in ('up0', 'a1', 'a2', 'up90')]
    bullets = ''.join(f'  · exact f5h_{tag(c)} u v H L s₁ s₂ hu hu1 hv hv1 hH1 hH2 hHl (by simpa using h1) '
                      f'(by simpa using h2) hmid\n' for c in order)
    with open(f'{LEAN}/Rect5Final.lean', 'w') as f:
        f.write('import Sqpack.Rect5HiCasesFF\nimport Sqpack.Rect5HiCasesFT\nimport Sqpack.Rect5HiCasesTF\n'
                'import Sqpack.Rect5HiCasesTT\n\n' + HDR + 'namespace SquarePacking\n\n'
                '/-- **`Final5Hi`**: the real inequality left by the chain for five squares in a strip of height\n'
                '`2√2 - 1 ≤ H < 2`.  Half-angle coordinates for both lines, then one case theorem per sign / branch /\n'
                'middle-square candidate (`Rect5HiCases*.lean`), each from polynomial leaves (`Rect5HiLeaves*.lean`,\n'
                '`Rect5HiTight.lean`). -/\n'
                'theorem final5_high : Final5Hi := by\n'
                '  intro ε₁ ε₂ C₁ S₁ C₂ S₂ H L s₁ s₂ x hC₁ hS₁ hCS₁ hC₂ hS₂ hCS₂ hH1 hH2 hHl ob od hx hmid\n'
                '  obtain ⟨u, hu0, hu1, rfl, rfl⟩ := half_angle hC₁ hS₁.le hCS₁\n'
                '  obtain ⟨v, hv0, hv1, rfl, rfl⟩ := half_angle hC₂ hS₂.le hCS₂\n'
                '  have hu : 0 < u := by\n    rcases eq_or_lt_of_le hu0 with h | h\n    · subst h; simp at hS₁\n    · exact h\n'
                '  have hv : 0 < v := by\n    rcases eq_or_lt_of_le hv0 with h | h\n    · subst h; simp at hS₂\n    · exact h\n'
                '  cases ε₁ <;> cases ε₂ <;> rcases ob with h1 | h1 <;> rcases od with h2 | h2 <;>\n'
                '    rcases hx with rfl | rfl | rfl | rfl\n' + bullets + '\nend SquarePacking\n')


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'emit':
    cmd_emit()


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'onecase':
    # emit one case theorem (leaves from the built Rect5HiLeaves*) into TmpHiCase.lean
    e1, e2, b1, b2, c = sys.argv[2] == 'T', sys.argv[3] == 'T', sys.argv[4], sys.argv[5], sys.argv[6]
    leaves, texts = collect()
    import subprocess
    txt = texts[(e1, e2, b1, b2, c)]
    # leaf names as emitted (canonical): reuse emit's renaming
    canon, rename = {}, {}
    tn = tight_names()
    for lname, P, Gh, Hint in leaves:
        key = (sp.srepr(P), tuple((g[0], sp.srepr(g[1])) for g in Gh))
        if lname in tn: continue
        canon.setdefault(key, lname); rename[lname] = canon[key]
    txt = re.sub(r'\bf5h_\w+', lambda m: rename.get(m.group(0), m.group(0)), txt)
    for nm in tn:
        txt = re.sub(rf'(have hP := {nm} u v H \(Real\.sqrt 2\)[^\n]*)', r'\1 (Real.sq_sqrt (by norm_num))', txt)
    txt = re.sub(r'^theorem f5h_', 'set_option maxHeartbeats 4000000 in\ntheorem f5h_', txt, flags=re.M)
    open(f'{LEAN}/TmpHiCase.lean', 'w').write('import Sqpack.Rect5HiGlue\nimport Sqpack.Rect5HiLeaves3\nimport Sqpack.Rect5HiTight\n\n'
        'namespace SquarePacking\n\nset_option linter.unusedVariables false\nset_option linter.unusedSimpArgs false\n\n' + txt + '\nend SquarePacking\n')
    print('written')
