#!/usr/bin/env python3
# SUPERSEDED 2026-10-25 by leg_emit.py (line lemmas + pair table); kept for history; leg_emit reuses band_emit leaf helpers and band_cert.
"""Emit the Lean leaves, cells and cores of the band theorem (2026-10-25).

Inputs: band_gen (pieces, leaves), band_cert cache (certificates).  Output (lean/Sqpack/):
  StripBandLeaves{k}.lean  leaf theorems: split the max atoms, then `linarith` with certificate products
  StripBandCells.lean      cell theorems (nested plx / plx_inf) and the five core theorems band_Sa, band_Sb,
                           band_C0, band_C1, band_N.
Generated; do not edit by hand.  Regenerate: python3 band_emit.py"""
import os, sys, pickle, itertools
import sympy as sp
import band_gen as G
import band_cert as BC

LEAN = os.path.join(os.path.dirname(__file__), '../../../lean/Sqpack')
C, S, H, r = G.C, G.S, G.H, G.r
NFILES = 4

def L(e):
    s = sp.sstr(sp.sympify(e)).replace('**', '^')
    return f'({s})'

BASE_ARGS = ('(C S H r : ℝ) (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hr2 : r ^ 2 = 2) '
             '(hrlo : (141421356 : ℝ) / 100000000 ≤ r) (hrhi : r ≤ (141421357 : ℝ) / 100000000) '
             '(h9 : (14 * r - 5) / 8 ≤ H) (h87 : H ≤ 187 / 100)')
BASE_CALL = 'C S H r hC hS hCS hr2 hrlo hrhi h9 h87'

GEN_PRELUDE = """  have gp1 : 0 ≤ C + S - 1 := by nlinarith [mul_nonneg hC hS]
  have gC1 : 0 ≤ 1 - C := by nlinarith
  have gS1 : 0 ≤ 1 - S := by nlinarith
  have g9 : 0 ≤ H - (14 * r - 5) / 8 := by linarith
  have g87 : 0 ≤ 187 / 100 - H := by linarith
  have grlo : 0 ≤ r - 141421356 / 100000000 := by linarith
  have grhi : 0 ≤ 141421357 / 100000000 - r := by linarith
  have hCS0 : C ^ 2 + S ^ 2 - 1 = 0 := by linarith
"""
BASE_GEN = {'C': 'hC', 'S': 'hS', 'p1': 'gp1', 'C1': 'gC1', 'S1': 'gS1', 'H9': 'g9', 'H87': 'g87',
            'rlo': 'grlo', 'rhi': 'grhi'}

def prod_term(names, gmap):
    ts = [gmap[n] if ' ' not in gmap[n] else f'({gmap[n]})' for n in names]
    t = ts[0]
    for u in ts[1:]: t = f'(mul_nonneg {t} {u})'
    return t

def mono(m):
    i, j, k, l = m
    fs = []
    for v, e in (('C', i), ('S', j), ('H', k), ('r', l)):
        if e == 1: fs.append(v)
        elif e > 1: fs.append(f'{v} ^ {e}')
    return ' * '.join(fs) if fs else '(1 : ℝ)'

def cert_text(plan, gmap, indent):
    d, used, eqs, gnames = plan[:4]
    hints = [prod_term(lab, gmap) for lab, w in used if lab]
    hints += [f'mul_eq_zero_of_right ({mono(m)}) hCS0' for m in eqs]
    return f'{indent}linarith [{", ".join(hints)}]\n'

def split_text(entry, gmap, indent):
    """Nested case splits on C at the cut values of the certificate (sign -1: C <= v, +1: v <= C)."""
    names, plans = entry
    signs = list(itertools.product((1, -1), repeat=len(names)))
    vals = []
    for k in range(len(names)):
        e = plans[0][4][k][1]                      # extra generator: sign * (C - v)
        sg = 1 if plans[0][4][k][0].endswith('+') else -1
        v = sp.nsimplify(sp.expand(C - sg * e))
        vals.append(f'({sp.sstr(v)} : ℝ)')
    def rec(k, cur, ind, first_on_bullet):
        if k == len(names):
            plan = plans[signs.index(tuple(cur))]
            gm = dict(gmap)
            for n, sg in enumerate(cur):
                gm[names[n] + ('+' if sg > 0 else '-')] = f'sub_nonneg.2 s{n}'
            return cert_text(plan, gm, ind)
        out = f'{ind}rcases le_total C {vals[k]} with s{k} | s{k}\n'
        out += f'{ind}·\n' + rec(k + 1, cur + [-1], ind + '  ', False)
        out += f'{ind}·\n' + rec(k + 1, cur + [1], ind + '  ', False)
        return out
    return rec(0, [], indent, False)

def leaf_hyp_terms(lf):
    """Lean proof terms for the leaf's hypothesis generators h0, h1, ... (kink parts), in band_gen order."""
    terms = []
    xk = lf.tag.startswith('xk'); yk = '_yk' in lf.tag
    if xk:
        terms += ['hCp.le', 'sub_nonneg.2 hk1']
        if any(k == 'ge' for k, _ in lf.hyps[2:3]) and len([h for h in lf.hyps[:3] if h[0] == 'ge']) == 2:
            terms += ['sub_nonneg.2 hk2']
    if yk:
        terms += ['hSp.le', 'sub_nonneg.2 hl1']
        nyge = len([h for h in lf.hyps if h[0] == 'ge']) - (len([t for t in terms if t.startswith('sub')]) - 1)
    # recompute exactly from lf.hyps to be safe
    terms = []
    seen_pos = []
    for k, h in lf.hyps:
        if k == 'pos':
            if sp.simplify(h - C) == 0: terms.append('hCp.le'); seen_pos.append('x')
            else: terms.append('hSp.le'); seen_pos.append('y')
        else:
            side = seen_pos[-1]
            if side == 'x':
                terms.append('sub_nonneg.2 hk1' if terms.count('sub_nonneg.2 hk1') == 0 else 'sub_nonneg.2 hk2')
            else:
                terms.append('sub_nonneg.2 hl1' if terms.count('sub_nonneg.2 hl1') == 0 else 'sub_nonneg.2 hl2')
    return terms

def emit_leaf(name, stmt, extra_args, lf, cache, maxargs):
    """stmt: Lean statement text; maxargs: list of max-argument texts in band_gen order (M1, M2)."""
    lines = [f'theorem {name} {BASE_ARGS}{extra_args} :\n    {stmt} := by\n', GEN_PRELUDE]
    kterms = leaf_hyp_terms(lf)
    tight = None
    def branch(choice, ind):
        lid = BC.leaf_id(lf, choice)
        gmap = dict(BASE_GEN)
        allterms = kterms + [('neg_nonneg.2 ' if ch == 0 else '') + f'm{k}' for k, ch in enumerate(choice)]
        for n, t in enumerate(allterms): gmap[f'h{n}'] = t
        if lid in BC.TIGHT:
            lemma = 'tight_P1' if lid.startswith(('Sa', 'C0')) else 'tight_P2'
            extra = ' (by linarith)' if lemma == 'tight_P2' else ''
            return f'{ind}nlinarith [{lemma} hC hS hCS hr2 (by linarith){extra}]\n'
        entry = cache.get(lid)
        if entry is None:
            return f'{ind}sorry -- MISSING {lid}\n'
        if entry[0] == 'split':
            return split_text(entry[1], gmap, ind)
        return cert_text(entry, gmap, ind)
    def rec(k, choice, ind):
        if k == len(maxargs):
            lines.append(branch(tuple(choice), ind)); return
        a = maxargs[k]
        lines.append(f'{ind}rcases le_total {a} 0 with m{k} | m{k}\n')
        lines.append(f'{ind}· rw [max_eq_right m{k}]\n')
        rec(k + 1, choice + [0], ind + '  ')
        lines.append(f'{ind}· rw [max_eq_left m{k}]\n')
        rec(k + 1, choice + [1], ind + '  ')
    rec(0, [], '  ')
    return ''.join(lines)

def pos_proof(s):
    return '(by norm_num)' if s.strip('()').replace(' ', '').isdigit() else '(by linarith)'

class Cell:
    def __init__(self, cond, i, j, xp, xs, yp, ys, const):
        self.cond, self.i, self.j = cond, i, j
        A, B = G.xform(xp, xs, const); A2, B2 = G.yform(yp, ys)
        _, (ln, ld), xh = xp; _, (lyn, lyd), yh = yp
        self.A, self.B, self.A2, self.B2 = L(A), L(B), L(A2), L(B2)
        self.ln, self.ld = L(ln), L(ld)
        self.hn, self.hd = (L(xh[0]), L(xh[1])) if xh else (None, None)
        self.lyn, self.lyd = L(lyn), L(lyd)
        self.yn, self.yd = (L(yh[0]), L(yh[1])) if yh else (None, None)
        self.xp, self.yp = xp, yp

def all_cells():
    out = []
    for cond, (xps, xs, yps, ys, const, _, _) in G.CONDS.items():
        for i, xp in enumerate(xps):
            for j, yp in enumerate(yps):
                out.append(Cell(cond, i, j, xp, xs, yp, ys, const))
    return out

def leaves_of_cell(cell):
    return [lf for lf in G.all_leaves() if lf.cond == cell.cond and lf.cell == (cell.i, cell.j)]

def emit_all():
    cache = pickle.load(open(BC.CACHE, 'rb'))
    leaf_texts = []
    cell_texts = []
    leafidx = {}
    allL = G.all_leaves()
    for cell in all_cells():
        c = cell
        Psi = f'({c.A2} + {c.B2} * y + max (S * y - C) 0)'
        lvs = {lf.tag: lf for lf in allL if lf.cond == c.cond and lf.cell == (c.i, c.j)}
        base = f'bandLeaf_{c.cond}_{c.i}{c.j}'
        # x ends: (tag, m text, E text, extra args for x-kink, binder text)
        xends = [('xlo', c.ld, f'(1 * ({c.A} * {c.ld} + {c.B} * {c.ln} + max (C * {c.ln} - S * {c.ld}) 0))', '', '', [f'(C * {c.ln} - S * {c.ld})'])]
        if c.hn:
            xends.append(('xhi', c.hd, f'(1 * ({c.A} * {c.hd} + {c.B} * {c.hn} + max (C * {c.hn} - S * {c.hd}) 0))', '', '', [f'(C * {c.hn} - S * {c.hd})']))
            xk_args = f' (hCp : 0 < C) (hk1 : {c.ln} * C ≤ S * {c.ld}) (hk2 : S * {c.hd} ≤ {c.hn} * C)'
            xk_call = ' hCp hk1 hk2'
        else:
            xk_args = f' (hCp : 0 < C) (hk1 : {c.ln} * C ≤ S * {c.ld})'
            xk_call = ' hCp hk1'
        xends.append(('xk', 'C', f'(1 * ({c.A} * C + {c.B} * S))', xk_args, xk_call, []))
        proof = []
        # outer application
        if c.hn:
            proof.append(f'  refine plx (m := 1) (E := {Psi}) (A := {c.A}) (B := {c.B}) (ln := {c.ln}) (ld := {c.ld}) '
                         f'(hn := {c.hn}) (hd := {c.hd}) hC hS one_pos {pos_proof(c.ld)} {pos_proof(c.hd)} hxl hxh ?_ ?_ ?_ (by ring)\n')
        else:
            proof.append(f'  refine plx_inf (m := 1) (E := {Psi}) (A := {c.A}) (B := {c.B}) (ln := {c.ln}) (ld := {c.ld}) '
                         f'hC hS one_pos {pos_proof(c.ld)} hxl ?_ ?_ ?_ (by ring)\n')
        order = ['xlo', 'xhi', 'xk'] if c.hn else ['xlo', 'xslope', 'xk']
        xmap = {e[0]: e for e in xends}
        for tag in order:
            if tag == 'xslope':
                lf = lvs['xslope']; nm = f'{base}_xslope'
                stmt = f'0 ≤ {c.B} + C'
                leaf_texts.append(emit_leaf(nm, stmt, '', lf, cache, []))
                proof.append(f'  · exact {nm} {BASE_CALL}\n'); continue
            _, m, E, xargs, xcall, xmax = xmap[tag]
            hm = 'hCp' if m == 'C' else pos_proof(m)
            if tag == 'xk':
                proof.append(f'  · intro{xcall}\n')
            else:
                proof.append('  ·\n')
            ind = '    '
            if c.yn:
                proof.append(f'{ind}refine plx (m := {m}) (E := {E}) (A := {c.A2}) (B := {c.B2}) (C := S) (S := C) (ln := {c.lyn}) (ld := {c.lyd}) '
                             f'(hn := {c.yn}) (hd := {c.yd}) hS hC {hm} {pos_proof(c.lyd)} {pos_proof(c.yd)} hyl hyh ?_ ?_ ?_ (by ring)\n')
                yorder = ['ylo', 'yhi', 'yk']
            else:
                proof.append(f'{ind}refine plx_inf (m := {m}) (E := {E}) (A := {c.A2}) (B := {c.B2}) (C := S) (S := C) (ln := {c.lyn}) (ld := {c.lyd}) '
                             f'hS hC {hm} {pos_proof(c.lyd)} hyl ?_ ?_ ?_ (by ring)\n')
                yorder = ['ylo', 'yslope', 'yk']
            for yt in yorder:
                lf = lvs[f'{tag}_{yt}']; nm = f'{base}_{tag}_{yt}'
                if yt == 'ylo':
                    stmt = f'0 ≤ {m} * ({c.A2} * {c.lyd} + {c.B2} * {c.lyn} + max (S * {c.lyn} - C * {c.lyd}) 0) + {E} * {c.lyd}'
                    maxargs = xmax + [f'(S * {c.lyn} - C * {c.lyd})']
                    leaf_texts.append(emit_leaf(nm, stmt, xargs, lf, cache, maxargs))
                    proof.append(f'{ind}· exact {nm} {BASE_CALL}{xcall}\n')
                elif yt == 'yhi':
                    stmt = f'0 ≤ {m} * ({c.A2} * {c.yd} + {c.B2} * {c.yn} + max (S * {c.yn} - C * {c.yd}) 0) + {E} * {c.yd}'
                    maxargs = xmax + [f'(S * {c.yn} - C * {c.yd})']
                    leaf_texts.append(emit_leaf(nm, stmt, xargs, lf, cache, maxargs))
                    proof.append(f'{ind}· exact {nm} {BASE_CALL}{xcall}\n')
                elif yt == 'yslope':
                    stmt = f'0 ≤ {c.B2} + S'
                    leaf_texts.append(emit_leaf(nm, stmt, '', lf, cache, []))
                    proof.append(f'{ind}· exact {nm} {BASE_CALL}\n')
                else:
                    if c.yn:
                        yargs = f' (hSp : 0 < S) (hl1 : {c.lyn} * S ≤ C * {c.lyd}) (hl2 : C * {c.yd} ≤ {c.yn} * S)'
                        ycall = ' hSp hl1 hl2'
                    else:
                        yargs = f' (hSp : 0 < S) (hl1 : {c.lyn} * S ≤ C * {c.lyd})'
                        ycall = ' hSp hl1'
                    stmt = f'0 ≤ {m} * ({c.A2} * S + {c.B2} * C) + {E} * S'
                    leaf_texts.append(emit_leaf(nm, stmt, xargs + yargs, lf, cache, xmax))
                    proof.append(f'{ind}· intro{ycall}; exact {nm} {BASE_CALL}{xcall}{ycall}\n')
        xh = f' (hxh : x * {c.hd} ≤ {c.hn})' if c.hn else ''
        yh = f' (hyh : y * {c.yd} ≤ {c.yn})' if c.yn else ''
        cell_texts.append(
            f'theorem bandCell_{c.cond}_{c.i}{c.j} {BASE_ARGS} (x y : ℝ) (hxl : {c.ln} ≤ x * {c.ld}){xh} '
            f'(hyl : {c.lyn} ≤ y * {c.lyd}){yh} :\n'
            f'    0 ≤ 1 * ({c.A} + {c.B} * x + max (C * x - S) 0) + {Psi} := by\n' + ''.join(proof) + '\n')
    hdr = ('import Sqpack.StripBand\n\n/-! Generated by `search/rectfront/chainn/band_emit.py`; do not edit by hand. -/\n\n'
           'set_option maxHeartbeats 2000000\n\nnamespace SquarePacking\n\n')
    chunk = (len(leaf_texts) + NFILES - 1) // NFILES
    for k in range(NFILES):
        imp = 'import Sqpack.StripBand\n'
        body = '\n'.join(leaf_texts[k * chunk:(k + 1) * chunk])
        open(f'{LEAN}/StripBandLeaves{k}.lean', 'w').write(
            imp + hdr.split('\n', 1)[1] + body + '\nend SquarePacking\n')
    open(f'{LEAN}/StripBandCells.lean', 'w').write(
        ''.join(f'import Sqpack.StripBandLeaves{k}\n' for k in range(NFILES)) + hdr.split('\n', 1)[1] + '\n'.join(cell_texts) + '\nend SquarePacking\n')
    print(len(leaf_texts), 'leaf theorems,', len(cell_texts), 'cells')

if __name__ == '__main__':
    emit_all()


# ---------------------------------------------------------------- core theorems
def LR(e):
    """Lean text with r = Real.sqrt 2 (for statements about the actual potentials)."""
    return L(e).replace('r', '(Real.sqrt 2)') if False else L(sp.sympify(e)).replace('*r', '*Real.sqrt 2').replace('(r', '(Real.sqrt 2').replace(' r', ' Real.sqrt 2').replace('-r', '-Real.sqrt 2')

PIECE_REWRITE = {
    # name: (lemma, needs nonneg var, number of piece hypotheses)
    'A0': 'u0_A0', 'Lz': 'u0_Lz', 'T': 'u0_T', 'O': 'u1_O', 'aT': 'u1_aT', '0': None,
}

def piece_cases(var, pieces, clipped_lo):
    """Yield (index, list of (hyp name, hyp text for the cell), rewrite command, prelude case-split lines) as a tree.
    pieces: [(name, (ln, ld), hi or None)].  Boundaries between consecutive pieces are the hi of the earlier."""
    return pieces

def emit_core(cond, cells):
    xps, xs, yps, ys, const, _, _ = G.CONDS[cond]
    names = {'Sa': ('bandU1 H y ≤ bandU0 H x + (kap C S H x (-y) - (2 * Real.sqrt 2 - H))', ''),
             'Sb': ('bandU0 H y ≤ bandU1 H x + (kap C S H x (-y) - (2 * Real.sqrt 2 - H))', ''),
             'C0': ('2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤ bandU0 H x + (kap C S H x (-y) - (2 * Real.sqrt 2 - H)) + bandU0 H y',
                    ' (hx2 : 1 / 2 ≤ x) (hy2 : 1 / 2 ≤ y)'),
             'C1': ('2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤ bandU1 H x + (kap C S H x (-y) - (2 * Real.sqrt 2 - H)) + bandU1 H y',
                    ' (hx2 : 1 / 2 ≤ x) (hy2 : 1 / 2 ≤ y)'),
             'N': ('1 ≤ kap C S H x (-y)', ' (hx2 : x ≤ 1 / 2)')}
    goal, extra = names[cond]
    out = [f'theorem band_{cond} {{C S H x y : ℝ}} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hb : InBand H)\n'
           f'    (hx : 0 ≤ x) (hy : 0 ≤ y){extra} :\n    {goal} := by\n',
           '  have r1 := sqrt2_lo\'; have r2 := sqrt2_hi\'\n',
           '  rw [kap_opp hC hS hx hy]\n']
    def lo_text(v, piece):
        _, (ln, ld), _ = piece
        return f'{LR(ln)} ≤ {v} * {LR(ld)}'
    def hi_text(v, piece):
        _, _, h = piece
        return f'{v} * {LR(h[1])} ≤ {LR(h[0])}'
    def rewrite(v, piece, hv):
        nm = piece[0]
        if nm == 'A0': return f'rw [u0_A0 hb {hv} (by linarith)]'
        if nm == 'Lz': return f'rw [u0_Lz hb {hv} (by linarith) (by linarith)]'
        if nm == 'T': return f'rw [u0_T hb (by linarith)]'
        if nm == 'O': return f'rw [u1_O hb {hv} (by linarith)]'
        if nm == 'aT': return f'rw [u1_aT hb (by linarith)]'
        return None
    def rec_x(i, ind):
        # case split x into piece i or later
        if i < len(xps) - 1:
            h = xps[i][2]
            out.append(f'{ind}rcases le_total (x * {LR(h[1])}) {LR(h[0])} with hx{i} | hx{i}\n')
            out.append(f'{ind}· '); rec_y_start(i, ind + '  ')
            if i + 1 < len(xps) - 1:
                out.append(f'{ind}·\n'); rec_x(i + 1, ind + '  ')
            else:
                out.append(f'{ind}· '); rec_y_start(i + 1, ind + '  ')
        else:
            rec_y_start(i, ind)
    def rec_y_start(i, ind):
        rw = rewrite('x', xps[i], 'hx')
        first = True
        if rw:
            out.append(f'{rw}\n'); first = False
        rec_y(i, 0, ind, first)
    def rec_y(i, j, ind, first):
        pre = '' if first else ind
        if j < len(yps) - 1:
            h = yps[j][2]
            out.append(f'{pre}rcases le_total (y * {LR(h[1])}) {LR(h[0])} with hy{j} | hy{j}\n')
            out.append(f'{ind}· '); finish(i, j, ind + '  ')
            out.append(f'{ind}· ')
            if j + 1 < len(yps) - 1: rec_y(i, j + 1, ind + '  ', True)
            else: finish(i, j + 1, ind + '  ')
        else:
            if not first or (i == 0 and len(xps) == 1 and xps[0][0] == '0'): out.append(ind)
            finish(i, j, ind)
    def finish(i, j, ind):
        rw = rewrite('y', yps[j], 'hy')
        lines = []
        if rw: lines.append(rw)
        xp, yp = xps[i], yps[j]
        args = [f'(by linarith)']                       # hxl
        if xp[2]: args.append('(by linarith)')          # hxh
        args.append('(by linarith)')                    # hyl
        if yp[2]: args.append('(by linarith)')          # hyh
        lines.append(f'have hc := bandCell_{cond}_{i}{j} C S H (Real.sqrt 2) hC hS hCS sqrt2_sq sqrt2_lo\' sqrt2_hi\' hb.1 hb.2 x y {" ".join(args)}')
        lines.append('try simp only [bandT, bandLz]')
        lines.append('linarith')
        out.append(lines[0] + '\n' + ''.join(f'{ind}{l}\n' for l in lines[1:]))
    rec_x(0, '  ')
    return ''.join(out)

def emit_cores():
    txt = ''.join(emit_core(c, None) + '\n' for c in G.CONDS)
    p = f'{LEAN}/StripBandCells.lean'
    s = open(p).read().replace('\nend SquarePacking\n', '\n' + txt + 'end SquarePacking\n')
    open(p, 'w').write(s)
