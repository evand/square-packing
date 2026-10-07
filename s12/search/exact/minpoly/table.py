"""The exact-forms table for the register records: EXACT_FORMS.md (readers) and exact_forms.json (machines).

Sources: minpoly/data/n-N.minpoly.json.gz (committed exact forms), minpoly/data/verify.txt (verify_exact.py replay),
minpoly/results.json (cross-checks), lean_batch/results.tsv (Lean `Packs n S*`), the band-lemma build list
(Lean local minima of the integer-side records) and the grade-A instances (N11L, N28L).

  python3 minpoly/table.py [--chain-ok FILE]     (run from search/exact)
"""
import argparse, gzip, json, os, re
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.dirname(HERE)


def poly_str(p, var='x'):
    """Integer coefficients, ascending -> readable descending string."""
    terms = []
    for e in range(len(p) - 1, -1, -1):
        c = p[e]
        if c == 0:
            continue
        sgn = '-' if c < 0 else '+'
        a = abs(c)
        mon = '' if e == 0 else (var if e == 1 else f'{var}^{e}')
        coef = str(a) if (a != 1 or e == 0) else ''
        terms.append((sgn, coef + mon))
    s = ''.join(f' {sg} {t}' for sg, t in terms).strip()
    return s[2:] if s.startswith('+ ') else '-' + s[2:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--chain-ok', default='', help='file with one n per line: band-lemma local minima built in Lean')
    A = ap.parse_args()
    R = {r['n']: r for r in json.load(open(os.path.join(HERE, 'results.json')))}
    ver = {}
    vf = os.path.join(HERE, 'data', 'verify.txt')
    if os.path.exists(vf):
        for line in open(vf):
            m = re.search(r'n-(\d+)\.minpoly\.json\.gz: (VALID)', line)
            if m:
                ver[int(m.group(1))] = True
    lean = {}
    tsv = os.path.join(EX, 'lean_batch', 'results.tsv')
    if os.path.exists(tsv):
        for line in open(tsv):
            n, deg, ch, tg, tt, st = line.rstrip('\n').split('\t')
            lean[int(n)] = st                                  # latest wins
    chain = set()
    if A.chain_ok and os.path.exists(A.chain_ok):
        chain = {int(x) for x in open(A.chain_ok).read().split()}
    gradeA = {11: 'N11L', 28: 'N28L'}
    rows, js = [], []
    for n in range(1, 325):
        f = os.path.join(HERE, 'data', f'n-{n}.minpoly.json.gz')
        r = R.get(n, {})
        if not os.path.exists(f):
            js.append({'n': n, 'status': 'open', 'reason': r.get('status', 'not run')[:200]})
            rows.append(f'| {n} | open | | | | | | | {r.get("status", "not run")[:60]} |')
            continue
        D = json.load(gzip.open(f, 'rt'))
        p = [int(F(c)) for c in D['S']['poly']]
        a, b = D['S']['interval']
        fdeg = len(D['field']['f']) - 1
        Sdec = r.get('S', '')
        lm = 'band lemma' if n in chain else (f'grade A ({gradeA[n]})' if n in gradeA else '')
        lp = 'yes' if lean.get(n) == 'ok' else ('' if n not in lean else 'failed')
        cc = r.get('crosscheck', '')
        js.append({'n': n, 'status': 'exact', 'S': Sdec, 'S_poly_ascending': p, 'S_interval': [a, b],
                   'S_degree': len(p) - 1, 'height_digits': max(len(str(abs(c))) for c in p),
                   'field_poly_ascending': D['field']['f'], 'field_degree': fdeg,
                   'verify_exact': bool(ver.get(n)), 'lean_packs': lp == 'yes', 'lean_local_min': lm,
                   'crosscheck': cc, 'data': f'minpoly/data/n-{n}.minpoly.json.gz'})
        ps = poly_str(p) if len(p) <= 9 and max(len(str(abs(c))) for c in p) <= 12 else f'degree {len(p) - 1}, {max(len(str(abs(c))) for c in p)}-digit coefficients (JSON)'
        rows.append(f'| {n} | {str(Sdec)[:22]} | {len(p) - 1} | `{ps}` | {fdeg} | {"✓" if ver.get(n) else ""} | '
                    f'{lp} | {lm} | {cc} |')
    ok = [j for j in js if j['status'] == 'exact']
    hdr = [
        '# Exact forms of the best-known packings, n ≤ 324',
        '',
        'For each register record: `S*`, the side of the packing, is the unique root of an irreducible integer polynomial '
        '`p` in a rational interval, and an exact configuration over a number field `K = ℚ[t]/f` realises it '
        '(`minpoly/data/n-N.minpoly.json.gz`).  Derived from the contact graph (`minpoly.py`), checked independently '
        'by `verify_exact.py` (stdlib only: contacts as identities mod `f`, every other pair and wall by rational '
        'intervals).  Machine-readable: `exact_forms.json`.',
        '',
        f'* exact forms: **{len(ok)}** of 324 (independent replay VALID: {sum(1 for j in ok if j["verify_exact"])}); '
        f'open: {", ".join(str(j["n"]) for j in js if j["status"] == "open")}',
        f'* Lean `Packs n S*` (kernel-checked, standard axioms; `lean_batch.py`): {sum(1 for j in ok if j["lean_packs"])}',
        f'* Lean local minimum (`IsLocalMinPacking`): {sum(1 for j in ok if j["lean_local_min"])} '
        '(band lemma: integer-side records; grade A: n = 11, 28)',
        '',
        '| n | S* | deg p | p | field deg | verify_exact | Lean Packs | Lean local min | cross-check |',
        '|---|---|---|---|---|---|---|---|---|']
    open(os.path.join(EX, 'EXACT_FORMS.md'), 'w').write('\n'.join(hdr + rows) + '\n')
    json.dump(js, open(os.path.join(EX, 'exact_forms.json'), 'w'), indent=1)
    print(f'{len(ok)} exact; table written')


if __name__ == '__main__':
    main()
