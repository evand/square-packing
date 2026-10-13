#!/usr/bin/env python3
"""Extend data/lower_bounds.json from n <= 100 to n <= 324 (Evan, 10-12: the site counts everything to 324).

Past 100 there are no per-n proofs in the literature; every floor is a general bound or a family:
  area-bound        s(n) >= sqrt(n)                                    (every n)
  trivial-square    s(k^2) = k                                         (perfect squares)
  Nagamochi2005     Thm 2(i) s(k^2-1) = s(k^2-2) = k, Thm 2(ii) general  status 'gap' (Lemma 1 is false; notes §2026-10-01)
  Karakus2026       Cor 1.2 s(k^2-1) = k (k >= 3); Cor 6.2 s(N) >= 1/2 + sqrt(N - floor(sqrt N) + 1/4), nonsquare N >= 8
  chelokot2026      s(k^2-2) = k, k >= 2 (Lean)
  EvanDaniel2026    s(k^2-3) = k, k >= 6 (k2m3);  s(k^2-4) = k, k >= 5 (k2m4; `awaiting` until jlevy/squares verifies it)
  wand125_2026      s(101) >= 257/25 (mixed_n101_L1028), and so for 102..105 by monotonicity, as jlevy/squares lists it
  Green2000         the survey's values as jlevy/squares carries them past 100 (DS7 unavoidable sets, monotonicity):
                    status 'claimed', as for n <= 100
`register: 'verified'` marks the entries jlevy/squares has replayed, read from its case files (packing/frontier/n-NNN.md,
`verified_lower_bound`, in the store's sparse clone search/packer/runs/jlevy).  `best` is the largest entry that is not
gap / claimed / reported, ties to the earliest (then the more specific), as for n <= 100.

Idempotent: rewrites n = 101..324 and leaves n <= 100 alone.  Run from site/:  python3 tools/lower_bounds_324.py
"""
import json, math, os, sys
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.dirname(HERE)
REG = os.path.join(SITE, '..', 'search', 'packer', 'runs', 'jlevy', 'packing', 'frontier')
P = os.path.join(SITE, 'data', 'lower_bounds.json')
NAGA = 'https://www.combinatorics.org/ojs/index.php/eljc/article/view/v12i1r37'
KARA = 'https://arxiv.org/abs/2609.37410'
CHEL = 'https://github.com/chelokot/square-packing-archive/blob/main/docs/nagamochi-compensation-proof.md'
GAP = ("published proof incomplete: it rests on Nagamochi's Lemma 1, which is false (chelokot 2026, Lean-checked "
       "counterexample; Karakus 2026, arXiv:2609.37410).  Not counted as a floor.")
AWAIT = ("jlevy/squares has not yet replayed the s(k^2-4) = k family (T-081); counted on the leaderboard once it is "
         "register-verified")


def register(n):
    t = open(os.path.join(REG, f'n-{n:03d}.md')).read()
    return yaml.safe_load(t.split('---')[1])['packing']


def entries(n, R):
    k = math.ceil(math.sqrt(n) - 1e-12); c = k * k - n; a = math.sqrt(n)
    H = [{'value': a, 'exact_form': f'sqrt({n})', 'source': 'area-bound', 'year': 1975, 'status': 'proved', 'url': None}]
    if c == 0:
        H.append({'value': float(k), 'exact_form': str(k), 'source': 'trivial-square', 'year': 1975, 'status': 'proved',
                  'url': None, 'note': f'n={n}=={k}^2, trivial grid packing meets area bound'})
        return H
    if c in (1, 2):
        H.append({'value': float(k), 'exact_form': str(k), 'source': 'Nagamochi2005', 'year': 2005, 'status': 'gap', 'url': NAGA,
                  'note': f'n={n}=={k}^2-{c}; Nagamochi 2005 Thm 2(i); ' + GAP})
    else:
        f = n - 2 * math.isqrt(n) + 1
        H.append({'value': math.sqrt(f) + 1, 'exact_form': f'sqrt({n} - 2*floor(sqrt({n})) + 1) + 1 = sqrt({f})+1',
                  'source': 'Nagamochi2005', 'year': 2005, 'status': 'gap', 'url': NAGA,
                  'note': 'general lower bound, Theorem 2(ii); ' + GAP})
    fl = math.isqrt(n)
    if c != 1:   # at k^2 - 1 Cor 6.2 reads exactly k: Cor 1.2 below is the same value, stated as such
      H.append({'value': min(float(k), 0.5 + math.sqrt(n - fl + 0.25)), 'exact_form': f'1/2 + sqrt({n} - {fl} + 1/4)',
              'source': 'Karakus2026', 'year': 2026, 'status': 'preprint', 'url': KARA, 'date': '2026-09-29',
              'note': 'general lower bound, Corollary 6.2: s(N) >= 1/2 + sqrt(N - floor(sqrt N) + 1/4) for nonsquare N >= 8 '
                      '(weaker than Nagamochi\'s stated bound; a conventional proof, unrefereed)'})
    if c == 1:
        H.append({'value': float(k), 'exact_form': str(k), 'source': 'Karakus2026', 'year': 2026, 'status': 'preprint',
                  'url': KARA, 'date': '2026-09-29',
                  'note': f'n = {k}^2 - 1: Corollary 1.2, s(k^2 - 1) = k for k >= 3, from a strip-measure rectangle bound '
                          'that does not use Nagamochi\'s Lemma 1; unrefereed'})
    if c == 2:
        H.append({'value': float(k), 'exact_form': str(k), 'source': 'chelokot2026', 'year': 2026, 'status': 'preprint',
                  'url': CHEL, 'date': '2026-09-04',
                  'note': f'n = {k}^2 - 2: kernel-checked Lean proof of s(k^2 - 2) = k for every k >= 2 that replaces '
                          'Nagamochi\'s false Lemma 1 by a compensation argument; unrefereed'})
    if c == 3:
        H.append({'value': float(k), 'exact_form': str(k), 'source': 'EvanDaniel2026', 'year': 2026, 'status': 'preprint',
                  'url': 'https://github.com/evand/square-packing/tree/main/certificates/k2m3', 'date': '2026-09-29',
                  'note': f'2026-09-29; s({n}) = {k} as the k = {k} case of s(k^2-3) = k for all k >= 6 (one fixed-profile '
                          'family of measures, its 7 x 7 box certified by an exact checker and re-checked by wand125\'s '
                          'separately written checker; the reduction to every k kernel-checked in Lean); unrefereed'})
    if c == 4:
        H.append({'value': float(k), 'exact_form': str(k), 'source': 'EvanDaniel2026', 'year': 2026, 'status': 'preprint',
                  'url': 'https://github.com/evand/square-packing/tree/main/certificates/k2m4', 'date': '2026-10-03',
                  'note': f'2026-10-03; s({n}) = {k} as the k = {k} case of s(k^2-4) = k for all k >= 5 (for k >= 8 one '
                          'fixed-profile family, its 9 x 9 box certified by one exact checker; the reduction to every k '
                          'kernel-checked in Lean, the tilted part of the 9 x 9 box the one unformalised step); unrefereed',
                  'awaiting': AWAIT})
    r = R.get('reported_lower_bound') or {}
    if 'wand125' in (r.get('proved_by') or []):
        H.append({'value': 257 / 25, 'exact_form': '257/25', 'source': 'wand125_2026', 'year': 2026, 'status': 'preprint',
                  'url': 'https://github.com/wand125/square-packing/tree/main/problems/square-lower-bounds', 'date': '2026-10-02',
                  'note': ('mixed point / segment / rectangle certificate mixed_n101_L1028 (total mass 10099999/100000)'
                           + ('' if n == 101 else f'; s({n}) >= s(101) by monotonicity, as jlevy/squares lists it'))})
    if 'Trevor Green' in (r.get('proved_by') or []):
        H.append({'value': float(r['value']), 'exact_form': r.get('exact_form'), 'source': 'Green2000', 'year': 2000,
                  'status': 'claimed', 'url': None,
                  'note': ('T. Green, private communication (Friedman DS7), as jlevy/squares carries it to this n ('
                           + str(r.get('kind')) + '); no proof to review, so claimed (notes 2026-09-30); not a floor')})
    # register-verified: the entry whose value jlevy/squares replayed
    v = R.get('verified_lower_bound') or {}
    ev = ' '.join(v.get('evidence') or [])
    who = ('Karakus2026' if 'karakus' in ev else 'chelokot2026' if 'chelokot' in ev else 'EvanDaniel2026' if 'evand' in ev
           else 'wand125_2026' if 'wand125' in ev else None)
    if who and v.get('value') is not None:
        vv = float(v['value'])
        for h in H:
            if h['source'] == who and abs(h['value'] - vv) < 1e-9 and not h.get('awaiting'):
                h['register'] = 'verified'
    H.sort(key=lambda h: (h.get('year') or 0, h.get('date') or '', h['value']))
    return H


def best_of(H):
    ok = [h for h in H if h['status'] in ('proved', 'preprint') and not h.get('reported')]
    top = max(h['value'] for h in ok)
    tied = [h for h in ok if h['value'] > top - 1e-12]
    spec = [h for h in tied if h['source'] not in ('area-bound',)] or tied
    return dict(min(spec, key=lambda h: (h.get('year') or 0, h.get('date') or '')))


def main():
    raw = open(P).read(); L = json.loads(raw)
    tail = {k: L.pop(k) for k in ('sources', 'meta')}
    for n in range(101, 325):
        H = entries(n, register(n))
        L[str(n)] = {'best': best_of(H), 'history': H}
    S = tail['sources']
    m = tail['meta']
    m['description'] = m['description'].replace('n=1..100', 'n=1..324')
    m['updated'] = '2026-10-12'
    m['extended'] = ('n = 101..324 added 2026-10-12 by tools/lower_bounds_324.py: general bounds and the families only '
                     '(no per-n proofs exist past 100 except wand125\'s s(101)); register flags from jlevy/squares\' case '
                     'files at the store\'s clone.  See notes/lower-bounds-notes.md, 2026-10-12.')
    L.update(tail)
    out = json.dumps(L, indent=1, ensure_ascii=False) + ('\n' if raw.endswith('\n') else '')
    open(P, 'w').write(out)
    cnt = {}
    for n in range(101, 325):
        b = L[str(n)]['best']; cnt[b['source']] = cnt.get(b['source'], 0) + 1
    print('lower bounds to 324; best past 100 by source:', cnt)


if __name__ == '__main__':
    sys.exit(main())
