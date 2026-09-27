# Literature check: s(32), the k^2 - 4 family, exact values (2026-09-26)

Bottom line: **s(32) = 6 is open.** No proof anywhere: journals, arXiv, the record catalogue, or
the 2026 GitHub certificate repos.  No member of `n = k^2 - 4` (k >= 4) has a known exact value
(s(12), s(21), s(32), s(45), ...).  A proof of s(32) = 6 would be the first exact k^2 - 4 value.
**But a competing repo is climbing on s(32) and s(45) right now** (wand125, rectangle-density
certificates: s(32) >= 5.95 and s(45) >= 6.945, both posted 2026-09-25/26; details in §3).

## 1. Proven exact values (n <= 50, non-square n)

| n | s(n) | proof | status |
|---|---|---|---|
| 2, 3 | 2 | Göbel 1979 | refereed |
| 5 | 2 + 1/√2 | Göbel 1979 | refereed |
| 6 (=3²-3) | 3 | Kearney–Shiu, EJC 9 (2002) #R14 | refereed |
| 7, 8 | 3 | Friedman, DS7 | survey |
| 10 | 3 + 1/√2 = 3.7071 | Stromquist, EJC 10 (2003) #R8 | refereed (yes, 10 is proved) |
| 13 (=4²-3) | 4 | Bentz, EJC 17 (2010) #R126 | refereed, but see caveat |
| 14, 15 | 4 | Friedman DS7 (14 is n²-2, also Nagamochi) | |
| 22, 33 (=5²-3, 6²-3) | 5, 6 | Bentz, arXiv:1606.03746 (v1 2016, PDF dated Oct 2018) | **preprint only, no journal ref**; Lean-formalized by chelokot (below) |
| 23, 24; 34, 35; 47, 48; 62, 63; ... | k | Nagamochi, EJC 12 (2005) #R37: s(k²-1) = s(k²-2) = k for all k | refereed, but see caveat |
| 46 (=7²-3) | 7 | Bentz 2010 | refereed |

- s(k²-3) = k: proved for k = 3..7 only (Bentz conjectures all k >= 3).  Open for k >= 8 (s(61), ...).
- s(k²-4): nothing proved.  (k = 3 fails: s(5) < 3.)  Open at k = 4 (s(12)), 5, 6, 7, ...
- Wikipedia's *Square packing* lists only k², k²-1, k²-2, 5, 6, 10, 13, 46 — omits 22, 33 (preprint).
- No 2023–2026 paper (arXiv or journal) proves a new exact value.  Recent arXiv work is asymptotic
  waste (Bui arXiv:2504.09489, 2508.04603; McClenagan arXiv:2602.01484) or Erdős #106.

Caveats found (not in our README):
- **Bentz s(13)**: chelokot (Lean) shows two printed auxiliary point sets (case R2 first subcase, R3)
  are *avoidable*; repaired with modified sets; theorem stands, formally checked.
  https://github.com/chelokot/square-packing-archive/blob/main/docs/bentz-13-formalization.md
- **Nagamochi 2005 Lemma 1** is false (explicit square of side 1.0001 in [0,4]² scoring 0.9775 < 1,
  Lean-checked); s(n²-2) = n re-proved by a replacement argument (square containers only).
  https://github.com/chelokot/square-packing-archive/blob/main/docs/nagamochi-score-counterexample.md
- **Stromquist's s(11) >= 2 + 4/√5**: jlevy (T-010) reports the printed Figure 14 unavoidable set has a
  strict counterexample; bound repaired with a new point set.  https://github.com/jlevy/squares

## 2. Friedman / Ellsworth "Squares in Squares"

Friedman's page https://erich-friedman.github.io/packing/squinsqu/ now just redirects to David
Ellsworth's https://kingbird.myphotos.cc/packing/squares_in_squares.html .  There:
- **n = 21, 32, 45 are not pictured**; the page states "For the n ≤ 324 not pictured, the trivial
  packing (with no tilted squares) is the best known packing."  So best known = 5, 6, 7; none marked proved.
- Marked "Proved": 2,3 (Göbel), 5 (Göbel), 6 (Kearney–Shiu 2001), 7,8 (Friedman), 10 (Stromquist 2003),
  13 (Bentz 2009), 14,15 (Friedman), 22 (Bentz Oct 2018), 23 (Nagamochi), 24 (Friedman), 33 (Bentz Oct 2018),
  34 (Nagamochi), 35 (Friedman), 46 (Bentz 2009), 47,48, 62,63, ... (Nagamochi).  The page carries no lower bounds.
- DS7 (last rev. 2009; copies at https://erich-friedman.github.io/papers/squares/squares.html and
  https://kingbird.myphotos.cc/packing/squares.html) Table 2: s(21) >= 4.7438 (Friedman); **no entry
  for 32 or 42–46**; nearest: 31 >= 5.6415 (Green), 40–41 >= 6.4061 (Green).

## 3. Lower bounds for s(21), s(32), s(45) (all upper bounds = trivial grid, no packing below it known)

| n | classical | Nagamochi general | 2026 (unrefereed, GitHub) | ours |
|---|---|---|---|---|
| 21 | 4.7438 (DS7) | 1+√14 = 4.7417 | jlevy 122/25 = 4.88 (09-23); **wand125 249/50 = 4.98 (09-26)** | **5000/1001 = 4.995005** (still best) |
| 32 | 5.6415 via s(31) | **1+√23 = 5.7958** | **wand125 119/20 = 5.95** (09-26; ladder 5.82 → 5.95 over 09-25/26) | — |
| 45 | 6.4061 via s(41) | **1+√34 = 6.8310** | **wand125 1389/200 = 6.945** (09-25/26) | — |

Nagamochi: s(N) >= min(⌈√N⌉, √(N - 2⌊√N⌋ + 1) + 1); at n = k²-4 this is 1 + √((k-3)(k+1)).

- wand125/square-packing-bounds https://github.com/wand125/square-packing-bounds — weighted point
  certificates plus "rectangle-density" certificates (tokoharu's format: D4-symmetrized uniform
  densities on axis-aligned rectangles, LP with row/column generation, total mass scaled to n - 1/100,
  exact verifier).  Very active (dozens of commits 09-24..09-26); expect further rungs.
- tokoharu/square-packing-density-bounds https://github.com/tokoharu/square-packing-density-bounds — the
  rectangle-density method (s(11) >= 3.81, s(26) >= 5.508, s(29) >= 5.71).
- jlevy/squares https://github.com/jlevy/squares — AI-agent project, survey of n <= 324 with
  reported-vs-verified bounds; its n-032 / n-045 records (reviewed 2026-08-24) list Nagamochi only;
  its n-021 record claims a method ceiling ⌈√n⌉·0.9977 = 4.9885 for *its* certificate format (our 4.995
  exceeds this, so that ceiling is format-specific).

Observation (ours, not from a source): an absolutely continuous density cannot certify at L = k exactly:
the k² grid tiles partition its mass, forcing mass >= k² > n.  Rectangle-density certificates are
therefore capped strictly below 6 for n = 32; exact s(32) = 6 needs point masses on tile boundaries /
closed semantics plus a limit argument, as in our zero-margin s(13) cover.  They can, however, get
arbitrarily close in principle, and are cheap — a "5.99" from them would take the shine off a 6.

## 4. Weighted points / unavoidable sets for lower bounds

- Classical unavoidable points: Göbel 1979; Friedman DS7 §5; Stromquist 2003; Kearney–Shiu 2002.
- Continuously varying unavoidable sets: Bentz arXiv:1606.03746 (s(22), s(33)).
- Pure unavoidable staggered lattices: Bentz 2010 s(46) uses a 45-point lattice (7 rows alternating
  6/7 points); chelokot generalizes (s(23): 3×4 + 2×5 = 22 pts; s(34): 3×5 + 3×6 = 33 pts).  For n = k²-4
  one needs k²-5 points (31 for n = 32); the row-spacing constraint (Δy² + 1/4 <= 1) keeps pure lattices
  at k²-4 points or more, which is presumably why k²-4 has resisted (our reading).
- Weighted resources (points, segments, areas): Nagamochi 2005 — closest classical precedent.
- Weighted/LP point covers, 2026 (all unrefereed, mostly AI-assisted): Burns (s(17) >= 4.4811), Massaccesi
  (4.5058), Fort, Mira, anabologyco-maker (4.5705), jlevy (T-017 s(12) >= 3.96 on 2026-09-04; s(11) >= 3.827),
  Guzhou0806 and Kleddamag (s(17) > 4.6200; **Kleddamag s(11) > 31/8 = 3.875**, Sept 2026, replayed by jlevy),
  tokoharu and wand125 (rectangle densities, §3).
- Lean: chelokot/square-packing-archive https://github.com/chelokot/square-packing-archive — kernel-checked
  s(n²-2) = n, s(6), s(10), s(13), s(22), s(33), s(46) (+47, 48, 23, 34).

## 5. Gaps vs. our README / status.md

- README "Credits" cites none of jlevy, tokoharu, wand125, Kleddamag, Guzhou0806, chelokot.
- **s(11) claim is superseded**: ours 3.8143 beats tokoharu's 3.81 but is below jlevy's 3.827 and
  Kleddamag's 3.875 (Sept 2026).  The README calls ours an improvement on Stromquist and should acknowledge these.
- s(12): "previous published bound 3.788854" is outdated — jlevy T-017 s(12) >= 99/25 = 3.96 (2026-09-04,
  unrefereed) predates our publication; ours (3.9686) is still the best.
- s(21): status.md says "previous best 4.7438 (DS7)"; actually jlevy 4.88 (09-23) and wand125 4.98 (09-26);
  ours 4.995 still best.
- Bentz s(13) printed-set errors and Nagamochi Lemma 1 counterexample (both chelokot) worth a line where we
  cite those proofs (proof-anatomy.md dissects Bentz).
- Bentz 22/33 is an arXiv preprint (no journal ref found) — cite as such.
