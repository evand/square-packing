# Lower bounds for s(n): sourcing notes and caveats

Data file: `data/lower_bounds.json` (copied to `www/data/` by `build.sh bounds`).

`s(n)` is the side length of the smallest square that can contain `n`
non-overlapping unit squares (arbitrary rotations allowed). This file documents
only **lower bounds** on `s(n)` (proofs that no packing fits in a smaller
square) — never best-known packings (upper bounds).

## JSON structure

Top-level object keyed by `"1"`..`"100"`, each holding `{"best": {...},
"history": [...]}` (history oldest-first), plus a `"sources"` dict of full
citations and a `"meta"` dict. Every `n` has at least an `area-bound` entry
(`s(n) >= sqrt(n)`) and, when `n` is not of the form `k^2`, `k^2-1`, `k^2-2`, a
`Nagamochi2005` general-bound entry. `"best"` is simply the entry with the
largest `value` (ties broken toward more specific/named results over the
generic area/formula bounds).

## Primary sources, in the order specified

1. **Erich Friedman's survey**, "Packing Unit Squares in Squares: A Survey and
   New Results." I fetched the maintained copy
   (https://erich-friedman.github.io/papers/squares/squares.html) directly
   with `curl` (WebFetch's HTML→markdown conversion was lossy for the table,
   so raw HTML was parsed instead) and transcribed **Table 2** verbatim
   (all rows, n=2 through 85, with the exact `n`, value, algebraic form, and
   author column). This maintained copy's content matches the EJC Dynamic
   Survey DS7 v5 (2009) "Lower Bounds" section (I independently fetched
   `ds7-2009.html` and cross-checked the theorems), but its reference list
   includes one entry newer than 2009 ("K. Morandi, 2010, private
   communication"), so it is a lightly-updated copy of the v5 text, not
   literally the 2009 PDF.
   - I also fetched the 1998 v1 PDF (`ds7v1-1998.pdf`) to establish *when*
     Friedman's own proofs of s(n) for n=2,3,5,8,15,24,35 first appeared
     (already in v1, 1998) versus n=7,14 (added later — present in the v5/2009
     text as Theorems 7–8 but not claimed as proved in the 1998 abstract).
     I did **not** separately re-verify the 2000 v2 PDF; the 1998 and
     "current" (2009+) fetches were judged sufficient to date every
     attribution used here. This is the main place completeness could be
     improved if someone wants to nail down exact intermediate-version dates.
   - Table 2 gives one value per **range** of `n` (e.g. "17-18", "26-27",
     "28-30"). I applied the listed value to every `n` in the range, per the
     table's own structure.

2. **Exact results (lower bound = upper bound).** Verified year/venue for
   each: Göbel 1979 (n=2,3,5 — in *Packing and Covering in Combinatorics*);
   Kearney & Shiu 2002, EJC 9 #R14 (n=6); El Moumni 1999, *Studia Sci. Math.
   Hungar.* 35 (n=7,8, predating Friedman's/Nagamochi's proofs of the same);
   Stromquist 2003, EJC 10 #R8 (n=10,11,12 — note 12 has since been
   *improved past* the exact-conjectured value 4 downward bound by nothing;
   Stromquist's 11-12 result is `2+4/sqrt(5)`, still short of 4, so 11 and 12
   remain open); Bentz 2010, EJC 17 #R126 (n=13,46, i.e. `m^2-3` for m=4,7);
   Bentz 2016, arXiv:1606.03746 (n=22,33, i.e. `m^2-3` for m=5,6 — I could
   not confirm a separate journal publication for this one beyond arXiv, but
   the proof itself is a conventional rigorous mathematical argument, not a
   computer search, so it is marked `"status": "proved"` here, not
   `"preprint"`); Nagamochi 2005, EJC 12 #R37, Theorem 2(i): `s(n^2) =
   s(n^2-1) = s(n^2-2) = n` for **every** integer `n >= 1` — I fetched the
   actual PDF and applied this exact theorem for every k=1..10 (covering
   n up to 100), which subsumes several of the historically-earlier specific
   results above (Göbel's n=2,3; Friedman's n=8,15,24,35; Kearney-Shiu's
   n=6 is *not* subsumed, since 6 isn't of the n²/n²-1/n²-2 form). Both the
   earlier specific-result citation and the later general Nagamochi citation
   are kept in `history` for these n, since both are independently valid
   proofs; `best` picks whichever is more specific (not `area-bound` or the
   generic per-n `Nagamochi2005` "general lower bound" formula entry) when
   values tie.

3. **General Nagamochi bound (Theorem 2(ii)).** Confirmed by reading the
   actual PDF (`ds7`... no — `combinatorics.org/.../v12i1r37/pdf`): for
   integer `N >= 4` not in `{n^2, n^2-1, n^2-2}`,
   `s(N) >= sqrt(N - 2*floor(sqrt(N)) + 1) + 1`, strictly stronger than the
   trivial `s(N) >= sqrt(N)`. Computed this exactly (via Python `math.isqrt`)
   for every `n` in 1..100 and included it as a `history` entry for every `n`
   not in the exact family, exactly as instructed.
   - **Finding worth flagging:** for `n` = 20, 30, 31, 41, 53, this generic
     formula (published 2005) actually gives a *larger* (better) lower bound
     than the specific ad-hoc value Friedman's Table 2 reports for those `n`
     (attributed to Green/Friedman, ~2000). This happens because Table 2
     reports one value per *range* of n (tuned to the range, e.g. worst case
     over "28-30"), so at the top edge of a range the generic per-n formula
     can overtake it. I used the maximum of the two in each case (i.e., the
     `Nagamochi2005` general-formula entry is `best` for n=20,30,31,41,53
     instead of the historically-reported Green/Friedman value). This is a
     mechanical, correct application of Nagamochi's published theorem, not
     new research on my part — but it does mean five of the `best` values in
     this table are technically *tighter* than the number one would read off
     Friedman's survey table directly for those five n. Flagging in case
     this surprises anyone diffing against the survey.

4. **2026 computer-assisted certificates (unrefereed).** All confirmed by
   fetching the actual pages/READMEs and, for the GitHub repos, the commit
   history via the GitHub API (to get exact dates and confirm the repos and
   claims are real):
   - **n=17**, chronological order of claims found:
     - Sam Burns, blog post, 2026-08-06, `s(17) >= 4.4811` (268-atom weighted
       certificate, produced with ChatGPT/GPT-5.6 Pro; explicitly
       unrefereed).
     - Stanislav Fort, github.com/stanislavfort/17squares, first commit
       2026-08-11, `s(17) > 4.456575` (16 rational witness points; the
       README states the author "doesn't really understand it => can't
       vouch for its correctness").
     - Mira-acc, github.com/Mira-acc/17squares, commit 2026-08-12, `s(17) >
       4.468292` ("exact lower bound" commit message).
     - Gustavo Massaccesi, blog post, 2026-08-21, `s(17) >= 4.5058` (LP over
       168 points in a 29x29 grid; explicitly built on and improves Burns'
       4.4811).
     - Mira-acc, commit 2026-09-08 (today, per system clock), two updates
       same day: `4.607028598640`, then a later commit the same day
       publishing a "reviewed" 1620-atom certificate at
       `s(17) > 4.6130286358861101094200443452...` = the exact form
       `sqrt(17650291964463886688094912400/829429719507765981945905041)`.
       This is the current best claimed value for n=17 and is what `best`
       uses, with `"status": "preprint"`. The repo itself says "external
       peer review remains pending."
     - All five values are recorded in `history` for n=17 in chronological
       order for the full trail.
     - **I did not propagate this n=17 result to n=18 via monotonicity**
       (`s` is non-decreasing, so `s(18) >= s(17)` would technically license
       raising n=18's lower bound too). I left n=18 at the proved
       Green/Friedman value `(40*sqrt(2)+19)/17 ≈ 4.4452` and noted the
       possible implication in that entry's `"note"` field, rather than
       silently inheriting an unrefereed number two hops removed from its
       source. Treat this as a judgment call, easy to revisit.
   - **n=12**: Evan Daniel (evand), github.com/evand/square-packing-12.
     Verified via the GitHub API commit log (all commits 2026-08-24 through
     2026-08-26): `3920/997 = 3.931795` (2026-08-24) → `980/247 = 3.967611`
     (2026-08-26) → `15680/3951 = 3.968616` (2026-08-26, current best,
     1736-point certificate). Repo explicitly states "nothing here has been
     peer reviewed." Recorded all three as `history`, `best` = 3.968616,
     `"status": "preprint"`.

## Friedman's Table 2, n = 82–85 row: closed form and decimal disagree

The published row reads `2√2+(288+12√3)/41 ≈ 9.2667`. That closed form
actually evaluates to `≈10.3598`, not `9.2667` — and 10.3598 is *impossible*
as a lower bound for n=82 since it exceeds the trivial upper bound
`ceil(sqrt(82)) = 10` (82 squares always fit in a 10×10 grid). I confirmed
this isn't a WebFetch/HTML-conversion artifact by re-downloading the raw
page with `curl` and grepping the literal bytes — the inconsistency is in
the source page itself, almost certainly a typo in the closed-form
expression (denominator or a coefficient), not in the decimal. Every other
row of the same "2√2 + .../..." family recomputes to match its stated
decimal exactly, so this one row is an outlier. I used the table's own
stated decimal, `9.2667`, as the `value` for n=82–85 and kept the (flagged
unreliable) closed form only for reference in `exact_form`.

## Coverage of the survey table beyond what Friedman's Table 2 lists

Table 2 (as published) stops at n=85. For n=86–97 there is no specific named
result available from the sources checked, so those fall back to the
Nagamochi 2005 general formula (which is itself a real, published, proved
bound — not a placeholder). For n=98,99,100 the exact Nagamochi Theorem
2(i) family applies directly (k=10), giving s=10 for all three, with 100
additionally trivial as a perfect square.

## Other minor notes / things not independently re-verified

- Green's results (T. Green, private communication, 2000) are cited only via
  Friedman's peer-reviewed survey, per the survey's own reference list — I
  did not attempt to independently trace or contact T. Green. Status
  `"proved"` reflects that the survey (a refereed EJC Dynamic Survey)
  vouches for them, not an independent re-derivation.
- I did not fetch the 2000 (v2) DS7 PDF separately; the 1998 (v1) and
  current maintained copy were sufficient to date every attribution used.
  If finer version-by-version dating of e.g. exactly when the n=13 value
  `3.8437` or n=21 value `4.7438` first appeared is wanted, v2/v3/v4 would
  need checking.
- Bentz's 2016 paper (arXiv:1606.03746) — I could not confirm a journal
  publication venue beyond arXiv; treated as a rigorous (non-preprint-in-the
  "unrefereed" sense) proof regardless, consistent with Bentz's track record
  (his 2010 paper on 13/46 was refereed and published in EJC).
- All `exact_form` fields for the 2026 preprints are `null` except Mira-acc's
  latest n=17 certificate, which explicitly publishes a closed nested-radical
  rational form; the others (Burns, Fort, Massaccesi, evand's n=12
  certificates) are recorded only as the decimal/fraction the source itself
  headlines.
- n = 11: the floor is `3040/797 = 3.814304` (EvanDaniel2026), re-added 2026-09-22 when
  `s12/certificates/s11_lower_3.8143.txt` became public in this repo; it had been reverted to
  Stromquist 2003 at the 2026-09-10 release because the certificate was not yet public.
- 2026-09-28 records pass (dates from each repository's commit history, read via the GitHub API):
  - n = 11: the floor is Kleddamag's `s(11) > 31/8 = 3.875` (Kleddamag/11-squares-certified-bound,
    2026-09-22; replayed in full by jlevy/squares).  The history adds jlevy's T-018 `3.81` (2026-09-04),
    T-026 `3.8264` (2026-09-09) and T-033 `3.8270` (2026-09-22), and tokoharu's `3.81` (2026-09-21).
    Ours (`3040/797`) is dated by its publication, 2026-09-22 (found 2026-08-26), so it is not a step
    on the chart: jlevy's `3.8264` was already public.
  - n = 12: ours (`15680/3951`, public 2026-08-26 in the old `evand/square-packing-12`) is still the
    floor; jlevy's T-017 `99/25 = 3.96` (2026-09-04) is added to the history.
  - n = 17: the floor is Guzhou0806's R067 `233009/50000 = 4.66018` (2026-09-28), a continuation of
    Kleddamag's `466001/100000 = 4.66001` (2026-09-27).  The history keeps the main steps since Mira's
    `4.613029`; the Kleddamag and Guzhou0806 repositories were updated daily that week, so the
    intermediate releases are left to them.
  - Others' certificates are linked as published; this project replays only its own.  Many
    `n = 18 … 91` floors from wand125's rectangle-density certificates (updated daily) are not yet
    in this file.
- 2026-09-30, the `s(k² − 3) = k` family (`EvanDaniel2026`, `s12/certificates/k2m3`, 2026-09-29; write-up `/k2m3/`):
  proved for every `k ≥ 6` by one fixed-profile family whose `k = 7` box is certified by a single exact checker, with
  the all-k reduction kernel-checked in Lean.  An infinite family does not fit a table keyed by `n`, so it enters as
  its concrete cases in range, `n = 33, 46, 61, 78, 97` (`k = 6..10`), and `meta.families` says so.  For 33, 46 and 61
  (already settled by Bentz 2016, Bentz 2010 and our `s(60) = 8`) it is a history entry and `best` keeps the earlier
  proof, per the tie rule above; for 78 and 97 it is `best` (previously Nagamochi's general bound).  The `n = 61`
  s(60)-corollary note no longer credits Bentz with `k = 3` (Kearney–Shiu 2002).
- 2026-09-30 records pass, from jlevy/squares at `5ddb1cd` (its per-n case records `packing/frontier/n-XXX.md`
  and results register), with dates from each upstream repository's commit history (GitHub API, UTC):
  - n = 11 is **settled**: `s(11) = T = 3.8770835900228141773…`, the side of Trump's 1979 packing (root of
    `s^8 - 20s^7 + 178s^6 - 842s^5 + 1923s^4 - 496s^3 - 6754s^2 + 12420s - 6865`), by Queuingtheorydotcom's
    Astra-assisted certificate proof (github.com/Queuingtheorydotcom/11SquaresOptimal, commit `f9e0de7`,
    2026-09-29 UTC), which builds on jlevy's T-026 and Kleddamag's `31/8`.  jlevy/squares replayed it
    independently (T-060, V4/C5; packet `packing/resources/web/n11-optimality-2026-09-29/`).  Computer-assisted
    and unrefereed, so `status` is `preprint` like the other 2026 results; the floor's `value` is the same double
    as the catalogue's upper bound, so the table reads it as settled.  Kleddamag's `31/8` stays in the history.
  - n = 17: the floor is Guzhou0806's R068 `116511/25000 = 4.66044` (2026-09-28, T-043, replayed).  R070
    (`4.6604427`, 2026-09-29) and R071 (`4.66044275`, 2026-09-30) are in that repository and not yet replayed
    anywhere; they are left to it.  jlevy's T-019 `459/100` (2026-09-04) is added to the history: it was the
    floor from 2026-09-04 until Mira's `4.6070` on 2026-09-08.
  - n = 18 … 72: the floors jlevy/squares has *replayed* (V3/V4) and that beat this file become `best`: jlevy's
    own T-030 (18), T-020 (19), T-021 (20); Tokoharu's T-047 (26, 29, 30); wand125's rectangle T-045 (27, 28, 31)
    and point T-044 (39, 40, 41, 52, 53, 55, 56, 68, 69, 70, 71, 72).  Several are one certificate serving
    another count, as jlevy records them: a certificate whose mass is below `m` holds for every `n ≥ m`
    (28 from the n = 27 file, 30 from 29, 52 from 53, 68 from 69), and 41 and 71 are carried from 39 and 70 by
    monotonicity.  These are jlevy's verified floors, taken as listed; no further propagation was done here
    (so, as before, nothing is carried forward from 17 to 18 or from an exact value to the next count).
    Earlier steps are added to the histories where they were real steps (T-019 for 17 and 18, T-020/T-021
    for 20 and 21, Tokoharu's `1377/250` for 27, 28 and `571/100` for 31).
  - **Reported values.**  wand125's rectangle-density bounds (T-046, 48 counts in 18 … 95) and his `s(50) ≥ 37/5`
    (T-048) are graded V0/C0 by jlevy: reported, with the author's checker run, but replayed by no one else.
    They go into `history` with `"reported": true` and a note saying so, and are **never `best`**, even where
    they are larger (at 45 counts they are); `best` stays the largest *non-reported* entry.  The four older
    wand125 rectangle entries at n = 21 and 45 were of the same kind and now carry the flag too.  This makes
    the `best = largest value` rule above read "largest value among entries not flagged `reported`".
  - Exact values kept: n = 21, 32, 45, 60, 61, 78, 97 stay on our own proofs.  n = 45 gains wand125's second,
    point-only certificate (T-054, V4/C3, 2026-09-28, a history entry under the tie rule); n = 21 gains his
    point-only certificate as a reported entry (T-055, V0/C0).  jlevy's register still shows the Nagamochi bound
    at 60, 61, 78 and 97, which it has not taken in yet; wand125's reported values there (7.94, 7.96, 8.955)
    are below our exact values and are not added.
  - Not taken: jlevy treats Green's 2000 values (Friedman's Table 2: 37, 38, 50, 51, 65–67, 82–85) as
    reported, since Green's proofs were a private communication and have not been recovered; this file keeps
    them as `proved` on the survey's authority, as before.  jlevy's record gives the n = 82–85 closed form as
    `94√2/41 + 247/41 = 9.2667`, which matches the decimal and would repair the flagged `exact_form` above;
    not changed here.  wand125's repository kept moving after jlevy's snapshot (e.g. `s(37) ≥ 161/25`,
    `s(87) ≥ 939/100`, a point-only `s(61) = 8`, all 2026-09-30) and those are left to it.
- 2026-09-30, **claimed** status (Evan's call: a bound is not proved without a source we can review).  T. Green's 2000
  values (Friedman's DS7 Theorems 9–10 and Table 2, cited as "[8] T. Green, 2000, private communication"; the survey
  gives no proof) and Friedman's own n = 13 (3.8437) and n = 21 (4.7438) table values (no figure, no proof) are
  `"status": "claimed"`: kept in the history, tagged on the Bounds page, never `best`, never a step of the floor.  `best`
  falls back to Nagamochi 2005's general bound for n = 37, 38, 50, 51, 65–67, 82–85 (jlevy/squares treats them the
  same way).  Friedman's own results with figures and arguments in the survey (n = 7, 8, 14, 15, 19–20, 24, 35) stay
  proved.  Open: DS7 Figure 34 draws Green's unavoidable set for n = 17–18, which our point checkers could verify.
- 2026-09-30, DS7 Figure 34 checked (`s12/search/GREEN_FIG34.md`).  n = 17–18: the 16-point set that reproduces Green's
  `(40√2+19)/17` exactly (derived from the bound; matches the figure's pixels) is **not** unavoidable at that side: in
  each band the second diagonal is `√((1−e)² + V²) ≈ 1.011 > 1`, and a closed unit square at ~30.4° misses all 16
  points (exact check; re-checked independently in mpmath).  Green may have used a different set or argument, so the
  value stays `claimed`, now with this note.  The same configuration re-tuned verifies `s(17) ≥ 111/25 = 4.44` (not
  entered: far below today's floors).  n = 19–20: Friedman's 18-point set verifies exactly at `1121/250 = 4.484`
  (just below `6√2 − 4 ≈ 4.48528`; the run at 4.485 hit its depth limit); stays `proved`.
- 2026-10-01, **Nagamochi 2005 has a proof gap** (status `"gap"`).  Nagamochi's Theorem 2 (the `k² − 1` and `k² − 2`
  families and the general bound `s(N) ≥ 1 + √(N − 2⌊√N⌋ + 1)`) all come from his rectangle bound (Theorem 1), whose
  proof rests on Lemma 1 (every square scores more than 1).  Lemma 1 is false: chelokot/square-packing-archive gives a
  Lean-checked counterexample (side 1.0001 in `[0,4]²`, score 0.9775, 2026-09-04), and H. Karakuş, "A counterexample to
  Nagamochi's scoring lemma and a new rectangle packing bound", arXiv:2609.37410 (2026-09-29; abstract and PDF read)
  gives a family of counterexamples near a corner of every rectangle with `a > 3`, `b > 2`.  Karakuş says the
  published proof of the rectangle bound "is incomplete, but [the counterexamples] do not disprove the bound itself",
  and that his own argument does not establish `s(k² − 2) = k` or Nagamochi's general bound.  So every `Nagamochi2005`
  entry (except at perfect squares) is now `"status": "gap"`, with a note: kept in the history, never `best`, and
  excluded from the floor on the Bounds page, like `claimed`.  Replacements:
  - `Karakus2026` (`preprint`): Corollary 1.2, `s(k² − 1) = k` for `k ≥ 3` (strip-measure rectangle bound
    `ν(a,b) < ab − Δ(a)`), at `n = 8, 15, 24, …, 99`; and Corollary 6.2, `s(N) ≥ 1/2 + √(N − ⌊√N⌋ + 1/4)` for every
    nonsquare `N ≥ 8`, at every such `n`.  It is weaker than Nagamochi's stated bound.
  - `chelokot2026` (`preprint`, kernel-checked Lean, unrefereed): `s(k² − 2) = k` for every `k ≥ 2` by a compensation
    argument that keeps Nagamochi's resource measure but not Lemma 1
    (`docs/nagamochi-compensation-proof.md` of that repository, commit of 2026-09-04), at `n = 2, 7, 14, …, 98`.
  - Our own route to `s(k² − 2) = k` (monotonicity from `s(k² − 3) = k`, `k ≥ 4`; `/k2m3/#k2m2`) is not entered: as
    before, nothing is carried here by monotonicity.  For `k ≥ 9` it would rest on the k2m3 single-checker certificate.
  - `best` changed where it was a Nagamochi entry: to the earlier specific proof where one exists (n = 2, 3 Göbel;
    8, 15, 24, 35 Friedman), to chelokot at 23, 34, 47, 62, 79, 98, to Karakuş's Corollary 1.2 at 48, 63, 80, 99, and to
    Karakuş's general bound at 37, 38, 42–44, 50, 51, 54, 57–59, 65–67, 73–77, 82–96 (a drop of 0.004–0.05 each).
    `n = 7, 14` keep Friedman.  This is a judgment call (Evan's to revisit): a refereed theorem whose proof is known to
    have a gap is treated like an unreviewable claim, although nobody has shown the statement false.
  - None of this project's own results uses Nagamochi's Lemma 1 or his rectangle bound (repository-wide check,
    2026-10-01).  The wall-strip chord lemma (his Lemma 7(i), also Stromquist's) is proved from scratch in
    `s12/lean/Sqpack/Chord.lean`.
- 2026-10-01, n = 17.  Lower bound: `best` stays Guzhou0806's R068 `116511/25000` (replayed by jlevy/squares, T-043;
  still its verified floor in `packing/frontier/n-017.md` on 2026-10-01).  Added to the history as `reported` (the
  author's own checkers and CI only): Kleddamag `46601/10000 = 4.6601` (2026-09-29, commit `17d18245`), Guzhou0806 R070
  `46604427/10000000` (2026-09-29, `8988d933`) and R071 `18641771/4000000 = 4.66044275` (2026-09-30, `8c11f696`;
  its "R071 exact replay" GitHub Actions run succeeded; its README says the old C027 partition records are missing).
  If R071 is to be `best`, drop its `reported` flag.  Upper bound: the catalogue's Bidwell value 4.67553009360455 is
  unchanged; jlevy/squares T-065 (PR #265, merged 2026-10-01) certifies the packing exactly at the root of its contact
  chart, `s(17) ≤ 4.6755300936045509516…` (a rational ceiling; Kleddamag's rational witness `4675530093604551/10^15`
  also replays).  The site's upper bounds come from the catalogue, so nothing changes in the data.
