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
