# t3-chord: prove `s(6) = 3` in the row language with chord rows — the one salvage attempt for the angle-space route  (2026-09-22)

**Why, and what is at stake.**  `notes/t3-existence.md` (2026-09-21) closed the go/no-go on the angle-space route as
formulated: the certificate side is complete and hypothesis-free, three strata are proved (axis-parallel, `>= 5`
axis-parallel, the far field `>= 25.8431°` by centre pigeonhole), and the surviving existence clause — "three squares
exactly axis-parallel and consecutively separated along one axis" — has an obstruction **X** (§4.4): it must separate
packings from configurations of margin `-0.36 t^2` / `-0.12 t^3` as `t -> 0`, and every mechanism tried (counting,
pigeonhole, unavoidable points, Menger/Dilworth) loses at *first* order in `t`.  Its §5.3 names the one thing untried:
**line-transversal (chord) rows**, the rows Kearney–Shiu's proof actually uses.  Evan's rule: we do not chase a route
for `s(12)` that cannot handle `s(6)`.  This is the single attempt to make it handle `s(6)`.  If it fails, the route is
closed and the note should say so in one sentence.

**Read first.**  `notes/t3-existence.md` §0, §1, §2 (the hypothesis-free identities: Lemma W′, the master formula,
Lemma Z″), §4 (the clause and X), §5.3 (chord rows), §8; `notes/t3-chain.md` §4 (Kearney–Shiu in this language — their
Lemmas 1–3, inequality (7), the two cases, and the dictionary §4.2); `notes/proof-anatomy.md` §5.2 (the K–S transcript);
`notes/chord-lemma.md` + `lean/Sqpack/Chord.lean` (the wall-strip chord lemma, already proved); `search/S6_SKELETON.md`
§3.1 (the disjunctive LP), §4.5 (chord sums); `notes/counting-ladder.md` §1 (chain counting, threshold order).
Instruments: `search/s6skel.py` (`chord_len(th, d)`, `Decider`, `Problem`), `search/t3_chain.py` (`all_rows`,
`dual_bound`), `search/t3_exist_*.py` (the 749 margin-0 packings: `runs/t3_exist_scope1.jsonl`, `t3_exist_row1.jsonl`,
`t3_exist_hunt1.jsonl`), `search/t3_chain_hform.py`, `search/ladder_hform.py`.

**The row.**  For a line `x = c` and squares `i, j` both meeting it, the chords `Q_i ∩ {x = c}`, `Q_j ∩ {x = c}` are
disjoint closed intervals, so `|y_i - y_j| >= (h_i + h_j)/2` with `h_i = chord_len(th_i, |x_i - c|)`; and `sum_i h_i <= T`
over all squares meeting the line.  **Be exact about linearity:** `h_i` is concave piecewise-linear in `|x_i - c|`
(`min(1/C, 1/S, (p - d)/(CS))`), so the row is linear in the centres only on a leaf that fixes a lower bound on `h_i`
— e.g. "`|x_i - c| <= d_0`" giving `h_i >= h(d_0)`, or "`Q_i` contains the point `(c, y_0)`" giving a chord through a
point.  That disjunction is where Kearney–Shiu's unavoidable points enter (a square covering a lattice point on `x = 1`
has a chord there of length `>= 2 sqrt2 - 2` by their Lemma 2 + (7)).  State every leaf hypothesis explicitly; a
"chord row" without its leaf hypothesis is not a row.

**Do.**  `T = 3`, `n = 6`, closed semantics (`notes/s13-casefree.md` §1).
1. **Reproduce Kearney–Shiu as a Farkas certificate.**  Take their two cases (`C` uncovered / `C` covered) and the
   incidence leaves they generate; write the rows each leaf satisfies (chord rows on `x = 1`, `x = 2`, `y = 2` with the
   constants `2 sqrt2 - 2`, `1/2`, `3/2`, `5/3`, `sqrt2 - 1/2`; wall rows; whatever separation rows they use), and
   exhibit the non-negative weights that give `delta <= 0` (their final `13/6 - sqrt2 < 2 sqrt2 - 2` is one such
   combination).  Verify each weight vector numerically (residual `< 1e-12`, `search/t3_chord_*.py`).  This is the
   calibration: if K–S cannot be written as certificates over leaves in this row system, say exactly which step resists.
2. **Replace the unavoidable set by pigeonhole where possible.**  Any tilt makes `sum_i width_i > 3T`... more precisely
   `sum_i (|cos th_i| + |sin th_i|) > n = 6` unless all six are axis-parallel, so some vertical or horizontal line meets
   `>= 3` squares with total chord `<= 3` — and with room to spare when tilts are large.  Quantify: as a function of the
   tilt vector, how many squares must some line `x = c` meet, and what is the best lower bound on `sum h_i(c)` over the
   choice of `c` (an integral over `c` of the chord function is the area, `= 1` per square — so the gain must come from
   *choosing* `c`, e.g. `c = 1, 2` as K–S do, or from the concavity of `h`).  Is there a `c` at which the chord row plus
   wall rows alone force `delta <= 0`?  Test on the 749 margin-0 packings and on the coherent-cone optima
   (`runs/t3_chain_scan1.jsonl`, families `coh*`, `unif*`), where X lives: does a chord certificate exist at
   `delta* = -0.36 t^2` points, and is its existence *provable* (i.e. does its leaf hypothesis follow from `delta >= 0`
   by counting, not by inspection)?
3. **Attack X directly.**  The clause to beat: a packing with `delta >= 0` and six tilts in `(0, 25.84°)`, one sign
   (the coherent cone), does not exist.  With chord rows, is there a certificate whose leaf hypotheses are implied by
   `delta >= 0` plus pigeonhole?  If yes, write the proof.  If no, give the configuration family (with `Decider` values)
   at which every chord certificate's leaf hypothesis fails while `delta*` is within `O(t^2)` of `0`.
4. **Assemble or refuse.**  Either a complete proof of `s(6) = 3` in this language (rows, leaves, weights, and the
   argument that the leaves are exhaustive), or a statement of the form "the proof needs [unavoidable set / lemma L]
   at step S, which is Kearney–Shiu's step and has no row-language substitute".
5. **Transfer to `T = 4`, in one section, honestly.**  Which parts of the chord argument are `T`-uniform (the rows,
   the chord lemma) and which are not (the accounting `n` squares vs `n+1` points).  Two concrete questions: (a) K–S's
   slack-1 accounting at `T = 4` needs a pure unavoidable set of **13** points in `[0,4]^2` (not 11: `t3-chain.md`
   §4.4 counted the bijection version); the fractional bound `COVER^closed(4) >= 12.2688` (`search/COVER4.md`) does not
   exclude 13, and the best known pure set is Bentz's 14 (`notes/s13-casefree.md` §1) — is 13 excluded by anything in
   the repo, and if not, what would it take to decide?  (b) Does the pigeonhole version of step 2 give anything at
   `T = 4` where the centre pigeonhole of `t3-existence.md` §3 does not fire?  Do not run `T = 4` computations beyond
   what answers (a)/(b) at the level of "excluded / open / decidable by X".

Label **[proved]** / **[measured]** / **[heuristic]** / **[guess]**.  Be adversarial with the repo and with this brief:
the last four briefs each had premises that were wrong (recorded in each note's last section).  Pen-and-paper first;
compute is for verifying weight vectors and testing candidate certificates on recorded packings.

**Deliverable.**  `notes/t3-chord.md` (verdict up front — one of: "`s(6) = 3` proved in row language", "K–S reproduced,
pigeonhole substitute fails at step S", "K–S does not fit the row system at step S"; then the proofs; the `T = 4`
transfer section), scripts `search/t3_chord_*.py` (import, never modify, existing ones), runs `runs/t3_chord_*`.
Do not edit `TODO.md`, `notes/status.md`, existing `notes/*.md`, existing `search/*`; do not commit.

**Working style.**  `<= 4` threads (a far-field job holds 8).  Anything over 10 minutes: `setsid nohup`, never `pkill`,
no `until … sleep` / `tail -f` watchers (a `pgrep -f <name>` loop matches itself).  Write the note section by section.
Stop when the verdict is clear; do not pad; a clean "no, because" is a full deliverable.
