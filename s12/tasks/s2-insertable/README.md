# S2: insertable box covers at the threshold (c = 3, k = 4; c = 4, k = 5)

Background: `search/FRIEDMAN.md` §2.1 (insertion lemma), §9.1, §9.4 (S2), §9.6 item 1.  LP machinery:
`search/line_cover.py` (the k = 5 LP cover `lc_m5germC_r12`, `search/LINE_COVER.md`, `S60_COVER.md` §4).

## Question
Does a valid closed cover measure μ of `[0,k]²` exist that is (i) **insertable** and (ii) saves more than c?
(i) means: for the centred slab `X = [t, t+1+√2] × [0,k]`, `t = (k−1−√2)/2`, μ restricted to X is 1-periodic in x
(`μ(A + e₁) = μ(A)` for Borel `A ⊂ [t, t+√2] × [0,k]`), and the same for the horizontal slab (D4 image).
If yes, the insertion lemma gives `s(k'² − c) = k'` for **every** k' ≥ k from that one box.

## Mass bookkeeping (derive and confirm; don't trust blindly)
`p_x = μ([t,t+1) × [0,k])`, `p_y` likewise horizontally, `q = μ([t,t+1)²)` (the crossing cell).  x-insertion adds
p_x; then y-insertion adds `p_y + q`.  Saving is preserved at the step iff `p_x + p_y + q ≤ 2k+1`, and at all
later steps iff additionally `q ≤ 1` (each step raises p_x, p_y by q).  Sanity: Lebesgue has `k + k + 1 = 2k+1`.
Check the bookkeeping on Lebesgue and on the k²−3 family box `search/qx2_data/L4_k02_box7.txt` (insertable by construction).

## Steps
1. **Baseline** at k = 4 and k = 5: the unrestricted box LP (D4, closed semantics) on the exact atom/row grid the
   S2 run will use, so the comparison is fair.  Existing values: k = 5 LP 20.7489; k = 4 cover LP bracketed
   `12.16 ≤ COVER^closed(4)`, best exact cover 12.956 (`certificates/rung2/`).
2. **Insertable LP**: same LP plus linear periodicity equalities on the atoms in each slab (split segments at the
   slab boundaries `t, t+√2, t+1, t+1+√2` so restrictions are well defined; the candidate atom set must be closed under
   ±e₁ within the slab) plus the mass constraints above.  Report total, saving, slack vs baseline.
3. **Two measure classes**, both for k = 4 and 5: (a) points + segments only (what the Lean `CovM` / ZMTree checks today,
   e.g. s(13) in 33 CPU-min); (b) additionally a central Lebesgue square `[a, k−a]²` with `a ≤ t` (periodic in both
   slabs automatically).  The (a)-vs-(b) gap is a data point for whether Lean polygon mass pays for itself.
4. If an insertable LP saves > c with real room (≳ 0.1 at k = 5; at k = 4 compare to the 0.044 margin of the rung-2
   cover), stop and report: the exact certificate (zm_mixed / zmx2 / zeromargin + an exact periodicity/mass check of
   the file) is the next task, not this one.  If no, measure the price `L(k) − S_ins(k)`; k = 6 for c = 4 is optional.

## Rules
- Compute: pin to physical cores 0–7 (`taskset -c 0-7`, ≤ 8 processes; logical i and i+16 are siblings), CPU budget
  ~10⁴ CPU-s total; ask before exceeding 3× that.  Nothing else is running.
- Write-up: `search/S2_INSERTABLE.md` (labels [proved]/[measured]/[heuristic]); runs in `runs/s2_*` (ignored).
  Commit code + write-up to the public repo `main` (not pushed).  Floats are fine; this is a go/no-go.
