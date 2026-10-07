# arch-review: adversarial review of `notes/proof-architecture.md`  (2026-09-20)

Try to break the sketch.  (1) Does S0+S1+A2..A5 actually imply Claim(T)?  Find every case of `theta` no lemma covers
(the note admits one hole; find the others).  (2) Check each **[new, unchecked]** item: the dilation argument in S0
under closed semantics; the Mirsky/chain counting in A4 (acyclicity, "every pair x- or y-separated" for tilts `< eps`,
the `h_x = T-1` case, what `eps` must be); `tan(t/2) >= 1/T`.  (3) Run the `T`-filter on every lemma: true at `T = 17`?
at `n = T^2-4`, `T = 3`?  Which lemma is false at `T = 17`, and is the note's guess (β) defensible?  (4) Is A3 as
stated even plausible — read `search/BENTZ.md` §0, §7 and `notes/review-2026-09-13b.md` ("which hypotheses help") and
say whether a near-axis mass cap cuts the smeared tiling.  Pure reading and reasoning; small numeric checks only
(`<= 4` threads).  Deliverable: `notes/proof-architecture-review.md` — numbered objections, each "fatal / repairable /
cosmetic", with the repair where you see one.  Do not edit any existing file.  Do not commit.
