# Corner-local covers: the corner deficit is exactly 0 (task corner-deficit, 2026-09-27)

Question: can one corner module (Lebesgue measure everywhere, modified only inside `C = [0,R]²`) save `D > 1`
(or `D > 3/4`) per corner, giving `s(k²−4) = k` (or `s(k²−3) = k`) for all large `k` from one certificate?
Formally `D*(R) = sup R² − μ(C)` over measures `μ ≥ 0` on `C` such that every closed unit square `S ⊂ Q = [0,∞)²`
captures `μ(S ∩ C) + area(S \ C) ≥ 1`, i.e. **`μ(S ∩ C) ≥ area(S ∩ C)`**.

Labels: **[proved]** (short exact argument, given here), **[measured]** (exact arithmetic on a finite object,
or a float LP/scan whose numbers are reported as-is), **[heuristic]** (interpretation).

Code: `search/corner_deficit.py` (modes `dual`, `lp`, `profile`).  Outputs: `runs/corner_dual.txt`,
`runs/corner_lp.txt`, `runs/corner_profile.txt`.  Wall time: minutes (the question closed analytically).

## 0. Answer

> **`D*(R) = 0` for every `R` [proved].  No corner module saves anything, so neither `D > 1` nor `D > 3/4` holds.**
> The same argument kills the 4-corner scheme in a `k`-box whenever `2⌈R⌉ < k` (§2), i.e. whenever the corners
> are not already the whole box.  The savings `≈ 4.1–4.3` of the `s(21)`, `s(32)`, `s(45)` covers are not
> corner-local.  A module can only save when the modified region cannot be almost tiled by pairwise-disjoint
> closed unit squares in the container, which rules out every region that stays clear of the walls on one side.

| `R` | LP `D` (rows: axis squares on the `1/q` lattice), `q = 10 / 20` | = `R² − (R − (⌈R⌉−1)/q)²` | off-lattice scan: min `μ/area` (q = 10 / 20) | honest `D` |
|---|---|---|---|---|
| 1.5 | 0.290 / 0.148 | yes / yes | 0.00 / 0.57 | **0** |
| 2   | 0.390 / 0.198 | yes / yes | 0.70 / 0.57 | **0** |
| 2.5 | 0.960 / 0.490 | yes / yes | 0.59 / 0.00 | **0** |
| 3   | 1.160 / 0.590 | yes / yes | 0.25 / 0.00 | **0** |
| 3.5 | 2.010 / 1.028 | yes / yes | 0.00 / 0.00 | **0** |
| 4   | 2.310 / 1.178 | yes / yes | 0.00 / 0.00 | **0** |

* **LP [measured].**  The LP value is exactly the dilated-grid dual of §1 with `ε = 1/q`.  It is a pure
  row-sampling artifact: it halves when `q` doubles.  The LP solutions are invalid off the lattice.  The axis-parallel
  scan at pitch `1/(4q)` finds squares capturing `0–70 %` of their area inside `C`, usually thin slivers at the
  far edges `x, y ≈ R` of `C` that the LP leaves empty.
* **Repair.**  The additive repair `μ + (1 − min ratio)·Leb(C)` gives `D = −0.8` to `−14.8` on the scanned poses.
  Scaling the modification cannot repair a `0` ratio.  Any mixture `(1−t)μ + t·Leb` is valid only at `t = 1`.  The
  honest value is therefore `0`, attained by Lebesgue itself, and §1 shows nothing beats it.
* **Growth with `R`:** none.  `D*(R) ≡ 0` (no plateau above 0, no sublinear growth).
* **Structure of the optimum:** Lebesgue.  The LP's apparent savings sit exactly in the `ε`-gaps of the dilated
  grid (width `1/q`, at `x, y ≈ 1, 2, …`), which is below the LP's resolution.

## 1. The proof: `D*(R) = 0`

`D*(R) ≥ 0`: `μ = Leb|_C` captures exactly 1 everywhere.

`D*(R) ≤ 0`: let `n = ⌈R⌉` and `ε > 0`, and take the `n²` squares `S_ij = [i(1+ε), i(1+ε)+1] × [j(1+ε), j(1+ε)+1]`,
`0 ≤ i, j < n`.
* They are closed unit squares in `Q`, pairwise disjoint (gaps `ε`).
* Each meets `C`, because `i(1+ε) ≤ n−1 + O(ε) < R` for small `ε`.
* Together they cover `C` except strips of total area `≤ 2(n−1)εR`, since the last column reaches
  `(n−1)(1+ε) + 1 ≥ R`.

Any feasible `μ` has `μ(S_ij ∩ C) ≥ area(S_ij ∩ C)`.  Summing over the disjoint family gives
`μ(C) ≥ area(C ∩ ⋃S_ij) ≥ R² − 2(n−1)εR`.  Let `ε → 0`.  ∎

`corner_deficit.py dual R --eps …` checks the family in exact arithmetic (`runs/corner_dual.txt`: e.g. `R = 4`,
16 squares, uncovered `2.4·10⁻⁵` at `ε = 10⁻⁶`).

The whole point is that squares crossing `x = R` or `y = R` may overlap each other freely outside `C` (there the
measure is Lebesgue and each square gets its full share).  So the dual is a packing that only has to be disjoint
inside `C`, and the dilated grid "exhales" its `O(ε)` excess through the open sides.  The walls `x = 0, y = 0`
push nothing back.  The same holds for any bounded modification region in `Q`, half-plane or wall strip that is
not bounded by walls on opposite sides.

## 2. Consequences for box covers (`[0,k]²`)

For a valid closed cover `μ` of `[0,k]²`, i.e. every closed unit square in the box has `μ(S) ≥ 1`:

1. **[proved] Localisation lemma.**  Suppose `μ = Leb` outside a set `M`, and some pairwise-disjoint closed unit
   squares in the box cover `M` up to area `η`.  Then `k² − μ([0,k]²) ≤ η`, because
   `μ(S ∩ M) ≥ area(S ∩ M)` for each `S`.
2. **[proved] Four corner modules** of size `R` (Lebesgue elsewhere) save nothing when `2⌈R⌉ < k`.  The four
   dilated corner families are then mutually disjoint and inside the box.  At `R = k/2` the four corners are the
   whole box, so the "corner" problem is just the D4 box problem again.
   * Under the natural non-interaction condition `k ≥ 2R + √2`, no unit square meets two corner regions.  A square
     in the box that meets `C_1` is then in `Q_1`, and a square in `Q_1` that meets `C_1` lies in the box, since it
     reaches at most `R + √2 ≤ k − R`.  So the box constraints are exactly the four single-corner constraints, and
     §1 applies corner by corner.
   * Either condition (non-interaction, or `2⌈R⌉ < k`) gives zero total saving.  The brief's intended route is
     closed, not just untested.
3. **[proved] Integer corners of any valid box cover carry no deficit.**  `μ([0,n]²) ≥ n²` for integer
   `n ≤ k−1`.  Same family: the slivers `(n, n+nε]` have mass `→ 0` as `ε → 0`.  Likewise
   `μ([w, k−w]²) ≥ (k−2w)²` for integer `w ≥ 1`.  So the whole saving `k² − μ` is at most the deficit of the
   width-1 L-strip `[0,k]² \ [0,k−1]²`.

## 3. Sanity check against the certified covers [measured, exact `Fraction`]

This is `corner_deficit.py profile` on `certificates/s21/s21_mixed_cover_5.txt`, `s32/s32_closed_cover_6.txt` and
`s45/s45_mixed_cover_7.txt` (`runs/corner_profile.txt`).  Each entry is the corner deficit `R² − μ([0,R]²)` of the
closed corner square, identical at the four corners (D4).

| `R` | k = 5 (savings 4.105) | k = 6 (4.286) | k = 7 (4.226) |
|---|---|---|---|
| 1   | −0.007 | −0.011 | −0.023 |
| 1.5 | **+0.803** | **+0.790** | **+0.812** |
| 2   | −0.028 | −0.046 | −0.090 |
| 2.5 | **+1.013** | **+1.069** | **+1.086** |
| 3   | −0.068 | −0.112 | −0.203 |
| 3.5 | +0.798 | +0.898 | +1.020 |
| 4   | −0.228 | −0.232 | −0.363 |

* Integer corners are negative, as §2.3 requires.  Centres `[1, k−1]²` and `[2, k−2]²` are negative too
  (`−2.56 / −0.71`, `−3.50 / −1.44`, `−4.53 / −2.14`).
* Half-integer corners show deficits `≈ 0.8–1.09`.  This is presumably where the "≈ 1 per corner" intuition came
  from.  **[heuristic]** It is a boundary artefact: the covers concentrate mass on the grid lines, and a cut at
  `x = n + ½` excludes the line `x = n+1` while counting half a cell's area.  These deficits are paid for by the
  mass just outside the cut.  They cannot be kept once the outside is replaced by Lebesgue measure (§1).
* The bottom strip `[0,k] × [0,1]` has deficit only `0.47 / 0.40 / 0.34`, falling with `k`.  The saving is spread
  around the whole boundary ring and is a global (rigidity) effect: `k` unit squares cannot sit disjointly in a row
  of length `k`.  It is not a sum of local corner contributions.

## 4. Caveats and what would still make sense

* The LP (§0 table) uses axis-parallel rows only.  That is enough here, because the exact argument uses only
  axis-parallel squares, so tilted rows could only lower `D` further.  Nothing heavier was run: the answer is
  analytic.  (A `q = 40` run was stopped after 10 min.  It is unnecessary, since the LP value is known in closed
  form.)
* Any certificate that works "for all `k`" must modify a region that the container's walls enclose on all sides
  (the localisation lemma), i.e. a family of box covers whose deficit lives on the whole wall ring.  By the
  Erdős–Graham-type bounds in the brief, a fixed per-unit-length wall profile has deficit `o(k)`.  Whether some
  wall-ring construction keeps the deficit above 4 uniformly in `k` is exactly as hard as the original problem.
  §1 does not decide it.
