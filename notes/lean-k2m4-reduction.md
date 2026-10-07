# Lean: `s(k² − 4) = k` for all `k ≥ 8`, conditional on one finite statement (2026-10-03)

Task `lean-k2m4-reduction`.  Files: `lean/Sqpack/BentzFam.lean` (the reduction, generic in the edge-zone width `R`),
`lean/Sqpack/Bentz4.lean` (the `k² − 4` instance and the end theorem), `lean/Sqpack/Bentz4Data.lean` (generated
data, committed), `lean/scripts/gen_bentzfam_data.py` (the generator).  Predecessor: `notes/lean-bentz-reduction.md`
(`k² − 3`, `Bentz.lean`, unchanged).

## The theorem

```lean
theorem SquarePacking.Bentz4.bentz4_of_valid9 (h : Valid9) :
    ∀ k : ℕ, 8 ≤ k → minSide (k ^ 2 - 4) = k
```

`#print axioms` → `[propext, Classical.choice, Quot.sound]` (also for `box9Cover_measure`, `famCover_total4`, and the
generic `BentzFam.fileCover_measure`, `mass_shift`, `famCover_total`, `minSide_eq`; all listed in `lean/Axioms.lean`).
No `sorry`, no `native_decide`.  `minSide`, `Packs`, `sq`, `box` are the definitions used by `bentz_of_valid7`.

**The hypothesis** (`Bentz4.lean` §1):

```lean
noncomputable def box9Cover : MixedCover Empty (Fin boxSegs.length) Unit :=
  fileCover 2000000000000 boxSegs boxPolyVerts boxPolyW
def Valid9 : Prop :=
  ∀ (c : ℝ × ℝ) (θ : ℝ), sq c θ 1 ⊆ box 9 → 1 ≤ box9Cover.measure (sq c θ 1)
```

`fileCover` (`BentzFam.lean` §2) is a FORMAT.md v1 box file verbatim, exactly as `box7Cover` was: no points; the
segments in file order, `(X0/5, Y0/5)–(X1/5, Y1/5)` with mass `w/den` spread uniformly by length (`segMeasure`); the
one polygon, the intersection of the left half-planes of its directed edges (`convPoly`; here `[14/5, 31/5]²`,
proved), mass `pw/den` spread uniformly by area.  The data are `search/qx2_data/K4_k008_box9.txt` (2076 segments,
mass denominator `2·10¹²`, polygon mass `23120000000000/(2·10¹²) = 11.56 = (17/5)²`).  That is the certificate of
record's statement (run `qx2_k4x_k008`, 2026-10-02/03, VERIFIED-D4 by `qx2_zm`, 0 uncertified; axis and germ scans
exact): every closed unit square in `[0,9]²`, every centre and angle, has mass `≥ 1`.

Data (copied verbatim from `runs/qx2_k4x_k008/`, which is gitignored):

| file | from | sha256 |
|---|---|---|
| `search/qx2_data/K4_k008_box9.txt` | `sol_exact_box9.txt` | `4151d7c4059d5dcf56130c9373e6b5b60635e1a4250f46562d7ecc4d64a27801` |
| `search/qx2_data/K4_k008_family.txt` | `sol_exact_family.txt` | `a66668d9d414ad1cf78b19826ce592c540c9033a0234a06b207eb776b3e8ff0f` |

Family: `R = w = 3`, pitch `1/5`, Lebesgue from `a = 14/5`, `σ = 0`,
`D = 214770225571/200000000000 ≈ 1.073851 > 1`.  `python3 lean/scripts/gen_bentzfam_data.py` (from `s12/`; defaults
are these two files) writes `Bentz4Data.lean` deterministically and runs a Python mirror of all the kernel checks
first.

## How it is proved

The `k² − 3` proof, written once for any `R` (`BentzFam.Fam` = `R`, Lebesgue corner `A/5`, mass denominator, two code
tables).  `Bentz.lean` is untouched and still builds.

1. **Codes** (`cellT`, `lineT`; `N = 5R`, tables `(N + 7) × (N + 7)`).  Cell of `[0,k]`: offset `0..N−1` from the
   nearer wall, `N + phase` in the band zone, `N + 6` outside.  Line: offset `0..N` from the nearer wall — **offset
   `N` is the band's closed end `R` / `k − R`**, `N + 1 + phase` in the open band zone, `N + 6` outside.  Table
   entries (generator): corner module `H y x0 x1 → (5x0, 5y)` (`V` = diagonal image, checked equal); profile
   `h y p0 p1 → (N + 5p0, 5y)`, `v p y0 y1 → (5y0, N + 1 + 5p)` and, for `p = 0`, also `(5y0, N)`; corner pieces
   on `y = R` (`H 3 …`) go to the second table.  `famCover F k` = both layers of all unit segments of `[0,k]²` + the
   Lebesgue square `[A/5, k − A/5]²` with mass = area.
2. **`famCover fam4 9` = the box file** (`fileCover_measure`, `box9Cover_measure`).  Kernel checks (`Bentz4.lean`
   §3): every file entry is a unit grid segment with positive mass equal to the family's in its layer (`segOK`;
   layer 2 iff the layer-1 mass does not match); the entries' keys (`keyN (fileIx e)`, generator-written `boxKeys`)
   merge-sorted (`msort`, fuel-based, proved a permutation) equal `boxGridKeys`, which is strictly increasing and
   equals the keys of the non-zero (segment, layer)s of `[0,9]²` in lexicographic order (`gridKeys`).  Hence the
   entries are in bijection with the non-zero (segment, layer)s, with equal masses.  OutZero (the tables vanish on the
   outside row/column), the polygon and its mass are also checked.
3. **Localisation** (`mass_shift`, `Good`): for `k ≥ 2R + 2` and a closed unit square in `[0,k]²`, per axis an
   integer shift `n`: `0` if `x₁ < min(R + 3, k − R)`; `k − (2R + 3)` if `x₀ > max(R, k − R − 3)`; else
   `[x₀, x₁] ⊂ (R, k − R)` and `n = ⌈x₀⌉ − R − 1` puts it in `(R, R + 3)`.  Codes, layers and Lebesgue membership
   agree on every cell/line the square meets, so `μ_k(Q) = μ_{2R+3}(Q − n)` exactly.
4. **Accounting** (`sum_cellT`, `sum_lineT`, `seg_total`, `famCover_total`): for `k = m + 2R + 1`, cells: edge codes
   twice, phases `m + 1` times, outside once; lines: edge codes `0..N` twice, phase 0 `m` times, phases 1–4 `m + 1`
   times; the segment mass is `2 Σ_layers (H00 + m (H01 + H10) + m² H11)` with the `H`s finite table sums.  For
   `fam4` (`decide +kernel`): layer 1 `H = (42732189097180, 7630891120660, 3569108879340, 0)`, layer 2
   `(12406391400, 0, 0, 0)`; segments `22400000000000 m + 85489190977160` (units `1/(2·10¹²)`), so
   `total(μ_k) = k² − 4D` for every `k ≥ 7` (`famCover_total4`; `H11 = 0` is `σ = 0`).
5. **End** (`minSide_eq`): `Valid9` ⇒ `μ_k(Q) ≥ 1` for all closed unit `Q ⊆ [0,k]²`, `k ≥ 8`;
   `μ_k([0,k]²) ≤ k² − 4D < k² − 4`; `not_packs_of_measure` ⇒ no packing in side `< k`; `packs_grid` ⇒
   `Packs (k² − 4) k`.

## The threshold

`k₀ = 8 = 2R + 2`, as the brief expected; it comes from localisation only.  A square's axis range has width
`w(θ) < 2`; with no good shift it would have to meet both `x ≤ R` and `x ≥ k − R`, impossible iff `k − 2R ≥ 2`.
`k = 8 < 9` works (shift `−1` for the far wall).  The accounting is valid from `k = 2R + 1 = 7`; `k = 7` is not
covered by this localisation: a square can meet both edge zones of `[0,7]` (they are one unit apart), and no integer
shift makes `[0,7]`'s codes agree with `[0,9]`'s on it.  (`s(32) = 6` is `S32Lower.lean`; `s(45) = 7` is not claimed here.)

## What differed from `k² − 3`

* **Corner-module mass on the band's end lines.**  `K4_k008_family.txt` has `H 3 4/5 1`, `H 3 1 6/5` (and `V 3 …`),
  i.e. corner pieces on `y = R` (`x = R`); `L4_k02_family.txt` had none on `y = 2`.  These lines are also the
  profile's phase-0 lines, so (a) the end line needs its own code (in `Bentz.lean` it was phase 0; here the end line
  and an interior phase-0 line carry different masses), and (b) the box file lists those unit segments twice, once per
  piece (16 duplicates: 2 cells × 2 orientations × 4 corners).  Handled by the second layer.
* **Mass denominator `2·10¹²`** (was `10¹²`); now a parameter.
* **Size**: 2076 segments (800 for `k² − 3`), 209 + 2 non-zero table entries (70).  The pairwise `Nodup` / `contains`
  checks of `Bentz.lean` do not scale (kernel ~100 µs per comparison: 400 entries 32 s, the full list > 5 min and
  26 GB RSS when stopped), hence the sorted-key check (§2 above).
* Lebesgue corner `a = 14/5 = R − 1/5` (was `9/5`), polygon `[14/5, 31/5]²`.
* Nothing wrong was found in the family: the Python mirror and the kernel agree that `famCover 9` is the box file
  exactly, and the accounting reproduces the file's `81 − 4D`.

## Generic

`BentzFam.lean` takes any `R`, `A ≤ 5R`, mass denominator and two tables; keys need `5K < 4096`.  `k² − 5`
(`R = 5`, box 13) should need only the data file (`gen_bentzfam_data.py --box … --family … --ns …`) and an instance
file like `Bentz4.lean` (five `decide +kernel` checks, eight table sums, the end theorem).  The `k² − 3` files also
pass the generator's mirror with an empty layer 2, so `bentz_of_valid7` could be re-derived from `BentzFam` (not done;
`Bentz.lean` is untouched).

## Build

In the default build (`Sqpack.lean` imports `Sqpack.BentzFam`, `Sqpack.Bentz4Data`, `Sqpack.Bentz4`).  Alone, on
cores 12–15 with `LEAN_NUM_THREADS=4`: `Bentz4Data` 5.8 s, `BentzFam` 10 s, `Bentz4` 22 s (kernel checks: `segOK`
4.6 s, keys 2.9 s, `gridKeys` 11 s, sort 3.2 s).
