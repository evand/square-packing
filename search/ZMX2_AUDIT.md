# Adversarial audit of `zmx2` (the Rust mixed-cover checker), 2026-09-27

Brief: `tasks/s21-finish/audit-zmx2.md`.  Object: `verify2/src/bin/zmx2.rs` (sha256 `cead27e1…7084`, identical
to `HEAD` and to the source recorded in `search/zmx2_manifest.txt`), `search/ZMX2.md`, `search/zmx2_tools.py`,
`search/zmx2_tests.sh`, `search/zmx2_manifest.txt`, `tasks/line-cover/FORMAT.md`.  `ZM_MIXED.md` was not used; zmx2
is judged on its own proofs.  `zmx2.rs` was not modified; a scratch copy was built outside the repo (cores 10-11).
Audit tooling (own exact evaluator written from `FORMAT.md` only, nested-box probe, random adversarial covers,
differential `cert` test, adversarial covers): `search/zmx2_audit/`, driver `search/zmx2_audit/run_audit.sh`.

## 0. Verdict

**No defect that touches the `s(21)` certificate.**  Every floating-point operation on the certification path is
outward-rounded or only ever *removes* mass from the bound; every lemma checked out against its code; the D4 and
unreduced root sets tile the pose space; ≈ 3.5 M adversarial boxes and 117 sharp differential `cert` runs found no
unsound bound and no missed violation.  The scratch build is **bit-identical** to the manifest's binary (sha256
`f89104bf…f25e`), and re-running both official sweeps gave exactly the manifest's census (`--d4`: 1,826,222 boxes,
0 uncertified; `--full`: 14,709,448 boxes, 0 uncertified, `VERIFIED`).

**One must-fix, at the tool level:** integer arithmetic on input values is unchecked `i128` in a release build (no
`overflow-checks`), and a crafted file gets **`VERIFIED-D4` for a false statement** (F1: demonstrated).  The
candidate is far inside the safe range, so the claim is unaffected, but the parser must bound its inputs before the
checker is released as a general tool.

| # | grade | finding |
|---|---|---|
| F1 | **must-fix** (tool; not the claim) | Unchecked `i128` overflow on input integers: `s_num·D` wrap makes zmx2 check `[0,5]²` while printing `s ≈ 4.25·10³⁷` and `VERIFIED-D4`; weight wrap makes it print total `20.8947` for a file whose total is `≈ 3.4·10²⁷`. |
| S1 | should-fix | Provenance: the input `runs/line-cover_m5_candidate_x1003.txt` and the logs `runs/zmx2/*.log` are git-ignored (`runs/`); the manifest's sha256s point at files not in the repo. |
| S2 | should-fix | The shipped soundness harness (T3) goes through `zmx2 boxes`, which bypasses the `cert`-only code (float point pre-filter, inherited point weight and candidate lists, root enumeration, pass-1 reflection).  Covered now by A7 below; add it to `zmx2_tests.sh`. |
| S3 | should-fix | `--log` resume trusts `ROOT … uncert 0` lines keyed by a 64-bit FNV-1a hash of the input plus settings; the binary is not in the header and FNV is not collision-resistant.  The official runs were fresh ("0 already done"), so no effect here. |
| N1–N7 | nit | §2 below. |

## 1. Must-fix and should-fix, in detail

**F1. Unchecked integer arithmetic on input values.**  `Cargo.toml`'s `[profile.release]` has no
`overflow-checks`, so `i128` overflow wraps silently.  The parser range-checks only signs and `D < 2^22`.
Demonstrated with the release binary (`search/zmx2_audit/mkcovers.py`, `run_audit.sh` F1):

* `wrap_s.txt` = the candidate with header `s_num = 5 + 2^125`, `s_den = 1`.  `s_num·D ≡ 5000 (mod 2^128)`, so
  `sx = 5000` and every check runs on `[0,5]²`, but the output reads
  `cover: s = 42535295865117307932921825928971026437/1 … total 522368729933/25000000000` and
  `VERIFIED-D4: every closed unit square in [0,s]^2 has mu >= 1`, i.e. "`s(21) ≥ 4·10³⁷`".
* `wrap_w.txt` = the candidate plus three points at `(0,0)` with weights `2^127−1, 2^127−1, 2`.  The aggregate
  wraps to 0 (so the checked measure is the candidate's) and so does the total: `info` and `cert` print
  `total = 522368729933/25000000000 = 20.894749197`, and `cert --d4` says `VERIFIED-D4`, for a file whose total
  mass is `≈ 3.4·10²⁷`.

Other unguarded places of the same kind: coordinates cast `i128 → i64` after the range check (truncates if
`sx ≥ 2^63`); the point-bucket allocation `nb² ∝ (10 s)²`; `F`, `target` and the bound sums
(`∝ W·Lc·2^30·mass`).  `ZMX2.md` §5(iv) ("`F < 2^76`, sums `< 2^80`") is true **for this input** (`W = 10^11`,
`Lc = 20`, per-line mass `< 5`: `|F|·2^30 < 2^74`), not asserted for general input.  (Lemma P is safe for every
accepted file: its offsets `A` are taken only for points within the `0.75` reach filters, so `|A| ≲ D·10·2^27 <
2^53` and all its products stay `< 2^112`; §2's "`|A|, E < 2^45`" is the `D = 1000` instance.)  Patch in §5.

**S1. Provenance.**  `s12/.gitignore` ignores `runs/`; `git ls-files` has neither the candidate nor the logs.  The
bundle must ship `runs/line-cover_m5_candidate_x1003.txt` (sha256 `8b415cee…fc23`) at a tracked path, and either
the two logs (sha256 `174d7d04…fdc4`, `2bc0562f…6af6`, verified) or a note that they are reproducible (both
reproduced here with identical censuses; the build is bit-reproducible with rustc 1.91.1).  At audit time the release
bundle in progress (uncommitted) already has `certificates/s21/s21_mixed_cover_5.txt` with the same sha256 and
`certificates/s21/zmx2_{d4,full}/` logs, which addresses this once committed.

**S2. Harness coverage of the `cert` path.**  `zmx2 boxes` re-tests every candidate point with the exact
`pt_certain_in`.  `cert` additionally (a) drops points by the float `pt_float_class`, (b) inherits the certain-in
weight and the undecided list to children, (c) builds roots and the reflected cover.  All three are sound by
inspection (§3), and are now exercised by A7 (differential `cert` with a known exact violation of `10^-5`), which
should join `zmx2_tests.sh`.

**S3. Resume.**  Suggest putting the sha256 of the input and of the binary in the log header, or having
`zmx2_run.sh` refuse to write a manifest when any root was resumed.

## 2. Nits

* **N1** Comment labels: `check_d4` cites "(sec 7, Lemma S)"; `g_range` cites "Lemma G" (it is part of Lemma M);
  `pair_bound` has its doc line twice; `reflect_y` says "theta in [45,90]" (pass 1 covers `[36.87°, 90°]`); a stray
  "lower bound of the segment mass …" comment sits above `pt_conds`.
* **N2** `Iv::rat` asserts `|n|, d < 2^53`, but `D < 2^22` with `cl ≤ 27` allows `2·den ≈ 2^54.3`: a reachable
  panic for `D > 2^21.4`.  Safe (a worker panic makes `h.join().unwrap()` panic in `main` before any verdict is
  printed; checked), but the documented limit is inconsistent: use `D < 2^20` or check `cl` against `D`.
* **N3** `ZMX2.md` §1 says `j, k ≤ 28`; the code stops at `cl, ul ≥ 28` (so `≤ 27`).
* **N4** `lcm(segment lengths) ≤ 2^24` rejected 114 of 800 random covers in A6 (completeness only; not an issue
  for the candidate, `Lc = 20`).
* **N5** `pt_float_class` uses an absolute slack `1e-9`: safe for coordinates up to `~10^5` (float error in `G`
  `~ 10^-16·s`), and it only ever drops points that are outside `Q` at every pose, or skips an exact test.
* **N6** `+5` is accepted as an integer (Rust's `parse`).  Harmless.
* **N7** Lemma H's proof uses `z = zlo = −∞`; in the code that is the clamp `−big`, which is correct by §5(iii)
  but worth saying in the proof.

## 3. What was checked, and how

### 3.1 Interval arithmetic and rounding (Lemma R)

* `Iv::add/sub/mul` widen each endpoint by one ulp after round-to-nearest (on IEEE-754, round-to-nearest error is
  `< 1 ulp`, also at overflow to `±∞`); `mul` takes the min/max of the four products; exact zero factors give exact
  0 (sound: the product is exactly 0).  `recip_pos` asserts `lo > 0` (no division by an interval containing 0;
  all call sites have `u > 0` or `1 − u² ≥ 3/4` or `1 + √disc ≥ 1`).  `sqrt` is correctly rounded, widened
  downward and clamped at 0.  NaN is asserted absent at every add and grid conversion.
* Conversions: `Iv::rat` converts `n`, `d < 2^53` exactly (asserted) and widens the quotient by 4 ulps.  For the
  candidate `den ≤ 10·2^27·1000 < 2^44`.  The grid values `glo/ghi` are `⌊next_down(fl(v·G))⌋`,
  `⌈next_up(fl(v·G))⌉` with `G = D·2^30` exact: correct directed bounds.  `f64 → i128` casts only happen after the
  clamp `|v| < 4s + 10`.
* Clamps (§5(iii)): the two "wrong-direction" clamps (`glo(v ≤ −big) = −big`, `ghi(v ≥ big) = big`) were checked
  at each use (`zlo`, `zhi`, `h1lo`, `l2hi`, and the pair function): `F` is constant outside `[0, s]`, atoms lie in
  `[0, s]`, shifts are `≤ u0 ≤ ½`, so no value changes.
* `u0g = ⌊u0·G⌋` (exact integer division) is `≤ u0·G ≤ u·G`: the direction the pair lemma needs.
* After the grid step everything is `i128` and exact (overflow margins: F1).
* **No float certifies.**  All float uses outside Lemma R: `initial_candidates` (bucket/reach filter: can only omit
  points), `pt_float_class` (drops points only when every `G_k` is `> 1e-9` over the box, i.e. outside `Q` at every
  pose, which is inherited correctly; class 1 only triggers the exact test), the `0.75` line/atom reach filters (omit
  only), `split` (choice of axis), and the diagnostics (`float_min_box`, `fmass`, `fscan`, `--tight`).

### 3.2 Lemmas against code

* **Lemma P**: re-derived `G0…G3` from `X = aC + bS`, `Y = −aS + bC`; the vertex test `4p2p0 − p1² < 0` with the
  interior condition in scaled integers is right; the same code serves atoms (per condition) and the rotated frame.
* **Lemma C/M**: `f1, f2` coefficients (`α = (d∓½)/2`, `β = ∓(d±½)/2`) match `line_ends`; affine in `d` for fixed `u`
  justifies evaluating at `d0, d1` only; `h_range`: interior critical point only when `αβ > 0` with certain signs,
  otherwise monotone or naive over `[u0, u1]`; the `u0 = 0` limits `±∞`/`0` are right; `g_range`: derivative
  `∝ −2d(1+u²) ± 2u`, critical point `2|d|/(1+√(1−4d²))`, none for `|d| > ½`, monotone cases match the sign
  conditions (`Iv::rat` of a non-zero numerator never contains 0, so the "straddles 0" branch is only the naive
  superset).  `T, R` increasing on `[0,1)` for the interval evaluation.
* **Lemma Z**: `f2(d−1,u) − f1(d,u) = u` verified symbolically.  `inf Φ`: `Φ` is continuous piecewise-linear in the
  *real* `z` with kinks only at `l_a, h_a, h1b−u0g, h2b−u0g, l_b−u0g`, the breakpoints of `a` and of `b` shifted,
  and atom positions, all grid integers, all in the candidate list; atoms: `A` left-continuous, `B`
  right-continuous, both one-sided combinations taken at every candidate and at the ends.
* **Lemma DP**: chains by `pos mod D`, pairs only at difference exactly `D`; atoms belong to one line; densities of
  crossing lines meet in a null set; ordinary points and atoms are disjoint.  Any partition is a lower bound.
* **Lemma W**: re-derived `q = (1 − u² + 2u³)/(2(1+u²))` from `f2((C+S)/2 − 1)`; the lower bound on `[u0,u1]` is
  valid since the numerator stays positive for `u ≤ ½`; the atom version uses the same `h1w`/`l1w`; at `u = 0`
  the modified conditions are exactly admissibility (`X ≤ ½ ⟺ c_x ≥ ½` for the line `x = 1`).
* **Lemma H** (`θ = 0` inside a `u0 = 0` box): checked the case analysis, including the atoms (conditional `a`-atoms
  at `z = −big` are in `Q` iff `d ≤ ½`), and the wall interaction: at an admissible `θ = 0` pose a pair `(0, 1)`
  with `b` the Lemma-W line has `d_a = c_x ≥ ½`, so the `z = −∞` branch never meets the modified `h1b`; the
  right-wall line `a = s − 1` has `d_a ≤ ½`, so the `z = +∞` branch never meets the modified `zhi`.
* **Lemma T**: `ρ(x,y) = (−y, x)` commutes with rotations; box, walls (`−s`, `0`), line positions (`−y`) and the
  atom frame (`pb`) match.
* **Lemma E**: `w' ∝ 2 − 4u − 2u²`, so the minimum of `w` on a `u`-interval in `[0,1]` is at an end; the
  enclosures are directed correctly.
* **Regions.**  `--d4`: generators `x ↦ s−x`, `x ↔ y` checked exactly on aggregated points and canonical line
  densities; reduction: the two axis reflections bring `c` into `[0, s/2]²` (each maps `θ ↦ 90° − θ`), the
  diagonal then fixes `θ ≤ 45°` (`u ≤ √2 − 1 < ½`) without leaving the quadrant; roots `25 × 25 × 4` closed
  boxes cover `[0, 2.5]² × [0, ½]`.  `--full`: pass 1 uses `y ↦ s − y`, `Q(c,θ) ↦ Q(ρc, 90° − θ)`, so
  `u ∈ [0, ½]` in both passes covers `θ ∈ [0°, 53.13°] ∪ [36.87°, 90°]`; roots `50 × 50 × 4 × 2`.  Verdict
  logic: missing roots → `INCOMPLETE`, capped roots count their stack as uncertified, restricted `--xlo/…/--bins`
  → no verdict, worker panic → no verdict.

### 3.3 Parser

`ERROR` (exit 2) confirmed for: version ≠ 1, `s_den ∤ s_num·D`, non-integer token, empty file, zero header
values, `W = 0`, `np` larger than the file, `D ≥ 2^22`, trailing tokens (including a plain file followed by
`ns npg`), negative coordinate, integer literal beyond `i128`, `mixed` alone, and the T1 set (negative weight,
outside point, degenerate/diagonal segment, polygon, short file).  Not refused: the wraps of F1.

## 4. New adversarial tests (all on cores 10-11; `run_audit.sh`)

| # | test | result |
|---|---|---|
| A1 | candidate × 0.9943 (`--d4`): exact `μ = 0.9999956` at `(0.573593, 1.440624, u = 0.0805601)` | refused; deepest uncertified box float-min `0.999991` at that pose (also refuses a `1.002` region at `38.6°`: completeness) |
| A2 | germ hole `g1` (asymmetric, `x = 3` on `[2.88, 2.92]`, `x = 2` on `[2.92, 2.96]`), `--full`, `c ∈ [2.4,2.6]×[2.3,2.7]`, bin 0 | refused; uncertified box `x,y ∈ [2.5, 2.5000122]`, `u ∈ [0, 7.6·10⁻⁶]`, float-min `0.99197`; exact `μ = 0.991819` at `(2.50000084, 2.5, u = 10⁻⁶)`, `1.698` at `θ = 0` |
| A3 | D4-symmetric germ hole `g2a` (one piece and its 7 images; germ `μ = 0.99405` at `θ → 0⁺`, `θ = 0` min `1.0037`), `--d4` | refused (4,619 boxes) |
| A4 | wall-germ hole `w1` (`x = 1` zeroed on `[2.94, 3.0]`), `--full`, `c ∈ [0.5,0.6]×[2.4,2.6]` | refused; float-min `0.998199` at `(0.500001, 2.500001, u = 2·10⁻⁷)` |
| A5 | nested-box probe at the violating poses of A1–A4 (boxes of every level `0…26` containing the pose, all placements when it is on a grid line / bin boundary, incl. `u = 0`, `c = 2.5` = D4 edge) | 0 fail; `μ − bound` down to `1.1·10⁻¹⁰` (the bounds are essentially exact in small boxes) |
| A6 | 800 random adversarial covers (`s ∈ {2, 12/5, 5/2, 3}`, `D ∈ {10, 20, 56, 100, 1000, 2^18}`, segments on grid, boundary, off-grid lines and partners at distance 1, overlapping; points on lines, crossings, free; 30 % reflected frame; `--pair-points` / `--no-atoms` mixed in), 60 poses each: germs (`c_x = ℓ ± ½ ± {0, 10⁻⁹…10⁻⁴}`), admissibility boundary `c = w/2`, points/segment ends/line crossings exactly on edges and corners of `Q`, `u ∈ {0, 10⁻⁹…10⁻³, 2⁻³⁰, 1/8…1/2}`; nested boxes to level 22 or 27 | 686 covers run (114 refused: `lcm > 2^24`), 41,160 poses, **3,503,598 boxes, 0 fail** (bound `≤` exact `μ` everywhere; `μ − bound = 0` attained) |
| A7 | differential `cert --full`: random line+point covers, float minimum located and refined, weights scaled so that exact `μ(P) = 1 − 10⁻⁵` at a rational pose `P` | 164 covers run (26 skipped: no usable minimum): **all refused**.  In the 117 sharp cases (the same cover × `1.003/μ_min` is `VERIFIED`, so `μ ≥ 0.997` everywhere) 15 minima are germs (`θ ≈ 0.003°`), 48 at `θ ≥ 45°` (pass 1); in 25 of them the location was checked: an uncertified box contains `P` every time |
| A8 | candidate `--full --tight 0.001`: the certified leaves with bound `< 1.001`, own exact evaluator at corners + random rational poses | 1,500 of the 797,100 dumped leaves (both passes), 28,322 admissible poses: **0 fail**; `μ − bound ≥ 0.0063` (median 0.0126) |
| — | the real covers (candidate, reflected candidate, `x9943`, `g1`, `g2a`, `s(13)` plain/reflected/`--no-atoms`, `s(32)` with `--pair-points`) at A6-style poses | 0 fail (≈ 64,000 boxes) |
| — | reruns of the official sweeps with the scratch build | identical censuses; `VERIFIED-D4` and `VERIFIED` |


## 5. Proposed patch for F1

Bound the inputs in `parse_cover` so that every later size claim of `ZMX2.md` §2/§5 holds for *all* accepted files
(`F·2^30 ≤ Σw·Lc·2^30 < 2^56+24+30 = 2^110`, point and atom weights likewise; `sx < 2^27` keeps the casts,
the bucket grid and the grid integers `≤ (4s+10)·D·2^30` small):

```rust
// after reading s_num, s_den, d, w:
const HDR: I = 1 << 40;
if s_num >= HDR || s_den >= HDR || w >= (1 << 50) {
    die("header value too large (s_num, s_den < 2^40, W < 2^50)");
}
if (s_num * d) % s_den != 0 { die("s_den must divide s_num*D"); }   // no overflow now: < 2^62
let sx = s_num * d / s_den;
if sx >= (1 << 27) { die("s*D too large for this checker (limit 2^27)"); }
// for every point and segment weight:
if wt >= (1 << 50) { die("weight too large (limit 2^50)"); }
// after the loops:
let tsum: I = raw_pts.iter().map(|p| p.2).chain(raw_segs.iter().map(|s| s.4)).sum();  // < 2^50 * #pieces
if tsum >= (1 << 56) { die("total weight too large (limit 2^56)"); }
```

(`#pieces` is bounded by the file; `2^50 · 2^20` pieces would still fit `i128`, the explicit `2^56` check keeps
`F·2^30·Lc < 2^110`.)  Also lower the `D` limit to `2^20` (N2).  An alternative is `overflow-checks = true` in
`[profile.release]`, but that edits `Cargo.toml` (shared with `zmcheck`, guarded by T0) and turns overflows into
panics rather than refusals.  Add `wrap_s.txt`/`wrap_w.txt` to T1.
