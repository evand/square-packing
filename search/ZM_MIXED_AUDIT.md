# Adversarial audit of `zm_mixed.py` and the `s(21) = 5` run (2026-09-27)

Brief: `tasks/s21-finish/audit.md`.  Auditor's code: `search/zm_mixed_audit.py` (new; independent file checks, an
independent exact mass oracle, component tests, hole finder).  All runs on cores 10–11; logs and holed covers in
`runs/zmm/audit/` (gitignored).  Audited: `zm_mixed.py` sha256 `ebf8bbc3…2a905` (= commit `d4b97f9` = HEAD), with
`zeromargin.py` `640fe453…086ab` (pinned) and `mixed_cover.py` `bb89de15…aae5` (= `b37ecf2`).

## Verdict

**No soundness defect found.**  Every lemma of `ZM_MIXED.md` §2 was checked line by line against its proof and its
code; I found no gap, and nothing a float can certify.  The run's 40,000 roots tile the D4 region exactly, the checked
file is the one whose total is claimed, and D4 invariance of the measure holds exactly.  New targeted tests (≈ 115 M
exact checks of every certified *component* against exact masses, plus rejection runs on three deliberately holed
versions of the real `m = 5` cover) found **0 violations**, and the checker was tight enough that several bounds hit
the exact mass with slack `0` (so the tests are sharp, not slack).

The claim `s(21) = 5` is **supported**, still conditional on the two items `ZM_MIXED.md` §7.5 already names: an
independent re-check (the xcheck agent) and the measure version of the FORMAT.md reduction in Lean (the lean agent;
the paper argument is correct, §1 below).  Before shipping, the should-fix items S1–S3 (provenance) matter; none of the
findings is must-fix.

| grade | # | finding |
|---|---|---|
| must-fix | — | none |
| should-fix | S1 | the certified cover, run log, jsonl and sha file live in gitignored `runs/`: nothing that was certified is committed |
| should-fix | S2 | the run log does not record `zm_mixed.py` / `mixed_cover.py` / input sha256 or all flags (SPLIT, theta-bias, chain-from); the sha file was written by hand |
| should-fix | S3 | `--resume` keys finished roots by `(label, root)` only: it would silently mix runs of different code / files / flags |
| should-fix | S4 | the "leaf stress tests, 0 exact violations" evidence (`ZM_MIXED.md` §4.4–4.5, §7.1, §7.5) cannot detect checker bugs on a cover with margin; replace/augment with component tests (§3 here) |
| should-fix | S5 | trusted surface larger than needed: Corollary T′ (points on germ lines) and the polygon code are never executed by the `s(21)` run (0 points on lines, 0 polygons) |
| nit | N1–N6 | dead code, docstring/proof wording, doc order (§5 below) |

## 1. Proofs (line by line)

Checked, with the closed semantics, `θ = 0` and `θ > 0`, clipping, endpoints, disjointness and inheritance:

* **Lemma A/B/C use, Fact 0/1, Lemma S.**  `cond_poly` builds exactly zeromargin's `_cond_poly(U, V)` with
  `U = 2(1+u²)p_x − Xn`, affine in `p`; the both-`R` branch divides by `1+u²`, which does not depend on `p`, so
  affinity survives.  `line_cond_iv` intersects five exact half-lines per bound choice; the hull over choices is
  certified because `C_k ∩ ℓ` is convex.  The slot table `((A,A),(B,B),(B,A),(A,B))` matches the monotonicity of
  `X = (aC+bS)/N`, `Y = (−aS+bC)/N` in the centre for `C, S ≥ 0` (θ ∈ [0°, 90°]).  All intervals closed. ✓
* **Lemma T** (both families, `θ = 0` with `T = ±∞`, the `u₀` shift for H).  Re-derived; correct.  `T + u₀ ≤ T + u`
  is the only use of `u ≥ u₀`. ✓  **Corollary T**: `f` is PL with breakpoints exactly as enumerated (`e ≤ τ↑` on the
  up line, `e − u₀` for `e ≥ τ↓` on the down line, `τ`'s themselves), constant on both tails; the code evaluates all
  breakpoints and one point per tail (`cand[0] − 1`, `cand[−1] + 1`).  The `'all'` case (`S = ℝ`) uses the whole `J₃`
  mass, as the proof's `τ = −∞` convention requires.  A bounded `S↑` is ignored (`τ↑ = +∞`), which is the safe side. ✓
  **Corollary T′**: steps `[min(T,τ↑) ≤ t_p]` (closed) match Lemma T's closed inequality; on open intervals STEP is
  constant; a crossing point is moved once (`idx in moved`); moved points get weight 0 in all three zeromargin arrays
  and are restored in a `finally`. ✓  (Not exercised by `s(21)`: the cover has no point on a segment line.)
* **Lemma L** (table re-derived for all 8 rows; typed/untyped classification; `r↑ ≥ a↑` from `α_k ∈ I_k ⊆ C_k`;
  untyped containment on `[b↓ − Δ↓, a↑ + Δ↑]` checked in `lemma_l_data`).  **General ends**: lower-hull edges of an
  increasing PL `g` with `g(0) = 0` have slope `≥ 0`, and a supporting line of the greatest convex minorant is `≤ g` on
  all of `[0, Δ]`; any `Δ > 0` is sound (the proof never uses `Δ ≥` the true movement), so the float-guided
  `move_bound` / `1e-12` floor are harmless. ✓  **L′ / hull form**: `x = min(r↑,B) ≥ A`, `y = max(r↓,A) ≤ B` hold; upper
  hull edges also have slope `≥ 0`; the `x < y` case gives a negative right-hand side. ✓
* **Corollary L.**  Both corner choices (`c₀` or `w/2`; `c₁` or `m − w/2`) give rectangles containing the admissible
  one at every `u`, so the float choice of split points and of the term per sub-bin is sound.  Concavity: every term is
  `min` of affine functions of `min_k t_k` (slope `≥ 0`) or of `−max_j t_j`.  The option `i*` chosen by float is
  corrected by the exact slack `σ ≥ 0`.  Bernstein-ratio bound requires all `β_i(Den) > 0`, checked;
  `S`-denominators refused when `u₀ = 0`. ✓
* **Lemma V.**  Sound, and more simply than stated: the Lemma L′ term is valid at *every* pose (it can only be
  negative), and `0` is valid everywhere, so the split is correct for any function `D`; soundness reduces to
  `region_phi`.  The float choice of `k, j` is therefore irrelevant to soundness. ✓ (nit N3)
* **Lemma R / SPLIT / `region_phi`.**  The region of a pose is `r = max{j : G_{q_j} ≤ 0}` (a prefix by the exact chain
  monotonicity).  `D_r`/`U_r` membership by binary search is monotone in `j` (also across kinds: `G_p ≤ G_q ≤ 0` puts
  `p` in `Q` whatever `p`'s kind).  `_gle0` is called with weight `−1`, outside its docstring's "positive" wording, but
  its corner argument (affine in the centre) holds for any real weights (nit N4).  `g_coeffs` equals zeromargin's
  `_gcoef` (re-derived from the geometry in my `g_exact`).  Region polygons: two parallel constraint lines, so every
  vertex is a rectangle corner or a constraint ∩ edge line; a true vertex at `u` is never "proven" excluded on its
  sub-bin; the constraint coefficients (`±2C`, `±2S`; for V `0` or `±N²/(SC)`) never vanish on a bin with `u₀ > 0`.
  `EMPTY` is therefore a proof of emptiness.  Points `T`, `D_r ∪ U_r` and pieces are disjoint; the phantom is excluded
  (`n = len(P) − 1`); the skip test with `Lbox` (possibly the parent's) is valid. ✓
* **Assembly, disjointness.**  (A) groups + non-group lines; (B) joint Lemma L over all lines + cores; within
  `lin_or_core` the three summands (`base_`, cores of eligible lines not in `datas`, `max(part_L, csd)`) cover disjoint
  line sets.  `rest = total − max(csd, part_L)` is exactly the pieces outside `datas`.  The float reach window drops
  whole segments, which lowers every bound (each lemma bounds the mass of *a sub-collection*). ✓
* **Lemma P, inheritance, stale state.**  zeromargin's primitives read `W`, `Wf`, `Wnum` live (no weight caches:
  `_adm_ctx_cache` is keyed by specs/bin and holds no weights), so the per-box phantom weight cannot go stale.  The
  phantom is set before every primitive call in a box; its `inh` bit is cleared in the mask handed to children
  (`kid[ph] = False`); `_last_inT` is reset to `None` per box and only ever records geometric membership, so a
  child cannot inherit a stale phantom weight or a weight-0 (moved) point's weight.  `P1`'s mask cannot select the
  phantom (`px + 1 ≥ m` and `cx₁ ≤ px + t` both false).  `L` is floored to `1/Wden` (`Wden | 10¹¹`, int64-safe);
  `L ≥ 0` always.  Children inherit `max(L, Lpar)` computed on the (clipped) parent box, of which they are sub-boxes;
  the with-points bound `L2` is never inherited (it would double count). ✓
* **Symmetry.**  D4: `θ ↦ θ` under the rotations, `θ ↦ 90° − θ` under reflections, so `[0, m/2]² × θ ∈ [0°, 45°]`
  suffices and `u ≤ ½` (θ ≤ 53.1°) covers it.  `--full` and the D2 default re-derived; correct (not used).
* **FORMAT.md reduction.**  A packing of 21 unit squares in side `s' < 5`, scaled by `5/s'` and shrunk
  concentrically, gives 21 pairwise disjoint closed unit squares in `[0,5]²`; `Σ μ(Q_i) ≥ 21 > total ≥ μ(∪Q_i)` for a
  finite positive measure.  Correct; `s(21) ≤ 5` is trivial.

## 2. The run

* **Coverage** (`zm_mixed_audit.py roots`, independent enumeration): `x1003_full.jsonl` has 40,000 distinct roots,
  exactly the grid `[0, 5/2]²` at pitch `1/20` × 16 `u`-bins of `[0, ½]`, no duplicates, volume `25/8` = the region;
  single label; boxes sum to 461,204; `UNCERT` and `unc` are 0 in every line.  The log has no `resume:` line: one
  uninterrupted invocation (06:44–08:43; `zm_mixed.py` last modified 04:46).
* **File** (`zm_mixed_audit.py file`, my own parser): 7,536 points, 1,872 segments (8 lines `x, y ∈ {1,2,3,4}`),
  all in `[0,5]²`, `w ≥ 0`; total `522368729933/25000000000 = 20.894749197320 < 21`, as logged.  Measure invariant
  under `x ↦ 5−x`, `x ↔ y`, `y ↦ 5−y`, rotation by 90°.  It is exactly the candidate × `1003/1000`.  sha256
  `8b415cee…c2a905` matches `x1003_sha.txt` and the checked file.
* **Oracle.**  `zm_mixed.exact_mass` agrees exactly with my independent oracle at 300 poses (θ = 0, `u ~ 10⁻⁶..10⁻¹²`,
  germ and wall poses).
* **Settings**: logged pitch/ubins/depth/CHAIN/T/L match `ZM_MIXED.md`; SPLIT, theta-bias and chain-from are not
  logged (defaults assumed) — S2.

## 3. Adversarial tests (new)

**Component tests** (`zm_mixed_audit.py leaves|boxes|tprime`).  Instead of "total ≥ 1" (which any correct cover
passes), each certified *component* is checked against the exact mass of that component, at adversarial rational
poses: box corners clamped to admissibility, `u₀ + 10⁻⁹, 10⁻¹²`, exact germ pivots `c = ξ ± 1/(2cos θ)` (and `± 10⁻¹²`),
a vertex of `Q` exactly on a line, the wall `c = w/2`, float minimisers of the piece mass and of the total mass:
`piece mass ≥ L` (this box) and `≥` the dumped phantom weight (incl. inherited `Lpar`); `point mass + phantom ≥ 1`
for ADM/CHAIN leaves; for SPLIT: every region's claimed points are in `Q`, region piece bound `≤` exact piece mass,
no pose in an `EMPTY` region, `points + region bound ≥ 1`.

| test | boxes / leaves | exact checks | violations | tightest slack |
|---|---|---|---|---|
| real m5 cover, 331 adversarial boxes (germs, wall band, vertex-on-line, random; sides to `1/2560`, `u`-bins to `2⁻¹⁶`) | 331 | 19.6 M | **0** | `0` (Lemma L / SPLIT region bound = exact mass) |
| real m5 leaves, hard cell `[1.2,1.3]×[1.5,1.6]` (SPLIT-heavy), 1,500 random of 21,201 | 1,500 | 43.7 M | **0** | `0` (region bound), `4·10⁻⁶` (points + phantom) |
| real m5 leaves, wall/dip area | 1,018 (all) | 0.78 M | **0** | `2.4·10⁻⁵` |
| real m5 leaves, germ `[1.45,1.55]²` | 744 (all) | 0.22 M | **0** | `1.1·10⁻⁴` |
| Corollary T′: grid cover with points **at the crossings** and along the lines | 150 | 7,403 (3,531 moves) | **0** | `0` |

Also rerun: the author's `zm_mixed_test.py selftest --n 100 --seed 41`: PASS (all 12 sections, 0 violations).

**Rejection tests on holed versions of the real cover** (the checker must leave every exact hole uncertified; the
certified leaves are then component- and total-tested — on a holed cover the total test is meaningful):

| cover | hole (exact, independent oracle) | run | result |
|---|---|---|---|
| R1: x1003 × 0.994 (tilted dip; a vertex of `Q` crosses `y = 2`: Lemma L′/V, SPLIT) | `0.9997246` at `(0.57364, 1.44062, 9.218°)` | 9 roots around it, depth 20 | NOT VERIFIED, 36 uncertified; the hole lies only in an UNCERT box; 826 certified leaves: 329k checks, 0 violations, min total in certified leaves `1.00024` |
| R3: x1003 minus the pieces `x = 1, y ∈ [1.50,1.52]` and `x = 2, y ∈ [1.48,1.50]` (D4 images) — a **germ-pivot** hole: `1.80` at `θ = 0`, `< 1` only at `θ > 0` | `0.9871818` at `(1.50124, 1.49923, 0.148°)` | the root containing it (`[1.5,1.55]×[1.45,1.5]`, `u`-bin 0), depth 14 | NOT VERIFIED, 1,703 uncertified; the hole lies only in an UNCERT box; 803 certified leaves (748 ADM, 55 CHAIN): 236k checks, 0 violations, min total in certified leaves `1.00048` |
| R2: x1003 × 0.9925 (anti-correlated `θ ≈ 33°` cell, SPLIT-heavy) | `0.9992940` at `(1.23383, 1.55544, 33.23°)` | 4 roots, `u`-bin 9, depth 20 (62,144 boxes, 1.4 CPU-h) | NOT VERIFIED, 10,385 uncertified; the hole lies only in an UNCERT box; 500 SPLIT + 800 ADM/CHAIN certified leaves sampled: 50.6 M checks, 0 violations, min total in certified leaves `1.000016`, min points + region bound − 1 = `1.0·10⁻⁶` |

## 4. Findings in detail

**S1 (should-fix).**  `runs/` is in `.gitignore`: `runs/line-cover_m5_candidate_x1003.txt`, `runs/zmm/m5/x1003_full.{log,jsonl}`
and `x1003_sha.txt` are not in the repository.  Bundle them in `certificates/s21/` (as for `s(32)`).

**S2 (should-fix).**  The `cert` header prints zeromargin's sha but not `zm_mixed.py`'s, `mixed_cover.py`'s or the input
file's, nor `use_split`, `theta_bias`, `chain_from`.  Patch (in `main()` after the zeromargin line):
```python
for p in (__file__, MC.__file__, a.path):
    print(f"sha256 {hashlib.sha256(open(p,'rb').read()).hexdigest()}  {os.path.relpath(p)}")
print('argv:', ' '.join(sys.argv))
```
**S3 (should-fix).**  `--resume` accepts any line with a matching `(label, root)`.  Write the header (shas + argv) as
the first jsonl line and refuse to resume on a mismatch.  (The audited run did not resume.)

**S4 (should-fix, evidence).**  `zm_mixed.py stress` re-checks only poses whose *total* float mass is `< 1 + 10⁻⁶`;
on a cover with `0.6 %` margin that set is empty, so "0 exact violations" on the m5 dumps is no evidence about the
checker.  The component tests above are; so are the rejection runs.  Suggest citing both in the write-up.

**S5 (should-fix, scope).**  For the certificate, run with Corollary T′ and polygons unreachable (or assert
`not cov.line_points and not cov.polys`), so a reviewer need not trust unexecuted code.

**Nits.**  N1 `min_density`, `cone_slope` are dead code.  N2 `ZM_MIXED.md`: §7 precedes §6; Lemma L is stated in its
`ρ_min` form, the code uses the (also proved) hull form — lead with the latter.  N3 Lemma V's proof can drop the role of
`D` (any split is valid).  N4 `_gle0` docstring says positive weights; `zm_mixed` passes `−1` (valid, see §1).
N5 `_group` appends to `moved` even when the group bound is finally not used (assembly B or `core` wins): sound (the
points are merely withheld), wasteful.  N6 the D4 check (multisets of segments) is sufficient, not necessary: a D4
measure cut into segments asymmetrically would be refused — fine for a certificate.

## 5. Reproduce

```
A=search/zm_mixed_audit.py; F=runs/line-cover_m5_candidate_x1003.txt
python3 $A file $F;  python3 $A roots runs/zmm/m5/x1003_full.jsonl;  python3 $A oracle $F --n 300
python3 $A boxes $F --n 400 --per 30 --seed 7
python3 $A leaves $F runs/zmm/m5/st_hard.leaves --per 20 --max 1500 --seed 5      # likewise st_wall, st_germ
python3 $A tprime x --n 150 --per 25 --seed 3
python3 search/zm_mixed_test.py ...   # holed covers: scale_cover(F, out, 994/1000), (9925/10000); germ hole: see runs/zmm/audit/
python3 search/zm_mixed.py cert runs/zmm/audit/m5_x1003_x0994.txt --d4 --disj --depth 20 --pitch 1/20 --ubins 16 \
    --nproc 1 --cx-lo 0.5 --cx-hi 0.65 --cy-lo 1.35 --cy-hi 1.5 --u-lo 0.0625 --u-hi 0.09375 --dump R1.leaves
python3 $A locate R1.leaves 356626453/621690037,756233892/524946169,71925427/892173648
python3 $A leaves runs/zmm/audit/m5_x1003_x0994.txt R1.leaves --per 30 --seed 21
```
