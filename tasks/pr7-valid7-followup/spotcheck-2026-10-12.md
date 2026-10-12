# PR #7 spot check (2026-10-12; wand125, "Valid7 follow-up", head 4eb59c4 on 7dfca4e)

Evan's call (10-12): spot check, then merge; budget 8 cores for an hour or several.  The work ran in a temporary
worktree of the PR head.  It shared main's `.lake/packages` and used a copy of main's `.lake/build`.

## 1. Static review
* 7 files: `ExTreeTilt.lean` and `LebMassBridge.lean` are new; `gen_face.py` is removed; `gen_main.py`,
  `Valid7Checker.lean` (imports + header), `VALID7.md` and `Axioms.lean` (one comment) change.  `Sqpack.lean`
  is unchanged, so the default target is untouched.  The PR changes no root module and no root generator
  (`emit_v7.py`, `runB.py`).
* The new `.lean` files have no forbidden tokens: no axiom, sorry, native_decide, attributes, set_option, instance,
  #-commands or macros.
* `Main` now gets `validTilt7` from `v7_tree` by `validTilt_of_tree`, then `valid7 := Bentz.valid7_of_tilt validTilt7`.
  That is main's own Lemma Z, `ValidSplit7.lean:87`, from `validAxis7`.  `bentz` is as before.

## 2. Builds (executed)
* `lake build Sqpack.Valid7Checker Sqpack.ExTreeTilt Sqpack.LebMassBridge` under a 60 GB cap and 800 % CPU: exit 0
  in 12 s wall.  Only the new modules and the checker were rebuilt; the only output was long-line lint warnings.
  `#print axioms` on `validTilt_of_tree`, `CovT.toLebMass`, `CovT.ofLebMass` and `volume_inter_le_mcoverP` gives
  `[propext, Classical.choice, Quot.sound]`.
* `Main` with the root tree as a hypothesis: the PR's `gen_main.py` generated `Main.lean` (9,800 roots; grid 35 × 35 × 8).
  The test turned `theorem v7_tree : CovT … := glue` into `variable (v7_tree : CovT …)` + `include v7_tree`, and
  dropped the root imports.  The rest is verbatim.  It compiles in 3 s.  `validTilt7`, `valid7` and `bentz` depend on
  `[propext, Classical.choice, Quot.sound]`, and `bentz : CovT 5 128 35 … 0 2240 0 2240 0 2361183241434822606848 →
  ∀ k, 6 ≤ k → minSide (k^2 - 3) = k`.

## 3. Sample root builds (executed; the item left open by the PR #6 spot check)
* Pickles: release `data/valid7-lean-certs-v1`.  The archive matches `SHA256SUMS`, and all 9,800 pickles match
  `v7cert-files.sha256`.
* Sample, chosen by pickle size (a cost proxy: 5.96 CPU-s per KB, from 6,650 CPU-h over 4.02 GB):
  * the heaviest root, R06407;
  * one root drawn at random (seed 20261012) from each size bin 50–75 %, 75–90 %, 90–97 % and 97–99.5 %: R06498,
    R07699, R06131, R02753;
  * one from the lower half (seed 20261013): R03361.
  R07699, R06131 and R02753 are in `v7_xroots.txt` and were emitted with `V7_X=1`; the others without.
* Emitted with `emit_v7.py --chunk 40` in Docker `v7emit:local` (`--network none --cap-drop ALL --security-opt
  no-new-privileges --user 65534 --read-only`; repo and pickles mounted read-only).  This gave 939 modules, 33 MB.
  `../pr6-valid7-kernel/scan_v7.py` found 0 problems: 30,876 declarations, 6 root theorems, imports only
  `V7.Data`, `LemmaEDec`, `LemmaESplit`, `CapK`, `ExTreeX`, `LemmaESplitX`.
* Kernel checks: `scripts/runB.py --list … --jobs 8 --mem 90` under a 100 GB cap.  All 939 modules passed, in
  2 h wall (08:22–10:22) and 12.1 CPU-h:

| root | X | pickle | CPU-h | peak GB |
|---|---|---|---|---|
| R06407 | no | 8.04 MB | 8.39 | 33.3 |
| R02753 | yes | 2.63 MB | 1.87 | 19.1 |
| R06131 | yes | 1.36 MB | 1.20 | 13.4 |
| R07699 | yes | 0.50 MB | 0.52 | 11.9 |
| R06498 | no | 0.26 MB | 0.10 | 13.1 |
| R03361 | no | 173 B | 0.00 | 5.9 |

  Measured cost tracks the size proxy: R06407 was estimated at 13.3 CPU-h and took 8.4.
* `#print axioms` on `root_02753 … root_07699` (all six) gives `[propext, Classical.choice, Quot.sound]`.
  `root_06407 : CovT 5 128 35 … 1408 1472 960 1024 0 295147905179352825856`, a cell of the grid that `v7_tree`
  glues.

## Not done
The other 9,794 roots, and `Main` with all of them (it needs 53–58 GB).  The sample is 0.3 % of the certificate bytes.
