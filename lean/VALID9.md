# ValidTilt9 in the kernel

`Sqpack/V9/Main.lean` (generated) proves

```lean
theorem validTilt9 : Bentz4.ValidTilt9                        -- the open part of ValidSplit9.lean
theorem valid9 : Bentz4.Valid9                                -- Bentz4.valid9_of_tilt validTilt9
theorem bentz4 : ∀ k : ℕ, 8 ≤ k → minSide (k ^ 2 - 4) = k     -- Bentz4.bentz4_of_validTilt9 validTilt9
```

in namespace `SquarePacking.LemmaELeaf`, with

```
'SquarePacking.LemmaELeaf.validTilt9' depends on axioms: [propext, Classical.choice, Quot.sound]
'SquarePacking.LemmaELeaf.valid9' depends on axioms: [propext, Classical.choice, Quot.sound]
'SquarePacking.LemmaELeaf.bentz4' depends on axioms: [propext, Classical.choice, Quot.sound]
```

No `sorry`, no `native_decide`, no new axioms. The certificate is the k2m4 run `qxzm_full`
(`certificates/k2m4/qx2_zm/run_k4x_k008_leaves.jsonl.gz`) with the cover `K4_k008_box9.txt`. It is
glued over the D4 root `[0, 9/2]² × (0, 1/2]` in 45 × 45 × 8 = 16,200 root cells.

This records one run as it was done. The run is a single kernel check (first system): it was not repeated,
and there is no independent second check.

## Files

- `Sqpack/V9/Data.lean`: the cover as the checker's data.
- `Sqpack/V9/Bridge.lean`: identifies `box9Cover` with the checker's cover.
- `Sqpack/V9/Tilt.lean`: `validTilt_of_tree`.
- `Sqpack/RowMinN.lean`: `rmChk`, the row-minimum check of `RowMinK.rmLoop` computed on natural-number
  digits, and `rmChk_sound : rmChk … = true → rmLoop … = true`.
- `Sqpack/LemmaEKM.lean`: calls `rmChk` instead of `rmLoop`, with `rmChk_sound` as one more proof step.
  No statement changes. This was done for kernel speed: heavy row modules ran 3 to 5 times faster.
  The Valid7 roots on main were built before this change; their statements are unaffected.
- `scripts/gen_exact.py`: adds `V7_TOUCH` (default `1e-12`, the old constant). This is the chord length
  above which a segment counts as meeting the square at a sample point.
- `scripts/subm.py`, `scripts/kmmirror.py`: the rows' minima sub-bins (`V7_SUBM=1`) and their Python mirror.
- `scripts/v9/`: the files of the run that differ from the Valid7 ones, as used.
  - `run_v7.py` and `emit_v7.py`: the k2m4 cover and leaves, `D = 5`, `W = 2·10¹²`, `MQ = 45`, `Sqpack.V9.*`
    modules. They also add `V7_SGK`, a centre grid `SGK` times finer, whose root theorem is moved back to
    `SG0 = 128` by `CovT.scaleS`. The default `1` gives the old output.
  - `gen_main9.py`: writes `Main.lean`.
  - `run9.py`: the driver. It runs one root after another, starts a `lean` only when enough memory is
    free, and resumes from per-task markers. Its paths are those of the machines it ran on.

  To run them, copy `scripts/v9/*.py` over `scripts/` in a separate checkout; they replace the Valid7
  `run_v7.py` and `emit_v7.py` there.

## Certificates and the record of the run

The release <https://github.com/wand125/square-packing/releases/tag/validtilt9-lean-v1> holds:

- `v9-certs.tar.zst` (109 MB): the 16,200 per-root certificates as used (Python pickles; check
  `CERTS.sha256` before loading). The 9 certificates made by hand are in `certs/manual/`, with
  `MANUAL.tsv`.
- `SOURCES.sha256.zst`: the sha256 of all 1,275,911 generated files (roots, leaves, row groups, `Main`).
- `v9-logs.tar.zst`: `tasks.jsonl` of both machines (one line per task, with CPU seconds and peak
  memory), every driver version and the trial directories of the hand-made certificates.

The bundle <https://github.com/wand125/square-packing/tree/main/problems/square-lower-bounds/certificates/k2m4_n77_L9/lean>
has the README of the run, `ENVIRONMENT.md` (machines, versions, OS settings), the built `Main.lean` and
`Axioms9.lean`, and the receipts.

## Certificates made by hand (9 of 16,200 roots)

In five roots, the automatic certificate step failed, and options that the generator already had closed
them: 2270, 4262, 11715, 11827 and 11397. Those options were then added to the automatic retries.

Four roots needed the two new options:

| root | options |
|---|---|
| 12100 | `V7_TOUCH=1e-9` |
| 11424 | `V7_SGK=8 V7_SPLIT=10` |
| 11768 | `V7_SGK=8 V7_SPLIT=18` |
| 8917 | `V7_SGK=8 V7_SPLIT=18` |

- 12100: a segment touches the square's corner at u = 0 and lies outside it for u > 0. With the old
  threshold it counted as meeting the square, and its case split failed a zero-margin vertex.
- 11424, 11768, 8917: splitting stopped at the 1/640 centre grid, which `V7_SGK` makes finer.

Every leaf, group and root of these roots went through the kernel like all the others.

## Rebuilding

From `lean/` on this branch:

1. `lake build Sqpack.V9.Tilt Sqpack.V9.Bridge Sqpack.CapK Sqpack.ExTreeX Sqpack.LemmaEDec Sqpack.LemmaESplit Sqpack.LemmaESplitX`.
   On 48 cores this took 21 min. `V9/Bridge` alone takes 830 s and 92 GB; the other modules take seconds.
2. Copy `scripts/v9/*.py` over `scripts/`. Unpack `v9-certs.tar.zst`, put all pickles into one certificate
   directory, and run the driver.
   - Per root, the driver runs `emit_v7.py --chunk 20` with `V7_X=1`, then one `lean -o` per leaf, row
     group and root module.
   - For 11424, 11768 and 8917, emit with `V7_SGK=8`.
   - Optionally compare the generated files with `SOURCES.sha256`.
3. `python3 scripts/gen_main9.py` (without `V7_SGK`). Then
   `lake env lean -o .lake/build/lib/lean/Sqpack/V9/Main.olean Sqpack/V9/Main.lean`, and `#print axioms`
   on the three theorems.
   - `Main` imports 1.27 M `.olean` files (340 GB).
   - Linux needs `vm.max_map_count` above that number. With the default 1,048,576, the import fails with
     `failed to read file` on an intact olean.

## The recorded run

The run used two 48-core Arm machines (376 GB each, Lean v4.33.1) for about five days.

| step | CPU hours |
|---|---|
| leaf modules (`lean -o`) | 2,694 (peak 22 GB per `lean`) |
| row groups | 4,665 (peak 20 GB) |
| roots | 31 |
| certificates (all attempts) | 177, plus 5.8 for the four rebuilt roots |

- The kernel checks total about 7,390 CPU hours.
- `Main` took 584 s (230 GB resident, mostly mapped `.olean` files), and the axioms 222 s.
