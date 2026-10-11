# PR #6 spot check (2026-10-10, after merge 060c1a1)

Evan's call (10-10): merge without replaying the ~7,000 process-hours of root checks ("optimistic scheduling"); run
the cheap checks.  Static review: `review-2026-10-10.md`.

## 1. Checker and Bridge build (executed)
`lake build Sqpack.Valid7Checker` (the new opt-in entry: all 30 modules of the PR, Bridge included), under a 60 GB
systemd cap: exit 0, 3 min 19 s wall, peak 19.5 GB; `V7.Bridge` 123 s.
`#print axioms SquarePacking.LemmaELeaf.box7Cover_eq_mcoverP` (main's `box7Cover` = the checker's cover):
`[propext, Classical.choice, Quot.sound]`.

## 2. All 9,800 root modules regenerated and scanned (executed; not built)
* Release `data/valid7-lean-certs-v1` of wand125/square-packing: `SHA256SUMS`, `v7cert-files.sha256`, `v7cert.tar.zst`
  match GitHub's asset digests and `SHA256SUMS`; the archive lists only `v7cert/R?????.pkl` + `v7cert.sha256` (no
  links, no absolute or `..` paths); all 9,800 pickles match `v7cert-files.sha256`.
* `lean/scripts/emit_v7.py --chunk 40` (merged code) in Docker `v7emit:local` (python 3.12, python-flint 0.9.0,
  sympy 1.13.3, numpy 2.4.2): `--network none --cap-drop ALL --security-opt no-new-privileges --user 65534
  --read-only`, repo and pickles mounted read-only.  Roots of `scripts/v7_xroots.txt` (2,612) with `V7_X=1`, the
  other 7,188 without, as in VALID7.md.  77 s + 35 s wall; 221,053 files, 10.62 GB.
* `scan_v7.py` (ours, text only): no forbidden token (axiom, sorry, native_decide, macro/elab/syntax, #-commands,
  attributes, instance, unsafe/extern/implemented_by, Lean/IO, debug, ...); `set_option` only `maxRecDepth 100000` /
  `maxHeartbeats 0`; one namespace (`SquarePacking.LemmaELeaf`), fixed header; 9,166,096 declarations (9,029,819
  theorems, 136,277 defs), every name a plain generated name (`root_NNNNN` / `rNNNNN…`, no dots), none repeated; all
  9,800 `root_NNNNN` present; imports outside the generated tree only `Sqpack.V7.Data`, `LemmaEDec`, `LemmaESplit`,
  `CapK`, `ExTreeX`, `LemmaESplitX` (all merged checker code).  Result: **0 problems**.
* So the generated sources cannot shadow `Bentz.Valid7` or weaken the kernel check; what remains is that each root
  compiles, which is wand125's build (not rebuilt here).  `V7/Main.lean` comes from `gen_main.py` and our own leaves
  file (no pickles); not generated or built here.

## 3. Not done (later, when cores are free)
A sample of root builds (5–10, heaviest included; ~40 min each on average) and `Main` (needs all roots).
Regenerating the sources takes ~2 min with the commands above.
