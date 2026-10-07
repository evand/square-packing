# k2m3-bundle: certificates/k2m3/ for s(k²−3) = k, k ≥ 6 (2026-09-30)

Paths relative to `~/math/square-packing/public/s12`.  Source of truth: `search/QUADRANT_EXACT.md` (read in full, incl.
§6 run V3), reviews `~/math/square-packing/private/s12/tasks/k2m3-review/*/REPORT.md` (read all six), Lean
`lean/Sqpack/Bentz.lean` + `notes/lean-bentz-reduction.md`.  Model bundles: `certificates/s60/` and `certificates/s21/`
(README structure, verify.sh, SHA256SUMS, checker/ copies, manifest style) — follow their structure closely.

1. **`search/qx2_records.py`** (new; independent of qx2_zm's control flow, like `search/mixed_records.py`): reads the
   V3 record `search/qx2_data/cert/runV3_6294052a_leaves.jsonl.gz` and checks, exactly: header shas = the bundle's
   checker/ files + box file; the roots are exactly the D4 grid over [0,7/2]² × u ∈ [0,1/2] (pitch 1/10, 8 u-bins); every
   leaf kind is a certifying kind; UNCERT 0; per root, **leaves plus the slabs removed by clip_bin tile the root exactly**
   (recompute the slabs: `clip_bin` is zm_mixed's; you may import it, or better re-derive the admissibility condition
   w(u) = cos θ + sin θ vs the centre range and check each uncovered part contains no admissible pose, exactly) —
   the reduction-coverage reviewer's `rerun_sample.py` / `cov_check.py` in the review dir show one way; census totals
   match the .out.  Keep it fast (minutes).  Exit 0 iff clean.
2. **`certificates/k2m3/`**: `README.md` (claim; status line below; proof in one paragraph; what is checked by what,
   incl. the Lean row — `bentz_of_valid7` reduces the all-k claim to `Valid7`, which is the Python certificate's claim;
   files; re-checking; what is not machine-verified; review record summarising the six reviews + guards + V3), the box
   cover and family and exact json (copies of `search/qx2_data/L4_k02_*`), `qx2_zm/` with `checker/` (the exact files
   that ran V3: qx2_zm.py 6294052a…, zm_mixed.py 1fd20346…, zeromargin.py, mixed_cover.py — verify shas against
   the V3 header), the V3 record (.jsonl.gz, .out), the Lemma Z output (run `qx2_zm.py axis` and save its output),
   `verify.sh` (quick: hashes, `qx2_family_check.py`, own-parser total/D4 check (mixed_records.py cover works on this
   file? check), Lemma Z, `qx2_records.py`, and optionally the Lean data regeneration check; `--full`: re-run qx2_zm.py
   V3 settings and compare census root for root), `SHA256SUMS`.
   **Status line (use verbatim near the top):** "Working in public: a single-implementation exact certificate,
   adversarially reviewed by six independent agents with no errors found; the all-k reduction is kernel-checked in Lean.
   Not yet independently re-implemented, externally reviewed, or fully formalised."
   **Small cases:** write the attribution sentence with a marker `<!-- LIT -->` — the literature check is running in
   parallel; current belief: k = 3 Kearney–Shiu 2002; k = 4, 7 Bentz 2010; k = 5, 6 Bentz arXiv:1606.03746 (preprint).
   Our result is k ≥ 6.  Apply the reviews' wording fixes in the README (not in QUADRANT_EXACT.md; list the doc fixes
   you'd make there at the end of your report instead).
3. Run `verify.sh` (quick) and report its output.  Do not run `--full`.

**Rules.**  Work in the main checkout on main; touch only `certificates/k2m3/`, `search/qx2_records.py`; don't commit.
Don't edit other files (list anything that should change).  Cores 0–7 only (taskset).  Report ≤ 400 words.
