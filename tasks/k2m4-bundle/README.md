# k2m4-bundle: certificates/k2m4/ for s(k²−4) = k, k ≥ 5 (2026-10-03)

Paths relative to `~/math/square-packing/public/s12`.  Model: `certificates/k2m3/` (README structure, verify.sh,
SHA256SUMS, qx2_zm/checker copies) and its brief `~/math/square-packing/private/s12/tasks/k2m3-bundle/README.md` — follow
them closely.  Sources: `search/K2M4_MARGIN.md` §4 (incl. the 2026-10-03 result), `runs/qx2_k4x_k008/` (gitignored: the run
of record `qxzm_full.{out,jsonl}`, `project.out`, `germscan.out`, `sol_exact.json`), `runs/qx2_k4x_k008_xloop.out` (Lemma Z
output), `search/qx2_data/K4_k008_{box9,family}.txt`, Lean `lean/Sqpack/Bentz4.lean` + `notes/lean-k2m4-reduction.md`.

1. **`search/qx2_records.py`**: generalise to the k² − 4 record (R = 3, box 9, roots over the D4 region at pitch 1/10 ×
   u-bins — read the run's header/argv for the exact root set; 16,200 roots).  The k2m3 invocation must behave exactly as
   before (k2m3's verify.sh must still pass unchanged).  Exact checks as for k2m3: header shas = shipped checker/ + box
   file; roots exactly the grid; every leaf kind certifying; UNCERT 0; per root leaves + clip_bin slabs tile the root
   exactly; census = the .out.  Fast (minutes).
2. **`certificates/k2m4/`**: README (claim; status line; proof in a paragraph; what is checked by what incl. Lean
   `bentz4_of_valid9`; files; re-checking; what is not machine-verified; review record — leave a `<!-- REVIEW -->` marker,
   the review runs in parallel), box + family + exact json, `qx2_zm/` with `checker/` (the files that ran: qx2_zm.py
   6294052a…, zm_mixed.py 1fd20346…, zeromargin.py 640fe453…, mixed_cover.py bb89de15… — verify against the run header;
   these equal k2m3's copies), the record (`.jsonl.gz`, `.out`), Lemma Z output (re-run `qx2_zm.py axis` on the shipped
   box and save it), germ scan output, `verify.sh` (quick: hashes, own-parser total + D4 check, `qx2_family_check.py` if it
   handles R = 3 (else say so), Lemma Z re-run compared, `qx2_records.py`, Lean data regeneration
   `lean/scripts/gen_bentzfam_data.py` compared byte for byte; `--full`: re-run qx2_zm with the run's settings — document
   cost ≈ 815k CPU-s, don't run it), `SHA256SUMS`.
   **Claim:** s(k² − 4) = k for every k ≥ 5 (Friedman's F₄ — check `search/FRIEDMAN.md` §0 and §5 for wording): k ≥ 8 by
   the family (Lean `bentz4_of_valid9` + Valid9 by this certificate); k = 5, 6, 7 by our bundles s21, s32, s45 (cite their
   status lines accurately; s(32) = 6 is hypothesis-free in Lean).  Prior literature: leave a `<!-- LIT -->` marker.
   **Status line (adapt from k2m3's, accurately):** working in public; single-implementation exact certificate of Valid9
   (same checker code as k2m3's, which six agents reviewed and which wand125's independent checker corroborated on Valid7);
   all-k reduction kernel-checked in Lean; not yet independently re-implemented for Valid9, externally reviewed, or fully
   formalised.
3. Run `verify.sh` (quick) and report its output.  Don't run `--full`.

Rules: main checkout, branch main; touch only `certificates/k2m4/`, `search/qx2_records.py`; don't commit.  List anything
else that should change.  Cores 0–3 only (taskset).  Report ≤ 400 words.
