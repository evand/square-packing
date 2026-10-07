# packer/data
* `cen4_rec.json`, `cen4_cand.json`: census push (10-04) from the s(110) record and the candidate (Couzo's packing): per trial
  σ, trial, soft-squeeze side, final side (slp2), first-order jam flag, full-line count, fingerprint.  Source of the numbers in
  `../s110-landscape.md`.
* `s110_minima/`: one configuration per distinct jammed sub-11 side (rounded 1e-7) from that census (packer text format).
  Caveat: slp2 false-jams at corner–corner touches, so some are not true local minima (exact solver: 5 of 7 tested were not).
  Re-check them with `../exact/` after the slp2 fix.
* **`cen7/` (10-05, the clean census; supersedes cen4 and `s110_minima/`):** 1000 trials (rec + cand, σ 0.01–0.1 × 125),
  slp2 after the false-jam fix, every sub-11 output solved by `../exact/exactsolve.py`.  `summary.txt` = `census_exact.py`
  output; `census_exact.json` per trial (side, jam, lines, fingerprint, exact class / status).  `minima/`: one exact
  configuration (`.exact.txt`, 70 digits) and rational certificate (`.cert`, check with `../exact/verify_cert.py`) per distinct
  certified sub-11 minimum (21 after the 10-05 exactsolve fixes), named by S_exact to 1e-20.
