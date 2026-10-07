# wand125's independent Valid7 checker: our check (2026-10-03)

Repo: https://github.com/wand125/valid7-independent-check @ 38dd31b3 (clone: `~/math/_untrusted-third-party/wand125-valid7-independent-check`,
provenance logged).  Reported in evand/square-packing#1 and jlevy/squares#296.  Evan approved (2026-10-03) running verify.sh's
record check only, sandboxed.

## Record check (executed)
Docker `v7check:local` (python:3.12-slim + python-flint 0.9.0), `--network none`, repo read-only, `--cap-drop ALL`, non-root,
cores 0–3, 12 GB cap.  Records of release `records-v1` downloaded outside the container; sha256 match `records.sha256`.
Result: 341 s wall, exit 0.  `roots 156800 leaf kinds {EMPTY 79927, TIERB2 2886043, CORE 6674090}`, `RECORD OK`; mutants M1
(wall, both tilt signs) and M3 (corner) `ok False`, M2 (Lebesgue) `COUNTEREXAMPLE` / NOT VERIFIED — all refused as expected.
Note: `--recheck 2000 --recheck-b 300` uses fixed seeds (random.seed(1), (2), check_record.py:115,122), so this re-certified
the same 2,300 leaves jlevy's run did.  A fresh-seed recheck would need a (one-line) edit to a copy: not approved yet.

## Method review (read-only agent; nothing executed)
Verdict: essentially sound, genuinely independent of ours (no D4, θ = 0 by closure, core polygon + outer hull, exact-in-u
symbolic execution with Sturm isolation).
* Cover byte-identical to `search/qx2_data/L4_k02_box7.txt` (sha c0a67507…3694b); parsing per FORMAT.md.
* Coverage/tiling: u ∈ [−½, ½] ⊃ a full 90° period; exact product grid of roots, exact bisection tiling, EMPTY rechecked.
* θ = 0 by closure: correct direction (closed squares, finite measure: limsup μ(Q_n) ≤ μ(Q)).
* Tier A, fixed-angle reduction (Brunn–Minkowski ⇒ Hessian has ≤ 1 positive eigenvalue ⇒ minima on face boundaries): sound.
* **Gap (logical, negligible here):** `tier_b2.py:279–298 nonneg_open` rejects odd-multiplicity roots then checks p(s) ≥ 0 at
  one rational s; passes wrongly if s is an even-multiplicity root of a p ≤ 0 elsewhere.  Fix: require p(s) > 0 or move s.
  Also used by `chain()`.
* Minor: line dedup / value caches keyed by Python `hash` (tier_b2.py:332, 432) — probabilistic, not exact; kind `TIERB`
  accepted syntactically only (never emitted); record header / cover hash not read by check_record.py; `mass_convex`
  collinear-polygon case unreachable here (tier_a.py:96–125).
* Bottom line: strong corroboration, not a fully re-verified second certificate — per-leaf assurance rests on their 626
  core-h run; verify.sh re-certifies a fixed 2,300-leaf sample.
