# arch-farfield: does a near-axis mass cap close the far field?  (lemma A3 of `notes/proof-architecture.md`)  (2026-09-20)

Question: at container side `t = T`, is the fractional-packing / cover LP value `< n = T^2-T` once the total mass on
poses of tilt `< eps` is capped at `T-1`?  Do `T = 3, n = 6` first if the pipeline allows a different container (the claim
is a theorem there), then `T = 4, n = 12`.  `eps` in `{0.5°, 1°, 2°, 5°, 10°}`.  Instruments: `search/bentz.py`
(tilt bands, orthogonal counts, `BENTZ.md` §7.2 — the `|theta| <= 1°` run is the complementary band),
`search/leaf_ceiling.py`, `search/closed4.py`; read `notes/status.md` first.  Report packing-side (LP lower bound on the
relaxation) and cover-side (certified upper bound) separately and never confuse them (`notes/review-2026-09-13b.md`
"Clarifications").  Also report the variant "no wall-to-wall unit strip holds `T` near-axis squares" if it is expressible.
A value `>= n` at every `eps` is a full answer: say so plainly and show the fractional packing that achieves it.
Compute: `<= 8` threads (another job holds 16 until ~23:00); anything over 10 min detached with `setsid nohup`, never
`pkill`, leave no `until ... sleep` / `tail -f` watchers.  Deliverable: `search/ARCH_FARFIELD.md` + new script(s) under
`search/arch_*.py`; runs in `runs/arch_farfield_*`.  Do not edit existing files.  Do not commit.
