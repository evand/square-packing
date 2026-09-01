# Loop speed: the exact verifier stops being the separation oracle (task K)

Code: `search/floatsep.py` (new), `search/branch.py` (`--sep`, `--exact-every`, `--sep-pitch`,
`--sep-dtheta`, `--stab`, `--seed-from`, `--corner-rows`), `verify/src/main.rs` (the sweep, and
`VERIFY_STATS`), `runs/` scripts named below.  **Every default reproduces the old behaviour**:
nothing changes unless one of the new flags is given.

**Soundness is untouched.**  The float oracle only decides *which rows the LP gets*.  A row is a
valid constraint by construction (see `floatsep.py`, "Rows are valid") whatever the scan computed,
the cutting-plane loop only ever *stops* on an exact verdict, and the file that is finally claimed
still goes through the unchanged Rust verifier and `xcheck.py`.  The sweep change in
`verify/src/main.rs` is guarded by bit-identity against the pre-change binary.

## 0. Where the time went

Measured on this box (32 logical / 16 physical cores) while the live leaves J12, J16 and J6 were
running with 6 verifier threads each; load average is quoted beside every timing and every run
below was pinned with `taskset` to 8 distinct physical cores (`8-15`) unless stated.

From the live logs the loop's shape was: J16 round 9 = 3.8 h of which the verifier 3.0 h and the
LP 40 min; J12 round 9 = 4.8 h of which the verifier 4.7 h and the LP 55 s.

TODO(fill): baseline table, after table, per-item speed-ups, validation.
