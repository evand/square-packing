# wand125's ValidTilt9 record check: our replay (2026-10-10)

Checker: wand125/valid7-independent-check tag `records-tilt9-v1` = c561dbb3fcea9178599db017f51f56b42ad9f91f (fetched into
`~/math/_untrusted-third-party/wand125-valid7-independent-check`, provenance logged; tarball of the commit byte-identical).
Records: release `records-tilt9-v1`, downloaded outside the container; all six sha256 lines of `records.sha256` match.
Cover `cover/K4_k008_box9.txt` byte-identical to `certificates/k2m4/K4_k008_box9.txt` (sha256 4151d7c4…7801).
Evan approved the run (10-10).

Sandbox: Docker `v7check:local` (python:3.12 + python-flint 0.9.0, as 10-03), `--network none --cap-drop ALL
--security-opt no-new-privileges --user 65534 --read-only`, source and records mounted read-only, 12 GB cap.

Command: `python src/check_record.py cover/K4_k008_box9.txt tilt9_{a,b,c}.jsonl --claim tilt --recheck 2000 --recheck-b 20 --seed 20261010`
(fresh seed: different leaves from wand125's seeds 1001–1015).

Result: exit 0, `RECORD OK`; `claim tilt roots 28350 leaf kinds {CORE 6565165, TIERB2 2940689, EMPTY 31489}` (as reported in #4).
7 min 16 s wall, one core.  Mutant covers not run this time.
