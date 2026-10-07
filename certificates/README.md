# Certificates

One directory per result, each with its own README (status, what was checked, by which checkers), `SHA256SUMS`
and `verify.sh`.  The loose files at this level are the s(11)/s(12) point certificates, checked by `../verify.sh`.
Files named in a `SHA256SUMS` are frozen.  The format is [`FORMAT.md`](FORMAT.md).

| Result | Where | Re-check |
|---|---|---|
| s(k² − 3) = k, k ≥ 6 | [`k2m3/`](k2m3/) | `certificates/k2m3/verify.sh [--full]` |
| s(k² − 4) = k, k ≥ 5 | [`k2m4/`](k2m4/) | `certificates/k2m4/verify.sh [--full]` |
| s(60) = s(61) = 8 | [`s60/`](s60/) | `certificates/s60/verify.sh [--full]` |
| s(45) = 7 | [`s45/`](s45/) | `certificates/s45/verify.sh [--full]` |
| s(32) = 6 | [`s32/`](s32/) | `certificates/s32/verify.sh [--full]` |
| s(21) = 5 | [`s21/`](s21/) | `certificates/s21/verify.sh [--full]` |
| s(13) = 4, case-free | [`rung2/`](rung2/) `s13_closed_cover_4.txt` | `./verify.sh` (band), `--full` (whole domain) |
| unavoidable point sets: p(3) = 7; no half-turn-symmetric 13-point set for [0,4]² | [`unavoid13/`](unavoid13/) | `./verify.sh` |

## Loose point certificates (weighted points; `verify/` checks them over every placement)

| File | Proves | Note |
|---|---|---|
| `s12_lower_3.9686.txt` (+ `.json`) | s(12) ≥ 15680/3951 = 3.968616 | 1,736 points; the headline s(12) bound, kernel-checked in Lean (`s12_ge_15680_3951`, opt-in build).  squarepacker's s(12) ≥ 3.9715 re-weights these points |
| `s12_lower_3.9676.txt` | s(12) ≥ 980/247 = 3.967611 | 764 points |
| `s12_lower_3.931795.txt`, `_sparse.txt` | s(12) ≥ 3920/997 = 3.931795 | the first bound (788 points) and a 224-point sparse version |
| `s12_56points_3.8.txt` | s(12) ≥ 19/5 | 56 points; also re-checked by `xcheck.py` |
| `s12_uniform_<k>of<m>_<s>.txt` | s(12) ≥ s | uniform certificates: every unit square holds k of m points (`search/uniform/UNIFORM.md`); `7of81` gives 35/9 |
| `s11_lower_3.8143.txt` (+ `.json`) | s(11) ≥ 3040/797 = 3.814304 | superseded: s(11) = 3.877084… is now proved (Queuingtheorydotcom) |
| `s12_boxclique_demo_3.9318_N2000.txt`, `s12_anchorclique_demo_3.9318.txt` | s(12) ≥ 3920/997 | demonstrations of the clique and anchor-clique atom types (`FORMAT.md`) |
| `branch/s12_t3.98_corner_k{0,1,2}.txt` | branch leaves at s = 3.98 | corner-occupancy branches k = 0, 1, 2 of a branch-and-bound attempt at s(12) ≥ 3.98 (k = 3, 4 open; `search/BRANCH.md`); `../verify_branch.sh` |

The `.json` files are the same certificates in a portable form; `../verify.sh` checks that each converts back to its
`.txt` byte for byte.
