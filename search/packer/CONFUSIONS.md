# Open confusions (packer)

Started 2026-10-11.  Things we observed and do not understand, contradictions between measurements, and places where a
tolerance or a definition could be hiding a wrong conclusion.  One line each: first seen (file:line or run), what would
resolve it, status.  Claims they touch: `CLAIMS.md`.  Close an entry by recording the answer and the evidence, not by
deleting it.

| ID | confusion | first seen | would resolve it | status |
|---|---|---|---|---|
| X01 | 225 n = 110 points stay "not jammed" (exactsolve corner-corner MILP dS < -1e-6) after both descend.py and kick + polish: a real descent neither method realises, or a tolerance-degenerate verdict? | `runs/census/n110/finish_backlog` (10-11) | dS distribution of those points; exactsolve at tighter / looser thresholds on a sample; inspect the MILP's descent direction | open |
| X02 | Expected 5.7 new minima per proc-h (L1/L2), live baseline 1.0 (insert110) | PACKER 10-10 (insert110) | census v2 with measured CPU and a frozen reference: saturation vs run length vs load (worker-wall unit) | open |
| X03 | Reports drift: they read the live `basin` / `minimum` tables (c110a: 30 vs 16 on 10-09, 33 vs 20 on 10-11) | REVIEW §3, PACKER 10-11 | snapshot the reference tables per report (pre-census-v2 item) | open |
| X04 | Discovery credit is order-dependent (first reach): matched not-min vs certified parents 29 vs 45 per 1000, but the minimum was expanded first in 236 / 308 pairs | `parent_retro.py` (10-11) | random parent order inside an arm, or score by per-minimum hit counts rather than first reach | open |
| X05 | Which minima count: same exact S with different contact graphs are merged (min_key 1e-17); near-twins like 10.9968869876 / ...869 | TODO (10-10), identity.py | contact-graph tiebreak for equal-S minima; re-solve near-twins keeping contacts | open |
| X06 | 40 / 760 minima at 110 are weak (PSD with zero modes / second order inconclusive): genuine degeneracy, or solver tolerance? | backfill (10-11) | re-solve the 40 with tighter tolerances; check whether the zero modes are realised flat motions | open |
| X07 | exactsolve: points with -1e-6 < MILP dS < 0 go on to KKT (comment assumes "MILP min >= 0") | REVIEW §1 (`exactsolve.py:843-848`) | read the branch; test on points with tiny negative dS | open |
| X08 | New vs old explorer recipe disagree by start group (L2: old faster on far groups, new on the near group) | PACKER 10-10 L2 | replicates (census v2 arms) | open |
| X09 | Crossover variant effect unresolvable (bestfit2 lowest in all three arms, 1-11 per cell) | PACKER 10-10 | more events, or drop the question (mate rule carries the effect) | open, low priority |
| X10 | tx237m vs tx237n / s: replicates differ by 3e-3 | PACKER 10-08 | replicates at 237 | open, low priority |
| X11 | What order parameter separates the sub-11 funnel from above-grid glass?  Q4 keeps packings off the grid but does not separate (u110a, bias110) | PACKER 10-09 | line-load bias, other coordinates; umbrella in a candidate q | open |
| X12 | The battery's clean-start filter leaks (compared against the 30 best sub-11 references only; the 0.069 "long jump") | PACKER 10-09 | refreeze with all minima as references (holdout work) | open |
| X13 | Is the 110 "second family" (d ~0.062, 1e-4 above the record; 10.9968870 most-picked) a separate funnel or a side valley? | PACKER 10-10 | far-start / holdout arm in census v2 | open |
| X14 | Does duplicate-biased parent picking hurt? (C16; C12 says copies expand like their minimum) | PACKER 10-10 | census v2 per-minimum picks arm | open |
| X15 | Why does the chain beat dedicated runs on 343 (archive memory?) | PACKER 10-10 Hunt 3 | compute-matched chain vs dedicated (TODO) | open |
| X16 | Do finisher settings, break-even and rates measured at 110 transfer to other n and packing families? | Evan 10-11 | out of scope for now; return after census v2 | parked |
| X17 | `5.7 / proc-h`, `48 / CPU-h` etc. before 10-11 are worker wall-clock on a machine at varying load: cross-run rate comparisons may be load artefacts | PACKER 10-11 (compute units) | measured CPU in manifests from now on; avoid cross-run comparisons of older rates | mitigated going forward |
