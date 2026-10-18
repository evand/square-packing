# Certification pipeline: one entry point for "is this packing certified, and how far" (10-17)

Trigger: the s(229) hit (jlevy#674).  exactsolve: "not first-order jammed" (MILP dS* ~ -1.5e-7
that no descender realises: descend.py, nlpdesc.py, near-axis snap, 540 harden runs all return to 15.6501632548).
`dilate_cert.py` certified the upper bound in seconds (both checkers VALID, +1.6e-11), but candidate.py never tries it.

## Pieces we have (each written for one investigation)
| tool | gives | written for |
|---|---|---|
| exact/exactsolve.py | exact KKT point, rational cert, second-order checks | exact batch |
| exact/descend.py | MILP corner-corner kick + slp2 | false jams |
| packer/nlpdesc.py | second-order fixed-line NLP descent | s(110) slow valleys |
| packer/weak{kkt,cone,farkas,modes,probe}.py, conelp.py | weak / degenerate minimum analysis | s(110) weak rows |
| packer/rigid.py | first-order LP jam / rigidity | general |
| exact/gradec.py | grade-C local-min certificate | exact forms |
| packer/dilate_cert.py | upper-bound cert for any feasible float packing | 10-10 |
| packer/analysis/pathcert.py | same-basin path certificate | same_minimum |

## Plan (Evan 10-17: test first or refactor first?)
1. Done now (small, known-good): candidate.py falls back to dilate_cert when exactsolve fails, verdict
   "upper bound certified, minimality unresolved".
2. Benchmark of stuck points, fixed before running anything: s(229) (#674), s(292) corner case (pre-fix input), the 6 left
   of s(110)'s 40 unresolved keys, a sample of the 225 stuck points, the 17 grade-C timeouts.  Expected class per point
   (false jam / weak multiplier / degenerate second order / not a minimum).
3. Run every existing tool on every point; table: tool x point -> verdict, time.
4. Only then refactor: one `certify.py` (classify -> route), pieces moved behind it, duplicates removed (three
   fixed-line / slp variants, two LP solvers).  Gaps found in 3 become new pieces.
