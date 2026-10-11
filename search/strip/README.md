# Float scripts for the s(12) strip leaf (10-11)

Results and reading: `notes/strip-split.md`; task list: `tasks/strip-leaf/README.md`.  All float, one-sided
(a fit found is real; "no fit" is search evidence only).  Copied from the session scratchpad as run; they expect
their inputs (`cfg_n8_H3.99.npy`, `census4*_H*.json`, ...) in the working directory, produced by `wbisect.py` and
`census4*.py`.  Core: `rectpack.py` (penalty + separating-axis gap), `census4*.py` (strip classes), `seam*.py`
(seam reach and region tests), `notch*.py` (notched region), `mkpage.py` (the shape-class page).
