import Sqpack.M2.Part0

set_option linter.style.longLine false

/-!
# `M2`: the mixed box tree covers the D4 fundamental region (generated; do not edit)

376 `Z` leaves, 664 `E` leaves, 600 clips, in 10 chunks, 1 files.
-/

namespace SquarePacking.M2

open BoxTree ZMTreeM

theorem cov_rootL : CovM 1000 4096 2000 4294967296 1000000 ptsL segsL 0 4096000 0 4096000 0 2147483648 :=
  (CovM.splitX 2048000 ok0 (CovM.splitX 2867200 ok1 (CovM.splitX 3276800 (CovM.splitY 2048000 ok2 (CovM.splitY 2867200 ok3 (CovM.splitY 3276800 (CovM.splitU 1073741824 ok4 (CovM.splitU 1610612736 ok5 (CovM.splitU 1879048192 ok6 ok7))) ok8))) ok9)))

theorem cov_root : CovM 1000 4096 2000 4294967296 1000000 pts.toList segs.toList 0 4096000 0 4096000 0 2147483648 := by
  rw [pts_toList, segs_toList]; exact cov_rootL

end SquarePacking.M2
