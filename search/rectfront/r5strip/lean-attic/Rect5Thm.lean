import Sqpack.Rect5Final
import Sqpack.Rect5Iff

/-! # R₅: exactly which rectangles hold five unit squares

`Fits5 w h` (see `Rect5Iff`) lists the cases: a row (one side ≥ 5), the wide part
(both sides ≥ 2 and one side ≥ 3 or `w + h ≥ 4 + √2`), and the zigzag strip
(`h < 2` and `w + 4h ≥ 4 + 6√2`, or the same with `w`, `h` swapped).  The strip lower
bound for `2√2 − 1 ≤ h < 2` is `final5_high` (machine certificates, `Rect5Final`). -/

namespace SquarePacking

theorem rect5_iff {w h : ℝ} : RectPacks 5 w h ↔ Fits5 w h :=
  rect5_iff_of final5_high

end SquarePacking
