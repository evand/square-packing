/-
Opt-in entry point (not in the default target) for the ValidTilt7 / Valid7 kernel checker of PR #6 (wand125, merged
2026-10-10; `lean/VALID7.md`).  `lake build Sqpack.Valid7Checker` builds the checker and `V7.Bridge` (Bridge peaks at
~19 GB).  The 9,800 generated root modules, `V7/Face.lean` and `V7/Main.lean` (`valid7 : Bentz.Valid7`,
`bentz : ∀ k ≥ 6, minSide (k^2 - 3) = k`) are not in the repository: `lean/VALID7.md` regenerates them from the
release `data/valid7-lean-certs-v1` (~7,000 process-hours).
-/
import Sqpack.BernsteinZ
import Sqpack.CapK
import Sqpack.ChordE
import Sqpack.CovMP
import Sqpack.ExTree
import Sqpack.ExTreeX
import Sqpack.KArith
import Sqpack.KBern
import Sqpack.LebCorner
import Sqpack.LebCredit
import Sqpack.LebE
import Sqpack.LebTan
import Sqpack.LemmaE
import Sqpack.LemmaEDec
import Sqpack.LemmaEK
import Sqpack.LemmaEKM
import Sqpack.LemmaELeaf
import Sqpack.LemmaELeafX
import Sqpack.LemmaEPoly
import Sqpack.LemmaESplit
import Sqpack.LemmaESplitX
import Sqpack.LemmaEVert
import Sqpack.LemmaZ
import Sqpack.PolyMin
import Sqpack.RatU
import Sqpack.RowMinK
import Sqpack.S5LowerLeb
import Sqpack.SquareArea
import Sqpack.V7.Bridge
import Sqpack.V7.Data
