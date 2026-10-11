import Sqpack.Exact.Packs.N282.Data

namespace UnitSquarePacking.EC.N282

set_option maxHeartbeats 0

theorem chunk_11 : ∀ i : Fin 282, 275 ≤ i.val → i.val < 282 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N282
