import Sqpack.Exact.Packs.N246.Data

namespace UnitSquarePacking.EC.N246

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 246, 0 ≤ i.val → i.val < 25 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N246
