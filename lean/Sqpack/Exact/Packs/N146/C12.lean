import Sqpack.Exact.Packs.N146.Data

namespace UnitSquarePacking.EC.N146

set_option maxHeartbeats 0

theorem chunk_12 : ∀ i : Fin 146, 120 ≤ i.val → i.val < 130 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N146
