import Sqpack.Exact.Packs.N146.Data

namespace UnitSquarePacking.EC.N146

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 146, 10 ≤ i.val → i.val < 20 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N146
