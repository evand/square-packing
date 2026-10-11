import Sqpack.Exact.Packs.N261.Data

namespace UnitSquarePacking.EC.N261

set_option maxHeartbeats 0

theorem chunk_10 : ∀ i : Fin 261, 250 ≤ i.val → i.val < 261 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N261
