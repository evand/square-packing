import Sqpack.Exact.Packs.N254.Data

namespace UnitSquarePacking.EC.N254

set_option maxHeartbeats 0

theorem chunk_9 : ∀ i : Fin 254, 225 ≤ i.val → i.val < 250 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N254
