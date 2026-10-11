import Sqpack.Exact.Packs.N135.Data

namespace UnitSquarePacking.EC.N135

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 135, 75 ≤ i.val → i.val < 100 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N135
