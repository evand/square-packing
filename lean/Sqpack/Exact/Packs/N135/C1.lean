import Sqpack.Exact.Packs.N135.Data

namespace UnitSquarePacking.EC.N135

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 135, 25 ≤ i.val → i.val < 50 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N135
