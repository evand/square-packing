import Sqpack.Exact.Packs.N318.Data

namespace UnitSquarePacking.EC.N318

set_option maxHeartbeats 0

theorem chunk_12 : ∀ i : Fin 318, 300 ≤ i.val → i.val < 318 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N318
