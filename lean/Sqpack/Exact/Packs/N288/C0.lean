import Sqpack.Exact.Packs.N288.Data

namespace UnitSquarePacking.EC.N288

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 288, 0 ≤ i.val → i.val < 25 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N288
