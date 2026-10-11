import Sqpack.Exact.Packs.N14.Data

namespace UnitSquarePacking.EC.N14

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 14, 0 ≤ i.val → i.val < 14 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N14
