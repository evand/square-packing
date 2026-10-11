import Sqpack.Exact.Packs.N16.Data

namespace UnitSquarePacking.EC.N16

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 16, 0 ≤ i.val → i.val < 16 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N16
