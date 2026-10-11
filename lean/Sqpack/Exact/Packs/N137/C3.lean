import Sqpack.Exact.Packs.N137.Data

namespace UnitSquarePacking.EC.N137

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 137, 75 ≤ i.val → i.val < 100 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N137
