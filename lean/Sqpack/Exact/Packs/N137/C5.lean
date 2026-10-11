import Sqpack.Exact.Packs.N137.Data

namespace UnitSquarePacking.EC.N137

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 137, 125 ≤ i.val → i.val < 137 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N137
