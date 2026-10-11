import Sqpack.Exact.Packs.N217.Data

namespace UnitSquarePacking.EC.N217

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 217, 100 ≤ i.val → i.val < 125 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N217
