import Sqpack.Exact.Packs.N245.Data

namespace UnitSquarePacking.EC.N245

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 245, 75 ≤ i.val → i.val < 100 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N245
