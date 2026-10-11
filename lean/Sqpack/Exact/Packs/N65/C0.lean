import Sqpack.Exact.Packs.N65.Data

namespace UnitSquarePacking.EC.N65

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 65, 0 ≤ i.val → i.val < 25 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N65
