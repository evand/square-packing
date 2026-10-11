import Sqpack.Exact.Packs.N90.Data

namespace UnitSquarePacking.EC.N90

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 90, 75 ≤ i.val → i.val < 90 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N90
