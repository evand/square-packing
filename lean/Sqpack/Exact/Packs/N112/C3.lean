import Sqpack.Exact.Packs.N112.Data

namespace UnitSquarePacking.EC.N112

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 112, 75 ≤ i.val → i.val < 100 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N112
