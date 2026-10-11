import Sqpack.Exact.Packs.N66.Data

namespace UnitSquarePacking.EC.N66

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 66, 50 ≤ i.val → i.val < 66 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N66
