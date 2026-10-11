import Sqpack.Exact.Packs.N198.Data

namespace UnitSquarePacking.EC.N198

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 198, 75 ≤ i.val → i.val < 100 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N198
