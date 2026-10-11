import Sqpack.Exact.Packs.N133.Data

namespace UnitSquarePacking.EC.N133

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 133, 100 ≤ i.val → i.val < 125 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N133
