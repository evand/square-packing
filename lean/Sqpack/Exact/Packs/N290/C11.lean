import Sqpack.Exact.Packs.N290.Data

namespace UnitSquarePacking.EC.N290

set_option maxHeartbeats 0

theorem chunk_11 : ∀ i : Fin 290, 275 ≤ i.val → i.val < 290 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N290
