import Sqpack.Exact.Packs.N285.Data

namespace UnitSquarePacking.EC.N285

set_option maxHeartbeats 0

theorem chunk_11 : ∀ i : Fin 285, 275 ≤ i.val → i.val < 285 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N285
