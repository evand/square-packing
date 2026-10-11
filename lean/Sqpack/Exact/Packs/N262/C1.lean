import Sqpack.Exact.Packs.N262.Data

namespace UnitSquarePacking.EC.N262

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 262, 25 ≤ i.val → i.val < 50 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N262
