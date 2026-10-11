import Sqpack.Exact.Packs.N262.Data

namespace UnitSquarePacking.EC.N262

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 262, 125 ≤ i.val → i.val < 150 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N262
