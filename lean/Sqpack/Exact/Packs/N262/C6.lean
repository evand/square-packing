import Sqpack.Exact.Packs.N262.Data

namespace UnitSquarePacking.EC.N262

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 262, 150 ≤ i.val → i.val < 175 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N262
