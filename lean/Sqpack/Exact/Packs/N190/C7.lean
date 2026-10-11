import Sqpack.Exact.Packs.N190.Data

namespace UnitSquarePacking.EC.N190

set_option maxHeartbeats 0

theorem chunk_7 : ∀ i : Fin 190, 175 ≤ i.val → i.val < 190 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N190
