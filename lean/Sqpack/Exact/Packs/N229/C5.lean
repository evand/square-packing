import Sqpack.Exact.Packs.N229.Data

namespace UnitSquarePacking.EC.N229

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 229, 125 ≤ i.val → i.val < 150 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N229
