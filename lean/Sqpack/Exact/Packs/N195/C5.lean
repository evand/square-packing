import Sqpack.Exact.Packs.N195.Data

namespace UnitSquarePacking.EC.N195

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 195, 125 ≤ i.val → i.val < 150 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N195
