import Sqpack.Exact.Packs.N191.Data

namespace UnitSquarePacking.EC.N191

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 191, 75 ≤ i.val → i.val < 100 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N191
