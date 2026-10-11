import Sqpack.Exact.Packs.N85.Data

namespace UnitSquarePacking.EC.N85

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 85, 75 ≤ i.val → i.val < 85 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N85
