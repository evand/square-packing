import Sqpack.Exact.Packs.N212.Data

namespace UnitSquarePacking.EC.N212

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 212, 75 ≤ i.val → i.val < 100 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N212
