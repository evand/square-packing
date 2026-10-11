import Sqpack.Exact.Packs.N212.Data

namespace UnitSquarePacking.EC.N212

set_option maxHeartbeats 0

theorem chunk_7 : ∀ i : Fin 212, 175 ≤ i.val → i.val < 200 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N212
