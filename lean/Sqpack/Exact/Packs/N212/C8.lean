import Sqpack.Exact.Packs.N212.Data

namespace UnitSquarePacking.EC.N212

set_option maxHeartbeats 0

theorem chunk_8 : ∀ i : Fin 212, 200 ≤ i.val → i.val < 212 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N212
