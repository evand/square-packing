import Sqpack.Exact.Packs.N280.Data

namespace UnitSquarePacking.EC.N280

set_option maxHeartbeats 0

theorem chunk_10 : ∀ i : Fin 280, 250 ≤ i.val → i.val < 275 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N280
