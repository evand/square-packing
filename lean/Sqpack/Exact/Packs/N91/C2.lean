import Sqpack.Exact.Packs.N91.Data

namespace UnitSquarePacking.EC.N91

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 91, 50 ≤ i.val → i.val < 75 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N91
