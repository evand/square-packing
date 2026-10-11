import Sqpack.Exact.Packs.N36.Data

namespace UnitSquarePacking.EC.N36

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 36, 25 ≤ i.val → i.val < 36 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N36
