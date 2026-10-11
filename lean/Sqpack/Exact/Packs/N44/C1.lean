import Sqpack.Exact.Packs.N44.Data

namespace UnitSquarePacking.EC.N44

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 44, 25 ≤ i.val → i.val < 44 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N44
