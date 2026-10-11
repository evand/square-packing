import Sqpack.Exact.Packs.N44.Data

namespace UnitSquarePacking.EC.N44

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 44, 0 ≤ i.val → i.val < 25 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N44
