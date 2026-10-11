import Sqpack.Exact.Packs.N43.Data

namespace UnitSquarePacking.EC.N43

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 43, 25 ≤ i.val → i.val < 43 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N43
