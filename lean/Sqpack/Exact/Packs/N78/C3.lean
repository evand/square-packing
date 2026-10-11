import Sqpack.Exact.Packs.N78.Data

namespace UnitSquarePacking.EC.N78

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 78, 75 ≤ i.val → i.val < 78 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N78
