import Sqpack.Exact.Packs.N117.Data

namespace UnitSquarePacking.EC.N117

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 117, 100 ≤ i.val → i.val < 117 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N117
