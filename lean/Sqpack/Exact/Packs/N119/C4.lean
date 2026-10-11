import Sqpack.Exact.Packs.N119.Data

namespace UnitSquarePacking.EC.N119

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 119, 100 ≤ i.val → i.val < 119 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N119
