import Sqpack.Exact.Packs.N233.Data

namespace UnitSquarePacking.EC.N233

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 233, 100 ≤ i.val → i.val < 125 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N233
