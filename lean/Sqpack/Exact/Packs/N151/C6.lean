import Sqpack.Exact.Packs.N151.Data

namespace UnitSquarePacking.EC.N151

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 151, 150 ≤ i.val → i.val < 151 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N151
