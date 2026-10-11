import Sqpack.Exact.Packs.N215.Data

namespace UnitSquarePacking.EC.N215

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 215, 100 ≤ i.val → i.val < 125 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N215
