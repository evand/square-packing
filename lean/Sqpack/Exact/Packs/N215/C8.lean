import Sqpack.Exact.Packs.N215.Data

namespace UnitSquarePacking.EC.N215

set_option maxHeartbeats 0

theorem chunk_8 : ∀ i : Fin 215, 200 ≤ i.val → i.val < 215 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N215
