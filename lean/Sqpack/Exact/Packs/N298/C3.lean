import Sqpack.Exact.Packs.N298.Data

namespace UnitSquarePacking.EC.N298

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 298, 75 ≤ i.val → i.val < 100 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N298
