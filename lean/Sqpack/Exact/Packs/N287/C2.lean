import Sqpack.Exact.Packs.N287.Data

namespace UnitSquarePacking.EC.N287

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 287, 50 ≤ i.val → i.val < 75 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N287
