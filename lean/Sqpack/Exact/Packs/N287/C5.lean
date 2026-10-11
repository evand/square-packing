import Sqpack.Exact.Packs.N287.Data

namespace UnitSquarePacking.EC.N287

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 287, 125 ≤ i.val → i.val < 150 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N287
