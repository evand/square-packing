import Sqpack.Exact.Packs.N287.Data

namespace UnitSquarePacking.EC.N287

set_option maxHeartbeats 0

theorem chunk_7 : ∀ i : Fin 287, 175 ≤ i.val → i.val < 200 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N287
