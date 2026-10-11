import Sqpack.Exact.Packs.N287.Data

namespace UnitSquarePacking.EC.N287

set_option maxHeartbeats 0

theorem chunk_11 : ∀ i : Fin 287, 275 ≤ i.val → i.val < 287 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N287
