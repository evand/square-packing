import Sqpack.Exact.Packs.N204.Data

namespace UnitSquarePacking.EC.N204

set_option maxHeartbeats 0

theorem chunk_8 : ∀ i : Fin 204, 200 ≤ i.val → i.val < 204 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N204
