import Sqpack.Exact.Packs.N308.Data

namespace UnitSquarePacking.EC.N308

set_option maxHeartbeats 0

theorem chunk_8 : ∀ i : Fin 308, 200 ≤ i.val → i.val < 225 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N308
