import Sqpack.Exact.Packs.N219.Data

namespace UnitSquarePacking.EC.N219

set_option maxHeartbeats 0

theorem chunk_8 : ∀ i : Fin 219, 200 ≤ i.val → i.val < 219 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N219
