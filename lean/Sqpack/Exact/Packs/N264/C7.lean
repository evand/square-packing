import Sqpack.Exact.Packs.N264.Data

namespace UnitSquarePacking.EC.N264

set_option maxHeartbeats 0

theorem chunk_7 : ∀ i : Fin 264, 175 ≤ i.val → i.val < 200 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N264
