import Sqpack.Exact.Packs.N310.Data

namespace UnitSquarePacking.EC.N310

set_option maxHeartbeats 0

theorem chunk_12 : ∀ i : Fin 310, 300 ≤ i.val → i.val < 310 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N310
