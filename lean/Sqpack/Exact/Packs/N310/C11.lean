import Sqpack.Exact.Packs.N310.Data

namespace UnitSquarePacking.EC.N310

set_option maxHeartbeats 0

theorem chunk_11 : ∀ i : Fin 310, 275 ≤ i.val → i.val < 300 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N310
