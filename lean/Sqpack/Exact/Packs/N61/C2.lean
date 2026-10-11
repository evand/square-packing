import Sqpack.Exact.Packs.N61.Data

namespace UnitSquarePacking.EC.N61

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 61, 50 ≤ i.val → i.val < 61 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N61
