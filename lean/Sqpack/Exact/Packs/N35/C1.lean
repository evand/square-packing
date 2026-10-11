import Sqpack.Exact.Packs.N35.Data

namespace UnitSquarePacking.EC.N35

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 35, 25 ≤ i.val → i.val < 35 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N35
