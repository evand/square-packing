import Sqpack.Exact.Packs.N58.Data

namespace UnitSquarePacking.EC.N58

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 58, 50 ≤ i.val → i.val < 58 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N58
