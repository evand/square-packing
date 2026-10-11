import Sqpack.Exact.Packs.N189.Data

namespace UnitSquarePacking.EC.N189

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 189, 25 ≤ i.val → i.val < 50 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N189
