import Sqpack.Exact.Packs.N189.Data

namespace UnitSquarePacking.EC.N189

set_option maxHeartbeats 0

theorem chunk_7 : ∀ i : Fin 189, 175 ≤ i.val → i.val < 189 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N189
