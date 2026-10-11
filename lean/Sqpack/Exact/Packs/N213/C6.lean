import Sqpack.Exact.Packs.N213.Data

namespace UnitSquarePacking.EC.N213

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 213, 150 ≤ i.val → i.val < 175 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N213
