import Sqpack.Exact.Packs.N147.Data

namespace UnitSquarePacking.EC.N147

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 147, 125 ≤ i.val → i.val < 147 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N147
