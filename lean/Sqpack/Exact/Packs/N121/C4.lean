import Sqpack.Exact.Packs.N121.Data

namespace UnitSquarePacking.EC.N121

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 121, 100 ≤ i.val → i.val < 121 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N121
