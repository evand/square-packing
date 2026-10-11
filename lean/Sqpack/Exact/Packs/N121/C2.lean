import Sqpack.Exact.Packs.N121.Data

namespace UnitSquarePacking.EC.N121

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 121, 50 ≤ i.val → i.val < 75 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N121
