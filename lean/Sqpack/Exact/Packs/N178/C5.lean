import Sqpack.Exact.Packs.N178.Data

namespace UnitSquarePacking.EC.N178

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 178, 125 ≤ i.val → i.val < 150 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N178
