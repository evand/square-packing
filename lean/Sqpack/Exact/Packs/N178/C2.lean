import Sqpack.Exact.Packs.N178.Data

namespace UnitSquarePacking.EC.N178

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 178, 50 ≤ i.val → i.val < 75 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N178
