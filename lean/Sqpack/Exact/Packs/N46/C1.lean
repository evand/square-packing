import Sqpack.Exact.Packs.N46.Data

namespace UnitSquarePacking.EC.N46

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 46, 25 ≤ i.val → i.val < 46 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N46
