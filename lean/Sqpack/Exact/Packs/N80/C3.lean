import Sqpack.Exact.Packs.N80.Data

namespace UnitSquarePacking.EC.N80

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 80, 75 ≤ i.val → i.val < 80 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N80
