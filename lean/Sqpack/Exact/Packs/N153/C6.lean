import Sqpack.Exact.Packs.N153.Data

namespace UnitSquarePacking.EC.N153

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 153, 150 ≤ i.val → i.val < 153 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N153
