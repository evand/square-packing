import Sqpack.Exact.Packs.N200.Data

namespace UnitSquarePacking.EC.N200

set_option maxHeartbeats 0

theorem chunk_5 : ∀ i : Fin 200, 125 ≤ i.val → i.val < 150 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N200
