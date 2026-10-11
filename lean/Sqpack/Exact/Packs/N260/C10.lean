import Sqpack.Exact.Packs.N260.Data

namespace UnitSquarePacking.EC.N260

set_option maxHeartbeats 0

theorem chunk_10 : ∀ i : Fin 260, 250 ≤ i.val → i.val < 260 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N260
