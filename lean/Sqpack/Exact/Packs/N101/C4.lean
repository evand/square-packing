import Sqpack.Exact.Packs.N101.Data

namespace UnitSquarePacking.EC.N101

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 101, 100 ≤ i.val → i.val < 101 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N101
