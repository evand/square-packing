import Sqpack.Exact.Packs.N53.Data

namespace UnitSquarePacking.EC.N53

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 53, 50 ≤ i.val → i.val < 53 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N53
