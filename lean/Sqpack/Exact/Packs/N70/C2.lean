import Sqpack.Exact.Packs.N70.Data

namespace UnitSquarePacking.EC.N70

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 70, 50 ≤ i.val → i.val < 70 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N70
