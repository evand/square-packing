import Sqpack.Exact.Packs.N73.Data

namespace UnitSquarePacking.EC.N73

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 73, 0 ≤ i.val → i.val < 25 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N73
