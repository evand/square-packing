import Sqpack.Exact.Packs.N73.Data

namespace UnitSquarePacking.EC.N73

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 73, 50 ≤ i.val → i.val < 73 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N73
