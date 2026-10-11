import Sqpack.Exact.Packs.N122.Data

namespace UnitSquarePacking.EC.N122

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 122, 75 ≤ i.val → i.val < 100 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N122
