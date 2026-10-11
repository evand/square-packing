import Sqpack.Exact.Packs.N122.Data

namespace UnitSquarePacking.EC.N122

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 122, 50 ≤ i.val → i.val < 75 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N122
