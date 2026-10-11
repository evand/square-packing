import Sqpack.Exact.Packs.N167.Data

namespace UnitSquarePacking.EC.N167

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 167, 0 ≤ i.val → i.val < 25 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N167
