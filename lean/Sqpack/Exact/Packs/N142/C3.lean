import Sqpack.Exact.Packs.N142.Data

namespace UnitSquarePacking.EC.N142

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 142, 75 ≤ i.val → i.val < 100 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N142
