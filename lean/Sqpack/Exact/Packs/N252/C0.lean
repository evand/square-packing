import Sqpack.Exact.Packs.N252.Data

namespace UnitSquarePacking.EC.N252

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 252, 0 ≤ i.val → i.val < 25 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N252
