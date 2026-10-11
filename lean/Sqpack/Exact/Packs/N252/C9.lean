import Sqpack.Exact.Packs.N252.Data

namespace UnitSquarePacking.EC.N252

set_option maxHeartbeats 0

theorem chunk_9 : ∀ i : Fin 252, 225 ≤ i.val → i.val < 250 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N252
