import Sqpack.Exact.Packs.N42.Data

namespace UnitSquarePacking.EC.N42

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 42, 25 ≤ i.val → i.val < 42 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N42
