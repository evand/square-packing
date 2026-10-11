import Sqpack.Exact.Packs.N19.Data

namespace UnitSquarePacking.EC.N19

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 19, 0 ≤ i.val → i.val < 19 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N19
