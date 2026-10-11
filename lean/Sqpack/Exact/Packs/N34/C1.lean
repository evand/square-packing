import Sqpack.Exact.Packs.N34.Data

namespace UnitSquarePacking.EC.N34

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 34, 25 ≤ i.val → i.val < 34 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N34
