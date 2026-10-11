import Sqpack.Exact.Packs.N26.Data

namespace UnitSquarePacking.EC.N26

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 26, 25 ≤ i.val → i.val < 26 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N26
