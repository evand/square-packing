import Sqpack.Exact.Packs.N9.Data

namespace UnitSquarePacking.EC.N9

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 9, 0 ≤ i.val → i.val < 9 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N9
