import Sqpack.Exact.Packs.N56.Data

namespace UnitSquarePacking.EC.N56

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 56, 50 ≤ i.val → i.val < 56 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N56
