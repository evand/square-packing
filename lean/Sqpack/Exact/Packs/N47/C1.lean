import Sqpack.Exact.Packs.N47.Data

namespace UnitSquarePacking.EC.N47

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 47, 25 ≤ i.val → i.val < 47 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N47
