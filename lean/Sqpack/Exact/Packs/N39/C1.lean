import Sqpack.Exact.Packs.N39.Data

namespace UnitSquarePacking.EC.N39

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 39, 25 ≤ i.val → i.val < 39 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N39
