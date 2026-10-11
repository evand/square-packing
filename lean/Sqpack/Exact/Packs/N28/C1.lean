import Sqpack.Exact.Packs.N28.Data

namespace UnitSquarePacking.EC.N28

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 28, 25 ≤ i.val → i.val < 28 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N28
