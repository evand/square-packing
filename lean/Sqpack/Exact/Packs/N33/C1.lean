import Sqpack.Exact.Packs.N33.Data

namespace UnitSquarePacking.EC.N33

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 33, 25 ≤ i.val → i.val < 33 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N33
