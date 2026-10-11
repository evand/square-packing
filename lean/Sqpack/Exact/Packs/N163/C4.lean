import Sqpack.Exact.Packs.N163.Data

namespace UnitSquarePacking.EC.N163

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 163, 100 ≤ i.val → i.val < 125 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N163
