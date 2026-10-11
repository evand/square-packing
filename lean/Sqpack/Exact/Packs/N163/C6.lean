import Sqpack.Exact.Packs.N163.Data

namespace UnitSquarePacking.EC.N163

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 163, 150 ≤ i.val → i.val < 163 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N163
