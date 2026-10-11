import Sqpack.Exact.Packs.N104.Data

namespace UnitSquarePacking.EC.N104

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 104, 100 ≤ i.val → i.val < 104 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N104
