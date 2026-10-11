import Sqpack.Exact.Packs.N157.Data

namespace UnitSquarePacking.EC.N157

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 157, 150 ≤ i.val → i.val < 157 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N157
