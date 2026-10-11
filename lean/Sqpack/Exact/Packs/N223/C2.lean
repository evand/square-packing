import Sqpack.Exact.Packs.N223.Data

namespace UnitSquarePacking.EC.N223

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 223, 50 ≤ i.val → i.val < 75 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N223
