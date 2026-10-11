import Sqpack.Exact.Packs.N187.Data

namespace UnitSquarePacking.EC.N187

set_option maxHeartbeats 0

theorem chunk_7 : ∀ i : Fin 187, 175 ≤ i.val → i.val < 187 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N187
