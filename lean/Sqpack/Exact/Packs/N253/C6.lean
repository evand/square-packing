import Sqpack.Exact.Packs.N253.Data

namespace UnitSquarePacking.EC.N253

set_option maxHeartbeats 0

theorem chunk_6 : ∀ i : Fin 253, 150 ≤ i.val → i.val < 175 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N253
