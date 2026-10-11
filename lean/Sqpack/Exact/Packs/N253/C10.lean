import Sqpack.Exact.Packs.N253.Data

namespace UnitSquarePacking.EC.N253

set_option maxHeartbeats 0

theorem chunk_10 : ∀ i : Fin 253, 250 ≤ i.val → i.val < 253 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N253
