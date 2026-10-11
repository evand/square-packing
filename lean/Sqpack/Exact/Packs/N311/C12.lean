import Sqpack.Exact.Packs.N311.Data

namespace UnitSquarePacking.EC.N311

set_option maxHeartbeats 0

theorem chunk_12 : ∀ i : Fin 311, 300 ≤ i.val → i.val < 311 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N311
