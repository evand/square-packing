import Sqpack.Exact.Packs.N77.Data

namespace UnitSquarePacking.EC.N77

set_option maxHeartbeats 0

theorem chunk_3 : ∀ i : Fin 77, 75 ≤ i.val → i.val < 77 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N77
