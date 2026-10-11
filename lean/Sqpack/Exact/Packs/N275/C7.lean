import Sqpack.Exact.Packs.N275.Data

namespace UnitSquarePacking.EC.N275

set_option maxHeartbeats 0

theorem chunk_7 : ∀ i : Fin 275, 175 ≤ i.val → i.val < 200 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N275
