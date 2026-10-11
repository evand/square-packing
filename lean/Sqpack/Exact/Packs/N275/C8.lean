import Sqpack.Exact.Packs.N275.Data

namespace UnitSquarePacking.EC.N275

set_option maxHeartbeats 0

theorem chunk_8 : ∀ i : Fin 275, 200 ≤ i.val → i.val < 225 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N275
