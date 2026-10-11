import Sqpack.Exact.Packs.N248.Data

namespace UnitSquarePacking.EC.N248

set_option maxHeartbeats 0

theorem chunk_9 : ∀ i : Fin 248, 225 ≤ i.val → i.val < 248 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N248
