import Sqpack.Exact.Packs.N299.Data

namespace UnitSquarePacking.EC.N299

set_option maxHeartbeats 0

theorem chunk_11 : ∀ i : Fin 299, 275 ≤ i.val → i.val < 299 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N299
