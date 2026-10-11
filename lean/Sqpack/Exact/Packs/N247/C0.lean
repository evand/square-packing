import Sqpack.Exact.Packs.N247.Data

namespace UnitSquarePacking.EC.N247

set_option maxHeartbeats 0

theorem chunk_0 : ∀ i : Fin 247, 0 ≤ i.val → i.val < 25 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N247
