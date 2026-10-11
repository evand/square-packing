import Sqpack.Exact.Packs.N247.Data

namespace UnitSquarePacking.EC.N247

set_option maxHeartbeats 0

theorem chunk_4 : ∀ i : Fin 247, 100 ≤ i.val → i.val < 125 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N247
