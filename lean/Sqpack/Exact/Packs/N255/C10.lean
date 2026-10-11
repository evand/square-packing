import Sqpack.Exact.Packs.N255.Data

namespace UnitSquarePacking.EC.N255

set_option maxHeartbeats 0

theorem chunk_10 : ∀ i : Fin 255, 250 ≤ i.val → i.val < 255 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N255
