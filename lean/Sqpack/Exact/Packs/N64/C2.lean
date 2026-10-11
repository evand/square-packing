import Sqpack.Exact.Packs.N64.Data

namespace UnitSquarePacking.EC.N64

set_option maxHeartbeats 0

theorem chunk_2 : ∀ i : Fin 64, 50 ≤ i.val → i.val < 64 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N64
