import Sqpack.Exact.Packs.N256.Data

namespace UnitSquarePacking.EC.N256

set_option maxHeartbeats 0

theorem chunk_10 : ∀ i : Fin 256, 250 ≤ i.val → i.val < 256 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N256
