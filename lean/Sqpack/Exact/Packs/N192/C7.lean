import Sqpack.Exact.Packs.N192.Data

namespace UnitSquarePacking.EC.N192

set_option maxHeartbeats 0

theorem chunk_7 : ∀ i : Fin 192, 175 ≤ i.val → i.val < 192 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N192
