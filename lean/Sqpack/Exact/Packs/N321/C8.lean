import Sqpack.Exact.Packs.N321.Data

namespace UnitSquarePacking.EC.N321

set_option maxHeartbeats 0

theorem chunk_8 : ∀ i : Fin 321, 200 ≤ i.val → i.val < 225 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N321
