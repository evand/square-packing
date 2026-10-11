import Sqpack.Exact.Packs.N48.Data

namespace UnitSquarePacking.EC.N48

set_option maxHeartbeats 0

theorem chunk_1 : ∀ i : Fin 48, 25 ≤ i.val → i.val < 48 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N48
