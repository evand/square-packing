import Sqpack.Exact.Packs.N323.Data

namespace UnitSquarePacking.EC.N323

set_option maxHeartbeats 0

theorem chunk_10 : ∀ i : Fin 323, 250 ≤ i.val → i.val < 275 → boxOK cert i = true ∧ rowOK cert i = true := by
  decide +kernel

end UnitSquarePacking.EC.N323
