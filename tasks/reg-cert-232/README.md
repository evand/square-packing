# n = 232: exact certificate of the regularized drawing, 3 angle groups (2026-10-17)

**Result.** The display regularization of the n = 232 record (Ostashev's packing, refined by PingYou Ltd and Polkovskyi;
`site/www/data/p/best-232.json`) is an exact packing at the record side, with exactly three angle groups:
172 squares at 0°, 36 at exactly 45°, and 24 at θ₁ = 2 atan(t*) ≈ −29.41455° (≡ 60.58545° mod 90).  Here t* is the
generator of the degree-7 field of the record's exact form.

- S* = the root near 15.7641743174866 of `x^7 − 84x^6 + 2984x^5 − 58024x^4 + 665867x^3 − 4501160x^2 + 16563960x − 25556432`.
  This is the same polynomial and the same element as `search/exact/cur/data/n-232.minpoly.json.gz`.
- Field L = K(√2), degree 14, primitive element w = t + √2,
  g(w) = `w^14 + 10w^13 − 3w^12 − 232w^11 − 259w^10 + 1706w^9 + 2141w^8 − 5856w^7 − 6749w^6 + 11006w^5 + 10543w^4 − 10488w^3 − 7329w^2 + 3470w + 1783`
  (irreducible; `field.py`; t and √2 as polynomials in w in `field.json`).
- Angle parameters: u = 0, u = t, u = √2 − 1 = tan 22.5°.  The file has exactly 3 distinct u values.
- Centres: exact in L, coefficients ≤ 37 digits.  Every square is in contact (the regularizer's coincidence
  placement).  2049 contact conditions (the 464 centre coordinates are fully determined, rank 464) all hold exactly in L.
  The drawing's centres agree with the exact ones to 3e-69; the site's float drawing agrees to 1.6e-15.

**Checks.**
- `search/exact/verify_exact.py` (stdlib only, unchanged): **VALID**, 5.5 s.  Field degree 14, 1518 touching
  incidences proved as identities, 5282 strict by intervals, 26024 pairs by discs.
- Lean: `lean/Sqpack/Exact/Packs/N232R3.lean` (`search/exact/lean_cert.py --split 25`; `lake build Sqpack.Exact.Packs.N232R3`, Data 25 s + 10 chunks in parallel + glue 9 s): `UnitSquarePacking.EC.N232R3.packs : ∃ s, pS s = 0 ∧ Sa ≤ s ≤ Sb ∧ Packs 232 s`, axioms propext, Classical.choice, Quot.sound only.  This is also the first Lean proof of the current record side (the older `N232.lean` proves the superseded 8 + 11√2/2).
- Separate check (10-17, python-flint, not in `build.py`): `f(t(w)) ≡ 0`, `r² ≡ 2`, `u_C² + 2u_C − 1 ≡ 0` (so 45° exactly)
  and `u_B ≡ t` mod g, and S is the old certificate's S element.

**How it was built** (`build.py`, ~10 s).
1. `reg.regularize(232, style='coincide')` gives the 80-digit drawing (3 groups after `merge_angles` / `merge_determined`).
2. The tight contacts: walls, plus face-normal gaps < 1e-40.
3. The same contact equations with exact coefficients in L, solved by sparse elimination (python-flint `fmpq_poly`
   mod g) on an independent subset chosen by pivoted QR.  Then every contact row is checked exactly.
4. Written in `minpoly.py`'s format, with the field polynomial f set to g.

**Why it needed √2.**  The 84 squares that carry no load (including the whole 45° block) are not in the record's
contact system.  So the degree-7 certificate leaves them at nearby rational angles, which is the "wonky" look, also
seen in jlevy's regularized view.  Putting the block at exactly 45° needs √2, which is not in K.

Files: `field.py`/`field.json` (L), `build.py`, `n-232-reg3.minpoly.json.gz` (the certificate), `build.json`, `reg232.json`/`.log` (the
regularizer run).  Rerun: `search/exact/.venv/bin/python3 build.py` (from this directory), then
`python3 -I search/exact/verify_exact.py tasks/reg-cert-232/n-232-reg3.minpoly.json.gz` (from the repo root).
