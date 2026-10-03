# farfield-strong: the near-axis mass cap in the strongest family at `eps >= 10 deg`  (2026-09-21)

**Why.**  `search/T4_CYCLES.md` §2.1: the `T = 4` far-field dual (`k <= 3` near-axis squares) is a *thicket* (`mu >= 2` at
104/126, median 13 pair rows on 11 squares) — not a chain, not a cycle, not any human-shaped certificate.  So the far field
is the region where a **machine certificate** is the right tool, and the repo's certified cover-LP pipeline is that tool.
`search/ARCH_FARFIELD.md` refuted the cap at small `eps` (`<= 2 deg`) in the points-only family, and left `eps >= 5 deg`
**unsettled because the pool was too weak** (its own uncapped value was `11.905 < 12`).  The strongest family (corner
`k = 4` + polygons + regions + chord, `search/BENTZ.md` §7.2, `LEAF_CEILING.md` §5.2, uncapped exactly `12`) was never
tested with a cap.  If the cap closes at some `eps_ff`, the far field `{<= 3 squares of tilt < eps_ff}` becomes the
**first proved region** of the theorem, and `eps_ff` is the width the hole lemma must cover.

**Read first.**  `search/ARCH_FARFIELD.md` (all of it, especially §1 direction discipline, §5 the pool caveat, §6),
`tasks/arch-farfield/README.md`, `search/BENTZ.md` §7, `search/LEAF_CEILING.md` §5, `notes/review-2026-09-13b.md`
"Clarifications", `notes/status.md` table.  Instruments: `search/arch_farfield.py` (`cg`, `exact`, `table`),
`search/packing_dual.py`, `search/dual_exact.py`, `search/bentz.py`, `search/leaf_ceiling.py`.  Certified measures:
`search/pgonly_corner_exact.txt`, `search/pgonly_pure_exact.txt`, `search/cover4_exact_support.txt`.

**Do.**  `T = 4`, `n = 12`, closed semantics.
1. **Build a pool whose uncapped value is `>= 12`** in the strongest family, seeded from the certified supports above (these
   sit at `12.0000` / `12.17` / `12.27`).  Report the uncapped value on that pool first; nothing below is meaningful until it
   is `>= 12` (`ARCH_FARFIELD.md` §5).
2. **Cap `mu(tilt < eps) <= 3`** for `eps in {5, 10, 15, 20, 30 deg}` (and `cap = 2, 1` at `10 deg` as a sanity ladder).
   Column generation with tilted poses allowed everywhere; then the exact pipeline (`exact`, `dual_exact.py check --full`)
   on any cell that reads `< 12`.  Report packing-side (LP lower bound on the relaxation) and cover-side (certified upper
   bound) **separately**.  A cell certified `< 12` cover-side is a theorem; a cell `>= 12` packing-side by an exact measure
   is a refutation at that `eps`; anything else is "pool-limited" and must be labelled so.
3. **Where does the cap bite?**  Give the value as a function of `eps`; identify `eps_ff` = the smallest `eps` at which the
   capped value is certified `< 12`, if any `<= 30 deg`.  Show the heaviest poses of the capped optimum (where does the LP
   move the mass — `ARCH_FARFIELD.md` §4 found "corner square at `1.2 deg`").
4. **`T = 3` control** at the same `eps` grid, cap `2`, `n = 6`: where does *it* close?  (`Claim(3)` is a theorem, so the
   cap must close somewhere if the mechanism is any good; `ARCH_FARFIELD.md` §3.4 says past `10 deg` on a thin pool.)
5. **Filter.**  At any `eps` where the cap `>= 4` is claimed `< 12`, that is a bug (a genuine margin-`0` packing with `>= 4`
   near-axis squares exists at every `eps <= 12.8 deg`, `proof-architecture.md` §0a item 2) — but note the closed-semantics
   caveat in `ARCH_FARFIELD.md` §5 before calling it one.

**Deliverable.**  `search/FARFIELD_STRONG.md` (verdict up front, one table: `eps` × {packing-side, cover-side, status}),
new scripts only under `search/farfield_*.py` (import, never modify, existing ones), runs `runs/farfield_strong_*`,
exact certificates as `search/farfield_strong_*_exact.txt` with `check` commands.  Do not edit `TODO.md`,
`notes/status.md`, existing `search/*.md`; do not commit.

**Working style.**  `<= 8` threads (three other agents share the machine).  Anything over 10 minutes detached with
`setsid nohup`, never `pkill`, no `until … sleep` / `tail -f` watchers (a `pgrep -f <name>` loop matches itself).
Write the note section by section as cells finish; a partial table with honest labels beats a complete one late.
