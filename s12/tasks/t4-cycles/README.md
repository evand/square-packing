# t4-cycles: which cycle certificates exist at `T = 4`, and what holds the chain-free cell shut  (2026-09-21)

**Why.**  The `T = 3` rehearsal (`notes/t3-chain.md`) found that the far-field certificate is a **cycle** — a closed loop of
pair separations winding once around, pinned by one wall row per wall, no rise term — exact at `45°` (`delta <= (93 - 66 sqrt2)/7`,
§3.2), and that chain H fails above `28.8°`.  The candidate lemma is now "chain with enough transverse links, **or** cycle"
(`TODO.md` top item; cycle lemma in `t3-chain.md` §5.3, `T = 4` guess §3.4: `delta <= 0` iff `k >= q sin(alpha/2)(T - u)`, i.e.
a cycle of **eight** at `45°`).  Nothing has been measured at `T = 4`.  Two questions: (1) is the chain-free `k = 8, 9` cell of
`search/BANDCUT_K.md` (margin `-0.1 eps^3`, the thinnest margin anywhere) held shut by a cycle — if so the chain route has been
chasing the wrong object there; (2) does the cycle lemma pass the filter at `T = 4`: it must **not** certify the `(4,1)` band
stack of `BANDCUT_SCAN.md` §5 (`runs/bandcut_scan_pT_T4.json`, a genuine margin-`0` point) below `0`, nor anything at `delta* = 0`.

**Read first.**  `notes/t3-chain.md` §0, §1.4, §3, §5.3; `search/BANDCUT_K.md` §0, §1, §3, §5; `search/S6_LOCAL.md` §2–3
(the `n = 12` pinwheel dual: staircase chain at `1/3` + transverse links at `t/3`, `t^2/3`).  Instruments: `search/t3_chain_shape.py`
(`support_graph`: dual support → graph, cyclomatic number, turns), `search/t3_chain.py`, `search/bandcut_k.py` (`cert` reads the
LP dual at a recorded optimum; `chain_certificate` is **buggy** — no common-normal check, `math.degrees` twice — do not use it),
`search/s6skel.py` (`Decider`).  Data: `runs/bandcut_k_T4*.jsonl` (the `(k, eps)` scan, chain-free optima included),
`runs/bandcut_scan_pT_T4.json`, `runs/bandcut_scan_bestbar_T4_p{3,4}.json`, `runs/s6local_n12*` if present.

**Do.**  `T = 4`, `n = 12`.
1. **Dual census.**  For every recorded `T = 4` optimum (all `(k, eps)` cells, chained and chain-free; the band stacks; uniform
   tilts `1°, 5°, 10°, 20°, 30°, 45°`; a few hundred random far-field angle vectors decided fresh), read the LP dual and classify
   its support: cyclomatic number, wall rows used (one per wall?), number of turns, length of the cycle, which squares are tilted
   and by how much.  Table by cell.  Same "chain / chain + rung / cycle / other" classes as `t3-chain.md` §3.1.
2. **The `eps^3` cell.**  At `k = 8, 9`, chain-free (`BANDCUT_K.md` §3, §5 says the dual is a `y`-chain of four through two
   tilted squares plus a transverse chain at a tenth of the weight): is that support a tree or does it contain a cycle?  Compute
   the best pure-cycle bound and the best H bound at those optima and compare with `delta*`.  Which object carries the `eps^3`?
3. **The cycle lemma at `T = 4`.**  Check `t3-chain.md` §5.3's formula against every cyclic dual found (weights `1` on links,
   `2 sin(alpha/2)` at turns; does it reproduce the LP value?).  At `45°`: does a cycle of eight exist among twelve squares, and
   is its bound `<= 0`?  If the true `45°` dual is not a single cycle, say what it is.
4. **Filter.**  On the `(4,1)` band stack (`delta* = 0`) and on the zero-margin points of `runs/bandcut_k_T4*`: every cycle and
   chain bound must be `>= 0` there.  Report the tightest.  Any bound `< 0` at a margin-`0` point is a bug in the bound —
   find it before anything else.
5. **Two-line verdict:** is the chain-or-cycle dichotomy true at every `T = 4` optimum sampled (`H <= 0` or cycle `<= 0` or
   `delta* = 0`)?  What fraction of the far field (`k <= 3`) is cyclic?

Label **[proved]** / **[measured]** / **[heuristic]** / **[guess]**.  Every number is a measurement on a feasible point unless
it is an exact dual evaluation.  Be adversarial: `t3-chain.md` found a repo bug and two wrong brief premises; expect more.

**Deliverable.**  `search/T4_CYCLES.md` (verdict up front), scripts `search/t4_cycles*.py`, runs `runs/t4_cycles_*`.
Do not edit `TODO.md`, `notes/status.md`, `notes/proof-architecture.md`, existing `search/*.py`; do not commit.

**Working style.**  Small steps; write the note section by section (append); keep derivations on disk.  Anything over
10 minutes: `setsid nohup`, `<= 16` threads, never `pkill`, no `until … sleep` / `tail -f` watchers (a `pgrep -f <name>` loop
matches itself).  The `n = 12` multistart is weak (`S6_LOCAL.md` §1): use recorded optima and structured starts, report hit rates.
