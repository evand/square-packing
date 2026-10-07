# counting-ladder: write the chain-counting lemma, the ladder closed form, and read both off `T = 11`  (2026-09-21)

**Why.**  Two small, self-contained pieces of the architecture are unwritten.  (1) The corrected chain-counting lemma
(`notes/proof-architecture.md` §0a item 1; `notes/proof-architecture-review.md` §1) is the **only proved existence
statement** the route has; nobody has written its proof, and writing it will show exactly what "x-separated" means for
near-axis squares and whether anything survives when the chain passes through a tilted square (which is what holds the
chain-free cells shut, `search/BANDCUT_K.md` §0(3)).  (2) The `H1` / three-level ladder of `search/T4_CYCLES.md` §3 has no
closed form; `notes/t3-chain.md` §1.3 has the two-level one (Lemma H).  With a closed form in `(T, eps)` the existence
clause the ladder needs can be stated exactly, and the `T = 11` packing tells us which clause fails there.

**Read first.**  `notes/proof-architecture.md` §0a (items 1–8, 14), `notes/proof-architecture-review.md` §1;
`notes/t3-chain.md` §1 (row system, Lemma H proof, Lemma W, §1.5 mixed tilts); `search/T4_CYCLES.md` §3 (the ladder's
rows and weights at `eps = 1, 5, 10 deg`), §7; `search/S6_LOCAL.md` §3 (`n = 12` pinwheel dual `1/3, t/3, t^2/3`);
`search/T11_CHAINS.md` (per-link inequality, the five tight chains of Cantrell's packing, "rise control fails");
`notes/bandcut-cost.md` §1 (MT).  Instruments: `search/t3_chain_hform.py` (closed forms + Farkas residual check — extend
it), `search/t3_chain.py`, `search/t4_cycles.py` (`cell89` rebuilds the ladder point), `search/t11_chains.py`.  Data:
`runs/t4_cycles_cell89.txt`, `runs/bandcut_scan_exact110.json` (`T = 11`), `runs/bandcut_scan_exact132.json` (`T = 12`).

**Do.**
1. **Chain-counting lemma, proved.**  State and prove: if `k > (T-1)^2` of the squares have tilt `< eps` with
   `eps < 1/(T-1)` (rad), the far squares arbitrary, then some `T` of the near-axis squares form a wall-to-wall chain in the
   x-DAG or the y-DAG at the configuration's own margin (define the DAG exactly as `BANDCUT_K.md` §1.2).  Handle: why
   "x-separated" is acyclic / transitive enough for Mirsky at that `eps` (ordering along the global axes); what "wall-to-wall"
   needs at the two ends; touching (closed semantics).  Then say precisely what breaks for `k <= (T-1)^2` and for a chain
   with one tilted member — is there any counting statement that forces a chain *through* tilted squares?  Check the lemma
   against `T = 17` (Cleemann: `199 > 256`? no — say what it gives there) and `n = T^2 - 4`, `T = 3`.  Write it Lean-ready
   (hypotheses as a list; no "clearly").
2. **Ladder closed form.**  Generalise Lemma H to `H1(T, eps)` / the three-level ladder: crossbar of `T` at weight `w_0`,
   `a + c` transverse legs at `w_0 tan eps`, one further main-direction link on a leg's foot at `w_0 tan^2 eps`, with the
   wall rows that close each level, at a common tilt (then with the two crossbar members tilted and the rest axis-parallel,
   as in the `eps^3` cell).  Derive the weights by the same square-by-square bookkeeping as `t3-chain.md` §1.3, verify the
   Farkas residual numerically (extend `t3_chain_hform.py`), and give the value and its `eps`-expansion.  It must reproduce
   `runs/t4_cycles_cell89.txt` (`-0.0989 eps^3` at `T = 4`, `k = 9`) and the `n = 12` pinwheel `1/3, t/3, t^2/3` of
   `S6_LOCAL.md` §3.  State the existence clause the ladder needs as a single sentence about the separation graph.
3. **Filter at `T = 11`.**  On Cantrell's packing (`delta = +2.75e-4`, exact): enumerate the ladders available on each of the
   five tight chains of `T11_CHAINS.md` §2 (transverse legs at the chain's ends through the band squares — "unpinned" was
   asserted there, not measured); evaluate the closed form; say which clause of the ladder lemma fails at `T = 11` — the
   inequality, chain existence, leg existence, or the common-normal hypothesis.  Same on `T = 12`.  One table.
4. **Also**: `L*(T, 0) = T - 2` and the ladder's `eps`-orders as a function of `T` — does the number of levels needed to
   reach `<= 0` depend on `T` in the closed form, or only on which rows exist?  (Cross-reference `tasks/break-h1`, which
   measures the depth at `T = 5`.)

**Deliverable.**  `notes/counting-ladder.md` (lemma + proof; closed form + verification; the `T = 11/12` table; verdict),
`search/t3_chain_hform.py` extended in place is allowed **only** as `search/ladder_hform.py` (new file; import the old
one), runs `runs/ladder_*`.  Do not edit `TODO.md`, `notes/status.md`, existing `notes/*.md`, existing `search/*`;
do not commit.

**Working style.**  Pen-and-paper first; `<= 1` thread of compute.  Label **[proved]** / **[measured]** /
**[heuristic]** / **[guess]**.  Be adversarial with the repo's own statements (item 1 of §0a was already corrected once).
No watchers (`until … sleep`, `tail -f`); nothing needs to run long here.
