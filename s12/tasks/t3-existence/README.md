# t3-existence: prove the existence clause at `T = 3`, or show the angle-space route is dead there  (2026-09-21)

**Why this task decides the route.**  `s(12) = 4` has been reduced (`notes/proof-architecture.md`, `notes/t3-chain.md`,
`search/T4_CYCLES.md`) to: for every angle vector `theta`, the centre-LP optimum carries a Farkas certificate of a known
shape (chain of `T` + transverse legs + at most one further main-direction row, "`H1`"; or a cycle).  Every *inequality* in
that reduction is proved and `T`-uniform.  Every *existence* clause — "a configuration with `delta >= 0` contains such a
support" — is unproved, at every `T`, and `notes/t3-chain.md` §5.4 concludes that at `T = 3` the existence clauses are "no
easier than Kearney–Shiu".  Evan's rule: **we will not chase a route for `s(12)` that cannot handle `s(6)`.**  So: make it
work at `T = 3`, or break it.

**Read first.**  `notes/t3-chain.md` in full (§1 the row system and Lemma H/W/Z, §2.5 the rung repair, §3 the cycle,
§4 Kearney–Shiu in this language, §5 the skeleton and the open clauses (O1)–(O3), §6 the filter); `search/T4_CYCLES.md` §0,
§7 (the `H1` statement); `search/S6_SKELETON.md` §3.1 (the disjunctive LP), §5 (pinwheel); `search/S6_LOCAL.md` §3, §5;
`notes/proof-architecture.md` §0a items 1, 5 (chain counting, `k > (T-1)^2`); `notes/proof-anatomy.md` §5.2 (K–S).
Instruments: `search/t3_chain*.py` (`all_rows`, `dual_bound`, `h_support`, `best_h`, `delta_star`; the repair and
existence scripts), `search/s6skel.py` (`Decider`), `search/s6local.py`, `search/chains.py`.  Data:
`runs/t3_chain_scan1.jsonl` (720 optima with duals), `runs/t3_chain_repair1.txt`, `runs/t3_chain_exist1.txt`.

**The statement to prove** (`T = 3`, `n = 6`, closed semantics as `notes/s13-casefree.md` §1).
> **(E3)**  Let six closed unit squares be packed in `[0,3]^2` with margin `delta >= 0` (every wall gap and every
> separating gap `>= delta`; `delta = 0` allowed, touching allowed).  Then the configuration's separation graph at level
> `delta` contains one of: (a) an H (`t3-chain.md` §1.2) with `L >= L*(3, t)` transverse links at a common tilt `t` and
> common link normals; (b) an H plus one further main-direction pair row whose Farkas combination is `<= 0`; (c) a cycle
> satisfying the hypotheses of `t3-chain.md` §5.3 with `k >= q sin(alpha/2)(3 - u)`.
Since each of (a)–(c) gives `delta <= 0` by a proved identity, (E3) with `delta > 0` is a contradiction, i.e. (E3) ⇒ `s(6) = 3`.
You may replace (a)–(c) by any other *proved-identity* support shape if that makes existence provable; say so.

**Do.**
1. **Carve off what is already provable.**  (i) `>= 5` near-axis squares (`tilt < 1/2` rad): chain counting gives a chain of
   three; does Lemma Z / Lemma H then close it, *including the tilts of the chain squares*?  Write that sub-case out fully or
   say exactly where it fails.  (ii) Uniform tilt: Lemma W reduces to "some certificate has `omega <= 1/(5 - u)`"; is that
   provable from `delta >= 0` alone?  (iii) The pinwheel stratum.  Record each sub-case as proved / open with the obstruction.
2. **Attack the remaining region honestly.**  Try at least three genuinely different mechanisms for existence: pigeonhole /
   strip counting through tilted squares (the hole's chains run *through* tilted squares, `BANDCUT_K.md` §0(3));
   K–S-style unavoidable points as the source of forced incidences (their dictionary is `t3-chain.md` §4.2 — can an
   unavoidable set *produce* the chain rather than replace it?); a topological / Menger-type argument (no wall-to-wall chain
   ⇒ a transverse "cut" of squares, and what a cut costs).  For each: either a proof, or the explicit configuration family
   that defeats it (verify with `Decider` that the family really has `delta >= 0` or `delta*` close to `0`).
3. **Test every candidate lemma at `delta >= 0` and not only at LP optima**: sample configurations of margin exactly `0`
   and slightly negative, including non-optimal ones (perturb optima; the `Z3` stratum of `t3_chain.md` §2.6 is
   under-sampled).  A candidate must also pass the `s(5) = 3` filter (`t3-chain.md` §6.3) and be silent at `T = 2`.
4. **Verdict, two lines, no hedging:** either "(E3) is proved [here]" with the proof written out, or "(E3) is open; the
   irreducible obstruction is X; the route at `T = 4` inherits it because Y" — and, if you can, your best judgement whether
   any *finite* support shape can have a provable existence clause, given that the `T = 4` far-field dual is a thicket
   (`T4_CYCLES.md` §2.1).

Label **[proved]** / **[measured]** / **[heuristic]** / **[guess]**.  Be adversarial with the repo: `t3-chain.md` and
`T4_CYCLES.md` each found errors in the brief and in earlier notes; expect more.  Partial proofs are welcome if their
hypotheses are stated exactly.  Do not spend the task re-measuring what `t3-chain.md` measured.

**Deliverable.**  `notes/t3-existence.md` (verdict up front; the proofs in full; a table of sub-regions proved / open),
scripts `search/t3_exist_*.py`, runs `runs/t3_exist_*`.  Do not edit `TODO.md`, `notes/status.md`, existing `notes/*.md`,
existing `search/*.py`; do not commit.

**Working style.**  Small steps; write the note section by section (append); keep derivations on disk.  Compute is
light here: `<= 2` threads.  Anything over 10 minutes: `setsid nohup`, never `pkill`, no `until … sleep` / `tail -f`
watchers (a `pgrep -f <name>` loop matches itself).  Stop when the verdict is clear; do not pad.
