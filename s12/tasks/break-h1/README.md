# break-h1: try to break `H1` at `T = 4`, and measure the ladder depth at `T = 5`  (2026-09-21)

**Why.**  `search/T4_CYCLES.md` §7: `H1` (chain of four + transverse legs + **one** further main-direction pair row) holds at
`727/727` sampled `T = 4` optima.  But (i) the `n = 12` multistart is weak, every `delta*` is a lower bound, and an
under-found `delta*` makes the `H`-tests *easier* (`t3-chain.md` §2.1 "Reliability") — successes may be artefacts, failures
are real; (ii) the pattern so far is `H` (fails) → `H1` (holds so far), and the `eps^3` cell needed exactly three ladder
levels because its margin is `eps^3`.  If some cell has margin `eps^4`, `H1` fails there and the lemma is "a ladder of
unbounded depth" — i.e. the LP dual in disguise, not a finite human statement.  Whether the ladder depth needed grows with
`T` is also exactly the kind of `T`-dependence the filter (`notes/proof-architecture.md` §5) demands.

**Read first.**  `search/T4_CYCLES.md` §0, §1, §3 (the ladder and its weights `1/5, tan(eps)/5, tan^2(eps)/5`), §7;
`search/BANDCUT_K.md` §0, §1, §3 (how the chain-free cells and the `eps^3` configuration were built), §5;
`search/S6_LOCAL.md` §1–3 (structured starts; `n = 12` pinwheel dual `1/3, t/3, t^2/3`); `search/ARCH_TLEDGER.md` §2
(permutation hole sets; `c_T`).  Instruments: `search/t4_cycles.py` (`H`, `H1`, `CYC` restricted duals; `cell89`,
`fresh`, `dich`), `search/bandcut_k.py` / `bandcut_khunt.py` (the `(k, eps)` scan and the chain-free hunt),
`search/t3_chain.py` (`all_rows`, `dual_bound`, `h_support`), `search/s6local.py`, `search/s6exact.py`.  Data:
`runs/bandcut_k_T4*.jsonl`, `runs/t4_cycles_*`.  Note `search/bandcut_k.py chain_certificate` is buggy (no common-normal
check); do not use it.

**Do.**
1. **Hunt for `H1 > 0` at `T = 4`** with *structured* starts, not random ones: tiling minus permutation hole sets with the
   hole squares tilted `eps` (both signs, coherent and mixed), central `2x2` / `3x3` pinwheels at tilt `eps`, the `(4,1)`
   band stack perturbed, and the chain-free cells `k = 4..9` with the chain cut *through* tilted squares as well as among
   near-axis ones.  At each optimum compute `H`, `H1`, and `H2` (= H + **two** further main-direction rows, best over pairs)
   and `delta*`.  Maximise `H1 - delta*` … no: maximise `H1` subject to `delta* < 0`; report the worst cells.  Re-solve every
   candidate failure at row tolerance `1e-10` and with a stronger multistart before calling it a failure.
2. **Margin exponent by cell.**  For the thinnest cells found, fit `delta* = -gamma eps^p` over `eps in {0.5, 1, 2, 5 deg}`;
   report `p` and the ladder depth (number of weight levels in ratio `tan eps`) of the dual.  Is depth `= p` always?
3. **`T = 5`, `n = 20`, top of the hole `k = 15, 16` chain-free** (the `T = 5` analogue of the `eps^3` cell): the chain-free
   maximum's exponent `p` and the dual's ladder depth, at `eps = 1, 2, 5 deg`.  Build from the tiling minus a permutation
   hole set with the central block a pinwheel, as `BANDCUT_K.md` §3 did at `T = 4`; the `T = 5` LP is `40` centre
   coordinates, fine for HiGHS.  If `p = 4` there, say so loudly: that is "ladder depth grows with `T`".
4. **Filter.**  Every restricted-dual bound is a Farkas bound and must be `>= 0` at every margin-`0` point; check on the
   `(4,1)` band stack and the `T = 5` tiling-minus-hole-set zeros.  Anything below `-1e-8` is a bug in your rows.
5. **Two-line verdict:** does `H1` survive at `T = 4` under adversarial search (with hit rates)?  What is the ladder depth
   needed at `T = 3, 4, 5`?

**Deliverable.**  `search/BREAK_H1.md` (verdict up front), scripts `search/break_h1*.py` (import, never modify, existing
ones), runs `runs/break_h1_*`.  Do not edit `TODO.md`, `notes/status.md`, existing `search/*.md`; do not commit.

**Working style.**  `<= 6` threads (three other agents share the machine).  Anything over 10 minutes detached with
`setsid nohup`, never `pkill`, no `until … sleep` / `tail -f` watchers (a `pgrep -f <name>` loop matches itself).
Write the note section by section; label **[proved]** / **[measured]** / **[heuristic]** / **[guess]**; report hit rates.
