# Move-kind bandit: timescale and per-n learning (2026-10-14, logs only)

Script: `analysis/bandit_timescale.py` (from `search/packer`: `python3 analysis/bandit_timescale.py --procs 8 [--json out.json]`,
about 5 min on 8 procs; `--no-scan` skips A and takes about 40 s).  It reads `runs/*/proposals.jsonl`, `archive.jsonl`, `manifest.json` and the store
(read-only).  No CPU was spent on search.

Data: 2,236 explorer rounds, 1.96 M proposals, 460 CPU-h.  Chain: night1014/nh (139 n, 593 rounds, 0.96 M proposals),
night1014/nb, hunt1012 h1-h5b, timetravel plain/sparse.  Single long runs at n = 110: census v2-v4 (32 runs, 60-240 min).
Pre-store, legacy reward only: hunt1, hunt1b, hunt2, hunt3/chain, lineage.  `sec` (the bandit's cost) equals measured CPU to 1-3 %.

## Numbers first

| question | answer |
|---|---|
| Do kind effects drift within a round or run? | **No drift detected.** The kind information in past rewards is largest with **no decay and carry 1** at every n bucket (A).  Between rounds at the same n, the spread of kind ratios beyond Poisson noise is 0.00 (DerSimonian-Laird). |
| What does the current setting (wall 600 s, x0.5 per round) cost? | It keeps 60-95 % of the available kind information: n < 70: 15.2 vs 15.8 nats/1000 rewards; 70-129: 9.9 vs 11.6; value reward at 70-129: 7.6 vs 12.9; census L: 8.6 vs 17.7. |
| How big are kind effects? | Small.  SD of the pooled log rate ratio between kinds is **0.13-0.14** (L/G) and 0.23 (V).  Best kinds (cross2, aswap, bigkick) are about 1.1x average; worst (mirror, rowslide) 0.75-0.8x. |
| Are they the same across n? | Mostly.  Between-n SD of a kind's log ratio is 0.16 (L), 0.14 (G) and 0 (V, too few rewards to resolve).  There are trends in n: cross2 falls by 0.42 per unit of ln n, bigkick by 0.23, kick rises by 0.25 and kicksym by 0.17. |
| What can a per-n bandit learn in its rounds? | night_nh: median 2 rounds per n (mean 4.3), about 230 L rewards per n but only **14 value rewards**.  Sampling SE of a per-(n, kind) log ratio after all rounds: **0.21 (L), 0.63 (V)**, larger than the effects it is trying to find (0.14).  Optimal per-n weight `tau^2 / (tau^2 + SE^2)` is about 0.37 for L and about 0 for V. |
| Pooling (C), predicting round j's kind split | Nats per 1000 rewards, L / G / V: per-n only 9.5 / 7.9 / 7.9; pooled other n of the run 9.1 / 7.5 / 8.1; **other runs only (fixed weights from history) 9.2 / 7.9 / 9.8**; hierarchical (pooled prior + per-n data, m = 10 pseudo-rewards) **12.1 / 10.2 / 10.4** (V best at m = 3, near-n prior: 12.2). |
| n >= 200 | Per-n kind statistics *hurt* prediction at every half-life (L -18, G -17 nats/1000 rewards): rewards there are too few and too clustered on parents. |
| Value of kind choice (E, held out) | Policy 1/3 uniform + 2/3 on the estimated best kind, vs uniform, n-balanced (rounds equal): **x1.17 (per-n), x1.17 (pooled), x1.17 (other runs), x1.20 (hierarchical)** for L; x1.11-1.18 for V.  The logged bandit's actual allocation got x1.07 vs uniform.  The in-sample oracle reaches x1.8, which is mostly noise. |
| Legacy vs global vs value ranking (D) | Chain: L and G agree almost exactly (Spearman 0.99); V agrees at 0.85.  Census at 110: L vs G 0.85, **L vs V 0.38**.  kicksym: chain L 0.86 / G 0.88 / V 1.02; census L 1.20 / **G 1.34** / V 1.04. |
| Base-rate drift (B) | Legacy reward per CPU-s falls within a chain round (Q4/Q1 = 0.80-0.83; census 0.31; lineage 0.25).  Global reward is flat within chain rounds (0.95-1.0).  The legacy decline is the run re-finding basins already found in earlier rounds.  It hits all kinds about equally, so it does not change which kind to pick. |

## Recommendation

1. **Unit: proposals or CPU-s, not wall.**  A wall half-life means 30x fewer samples at n > 200 (250 proposals per 4-min round
   at 1.9 s each) than at n < 70 (4,100 at 0.1 s).  Within one n, proposals and CPU-s are equivalent.  If a decay is kept at
   all, use **h >= 3,000 proposals (or >= 1,000 CPU-s)**: the gain curve is flat from there to infinity at every n bucket.
   Half-lives below 300 proposals lose most of the information, and at large n the current 600 s wall is one of them.
2. **Half-life: effectively infinite; carry 1, not x0.5.**  Nothing in the logs supports forgetting on the scale of rounds.
   Carry 0 (reset every round) is clearly worse (n < 70: 9.8 vs 15.4).
3. **Pool across n.**  Kind effects are mostly global, with mild log-n trends.  Use a hierarchical estimate: a prior from
   all n, or from nearby n and history, plus per-n counts at about 10 pseudo-rewards of prior weight.  With the value reward,
   per-n data adds almost nothing (14 rewards per n).
4. **`--adapt` in chain: drop the per-n bandit.**  Replace it with fixed kind weights from the pooled table (refit per
   n-band, e.g. favour cross2 / bigkick at small n and kick / kicksym at large n), or with one pooled bandit with no decay.
   The ceiling over uniform is about x1.2 (assuming no depletion), and the current bandit realises x1.07.  This is a
   second-order knob; the reward definition matters more.
5. **Floor: keep about 1/3 uniform (0.03 per kind).**  With effects of +-15 % and no drift, it costs at most about 7 % (1/3
   x (best - average)).  It also guards against the depletion that this observational data cannot measure.

## Methods (A-E as in the script docstring)

- Rewards recomputed for every logged proposal:
  - **L** (legacy): new archive entry below k = ceil(sqrt n).
  - **G** (global): L, and not within 2e-9 of any store packing at n added before the run's start, nor of any archive
    entry from an earlier round of the same chain run at n.
  - **V** (value): new and globally novel (no k cut), weighted by min(1, 1e-3 / gap), where gap is measured to the best
    known before the run (store plus earlier rounds).  The robustness factor of `--reward value` is omitted.
- **A** (online predictive log-likelihood):
  - Base expectation E_i: an all-kind decayed rate (half-life 200 proposals, x0.5 at round starts) times sec_i.
  - Kind multiplier: rho_k = (R_k + 2) / (E_k + 2), decayed with half-life h in proposals / CPU-s / round wall time, with
    carry {1, 0.5, 0} at chain round boundaries.
  - Gain: sum over proposals of y log rho - (rho - 1) E, in nats per 1000 rewards.  The argmax over h is the timescale
    on which kind information goes stale.
- **C**: multinomial log-likelihood of round j's reward split among kinds, given round j's own CPU per kind, under
  estimates from data before round j.  Rates are indirectly standardised (expected = kind CPU x round reward rate).
- **D**: pooled indirectly standardised ratios.  Between-n and between-round heterogeneity by DerSimonian-Laird on the
  log ratios.
- **E**: expected reward of a policy's kind weights on round j's empirical per-kind rates (held out), divided by uniform.

## Caveats

- **Adaptive logging policy.**  The kinds were chosen by the per-n bandit (floor about 1/3 uniform) and parents by elite
  share 0.3 / frontier 0.7.  Per-kind rates are therefore measured at policy-dependent allocations.
- **Depletion is not identified.**  A kind used more finds its own basins first, so its marginal rate drops.  The E ratios
  assume rates do not depend on allocation and are optimistic for concentrated policies.  The realized x1.07 is in-sample
  for the round.
- **Global novelty is weak outside n = 110.**  The store holds few explorer archives, so in chain G is about 86 % of L.  G
  then means "new to the chain run", not "new to the world".  At 110 (150 k store packings) G is 8 % of L, and there G
  ranks kicksym first, so the 10-08 ex3 finding (legacy over-rewards kicksym) does not reproduce with this reference.
  Value grading does demote kicksym at 110.
- **No certified classes.**  No `exact.json` exists in these runs.  Novelty is side-based (2e-9); certification could merge
  or split basins.
- **Units and coverage.**  Wall time is round-relative (rounds concatenated).  Pre-store logs carry L only.  The values
  rewards at n >= 200 are a few hundred events: those rows are noisy.
