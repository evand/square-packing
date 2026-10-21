# Allocator: standing generators as a bandit (design, 10-20)

Evan 10-20: "All generators should be standing generators.  If that has downsides, let's solve them in a disciplined and
general way.  This is a bandit arm / bayesian inference / value of information math problem, not a judgment call."
Also: cluster weighting in start selection, "moderately", with weighting diverse enough to treat it as a math problem;
online or retrospective evaluation both fine.  Store policy: inclusive intake, conservative dedupe (done 10-20:
`pk/policy.yaml`, `pk/sink.py`).

Status: design for review.  Nothing below is built except where marked (built).

## 1. Objects

* **Arm** a = (g, n, θ): generator g (explorer-from-pool, dgen band sweep, hborder, lcut / L-cut, bulbgen, mondrian
  compose, EV2, junction refill, near-miss finisher, ...), target n, parameters θ (start-selection weights α, τ below;
  generator knobs).  Hierarchy g ⊃ (g, n-band) ⊃ (g, n) ⊃ θ for pooling.
* **Job**: one launch of an arm with a CPU budget c (measured CPU, `cpu_limit`, not wall; TODO "CPU budgets" is a
  prerequisite).  Its outputs reach the store through the sink, attributed by run dir.
* **Pool** P_n: the store's distinct packings at n (inclusive intake), clustered (pk/cluster.py).

## 2. Value: one potential over the pool

Records are rare (SW + SW2: 3 in 146 CPU-h), and both SW2 records came from stored non-record packings two
generations later.  So reward must credit pool additions, not only records, and must not double count.

Define the pool potential Φ(P) = Σ_n V_n · E[records at n from the standing policy expanding P over a horizon H].
A job's reward is its realised value, records plus the change in potential:

    R(job) = Σ_{records r} V_n(r)  +  [Φ(P ∪ Y) − Φ(P)]          (Y = its new distinct packings)

so a generator that only feeds the pool is paid at creation (the expected value of what it added), and an expander that
turns a pool packing into a record is paid the record minus nothing it did not cause (the parent's potential was
already in Φ(P); records reset the frontier and Φ is recomputed).  This is potential-based shaping: the total over
time telescopes to records, so no allocation can game it.

Approximation (submodular, cluster-aware; this is where diversity enters the value):

    Φ_n(P) ≈ V_n · Σ_{clusters C} ρ_C(P),   ρ_C = 1 − Π_{p ∈ C, distinct minima} (1 − h(p))^{w_dup(p)}

h(p) = P(a record descends from p within H | φ(p)); φ = log gap to the frontier, screened / polished / certified,
times expanded (CPU already spent below p), children, n band, ntilted, generator of origin.  A second copy of a
minimum adds 0 (it is the same row); a new minimum in a crowded cluster adds less than one in an empty cluster.
h is the explorer's existing parent value model (Poisson rate with saturation terms, PACKER "TT3 round 3"; PF49:
2.81x [2.22, 3.45] over G) lifted from within-run to the store: label = records in p's lineage subtree, exposure = CPU
spent in that subtree (a hazard per CPU, so recent packings are censored correctly: little CPU, little evidence).

V_n is the one true input (a preference; Evan 10-20, as guidance: 1 for a record at n ≤ 324, 1/3 at n > 324, a 274
record the campaign's focus).  It is not a rule either: report how the allocation changes over a range of V (§9) and
only ask when the decision is sensitive to it.  V moves compute only between arms whose outputs land in different n
classes (TS compares Σ_n λ_{a,n} V_n).

Scope is derived, not stated.  Many generators fix an input and produce whatever n results, so an arm's value is over
its output distribution, value rate(a) = Σ_n λ_{a,n} V_n + its potential gain (incl. stepping stones toward 274 via
cross-n lineage).  "Only arms that can produce a 274 outcome" is then a consequence: an arm with no 274 output and no
274 potential scores low on a 274-weighted V, not zero; nothing is excluded by fiat.

## 3. Posterior per arm

Per job j of arm a: realised reward R_j, CPU c_j.  R is zero-inflated and heavy-tailed, so model its two parts:

* records: K_j ~ Poisson(λ_a c_j)
* pool value: ΔΦ_j ~ compound Poisson, rate μ_a c_j, marks from a pooled (per g) log-normal

Hierarchical Gamma priors λ_a ~ Gamma(k, k / λ̄_parent(a)) (same for μ), λ̄ from the next level up (g, n-band), k
fitted by empirical Bayes.  A new arm starts at its family's mean, not at an optimistic or zero guess.

**Depletion** (why arms "play out", SR 10-19): rates decay with CPU spent since the arm's pool last changed,
λ_a(t) = λ_a0 · (1 + C_a(t)/κ_g)^(−β_g), the same form as the explorer's (1 + expanded)^−0.78.  New pool packings from
other arms reset C_a in proportion to their potential: this is how generators raise expanders' rates.

Everything is recomputed from logs (stateless): when h is refitted, past rewards are rescored.

## 4. Decision rule: Thompson sampling, batched

When cores free up: draw (λ̃, μ̃) for every arm from its posterior (with depletion at the current C_a), score
(λ̃ V + μ̃ m̄) per CPU-h, launch the argmax with quantum c (≈ 1-2 CPU-h; smaller for new arms).  Repeat per free slot,
with the posterior updated by a pending-job fantasy (the job counts as spent CPU with its expected reward) so one
batch does not pile onto one arm.

Why TS: the objective is records over an open horizon (regret), not best-arm identification; TS allocates to each arm
in proportion to the posterior probability that it is best, so **every arm keeps a standing share that shrinks only as
evidence against it accumulates**.  That is "standing generators" made exact.  Knowledge-gradient / VOI matters near
the end of a fixed horizon; flag for the 274 campaign if it gets a deadline.

**Information budget** (why the proxy is not optional): at ~1 record / 49 CPU-h, telling two arms 2x apart from records
alone needs ~15 events each, ~700 CPU-h per arm.  Pool value events (new near-frontier minima) are 10^2–10^3x more
frequent, so arms are ranked by ΔΦ in practice, and records calibrate h (and so ΔΦ) across all arms together.

## 5. Downsides of standing generators, and the general fix for each

| Downside | Fix in the model |
|---|---|
| dilution over many arms | hierarchical prior + TS: weak arms' shares shrink with evidence, never by fiat |
| arms play out | depletion term (§3), fitted, not a stop rule |
| delayed credit (records 2 generations later) | potential-based reward (§2); h refitted from lineage |
| redundant output (same region again) | submodular Φ: copies add 0, crowded clusters add little |
| arm proliferation / multiple comparisons | new arms enter at the family mean; shrinkage |
| noisy cost (load, wall-clock limits) | CPU budgets only (prerequisite) |
| store growth | inclusive intake; measured growth per CPU-h in the sink log; revisit only if disk forces it |

## 6. Start selection: cluster weighting as a parameter family

Within a job, parents are drawn with

    π_{α,τ}(p) ∝ exp(−gap(p)/τ) · |C(p)|^(−α) · (1 + cpu(p)/CPU0)^(−1)

|C| = the number of distinct minima in p's cluster.  α = 0 is row-uniform (funnel bias toward crowded basins), α = 1 is
cluster-uniform (pk/parents 'cluster'); "moderately" = α = 0.5 as the default.  τ sets depth (gap scale).

To treat α as a math problem rather than a setting:

* **Logging policy = a mixture.**  Each job draws (α, τ) from the grid α ∈ {0, 0.25, 0.5, 0.75, 1} × τ ∈ {3 values}
  (uniform at first, then TS over the grid as θ-arms inside each generator).  Every pick logs its probability under
  the policy that made it (explore.py already logs `pp` for the model pick).
* **Off-policy evaluation for any α (online or retrospective).**  With logged propensities, the doubly-robust estimate
  V̂(π) = mean[ π(p)/π_log(p) · (y − ĥ(p)) ] + E_{p∼π}[ĥ(p)] scores any policy in the family, including untried α.
  The mixture keeps importance ratios bounded (≤ the number of mixture components per logged pick), which is the
  precise sense of "diverse enough weighting": overlap of support, so the estimator's variance is finite and small.
* Outcome y per pick: the pick's realised ΔΦ (+ records).

## 7. Retrospective checks before trusting it

1. h calibration: fit on lineage before a cutoff, predict records / near-frontier minima per CPU after it (by n,
   generator); reliability plot.  Cutoff 10-17, test 10-17..10-20.
2. Posterior predictive: arm-level predicted vs realised ΔΦ on the test window.
3. Depletion: SR's data (records as seeds play out) and SW2's per-n trajectories give β, κ.
4. Pre-register the first TS window in PACKER (arms, quantum, V_n, metric: records + ΔΦ per CPU-h vs the previous
   week's hand allocation), like every other run.

## 8. Build order

0. (built 10-20) sink: every run's results in the store, inclusive policy; hits + near misses on import.
1. CPU budgets for explore / dgen / hborder / lcut / bulbgen (`cpu_limit`), so cost is measured.
2. Arm registry + launcher wrapper over `jobs.py queue`: arm id, θ, propensity, budget into the manifest; sink
   attributes outputs to the job.
3. h from store lineage (extend the explorer's Poisson model), backtest (§7.1).
4. π_{α,τ} in pk/parents + explore.py with per-pick propensity logging; the (α, τ) mixture.
5. Posterior + TS scheduler (feeds the queue); reports (allocation, posterior per arm, OPE of α).

Every arm is an arm: the finisher (few samples: its posterior stays near its prior, which is the honest state),
arms with bad history (CS274: 1,662 rows at 274, all grid, is a likelihood, not a verdict), arms never run at 274.
candidate.py is a checker, not an arm.  Pre-registered predictions stay: they are elicited priors, and scoring them
(Brier / log score over PACKER's history) calibrates how much weight elicited priors get (§9).

## 9. What the framing resolves, and what is left open

Resolved by the math (no judgment call): arm inclusion / retirement (posterior shares); thresholds like SW2's
"1 / 40 vs 1 / 100" (replaced by the posterior); the default α (a prior over the α grid, then learned); horizon H
(the empirical lineage delay distribution, records vs generations / CPU since the ancestor); the quantum (setup cost vs
adaptivity; per-arm, from measured overhead); pipeline credit (potential Φ).

Left open, by kind:
1. Utility V_n: a preference.  Handled by sensitivity analysis; asked only where allocation flips.
2. Model specification (hierarchy groupings, h's features, depletion form, stationarity): assumptions, but testable:
   chosen by held-out predictive score (§7), not taste.
3. Priors where data cannot speak (new arms, the finisher): from the family level, and from the calibration of past
   elicited predictions.
4. Building new arms (junction refill, bulbgen at 274, ...): a VOI question in developer time, not CPU; the same
   prior-from-family gives an expected value per arm, but the cost side (our time) is outside the CPU bandit.
5. CPU outside the campaign (other projects, census, exact forms): the same machinery if their outputs get a V; until
   then a fixed split.

## 10. Target choice is the same calculation (Evan 10-20: "even the 274 record hunt is a math problem ... where value
and likelihood align")

Events have values, not only n: V(first packing below the grid at n) ≫ V(improving an existing sub-grid record).
The campaign target is argmax over (n, event) of V × P(event | effort), so it is derived the same way as arm shares.

Likelihood data: the sub-grid onset m(k) = max{n : s(n) < k known} (frontier 10-20):

| k | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| m(k) − (k² − k) | −1 | −1 | −1 | 0 | 0 | 0 | +1 | +1 | +1 | **+1** | +2 |
| slack k − s(m) | .054 | .056 | .050 | .0032 | .0143 | .0208 | .0004 | .0065 | .0129 | **.0161** | .0029 |

Plateaus of 3-4 k; slack grows inside a plateau, then the offset steps up and the slack resets.  274 = 17² − 17 + 2
would be the step at k = 17 (k = 18 has already stepped: 308 = 18² − 18 + 2).  Crude reading: the step at 13 → 14
came after slack 0.0208; 17's slack is 0.0161, so the extrapolation says "close call", and our best non-grid 274
(17.000478, +4.8e-4) is already nearer than that.  Two likelihood signals per n disagree elsewhere: 157 (k = 13,
offset +1) has the largest unstepped slack (0.0208) but our best non-grid 157 is +1.1e-2, and smaller n have had more
outside search (evidence against an easy sub-grid packing).  A model of P(first sub-grid at n | slack, our closest,
outside effort) reconciles these and ranks the next targets (candidates: 157, 242, 133, 309); 274 tops both signals.
