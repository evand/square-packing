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

V_n (decided, Evan 10-20): **V = 1 for a record at n ≤ 324, 1/3 at n > 324** ("we'll see what we get").  What the
ratio does: TS compares rate × V, so V only moves compute between arms whose output lands in different n classes;
an n > 324 record must come 3x cheaper to draw compute away from n ≤ 324.

**Campaign scope (274, Evan 10-20): only arms that can produce a 274 outcome are launched**, but many generators fix
an input and produce whatever n results (deletion / insertion, L-cut, refill, compose).  So an arm's value is over its
output distribution:

    value rate(a) = Σ_n λ_{a,n} V_n,     eligible iff P_a(output at 274, or a 274 pool packing) > 0

and a record at any n counts at its V_n.  λ_{a,n} pools across n through the hierarchy (an arm's rate at 273 informs
its rate at 274).  The campaign's own target enters through h: pool packings at 272-276 that are stepping stones to 274
carry potential at 274 through cross-n lineage.

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

The near-miss finisher is an arm like any other (Evan 10-20: probably high value, small sample, hard to test, so it runs
mostly on its prior; the hierarchical prior is exactly that, and its posterior moves as items finish).  candidate.py is
a checker, not an arm.  Still open: horizon H, the quantum.
