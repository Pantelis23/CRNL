# CRNL — findings

Measured results, with the numbers, the caveats, and what is still open. Raw data
for every table is committed under `results/`; each figure is regenerable by the
named experiment.

Speculation lives next door in [`THEORIES.md`](THEORIES.md) — conjectures that make
falsifiable predictions, questions with no theory yet, and the four confident wrong
results this project has produced, kept with what killed them.

Everything here is Approximate Majority (AM) or its n-winner generalization, all
rate constants k = 1, simulated with the exact Gillespie SSA against the
mass-action ODE. Ω is the molecule count at concentration 1.

---

## 1. The restoration wall — `restoration_wall.py`

Error probability from a biased start falls exponentially in population:

    P(error) ~ exp(−c(ε)·Ω)

Measured **c(0.10) ≈ 0.018, R² ≈ 0.94**. The deterministic ODE from the same start
reaches the correct rail at every Ω — error exactly 0. That gap between the two
engines is the entire subject of the project.

The all-blank outcome, which the ODE calls an unreachable repeller, genuinely
occurs at small Ω (≈5% at Ω=6, gone by Ω=16) — a pure finite-count effect.

**Note on the bias.** The design doc's illustrative 51/49 (ε=0.02) has an
intrinsically tiny barrier: c·Ω stays below ~1 across the whole observable window,
so the wall never clears the algebraic-prefactor crossover until Ω ~ thousands,
where errors fall to ~e⁻²⁵ and read as exactly zero. The default is ε=0.10, which
puts a clean exponential in the accessible window; 51/49 remains one flag away.

---

## 2. Predicting the barrier from the saddle — `quasipotential.py`

The one result here that is *theory*, not measurement. Derived in `docs/design.md`
§9 from the saddle geometry and finite-count noise:

    c(ε) = (3/2)·ε²          [quadratic in bias, prefactor κ = 3/2]

Measured across ε ∈ [0.04, 0.20] (`results/quasipotential_ceps.json`):

| ε | c(ε) | c/ε² | R² |
|---|------|------|-----|
| 0.04 | 0.00254 | **1.586** | 0.999 |
| 0.06 | 0.00590 | 1.638 | 0.991 |
| 0.08 | 0.01099 | 1.716 | 0.976 |
| 0.10 | 0.01720 | 1.720 | 0.977 |
| 0.14 | 0.03435 | 1.752 | 0.895 |
| 0.20 | 0.07237 | 1.809 | 0.899 |

Power-law fit: **p = 2.08** — quadratic confirmed. The prefactor drifts *down*
toward the predicted 1.5 as ε → 0, reaching **1.586 at the smallest and cleanest
point** (within 6% of a first-principles prediction with no fitted parameters).

The residual upward drift at larger ε is expected: the derivation linearizes about
the saddle, so it is a small-bias expansion. Part of the gap is also real physics —
the bias amplifies during the initial B-buildup transient (dδ/dt = δ·b > 0), so the
effective δ₀ entering the saddle exceeds ε.

---

## 3. Radix vs. margin — `radix_scaling.py`, `radix_wall.py`

n committed symbols, champion leading **each** rival by a fixed pairwise margin
δ = 0.10 (= 55/45 at n=2). Error = champion does not win.

| n | c(n) | Ω for ≤5% loss | R² |
|---|------|-----|-----|
| 2 | 0.01567 | 96 | 0.78 |
| 3 | 0.00659 | 341 | 0.99 |
| 4 | 0.00455 | 530 | 1.00 |
| 6 | 0.00349 | 767 | 1.00 |
| 8 | 0.00297 | 905 | 1.00 |
| 12 | 0.00265 | 1058 | 0.99 |
| 16 | 0.00235 | 1179 | 1.00 |
| 24 | 0.00230 | 1255 | 1.00 |
| 32 | 0.00227 | 1293 | 1.00 |
| 64 | 0.00224 | — | 0.99 |

**The barrier saturates.** c(n) falls ~7× from n=2 to n=32 and then stops: 24→32 is
a 1.3% change, and n=64 (measured separately, band Ω∈[600,1100]) lands on the same
floor ≈ 0.0022. It is **not** a clean power law — an `n^−0.60` fit passes through
the middle and misses both ends.

**Why it saturates, and the scope of the claim.** With a *fixed pairwise margin*,
the champion's share of the whole population converges to δ as n grows (each of the
n−1 rivals shrinks toward zero). Past n ≈ 16 the contest stops changing — it is
asymptotically "one species at 10% against a 90% field split infinitely finely" —
so adding symbols no longer crowds the champion further. **Under this convention
the radix penalty on the margin is bounded.** A different convention (fixed
champion *share*, or symmetric plurality) could well be unbounded; that is untested.

The cost that grows monotonically is **population**: Ω_required rises ~13× (96 →
1293) over the same range, itself flattening.

**Symmetric start** (`radix_discovery.py`): from an unbiased start the system
resolves to a single winner with probability ≈1.0 at every n tested (2–48), under
both fixed-total-Ω and fixed-density conventions. The hypothesised high-n
"all-blank collapse" **does not occur** — an earlier prediction the simulation
overturned. The cost of radix appears instead as rising consensus *time*
(fixed-total Ω=120: 16.1 → 26; fixed-density: 13.2 → ~35 out to n=48).

### 3.1 The penalty is convention-dependent — and fixed margin is the right one

§3 flagged its own scope: "A different convention (fixed champion *share*, or
symmetric plurality) could well be unbounded; that is untested." Now tested
(`radix_convention.py`), with both conventions anchored to the same 0.55/0.45 start
at n=2, Ω=120, 3000 trials:

| n | 2 | 3 | 4 | 8 | 16 | 24 |
|---|---|---|---|---|---|---|
| P(win), fixed margin | 0.971 | 0.802 | 0.685 | 0.522 | 0.400 | **0.370** |
| P(win), fixed share | 0.971 | 1.000 | 1.000 | 1.000 | 1.000 | **1.000** |

Under fixed share the penalty does not weaken — it **vanishes entirely**. But the two
conventions are not two readings of one experiment, because only one holds the
contest fixed:

| n | 2 | 3 | 8 | 24 |
|---|---|---|---|---|
| pairwise lead, fixed margin | 0.100 | 0.100 | 0.100 | 0.100 |
| pairwise lead, fixed share | 0.100 | 0.325 | 0.483 | **0.525** |

Fixed share hands the champion a **5.2× larger lead** by n=24. It is asking an easier
question at every n, not measuring the same one differently. **So §3's choice of a
fixed pairwise margin is the convention that isolates alphabet size, and its finding
stands.**

Nor is *share* the governing variable: at a fixed share of 0.50, P(win) is 0.606 at
n=2 and 0.997 at n=3, because the same share is a pairwise lead of ~0 at n=2 and 0.25
at n=3. Neither share nor symbol count alone sets the difficulty.

**A prediction that failed, recorded.** `THEORIES.md` T3 predicted fixed share would
give an *unbounded* penalty, "because the champion's per-rival lead keeps shrinking."
The lead is `s − (1−s)/(n−1)`, which **grows** toward `s`. The reasoning was wrong
before any code ran — caught by writing the prediction down before the sweep instead
of after.

**Untested:** symmetric plurality, and whether §3's saturation floor itself moves
under any convention that keeps the pairwise lead fixed by other means.

---

## 4. Chemical freeze-out in an expanding volume — `expansion.py`

Volume expands as Ω(t) = Ω₀·e^{Ht}; bimolecular propensities dilute as 1/Ω, so the
restoring reactions slow while the decision is still being made. The SSA is exact
for this case (closed-form waiting time; freeze-out is the event that the remaining
integrated propensity a₀/λ is finite and never reaches the target — see `expanding.py`).

**There is a critical expansion rate.** Below it consensus completes and the system
reaches a clean rail; above it the reaction freezes mid-decision, locking in a relic
minority abundance that grows toward 0.5. This is the chemical analogue of
cosmological freeze-out (the Γ-vs-H competition that set the relic dark-matter
abundance and primordial helium) — shared mathematics, not a claim about cosmology.

At the transition the frozen composition is maximally broad (critical fluctuations);
below it the distribution is sharply bimodal at the rails, above it a single spike
at the initial 50/50.

---

## 5. Freeze-out is a genuine transition — `freezeout_scaling.py`

Six system sizes spanning ×32 collapse onto a single master curve under

    D(H, Ω) = F((H − Hc)·Ω^a)

with **Hc ≈ 0.055** and **a ≈ 0.38** (`results/freezeout_fss.json`). Finite-size
crossings drift toward Hc as predicted:

| Ω | 40 | 80 | 160 | 320 | 640 | 1280 |
|---|----|----|-----|-----|-----|------|
| H*(D=0.5) | 0.1658 | 0.1409 | 0.1219 | 0.1078 | 0.0972 | 0.0879 |

A crossover would not collapse; this does. So the sharpening seen in §4 is genuine
finite-size scaling about a critical point, and in the Ω→∞ limit the transition is
sharp: expand slower than Hc and you always decide, faster and you never do.

**Open — do not over-read the exponent.** a ≈ 0.38 comes from a two-parameter grid
search with no error bars, and sits between 1/3 (0.333) and 2/5 (0.400). The
universality class is **not** identified. Connecting it to the quasipotential of §2
is the obvious next theoretical step.

---

## 6. Bigger alphabets freeze more easily — `expansion_radix.py`

Running n-winner AM under expansion ties §3 and §4 together. The critical rate falls
monotonically with alphabet size (Ω=160, D=0.5 crossing):

| n | 2 | 3 | 4 | 6 | 8 | 12 | 16 |
|---|---|---|---|---|---|----|----|
| H*(n) | 0.121 | 0.110 | 0.101 | 0.089 | 0.079 | 0.074 | 0.071 |

H*(n) also saturates toward a floor, mirroring c(n). The frozen relic is richer:
fast expansion leaves ≈n coexisting species, versus 2 for AM.

**So a larger alphabet is penalised twice** — a lower restoration margin *and* a
lower tolerance for expansion — and both penalties bottom out rather than diverging.

---

## 7. Why restoration matters at all — `cascade.py`

The founding claim, finally built. A one-bit signal passes through D stages, each
injecting the same channel noise (σ=0.35); the two modes differ only in what the
stage does with it — restoring (AM snaps the noisy signal back to a full-magnitude
rail) versus non-restoring (analog passthrough).

| depth | 1 | 10 | 22 | 45 |
|---|---|----|----|----|
| non-restoring | 0.80 | 0.57 | 0.51 | 0.50 |
| restoring (Ω=80) | 0.80 | 0.79 | 0.77 | 0.74 |

The non-restoring cascade reaches the coin flip by depth ~22 — the signal is gone,
deep computation impossible. The restoring cascade still carries the bit at depth 45.

**Stated precisely:** restoration does not make per-stage error zero. It makes it
exponentially small in Ω, so survivable depth scales like e^{cΩ} instead of being
set by the noise variance. Both curves decay; the rates differ by orders of
magnitude. The residual slope of the restoring curve is exactly the finite-Ω
reliability of §1 showing through, and it flattens as Ω grows.

---

## 8. Cross-project: does this predict a real planner? — EIR bridge

[EIR](https://github.com/Pantelis23/eir) relaxes discrete logic into a continuous
energy field, anneals noise, and decodes by argmax — the same "logic in, landscape
out" move CRNL studies as physics. Mapping: EIR's annealing budget (`steps`) ↔ 1/H,
option count ↔ alphabet size n.

An n-option decision (one `DiscreteVar`, `OneHotFactor` + `ActionBiasFactor` with
margin 0.1), sweeping budget × n (`results/eir_bridge.json`):

- **Freeze-out appears.** Decision quality rises with annealing budget and
  saturates; a fast quench leaves the planner frozen near random.
- **Decode reliability falls with option count**, staying above the 1/n random
  baseline but by a shrinking margin.

**The radix half was then tested under control, and it holds.** The first pass used
*fixed* hyperparameters at every n, which conflates a real effect with settings that
merely suited small n. Re-run with **each n given its own best `(lr, noise_start)`**
over a 3×3 grid at a large fixed budget (`eir/tools/crnl_bridge.py --mode radix`):

| n | 2 | 3 | 4 | 6 | 8 | 12 | 16 | 24 | 32 |
|---|---|---|---|---|---|----|----|----|----|
| P tuned | 0.957 | 0.847 | 0.797 | 0.687 | 0.590 | 0.507 | 0.423 | 0.330 | 0.267 |
| P untuned | 0.867 | 0.797 | 0.747 | 0.657 | 0.557 | 0.480 | 0.397 | 0.300 | 0.233 |
| tuning gain | +0.09 | +0.05 | +0.05 | +0.03 | +0.03 | +0.03 | +0.03 | +0.03 | +0.03 |

Per-n tuning lifts every point by a near-constant +0.03…+0.09 — it shifts the curve
without flattening it. **The falloff is structural, not a hyperparameter artifact.**

**Remaining caveat.** This is one test problem (a single decision variable with a
uniform cost gap). It shows the phenomenon exists in EIR; it does not establish the
same *mechanism* as CRNL's basin crowding — EIR's `FieldBackend` also stops early
once it ever hard-decodes the true optimum, so part of the n-dependence is a search
effect (fewer chances to hit a smaller target) rather than margin erosion. Telling
those two apart would need a variant that disables early stopping.

---

## 9. What restoration costs in free energy — `reversible_landscape.py`, `dissipation_decision.py`, `dissipation_memory.py`

Every result above takes the project's founding thermodynamic claim on faith.
Irreversible AM has *formally infinite* dissipation — there is no number to report —
so measuring the cost required rebuilding AM as a proper thermodynamic CRN: every
reaction reversible, all reverse rates scaled by one parameter γ.

    X + Y ⇌ 2B        B + X ⇌ 2X        B + Y ⇌ 2Y      (reverse rate = γ · forward)

Detailed balance requires the Wegscheider condition `k_f1k_f2k_f3 = k_r1k_r2k_r3`,
i.e. `γ³ = 1`, so **every γ < 1 is genuinely driven**. `rank(S) = 2` makes the cycle
space exactly one-dimensional, so the whole drive is a single number:

    A(γ) = −3·ln γ          (dimensionless; the free energy is k_B T · ΔS)

γ→1 is equilibrium; γ→0 recovers today's irreversible AM. All measurements here are
**exact** — the chemical master equation solved by sparse linear algebra on the
conserved simplex (7381 states at Ω=120, 0.20 s), not sampled.

**Why exact, stated correctly (a corrected claim).** An earlier version of this
section said a direct SSA measurement at Ω=120 "would take hundreds of hours." That
is wrong wherever the exact solve is trustworthy: at measured throughput (84,000
steps/s, ≈0.4·Ω steps per unit time) one flip at Ω=120, γ=0.35 costs **5.5 minutes**,
and at Ω=60, γ=0.30 about **one minute**. The hundreds-of-hours regime is real but
sits at γ ≤ 0.30 at Ω=120 (58 h per flip at γ=0.30, 3·10⁵ h at γ=0.25) — which is
**exactly where the CME's own first-passage solve fails its validity guard**. So the
honest statement is not "the CME is cheaper": it is that the CME is *exact, gives the
whole MFPT field from one solve, and carries no sampling error*, while neither
instrument currently reaches the strong-drive corner at large Ω. Rescuing that corner
(the direct solve loses precision; a quasi-stationary eigenvalue would not) is open.

### 9.1 A landscape has a minimum price: A > 3 ln 2

The symmetric point sits at `(⅓,⅓,⅓)` for **every** γ, and the decision mode there has

    λ(γ) = (1 − 2γ)/3

which is `+1/3` at γ=0 — exactly the irreversible AM saddle eigenvalue of
`design.md` §2.3, so the reversible family contains today's AM as its γ→0 limit — and
vanishes at

    **γ_c = 1/2      A(γ_c) = 3·ln 2 = 2.0794**

Above γ_c the two rails have merged into the symmetric point (a pitchfork, with
`δ* ∝ √(γ_c−γ)` approaching `4√2/3`): a single minimum, and **no population size Ω
can restore, because there is nothing to restore toward.** The attractors also move
inward as the drive weakens, with `b* = γ/(1+γ)` exactly — the thermal population is
mainly *blank*, while the losing symbol stays quadratically small.

Verified four independent ways: closed-form algebra, the engine's own numeric
Jacobian (agreement **5.6e-17**), a `fsolve` enumeration from 1830 simplex starts
(exactly 3 fixed points below γ_c, 1 above), and two adversarial reviews.

> `A(γ_c) = 3 ln 2` is 3 reactions × ln 2. The resemblance to Landauer's `k_B T ln 2`
> is arithmetic coincidence, not physics.

### 9.2 The cost of deciding, and why the exchange rate is not constant

Entropy production is exact per jump, `ΔS = ln[a_ρ(n)/a_{−ρ}(n′)]` with the reverse
propensity evaluated *post*-jump. For this network it also has a closed form that
holds at **every** state and every γ (verified to 9.3e-15 over 6648 jumps):

    ΔS = ln W(n′) − ln W(n) + s_ρ·ln(1/γ),      W = multinomial coefficient

so along a trajectory `⟨ΔS⟩ = ln[W(n_stop)/W(n₀)] + (A/3)·⟨M_fwd − M_rev⟩` — a
boundary term and a cycle term, both independently meaningful, no cancellation.

At Ω=120, holding the decision threshold at a constant fraction of the landscape
(θ/δ* = 0.70) and the bias at 0.20 δ*:

| γ | A | P(error) | boundary | cycle | total ⟨ΔS⟩ |
|---|---|----------|----------|-------|------------|
| 0.05 | 8.99 | 0.0006 | +4.5 | 562 | **567** |
| 0.15 | 5.69 | 0.0092 | +13.2 | 417 | 430 |
| 0.25 | 4.16 | 0.066 | +20.2 | 421 | 441 |
| 0.35 | 3.15 | 0.167 | +28.7 | 431 | 460 |
| 0.49 | 2.14 | 0.398 | +44.2 | 88 | **133** |

**4.3× more free energy buys 664× lower error** — but the exchange rate is very far
from constant. Across γ ∈ [0.15, 0.40] the cost sits flat at **430–470 k_BT while
the error varies 25×** (0.009 → 0.25). The mechanism: raising γ lowers the affinity
per cycle (each firing is cheaper) but requires more firings and more time, and the
two partly cancel. So "restoration costs dissipation" is true, and the naive
monotone reading of it is not.

Cost is extensive in Ω and grows as `ln(1/γ)` — roughly `1.9·Ω·ln(1/γ) k_BT`, i.e.
~2 k_BT per molecule per unit of affinity. Any "orders of magnitude above
`k_B T ln 2`" statement is an instance at a stated (Ω, γ), not a bound; Landauer does
not apply to this protocol.

**A protocol trap, recorded because it produced a convincing false result.** Two
distinct artifacts bite here and they need *opposite* fixes: the initial bias must
not jitter on the integer lattice (one molecule ≈ 20 k_BT), while the threshold must
scale with the shrinking landscape. Holding *both* fixed put θ outside the landscape
above γ≈0.42 (θ/δ* = 1.88 at γ=0.49), which turned "deciding" into "fluctuating past
the attractor" and produced a clean U-shaped dissipation curve with a minimum near
γ≈0.3 — a plausible-looking dissipation optimum that was pure artifact. Holding
θ/δ* constant instead removes it. Both figures were generated; only the corrected
protocol is reported.

### 9.3 The cost of remembering — retention is bought by drive, not by power

At finite γ a decided state is only metastable: the reverse reactions regenerate
blank and let the loser take over. Exact mean first-passage lifetimes τ and
stationary dissipation rates σ:

| Ω | γ=0.30 | γ=0.40 | γ=0.49 |
|---|--------|--------|--------|
| 30 | τ=4.6e3, σ=0.82 | τ=290, σ=1.90 | τ=42, σ=2.19 |
| 60 | τ=1.9e5, σ=1.54 | τ=853, σ=4.00 | τ=65, σ=5.18 |
| 120 | **dropped** (σ=3.01) | τ=6.1e3, σ=7.67 | τ=102, σ=11.5 |

> **Corrected cell.** The (Ω=120, γ=0.30) entry previously read "τ=5.5e5, σ=4.93".
> That is the **γ=0.35** row (τ=5.481e5, σ=4.927) published under the wrong label; the
> real (120, 0.30) solve returns τ=3.5e8 with residual 2.5e-6 and is `valid=False`, so
> it is one of the 11 drops counted in the caveat below. σ is unaffected — it comes
> from the stationary solve, not the first-passage solve — and the true σ(120, 0.30) =
> 3.01 leaves the σ story unchanged (still rising monotonically through γ=0.49 at this
> Ω, just more steeply). Every other cell in the table was re-checked against
> `results/dissipation_memory.json` and is correct, as is the 17,800× figure below
> (τ = 41.7 → 7.43e5 at Ω=30, γ = 0.49 → 0.20, both valid).

**Retention is exponentially sensitive to drive**: at Ω=30, raising A by 2.3×
(γ 0.49→0.20) buys **17,800× longer memory**.

Three facts break the naive "restoration costs dissipation" reading as a statement
about *power*:

1. **σ → 0 in *both* limits** — at γ=1 by detailed balance, and as γ→0 because
   `σ = A·J` and the cycle flux collapses faster than `A = −3lnγ` grows.
2. **σ and τ move in opposite directions** across the bistable range: more steady
   power, *less* retention. But σ is a rate and τ a time — that pairing is not a
   correlation, and the dimensionless `σ·τ` (total dissipation per lifetime) is
   monotone and reads the other way.
3. **The zero-power memory limit is textbook, not news.** γ→0 is a *singular* limit
   in which the memory states become absorbing (a halted dynamics dissipates
   nothing — so does a rock), and γ→0 *is* A→∞, so infinite affinity is being
   treated as free. This is the ideal zero-leak ratchet.

**The corrected claim.** `design.md` §2.3's landscape statement — "a system at
equilibrium sits in a single free-energy minimum: no threshold, no restoration" — is
**confirmed and quantified** by γ_c = 1/2. What needs narrowing is only the cost
clause: *deciding* costs `O(Ω)·A`, and each cascade stage pays again; *holding* a
decided state costs no power in the zero-leak limit. Restoration requires a minimum
**affinity**, not a minimum dissipation rate.

**Caveats.** Single test problem, symmetric rate constants, one threshold
convention. σ's peak location is Ω-dependent and can sit *outside* the bistable
window (at Ω=30 it peaks at γ=0.45; at Ω=60 and 120 it is still rising at γ=0.49).
11 of 27 (Ω,γ) points were **dropped** by the solver's validity guard at small γ and
large Ω, where the MFPT linear solve returns a negative time (−5.9e15 at Ω=60,
γ=0.10) — those are reported, never fitted.

## 10. Sampling confirms the exact numbers; and what a restoring stage costs — `dissipation_decision.py --ssa-trials`, `dissipation_memory.py --ssa`, `dissipation_cascade.py`

§9 is entirely exact linear algebra. This section does two things: checks it by
sampling, and finally prices a *restoring cascade stage* — the measurement
`design.md` §8 asked for and §7 left unpriced.

### 10.1 The exact solve and the SSA agree

An instrumented Gillespie loop (`crnl/thermo.py`) accumulates entropy production as
an **integer counter** — no logarithms in the hot loop — using §9.2's closed form. It
is pinned bit-for-bit against the verified `gillespie_fast` on the same seed.

Part A, Ω=120, 20 000 trials per γ:

| γ | ⟨M⟩ exact | ⟨M⟩ sampled | agreement | P(err) exact | P(err) sampled |
|---|-----------|-------------|-----------|--------------|----------------|
| 0.30 | 339.1 | 337.9 ± 1.40 | 0.86 SEM | 0.0971 | 0.0967 ± 0.0021 |
| 0.45 | 445.7 | 445.7 ± 2.57 | 0.00 SEM | 0.3457 | 0.3509 ± 0.0034 |

Part B compares the lifetime τ against a hysteretic flip counter. Seed-averaged over
8 trajectories, the ratio τ_SSA/τ_CME is **0.868–1.043** across six (Ω, γ) points.

Two things this cross-check caught, both of which would have been published:

- **The flip-rate convention was inverted.** A one-way Schmitt counter gives
  `flips/T → 1/τ`; `1/(2τ)` is the *round-trip* rate. The wrong convention produced
  ratios of 0.37–0.58 with 52–82 flips — a confident false alarm, and the original
  note told the reader to "fix" the correct code.
- **A single trajectory is not a measurement.** One run has sd ≈ 0.26 in the ratio;
  the first sweep produced 1.50 and 1.58, which over 8 seeds became 0.950 ± 0.091 and
  1.047 ± 0.112. A predicted arm-vs-attractor offset of ~0.8 also **failed to appear**
  (mean 0.97) — it confused an MFPT-from-the-arm with the mean time between crossings
  of a long trajectory, which is set by the full dwell.

**The ⟨M⟩ column carries the power**, not P(error): a 20% protocol error is 15.7 SEM
in ⟨M⟩ at 2000 trials but only 2.9 SEM in P(error).

### 10.2 A restoring stage costs more exactly where it works less

A stage seeds a fresh vessel (B=0) from the previous stage's transmitted output, runs
for a **fixed time**, and emits the δ the chemistry actually reached — no threshold,
no `sign()`, no renormalization. Channel noise is in landscape units,
`σ_ch = 0.35·δ*(γ)`; since **δ*(0) = 1 exactly**, §7 is the γ→0 member of this family.
Everything is exact (an augmented generator, no sampling, no quadrature).

Depth 30, both control conventions reported:

| γ | A | Ω=30 ΔS/stage → fidelity | Ω=60 | Ω=120 |
|---|---|---|---|---|
| 0.05 | 8.99 | 21.7 → 0.886 | 44.1 → 0.910 | 89.5 → **0.921** |
| 0.15 | 5.69 | 20.0 → 0.807 | 41.1 → 0.862 | 83.3 → 0.890 |
| 0.30 | 3.61 | 26.7 → 0.562 | 57.0 → 0.632 | 116.9 → 0.698 |
| 0.45 | 2.40 | 30.5 → 0.500 | 69.3 → 0.501 | 149.1 → **0.502** |

**At Ω=120, 1.67× the free energy per stage buys a total loss of function** (89.5 k_BT
→ 0.921 versus 149.1 k_BT → 0.502, a coin flip). Cost rises toward γ_c because the
shrinking landscape demands more cycling to hold a bit it can no longer hold. So
restoration does not degrade gracefully into cheapness — **it degrades into paying
more for nothing.** Cost is extensive in Ω throughout, consistent with §9.2.

Robust across both controls and every Ω: fidelity falls monotonically as γ→γ_c, and
rises monotonically with Ω. Near γ_c a "restoring" stage is genuinely worse than a
passive channel — at γ=0.45 it sits at the coin flip for Ω = 30, 60, 120 **and 240**,
at every stage time from 0.5 to 150, while its single-stage flip rate falls 340×.
**8× the population buys 0.0045 of depth-30 fidelity there.**

### 10.3 A withdrawn claim, and the two designs that produced it

An earlier version of this section claimed **"restoration requires a minimum Ω as well
as a minimum affinity."** *It is withdrawn.* It was an artifact of the control.

The control walked a lattice hard-limited at ±1 — an **absolute** dynamic range —
while its noise was scaled by δ*(γ). That mismatch grows with γ, in exactly the
direction of the claimed effect. Rail the control to the chemistry's own ±δ*(γ) and it
becomes γ-independent (0.5245/0.5226/0.5257/0.5269 — spread 0.003, versus 0.154 for
the ±1 control), and the crossover disappears:

| γ=0.30 | chemistry | control ±1 | control ±δ* |
|---|---|---|---|
| Ω=30 | 0.5617 | 0.5945 → **loses** | 0.5196 → **wins** |
| Ω=60 | 0.6321 | 0.5921 → wins | 0.5257 → wins |

The chemistry arm is *identical* in both columns. Under §7's absolute-σ convention the
crossover exists but relocates to Ω ∈ (60, 120]. A crossover whose location depends on
the comparator is a property of the comparator. `dissipation_cascade.py` now prints
both columns on every row and flags disagreement; **1 of 12 cells disagrees, and it is
exactly this one.**

The design before that failed differently and worse: it stopped a stage at
`0.7·δ*(γ)` and emitted ±1. Since ±1 exceeds δ*, the stop predicate fired on the
*initial* state — **83–96% of stages ran zero reactions**, the harness performed the
restoration for free, and the survival curve was a bare `sign()` limiter. The reported
cost fell with γ mostly because the *no-op fraction* rose (5.0× of the measured 7.7×
fall). Both failures are the §9.2 artifact class: a landscape that shrinks with γ,
hidden once in a duty cycle and once in a comparator's dynamic range.

**Caveats.** Single test problem, symmetric rate constants. `σ_ch = 0.35·δ*` is a
*choice* (δ*(0)=1 makes §7 its γ→0 member, but §7 starts from a weak 0.3 signal and so
opens at 0.80 where a rail start opens near 1.0 — the comparison is directional, not
quantitative). `ΔS/stage` averages a non-stationary sequence and is meaningless without
its depth. Ω ≤ 120 by cost. `t_stage = 8` is one point on an axis, not a canonical
value.

## 11. The cost of a bit, with no comparator — `bit_cost.py`

Every verdict in §10 needed a passive control, and §10.3 withdrew a claim because
that control's dynamic range was a free parameter. This removes the comparator.

Send `b = ±1` equiprobably as `δ₀ = ±δ*(γ)`. After D stages the two output
distributions are mirror images, so `I(b;X_D) = H(mixture) − H(p₊)` is exactly what
survived. Divide cumulative dissipation by it: **k_BT per bit delivered**. No rails
to choose, no tie band, and comparable to `k_B T ln 2` in the units it is stated in.

**Depth is part of the question, not a nuisance parameter** — and this is the third
time this project has hit the same trap. At **depth 1 the measure is degenerate**: a
stage with `t_stage → 0` does nothing, costs nothing, and still scores well because
one channel application barely damages a bit sitting on a rail. Measured at γ=0.15,
Ω=30: **0.89 k_BT/bit at t=0.05 versus 20.2 at t=16** — the cheapest "restoration
event" is the one that does not restore. (The earlier two forms were a stop predicate
firing on the initial state, and a control free to have different rails.) At depth 30
the ordering **reverses** — t=0.05 costs 5493 and t=1 costs 155 984 — because a
passive channel loses the bit outright there. `cost_per_bit` therefore requires a
depth and the experiment refuses `--depth < 5`.

### Cheapest bit at depth 30 (exact, per (γ, Ω) over a t_stage grid)

| Ω | γ=0.05 | γ=0.15 | γ=0.30 | γ=0.45 |
|---|--------|--------|--------|--------|
| 30 | **1239** | 1649 | 52 337 | ~6·10⁸ |
| 60 | 2195 | 2421 | 14 257 | ~3·10⁸ |
| 120 | 4191 | 4244 | 12 345 | ~3·10⁸ |

(k_BT per bit; best `t_stage` per cell. Full grid in `results/bit_cost.json`.)

**Cheapest bit measured: 1239 k_BT at γ=0.05, Ω=30, t_stage=16** — carrying 0.52
bits of the original 1. That is **1787× `k_B T ln 2`**. Landauer bounds *erasure*,
not transmission, so this is a **scale comparison, not a claim that the bound is
approached**; it is nonetheless the first number in this project that can be put
beside `k_B T ln 2` without inventing a protocol first.

Three things follow, and two invert the naive reading:

1. **Weak drive is not cheap — it delivers nothing.** Cost per bit rises ~40× from
   γ=0.05 to γ=0.30 and diverges by γ=0.45. §9.3 showed a weak drive dissipates less
   *power*; per bit actually delivered it is catastrophically more expensive.
2. **Reliability is bought superlinearly.** Cost per bit *rises* with population —
   1239 → 2195 → 4191 at Ω = 30/60/120 — while information only creeps up (0.52 →
   0.60 → 0.63 bits). Quadrupling Ω buys 21% more information at 3.4× the price. Big
   populations are more *reliable* and less *efficient*.
3. **`t_stage` has an interior optimum** (~16 here): too short fails to restore, too
   long pays for idle cycling once the state has relaxed.

### 11.1 There is no efficiency optimum in Ω — and why that is the wrong question

Extending the grid down to Ω=4 finds **no turnaround**: cost per bit falls
monotonically as the system shrinks (4191 → 1239 → 915 → 658 → 520 k_BT/bit at
Ω = 120 / 30 / 20 / 10 / 4). It is not a quantization artifact — the start snaps to
the full rail (+5.0%) for every Ω ≤ 20 but to −2.0% at Ω=30 and +0.4% at Ω=45, and
the trend runs smoothly across that sign change.

The reason is a second degeneracy, pointing the opposite way to the depth-1 one:
**the ratio is minimized by a system that barely transmits.** At Ω=4 the "cheapest"
cell carries **0.12 bits**. So unconstrained cost-per-bit does not identify an
operating point either, and the well-posed question fixes the information first.

**Efficient frontier** — cheapest total ΔS for each level of information actually
delivered to depth 30 (every point at γ=0.05, t\*=16):

| I (bits) | 0.12 | 0.25 | 0.38 | 0.46 | 0.52 | 0.57 | 0.60 | 0.62 | 0.63 |
|---|---|---|---|---|---|---|---|---|---|
| Ω | 4 | 8 | 14 | 20 | 30 | 45 | 60 | 90 | 120 |
| total ΔS | 62 | 153 | 287 | 420 | 646 | 978 | 1310 | 1982 | 2655 |
| marginal k_BT/bit | — | 712 | 992 | 1741 | 3615 | 6655 | 13014 | 27548 | **54689** |

**The marginal cost of information rises 77×** along that ladder. Half a bit costs
646 k_BT; the next 0.11 bits cost 2000 more. Fidelity is not merely expensive — its
price accelerates, and the last fraction of a bit is unaffordable at any population
this method can reach. (Steps at constant Ω are excluded from that ratio: those are
stage-time tuning, a free lunch, not purchased fidelity.)

**Caveats.** Single test problem, symmetric rates, `σ_ch = 0.35·δ*` still a choice
(δ*(0)=1 makes §7 the γ→0 member). Ω ≤ 120 by cost, so the frontier's top end is
grid-limited, not physical. The measure is comparator-free but **not**
convention-free: `noise_frac`, the depth D, and the 0.05-bit floor on the frontier
remain stated inputs. The frontier is a Pareto set over a discrete grid, so its
points are the best *tested* cells, not proven optima.

## 12. Two regimes, one formula — `channel_wall.py`

§11.1 left a puzzle: cost per bit *rises* with population and the frontier's
marginal cost explodes. The reason is a crossover between two things this project
had measured separately and never connected.

**The restoration wall** (§1–2): finite-count error falls like `exp(−κ·ε²·Ω)`, so
molecules buy exponentially better reliability. **The channel floor** (§11):
injected noise flips the sign outright with an Ω-*independent* probability, so
molecules buy nothing. A stage sits on one side or the other.

A saddle point over where the flip happens gives **one expression for both**. A flip
needs the channel to displace the state from the rail to some δ, then finite-count
noise to finish it; the costs add in the exponent, and minimizing over δ gives

    −ln p  ≈  κ·Ω·δ*² / (1 + 2·κ·Ω·σ²)          κ(γ) = (3/2)(1 − 2γ)

with **no fitted parameter**. `κ` is `design.md` §9's `3/2` corrected for the
landscape's restoring gain — `κ = (9/2)·λ_antisym(γ)`, so it vanishes at γ_c along
with the ability to restore. The limits are the two known results:
`2κΩσ² ≪ 1` → `κΩδ*²`, the wall; `2κΩσ² ≫ 1` → `δ*²/(2σ²)`, the floor. The crossover
is at **Ω× = 1/(2κσ²)**.

**Measured, 216 cells** (γ × σ_ch/δ* × Ω, exact):

| subset | slope | R² | n |
|---|---|---|---|
| pooled | 0.742 | **0.9329** | 216 |
| γ=0.05 | 0.795 | 0.9894 | 54 |
| γ=0.15 | 0.626 | 0.9564 | 54 |
| γ=0.30 | 0.419 | 0.8944 | 54 |
| γ=0.45 | 0.497 | 0.9688 | 54 |

The γ-correction is what makes this work: with `κ` fixed at 3/2 the pooled fit is
**R² = 0.69** with per-γ slopes running 0.79 → 0.08. Slope is not 1 because the
saddle point drops the prefactor and the Gaussian-tail correction — **the claim is
that one expression collapses both regimes, not that it is exact.**

**The population-limited side is dramatic and was unexplored.** At γ=0.05,
σ_ch/δ* = 0.10, the per-stage flip probability falls **eleven orders of magnitude**
(2.0e-3 at Ω=4 → 1.8e-14 at Ω=96), a clean exponential with R² = 1.000. At
σ_ch/δ* = 0.45 the same population change moves it 2.2×. Every cascade result in
§10–§11 used 0.35, which sits on the channel-limited side — **which is precisely why
§11.1 found no efficiency optimum in Ω.** The frontier saturates because the protocol
was operating where molecules cannot help.

**Caveat that matters more than the R².** At γ=0.45 the formula fits well (0.969)
**for the wrong reason**: p there is nearly independent of the channel, moving only
1.36× (0.089 → 0.121) across a 4.5× change in σ_ch at Ω=96. The channel term is not
doing the work — the shallow landscape fails to hold the state on its own, with
`t_stage = 16` under one relaxation time (`1/λ_antisym = 30`). A good fit in that
block is not evidence for the mechanism, and is asserted as such in the tests.

Other caveats: one γ-family, one `t_stage`, and `p` is inferred from the decay of
`I(D)` rather than counted directly, so cells below the 1e-15 resolution of that fit
are dropped and never fitted.

### 12.1 A depth ceiling no population can raise

§11.1 found no efficiency optimum in Ω and blamed the channel-limited regime.
Re-running the frontier on the **population-limited** side (σ_ch/δ* = 0.15) shows
that diagnosis was wrong: cost per bit is still monotone (28 → 887 k_BT from Ω=4 to
96). The actual reason is simpler and more general — **information is bounded by one
bit while cost is linear in Ω**, so at depth 30 even Ω=4 already delivers 0.71 bits
and the numerator can improve 1.4× while the denominator rises 44×.

An interior optimum therefore requires a cascade **deep enough that small systems
fail**, and it does appear:

| depth | 30 | 100 | 300 | 1000 |
|---|---|---|---|---|
| optimal Ω | 4 (edge) | 4 (edge) | **10** | **12** |
| cost there | 28 | 181 | 1097 | 5061 k_BT/bit |

**And the same formula predicts a hard ceiling.** At large Ω the exponent saturates
at `δ*²/(2σ²)`, so `p` has an Ω-*independent* floor and the bit dies at a depth no
population can extend:

    D_max ~ exp(δ*² / 2σ²) / 4

Measured, as the depth at which `I` falls through 0.5:

| σ_ch/δ* | 0.45 | 0.35 | 0.28 | 0.22 |
|---|---|---|---|---|
| predicted D_max | 3.0 | 14.8 | 147 | 7664 |
| measured, Ω=64 | 9 | 44 | 355 | >4000 |
| measured, Ω=128 | 9 | 50 | 489 | >4000 |

The predicted scaling holds over the tested range (a 50× span) with a constant
prefactor ≈3 that the saddle point drops — the same missing prefactor as §12's slope
of 0.74.

**A correction, because the first version of this section overstated it.** It read
"at σ/δ* = 0.45 the two populations die at *exactly* the same depth." They do — but
that is **integer rounding across two adjacent points**, not an exact invariance, and
it breaks at Ω=256, where the integer crossing reads **10**. Interpolating the
`I = 0.5` crossing shows what is really happening:

| Ω | 16 | 32 | 64 | 128 | 256 |
|---|---|---|---|---|---|
| exponent, % of ceiling | 88.8 | 94.1 | 96.9 | 98.5 | **99.2** |
| interpolated death depth | 6.53 | 7.44 | 8.27 | 8.86 | **9.14** |
| increment per doubling | — | +0.91 | +0.83 | +0.59 | **+0.28** |

The increments shrink by roughly half each doubling (ratios 0.91, 0.71, 0.47), giving
a geometric limit of **D∞ ≈ 9.4**. The two "9"s were never the same number — at depth
9 the information is **0.474796 at Ω=64 and 0.495227 at Ω=128**, and Ω=128 sits
0.005 away from not crossing there at all, which is the entire margin by which the
integer agreed. `information.depth_at_information` now returns the float crossing so
the continuous quantity is the default and rounding happens only at display. So the ceiling is real — a **16× population change
buys 1.4× depth, and convergent** — but it is a limit approached, not a constant hit.
For contrast, in the wall regime (σ/δ* = 0.10) `p` falls **eleven orders of
magnitude** over the same population range, so there the depth grows without bound.
That contrast, not the equality of two integers, is the finding.

**How tight the ceiling is depends on the noise, and σ/δ* = 0.45 is its tightest
case.** Quoting only that row would overstate the saturation. Across an 8× population
sweep (Ω = 16 → 128): depth grows **5.4×** at σ/δ* = 0.28 and **2.2×** at 0.35,
against 1.4× (for 16×) at 0.45. Molecules do buy depth — with sharply diminishing
returns that set in earlier the noisier the channel, because `Ω× = 1/(2κσ²)` falls as
σ rises.

**This retroactively explains §10–§11.** Those experiments ran depth 30 at
σ_ch/δ* = 0.35, where the ceiling is ≈44–50 stages. They were operating near a limit
that no amount of population could move, which is why fidelity plateaued at 0.63 and
why 8× the molecules bought 0.0045. The limit was never about the chemistry's
landscape; it was about the ratio of channel noise to landscape width.

**Caveats.** `D_max` is defined by an arbitrary `I = 0.5` crossing; a different
threshold shifts the prefactor, not the scaling. One γ-family, one `t_stage`. The
σ/δ*=0.22 row is a lower bound only — the ceiling there exceeds the tested depth.

## 13. What an n-symbol landscape costs in drive — `n_winner_affinity.py`

§9.1 priced the two-symbol landscape: it dies at γ_c = 1/2, so it costs a minimum
cycle affinity **A_c = 3 ln 2**. This sweeps n.

The elementary cycle is **three reactions for every n** — fire `X_i+X_j→2B`, then
`B+X_i→2X_i`, then `B+X_j→2X_j`, and every count returns to its start — so the
affinity per cycle stays `A = −3 ln γ` and the whole question is where γ_c(n) sits.
Measured exactly, from the analytic Jacobian at the symmetric point:

| n | 2 | 4 | 8 | 32 | 128 | 256 |
|---|---|---|---|---|---|---|
| γ_c | 0.5 | 6.81e-2 | 3.94e-3 | 3.58e-5 | 4.96e-7 | **6.08e-8** |
| A_c | 2.079 | 8.061 | 16.61 | 30.71 | 43.55 | **49.85** |
| A_c / ln n | 3.000 | 5.815 | 7.987 | 8.861 | 8.976 | **8.989** |
| local exponent | — | −3.79 | −3.93 | −3.19 | −3.05 | **−3.02** |

**γ_c(n) → n⁻³, hence A_c(n) → 9 ln n.** The local exponent peaks near −4.2 around
n=6 and converges monotonically to −3.02 by n=256; the fit over n ≥ 64 gives
`γ_c ~ n^−3.042`, i.e. `A_c ≈ 9.13 ln n`.

Since a symbol carries `log₂ n` bits, that is **9 ln 2 ≈ 6.24 k_BT of drive per bit
of alphabet** — exactly 9× Landauer's `k_B T ln 2` per bit. The 9 is `3 × 3`: three
reactions in the cycle, and a critical γ suppressed by the cube.

**n = 2 is not on the asymptote**, and that matters for reading §9.1. A_c(2) = 3 ln 2
= 2.079 where 9 ln 2 = 6.238 — the ratio is 3, not 9. The famous case is the special
one, and the law is approached from below over decades in n.

**Two predictions this killed.** THEORIES Q3 guessed `A_c = n ln 2`; before running,
I predicted `γ_c = 1/n` giving `A_c = 3 ln n`, because it reproduces the exact n=2
value. Measured γ_c falls ~870× below 1/n by n=32.

**Three ways n ≥ 3 is structurally unlike n = 2**, all verified:
- The **cycle space is no longer one-dimensional**: counting each reversible pair as
  one edge (the counting under which §9.1 called AM's cycle space 1-D), the dimension
  is exactly **C(n,2)** — 1, 3, 6, 15 at n = 2, 3, 4, 6.
- The **symmetric point moves.** At n=2 it is pinned at (⅓,⅓,⅓) for every γ; for n≥3
  it is not at 1/(n+1) and it depends on γ (n=3: x = 0.2050 → 0.2431 over
  γ = 0.02 → 0.6).
- λ_breaking still reduces to §9.1's `(1−2γ)/3` at n=2, exactly.

**Verification.** The symmetric state is an exact fixed point (|rhs| < 1e-12); the
breaking vector is a true eigenvector, not merely a Rayleigh bound (relative residual
1e-16); and the eigenvalue crossing is a **real bifurcation** — integrating the ODE
from a 1e-6 perturbation gives a surviving broken state below γ_c (winner share 0.858
/ 0.953 / 0.997 at n = 3 / 4 / 8) and decay to 1e-16 above it.

**Caveats.** Deterministic and exact, so this is the Ω→∞ landscape question only —
it says nothing about the population needed to *resolve* the landscape at finite Ω,
which §1 and §12 govern. Uniform γ throughout; asymmetric rates are untested (Q4).
The n⁻³ exponent is measured to −3.02 and still drifting, so it is an approach, not a
proof.

## 14. Why the radix barrier saturates — `n_winner_reversible.py`

§3 measured that the n-winner barrier `c(n)` falls ~7× and then **saturates** at
≈0.0022, and could not say why. §13's machinery plus §2's quasipotential now
explain the saturation, though not its size.

Two exact closed forms, both verified against the engine to 7 decimals over
n = 2..64:

    λ(n)  = 1/(2n−1)                  symmetry-breaking eigenvalue at γ=0
    D₀(n) = (2n−3)/(2n−1)²            van Kampen diffusion in the same direction

`D₀(2) = 1/9` is exactly `design.md` §9's `D = 1/(9Ω)` for irreversible AM, so this
is the same reduction, generalised.

**The mechanism.** A quasipotential barrier goes like `c ∝ λ/D₀`, and here

    λ/D₀ = (2n−1)/(2n−3) → 1

**λ and D₀ vanish at the same rate**, so their ratio saturates instead of
diverging. That is why `c(n)` has a floor. The instability rate does die like
1/(2n) — but so does the noise that would exploit it.

So the prediction is `c(n) = δ²(2n−1)/(2(2n−3)) → δ²/2`:

| n | 2 | 4 | 8 | 16 | 32 | 64 |
|---|---|---|---|---|---|---|
| predicted c | 0.01500 | 0.00700 | 0.00577 | 0.00534 | 0.00516 | 0.00508 |
| measured c (§3) | 0.01567 | 0.00455 | 0.00297 | 0.00235 | 0.00227 | 0.00224 |
| ratio | **0.957** | 1.538 | 1.943 | 2.274 | 2.275 | **2.268** |

**What this does and does not settle.** At n=2 it is right to 4% — it *is* §2's
`c(ε) = (3/2)ε²`, recovered as the n=2 member. The saturation is genuinely derived.
But the predicted floor is `δ²/2 = 0.0050` against a measured 0.0022, and the ratio
climbs to **a constant 2.27** and stays there from n=16 to n=64. Predicted fall
2.95×, measured 7.0×.

A constant offset that itself saturates is a **prefactor**, not a wrong shape — the
same species as §12's slope of 0.74 and §12.1's factor ≈3, and logged with them
under THEORIES Q7. Candidates not yet separated: the n−1 competing escape
directions, which multiply the escape rate without changing the exponent; and the
transient bias amplification §2 already flags at n=2.

**Caveats.** The closed forms are for γ=0 (irreversible), which is what §3 measured;
`breaking_diffusion` computes the general-γ case numerically but it is untested
against a barrier measurement. The comparison inherits §3's fixed-margin convention
(§3.1) and its δ=0.10. The barrier is read from a 1-D reduction along one mode,
which is exact for the *rate* by symmetry but is an approximation for the escape.

## Open questions

1. **Universality class of the freeze-out transition** (§5). Is a = 0.38 really 1/3
   or 2/5, and can the quasipotential of §2 predict it? Needs error bars on the
   collapse fit. Now better motivated: §14 shows the same λ/D reduction predicts a
   saturation correctly in shape, so the technique is worth pointing at §5.
2. ~~Is the radix saturation convention-dependent?~~ **Answered** (§3.1). Yes, the
   penalty's *existence* is — but for a mundane reason that vindicates §3's choice
   rather than undermining it. Still open: symmetric plurality, which was not tested.
3. ~~A controlled EIR radix test with per-n hyperparameter tuning.~~ **Done** — the
   penalty survives tuning (§8). Still open: whether EIR's falloff shares CRNL's
   *mechanism* (basin crowding) or is partly a search effect, which needs a variant
   with early stopping disabled.
4. ~~Free-energy accounting.~~ **Done** (§9): reversible AM, exact CME. A landscape
   costs `A > 3 ln 2`; deciding costs `O(Ω)·A`; holding costs no *power*. Still open:
   whether the flat 430–470 k_BT middle range of §9.2 has a clean analytic form, and
   whether the σ peak crossing γ_c as Ω grows (§9.3) means anything.
5. **Structured (asymmetric) landscapes** — unequal rate constants, deformed basins;
   the honest bridge toward real chemistry rather than a flat symmetric democracy.
   Now unblocked: §9's reversible infrastructure is what it needs. Note the closed-form
   EP identity of §9.2 fails there, which is why `thermo.entropy_step` exists as the
   general primitive.
6. ~~The SSA half of the dissipation work.~~ **Done** (§10). Sampling confirms §9
   (⟨M⟩ to 0.86 SEM at Ω=120), and a restoring stage is priced. Still open: whether
   there is *any* population at which a stage restores near γ_c — 8× buys 0.0045 of
   fidelity at γ=0.45, so if a threshold exists it is far outside reach, and the
   honest answer may be that there is none.
7. ~~A comparator that needs no convention.~~ **Done** (§11): k_BT per bit
   delivered needs no control at all, and §11.1 shows there is no efficiency optimum
   in Ω — the ratio is minimized by a system that barely transmits, so the
   information must be fixed first. Still open: the frontier's marginal cost rises
   77× over the tested range and is still climbing at Ω=120, so **where it ends is
   unknown**; and the measure remains convention-dependent through `noise_frac`, D,
   and the 0.05-bit floor.
8. ~~Is `t_stage` hiding anything?~~ **Mostly answered** (THEORIES T6): the optimal
   `t_stage` is flat in absolute time across γ ≤ 0.30, not proportional to the
   relaxation time, so the fixed value used in §10–§12 does not smuggle in a
   γ-dependent effort. Still true that cost grows linearly in t while fidelity
   saturates, so "cost per stage" must always be quoted with its t.
