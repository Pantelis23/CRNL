# CRNL — findings

Measured results, with the numbers, the caveats, and what is still open. Raw data
for every table is committed under `results/`; each figure is regenerable by the
named experiment.

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
conserved simplex (7381 states at Ω=120, 0.20 s), not sampled. A direct stochastic
measurement of one rare flip at Ω=120 would take hundreds of hours.

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
| 120 | τ=5.5e5, σ=4.93 | τ=6.1e3, σ=7.67 | τ=102, σ=11.5 |

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

## Open questions

1. **Universality class of the freeze-out transition** (§5). Is a = 0.38 really 1/3
   or 2/5, and can the quasipotential of §2 predict it? Needs error bars on the
   collapse fit.
2. **Is the radix saturation convention-dependent?** (§3) Predicted yes — fixed
   champion-share or symmetric-plurality framings should behave differently.
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
6. **The SSA half of the dissipation work** (Plan 2): instrument the Gillespie loop
   with an integer cycle counter and cross-check §9's exact numbers by sampling, plus
   the cost of a *restoring cascade stage* against the non-restoring baseline of §7 —
   the one measurement that prices restoration as this project defines it.
