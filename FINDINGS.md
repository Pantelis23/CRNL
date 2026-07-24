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
- **Decode reliability falls with option count:** large-budget ceiling 0.90 (n=2) →
  0.535 (n=8) → 0.223 (n=32), staying above the 1/n random baseline but by a
  shrinking margin.

**Caveat — the radix half is suggestive, not established.** EIR was run with
*fixed* hyperparameters (`noise_start=0.15`, `lr=0.1`, `spsa_samples=6`) at every n.
A fair test would tune per n; some of that falloff is untuned hyperparameters rather
than intrinsic basin crowding. The freeze-out half (quality vs. budget at fixed n)
is not subject to that objection.

---

## Open questions

1. **Universality class of the freeze-out transition** (§5). Is a = 0.38 really 1/3
   or 2/5, and can the quasipotential of §2 predict it? Needs error bars on the
   collapse fit.
2. **Is the radix saturation convention-dependent?** (§3) Predicted yes — fixed
   champion-share or symmetric-plurality framings should behave differently.
3. **A controlled EIR radix test** (§8) with per-n hyperparameter tuning.
4. **Free-energy accounting.** Not started. AM as written is irreversible, so its
   dissipation is formally infinite; measuring the cost of restoration requires
   rebuilding it as a reversible CRN with finite ΔG. This is the design doc's §8
   summit and the natural next phase.
5. **Structured (asymmetric) landscapes** — unequal rate constants, deformed basins;
   the honest bridge toward real chemistry rather than a flat symmetric democracy.
