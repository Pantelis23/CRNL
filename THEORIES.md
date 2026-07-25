# CRNL — theories and questions

The speculative companion to [`FINDINGS.md`](FINDINGS.md). Findings holds what has
been *measured*; this holds what we think might be true, what we don't understand,
and what we believed and were wrong about.

**Rules, so this stays useful rather than becoming a wish list:**

1. A **theory** must make a prediction that could fail. "Restoration is
   interesting" is not a theory. "`c(γ) ∝ (γ_c − γ)` near the bifurcation" is.
2. Every theory carries a **status** and a **how to kill it**.
3. When a theory is measured, it moves to `FINDINGS.md` and leaves a one-line
   stub here pointing at the section.
4. **Disproven theories stay**, with what killed them. This project has produced
   four confident wrong results so far; the record of *how* they were wrong is
   worth more than the theories that survived.
5. Numbers quoted here must come from a run. This file is the place for
   speculation about *mechanism*, not for unverified arithmetic — every wrong
   claim in §4 below began as a number written without running it.

---

## 1. Live theories

### T1. The cascade is channel-limited, not population-limited
**Status: partly confirmed, mechanism identified, crossover not yet mapped.**

`FINDINGS.md` §11.1 found that cost per bit *rises* with Ω and that the marginal
cost of information explodes. The proposed reason: once Ω is large enough that
finite-count noise is subdominant to the injected channel noise `σ_ch = 0.35·δ*`,
**more molecules buy essentially nothing**, while cost stays extensive in Ω.

First check (γ=0.05, t=16, per-stage flip probability from the decay of I(D)):

| σ_ch/δ* | Ω=30 | Ω=60 | Ω=120 | p(120)/p(30) |
|---|---|---|---|---|
| 0.20 | 5.1e-5 | 1.3e-5 | 4.3e-6 | **0.08** |
| 0.28 | 1.12e-3 | 7.0e-4 | 5.2e-4 | 0.47 |
| 0.35 | 4.93e-3 | 3.96e-3 | 3.48e-3 | 0.71 |
| 0.45 | 1.78e-2 | 1.61e-2 | 1.53e-2 | **0.86** |

Across a row (noise) p moves ~350×; down a column (population) it moves 1.2–2×,
and the population dependence *weakens* as noise grows. So there is a **crossover
from population-limited to channel-limited**, and the default protocol sits mostly
on the channel-limited side — which is why §11.1's frontier saturates.

**Prediction:** the crossover sits where finite-count noise `~1/√Ω` matches
`σ_ch/δ*`, i.e. at `Ω_× ~ (δ*/σ_ch)²`. At `σ_ch/δ* = 0.35` that is `Ω_× ≈ 8`;
at 0.20 it is `Ω_× ≈ 25`, which is consistent with the 0.20 row still showing
strong Ω dependence at Ω=30 while the 0.45 row does not.
**How to kill it:** map p(Ω) finely at fixed `σ_ch/δ*` and show the knee is *not*
near `(δ*/σ_ch)²`, or that p keeps falling with Ω on the supposedly saturated side.

### T2. The freeze-out exponent is predictable from the quasipotential
**Status: untested. The most valuable open theoretical target.**

§5 measured `a ≈ 0.38` for the collapse `D(H,Ω) = F((H−H_c)·Ω^a)`, sitting
awkwardly between 1/3 and 2/5 with no error bars and no identified universality
class. §2 derived `c(ε) = (3/2)ε²` from the saddle geometry — λ=⅓ and diffusion
`D = 1/(9Ω)`.

**Conjecture:** freeze-out is the competition between the decision rate near the
saddle and the dilution rate H, so the same two ingredients should fix `a`.
Dimensionally, if the escape time from the saddle goes like `Ω^{1/2}` in the
critical region and the dilution acts on `1/H`, an exponent of 1/3 or 2/5 falls
out of matching them — but *which* depends on how the barrier scales with the
frozen composition, which has not been worked out.
**How to kill it:** derive a specific exponent and measure a ≠ that, with error
bars. First step is error bars on the collapse fit, which §5 flags as missing.

### T3. Radix saturation is an artifact of the fixed-margin convention
**Status: predicted, untested. Cheap to run.**

§3 found `c(n)` saturates at ≈0.0022 past n≈16 under a **fixed pairwise margin**
δ. The stated reason: the champion's share converges to δ as n grows, so the
contest stops changing.

**Prediction:** under a **fixed champion share** (champion holds a constant
fraction of the whole population, rivals split the rest) the penalty should be
**unbounded** — `c(n) → 0` with no floor, because the champion's per-rival lead
keeps shrinking. Under symmetric plurality it should differ again.
**How to kill it:** run `radix_wall.py` with a share-based convention and find the
same floor ≈0.0022.

### T4. The restoration barrier vanishes quadratically in the landscape width
**Status: predicted, untested.**

§2 has `c(ε) = (3/2)ε²` where ε is the bias. §9.1 has the landscape width
`δ*(γ) ∝ √(γ_c − γ)`. If the barrier is governed by the bias *in units of the
landscape*, then at fixed relative bias

    c(γ) ∝ δ*(γ)² ∝ (γ_c − γ)

**Prediction:** the effective restoration barrier falls **linearly** in `(γ_c − γ)`,
so the population needed for a given reliability diverges as `1/(γ_c − γ)`. This
would explain §10.2's γ=0.45 result (Ω=240 buys 0.0045 of fidelity) without any
"minimum Ω" — there is no threshold, just a diverging requirement.
**How to kill it:** measure c(γ) at several γ near γ_c and find a different power,
or a genuine threshold in Ω.

### T5. The flat 430–470 k_BT middle range has an analytic form
**Status: open, no candidate expression.**

§9.2: across γ ∈ [0.15, 0.40] the cost of a decision sits flat at 430–470 k_BT
while the error varies 25×. The stated mechanism is that affinity per cycle
(`ln(1/γ)`, falling) and the number of cycles required (rising) nearly cancel.
"Nearly" is doing a lot of work — a near-cancellation over a 2.7× range of γ is
suspicious enough to have a reason.
**How to kill it:** show the flatness is coincidental by finding a protocol
(different θ/δ* or bias fraction) where it disappears.

### T6. The optimal stage time is set by the relaxation time
**Status: suggestive, one data point.**

§11's frontier selects `t* = 16` at essentially every Ω and γ tested. The natural
scale is the landscape's own relaxation time `1/|λ_antisym(γ)| = 3/(1−2γ)`, which
is 3.33 at γ=0.05 — so `t* ≈ 4.8` relaxation times.
**Prediction:** `t*` should track `1/|λ_antisym(γ)|` across γ, i.e. grow and
diverge as γ→γ_c, rather than staying pinned at 16.
**How to kill it:** sweep t_stage per γ and find t* flat in absolute time.

---

## 2. Open questions with no theory yet

- **Q1. Where does the efficiency frontier end?** §11.1's marginal cost rises 77×
  and is *still climbing* at Ω=120. Is there an asymptote, a divergence at finite
  information, or does it continue indefinitely? Blocked on Ω>120 cost, which the
  chunked augmented generator partly relieves.
- **Q2. Does σ's peak crossing γ_c mean anything?** §9.3: at Ω=30 the stationary
  dissipation rate peaks at γ=0.45; at Ω=60 and 120 it is still rising at γ=0.49.
  An Ω-dependent peak location that crosses the bifurcation point is either a real
  finite-size effect or a coincidence, and we cannot currently tell.
- **Q3. n-winner reversible thermodynamics.** The engine supports it for free and
  no experiment uses it. Does the affinity floor `A > 3 ln 2` become `A > n ln 2`,
  or something else? This is the cheapest genuinely new territory available.
- **Q4. Structured (asymmetric) landscapes.** Unequal rate constants, deformed
  basins — the honest bridge toward real chemistry. Note §9.2's closed-form EP
  identity **fails** there, which is exactly why `thermo.entropy_step` exists as
  the general primitive.
- **Q5. Does EIR's decode falloff share CRNL's mechanism?** §8 showed the radix
  penalty survives per-n hyperparameter tuning, so it is structural. Whether it is
  *basin crowding* (CRNL's mechanism) or partly a search effect needs a variant
  with early stopping disabled.
- **Q6. Is the 0.35 noise fraction hiding a regime?** Every cascade result uses
  `σ_ch = 0.35·δ*`. T1 says this sits on the channel-limited side. The
  population-limited side (σ_ch/δ* ≲ 0.2) is essentially unexplored and is where
  the restoration wall of §1 should actually govern.

---

## 3. Retired — measured and moved

- ~~Does a landscape exist, and at what drive?~~ → §9.1, `γ_c = 1/2`, `A = 3 ln 2`.
- ~~What does a decision cost?~~ → §9.2.
- ~~What drive does remembering require?~~ → §9.3.
- ~~Can restoration be priced per bit without a comparator?~~ → §11.

---

## 4. Disproven — kept deliberately

Four confident, plausible, **wrong** results. Every one of them looked like
physics, and every one was a protocol artifact of the same family: *something that
scales with the landscape was held fixed while the landscape shrank.*

| Claim | What killed it |
|---|---|
| **"Restoration requires a minimum Ω as well as a minimum affinity."** | The passive control's dynamic range was absolute (±1) while its noise was scaled by δ*(γ). Rail the control to ±δ*(γ) and it becomes γ-independent (spread 0.003 vs 0.154); the crossover vanishes. Identical chemistry arm in both. §10.3. |
| **"Dissipation has a minimum near γ≈0.3"** — a clean U-shaped curve. | The decision threshold was held fixed while δ*(γ) shrank, so above γ≈0.42 the threshold sat *outside* the landscape (θ/δ* = 1.88 at γ=0.49) and "deciding" became "fluctuating past the attractor". Scaling θ with δ* makes the curve monotone. §9.2. |
| **"High-n AM collapses to all-blank."** | Simply ran it: single-winner probability ≈1.0 at every n from 2 to 48, under two conventions. The cost of radix is consensus *time*, not collapse. §3. |
| **"A direct SSA measurement at Ω=120 would take hundreds of hours."** | Measured throughput: 84,000 steps/s and ~0.4·Ω steps per unit time make one flip at γ=0.35 cost **5.5 minutes**. The hundreds-of-hours regime is γ ≤ 0.30 — exactly where the CME's own solve is rejected, so neither instrument reaches it. §9. |
| **"τ_SSA/τ_CME ≈ 0.8 from an arm-vs-attractor offset."** | Seed-averaged measurement gives 0.97 with no offset visible. The prediction conflated an MFPT-from-the-arm with the mean time between crossings of a long stationary trajectory, which is dominated by the dwell near the attractor. §10.1. |
| **"The cascade's decay length ξ grows like the restoration wall e^{cΩ}."** | Fit directly: ln ξ vs Ω has R²=0.69 and ξ saturates (24→76 for 15× Ω). Superseded by **T1**, which explains why — the cascade is channel-limited, so ξ cannot inherit the finite-count scaling. |

**The pattern worth remembering.** In every case the error was invisible to the
guard that was supposed to catch it, because the guard watched the wrong quantity:
`θ/δ*` was correctly constant while the *input* relative to θ was not; the control's
rails were never checked against the chemistry's. And in every case the claim was
first written down as a number that had not been run.
