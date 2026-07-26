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
**Status: CONFIRMED and superseded by a formula → `FINDINGS.md` §12.**

Confirmed, and then sharpened into something better than the original claim. A
saddle point over where a flip happens gives one parameter-free expression for
both regimes:

    −ln p ≈ κ·Ω·δ*² / (1 + 2κΩσ²),    κ(γ) = (3/2)(1−2γ) = (9/2)·λ_antisym(γ)

Pooled collapse over 216 cells: **R² = 0.933**. The population-limited side, which
nothing had explored, is dramatic — 11 orders of magnitude in p across Ω at
σ_ch/δ* = 0.10. §11's protocol sat on the channel-limited side, which is *why* it
found no efficiency optimum in Ω.

**What is still open here (T1a):** the fitted slope is 0.74 pooled and drifts with
γ (0.80 / 0.63 / 0.42 / 0.50), where an exact saddle point would give 1. The
missing piece is the prefactor and the Gaussian-tail correction, neither of which
was worked out. A derivation that predicts the slope — or shows why it should
drift with γ — would turn a collapse into a law.
**How to kill it:** find a (γ, σ, Ω) region where p departs from the formula by
more than the prefactor could explain, i.e. where the *shape* in Ω is wrong rather
than the scale.

### T2. The freeze-out exponent is predictable from the quasipotential
**Status: ANSWERED — and the answer dissolves the question → `FINDINGS.md` §5.1.**

T2 conjectured that §2's saddle ingredients (λ = ⅓, `D_δ = 1/(9Ω)`) should fix
§5's collapse exponent `a ≈ 0.38`, and called this "the most valuable open
theoretical target". The conjecture's *premise* was right — freeze-out is exactly
the competition between the saddle escape and the dilution — and its *form* was
wrong. Worked out, those two ingredients say there is no exponent:

- The expanding SSA is **exactly** ordinary SSA stopped at internal time `τ = 1/H`
  (`crnl/freezeout.py`, verified bit-for-bit). Freeze-out is not its own
  dynamics; `H*(Ω) = 1/τ*(Ω)` with `τ*` the AM consensus time.
- `design.md` §9 already gives the effective seed at the saddle,
  `σ² = D_δ/λ = 1/(3Ω)`. From an **exactly symmetric** start that shot noise is
  the only seed, so `τ* = (1/λ)·ln(1/σ) + O(1) = (3/2)·ln Ω + O(1)`.
- Hence `Hc = 0`. No critical point, no universality class, and `a` is a
  parameter of the wrong functional form.

T2 asked whether `a` is 1/3 or 2/5. It is neither. FINDINGS open question 1 is
**void**, not answered — there was nothing there to have a universality class.

**T2a, what is actually left:** the intercept `B` in `1/H* = (3/2)lnΩ + B`. The
slope is parameter-free and measured; `B` is not derived. Estimating it from §9's
`σ` plus "decide when δ = ⅓" lands ~1 unit of τ away, because the decision level
is reached in the *nonlinear* phase where the threshold picture fails. Small, and
it needs the shape of the deterministic sigmoid rather than a threshold.
**How to kill T2a:** derive `B` and measure a different one.

### T3. Radix saturation is an artifact of the fixed-margin convention
**Status: TESTED → `FINDINGS.md` §3.1. The prediction was wrong; the conclusion is
the opposite of what T3 expected.**

T3 predicted that a fixed champion *share* would give an **unbounded** penalty,
"because the champion's per-rival lead keeps shrinking." **The lead is
`s − (1−s)/(n−1)`, which grows toward `s`** — 0.100 → 0.525 over n=2..24, a 5.2×
increase. The reasoning was wrong on paper, before any code ran.

Measured, the penalty under fixed share does not become unbounded; it **vanishes**
(P(win) = 1.0000 at every n ≥ 3). But that is because fixed share is asking an
easier question at every n, not a different reading of the same one. **So §3's fixed
pairwise margin is the convention that isolates alphabet size, and §3 stands.**

Nor is share the governing variable: at a fixed share of 0.50, P(win) is 0.606 at
n=2 and 0.997 at n=3.

**Still open:** symmetric plurality, untested; and whether the saturation *floor* of
§3 moves under any other convention that genuinely holds the pairwise lead fixed.

### T4. The restoration barrier vanishes quadratically in the landscape width
**Status: PARTLY ANSWERED by §12 — and the answer was not what T4 guessed.**

§12 measures the barrier coefficient as `κ(γ)·δ*²` with `κ(γ) = (3/2)(1−2γ)`. So the
barrier carries **two** vanishing factors as γ→γ_c, not one: `δ*² ∝ (γ_c−γ)` as T4
guessed, *times* `κ ∝ (1−2γ) ∝ (γ_c−γ)`. The barrier therefore falls like
**(γ_c−γ)²**, and the population needed for fixed reliability diverges like
`1/(γ_c−γ)²` — faster than T4 predicted, and enough to explain §10.2's γ=0.45 result
(Ω=240 buying 0.0045 of fidelity) with no threshold anywhere. Still untested
directly: the quadratic form has been inferred from the collapse, not measured by
sweeping γ finely near γ_c.

**How to kill it:** sweep γ finely in [0.40, 0.499] and fit the barrier's exponent
in `(γ_c − γ)`. A measured exponent of 1 would restore T4's original guess and
falsify the κ correction; anything other than 2 means `κ ∝ λ_antisym` is wrong even
though it improves the §12 collapse.

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
**Status: DEAD. Probed and falsified.**

T6 predicted `t*` would track `1/λ_antisym(γ) = 3/(1−2γ)` and so grow ~9× across
γ = 0.05 → 0.45. Measured (Ω=30, depth 30, σ/δ*=0.35, minimising k_BT per bit):

| γ | 0.05 | 0.15 | 0.30 | 0.45 |
|---|---|---|---|---|
| 1/λ_antisym | 3.33 | 4.29 | 7.50 | 30.0 |
| t* | 16 | 16 | 16 | 1 |
| t*/relax | 4.80 | 3.73 | 2.13 | 0.03 |

`t*` is **flat in absolute time** over γ ≤ 0.30, not proportional to the relaxation
time. (The collapse to t*=1 at γ=0.45 is the do-nothing degeneracy of §11 — where
the chemistry cannot hold the bit, the cheapest stage is one that barely runs — not
a relaxation effect.) So the fixed `t_stage` used across §10–§12 is defensible, and
FINDINGS open question 8 is answered in the same stroke for γ ≤ 0.30.

### T7. The n-winner barrier follows the symmetry-breaking eigenvalue
**Status: DEAD in the simple form, and the failure sharpens §3.**

§12 found `κ(γ) = (9/2)·λ_antisym(γ)` — the restoration barrier proportional to the
symmetry-breaking eigenvalue. §13 computes that eigenvalue for any n, so the hope
was that it would predict §3's measured `c(n)` and explain its unexplained
saturation. It does not:

| n | 2 | 4 | 8 | 16 | 64 |
|---|---|---|---|---|---|
| λ_breaking(n, 0) | 1/3 | 1/7 | 1/15 | 1/31 | 1/127 |
| λ / λ(2) | 1.000 | 0.429 | 0.200 | 0.097 | **0.024** |
| c(n) / c(2) | 1.000 | 0.290 | 0.190 | 0.150 | **0.143** |

The ratio of the two drifts 1.0 → 6.05. **λ_breaking vanishes like 1/(2n−1) while
c(n) saturates**, so §3's saturation is *not* a linear-stability effect.

Worth keeping from the attempt: **λ_breaking(n, γ=0) = 1/(2n−1) exactly** (checked
against ten values of n), and the symmetric state at γ=0 is
`x = 1/(2n−1)`, `b = (n−1)/(2n−1)`.

**T7b: λ_breaking DOES set a barrier — just not this one → `FINDINGS.md` §6.1.**
T7 failed because §3's `c(n)` comes from a *biased* start with a fixed pairwise
margin, where the seed does not shrink and λ is not the bottleneck. From a
**symmetric** start it is: the freeze-out decision time obeys
`dτ*/dlnΩ = 1/(2λ(n)) = (2n−1)/2`, measured 1.492 ± 0.034 / 2.471 ± 0.040 /
3.373 ± 0.041 / 4.921 ± 0.118 at n = 2 / 3 / 4 / 6. So the same `λ(n) = 1/(2n−1)`
that fails to explain `c(n)` exactly governs the freeze-out penalty, and that
penalty is therefore **unbounded in n** where `c(n)` saturates. Two different
questions, and which one λ answers depends entirely on whether the initial
asymmetry shrinks with Ω.

**T7a: CONFIRMED for the shape → `FINDINGS.md` §14.** `D₀(n) = (2n−3)/(2n−1)²`
exactly, so `λ/D₀ = (2n−1)/(2n−3) → 1` — λ and the diffusion vanish at the same
rate and their ratio saturates. **That derives §3's saturation**, unexplained since
it was measured. `D₀(2) = 1/9` recovers `design.md` §9, and the n=2 barrier comes
out at `1.5 δ²`, i.e. §2's result, to 4%.

**Not settled: the size — and it is a shape problem, not a size problem.**
Predicted floor `δ²/2 = 0.0050` against a measured 0.0022. This paragraph used to
read: *"the ratio climbs to a constant 2.27 and holds from n=16 to n=64. A constant
offset is a prefactor, so this folds into Q7."* **Both halves were wrong.** The
ratio is constant only along n; measured at two further δ it falls ~40%
(2.274 → 1.364 at n=16), because the measured exponent is `c ∝ δ^~2.5` rather than
`δ²` (§14.1). So there is no prefactor to name, Q7 is dissolved, and the open
question is why the exponent in δ exceeds 2 — that is **Q9a**. Two candidates for
the *shape* remain unseparated: the n−1 competing escape directions, and the
transient bias amplification §2 already flags at n=2.

---

## 2. Open questions with no theory yet

- **Q1. Where does the efficiency frontier end?** §11.1's marginal cost rises 77×
  and is *still climbing* at Ω=120. §12.1 partly answers it — the frontier is
  bounded above by the depth ceiling, since past `D_max` no Ω delivers the bit at
  all — but the shape of the approach is unmeasured.
- ~~**Q7. THE PREFACTOR. Three measured-but-underived numbers are one question.**~~
  **DISSOLVED. Every member has now been checked along a second axis and none of
  them is a prefactor.** This entry was once "the hub of the whole open list"; it
  consolidated by accident, and it was wrong about all three members.

  | where | saddle prediction | measured | what the residual actually is |
  |---|---|---|---|
  | §12 collapse | exponent × 1 | exponent × **0.74** | exponent scale, −26% |
  | §12.1 depth ceiling | `exp(1·δ*²/2σ²)/4` | `exp(`**`1.0695`**`·δ*²/2σ²)·0.663` | exponent scale, +7% |
  | §14 radix floor | `c ∝ δ²` | `c ∝ δ^~2.5` | exponent in δ (Q9a) |

  **None of the three is an amplitude.** All three are multiplicative errors inside
  an *exponent*, which is exactly why each looked constant along whichever axis was
  swept first and then drifted along a second one: forcing `k = 1` turns an
  exponent error into an apparent prefactor `exp((k−1)·x)`, which is flat only over
  a short range in `x`. §14's "constant" 2.27 ran to 1.36 across δ; §12.1's "≈3"
  runs 3.07 → 4.05 across σ, and freeing its exponent fits at R² = 0.9998 with 2.6%
  residuals. A Laplace correction — the one technique this entry was built around —
  supplies an amplitude and would not have fixed any of them.

  **The obvious next move, and why it is not being made.** It is tempting to say
  "then they are still one question: the saddle-point *exponent* is inexact." That
  may be true. But the three differ in sign and in size (−26%, +7%, and a shift in
  a different variable's exponent), and merging them on that basis would be the
  **third** consolidation in this file resting on the same thin evidence — after
  §14's ejection and now §12.1's. It is logged as a conjecture with three points,
  not a finding, and it needs a derivation rather than another table.

  **The original reasoning, kept because the error is instructive:**

  > All three come from the same move: a **saddle point keeps only the exponent**
  > and throws away the Gaussian fluctuations around it. In every case the *shape*
  > is right — the collapse holds at R²=0.93, the ceiling scales correctly over 50×,
  > the saturation is derived — and only the amplitude is missing. §14's is the
  > cleanest target because its offset is *constant* (2.274 / 2.271 / 2.275 / 2.268
  > at n=16/24/32/64), so it is a pure number waiting to be named rather than a
  > drifting discrepancy.
  >
  > **Why it is worth doing:** a Laplace correction is one technique, and if it
  > supplies either of the two it plausibly supplies both.

  The premise "the shape is right and only the amplitude is missing" was the whole
  argument, and it was never tested — it was *inferred* from the numbers looking
  constant. Note also that the entry nominated as its "cleanest target" the member
  that turned out to be wrongest.

  > **This entry used to claim a third beneficiary and it is gone.** It read
  > "the same machinery is what **T2** (the freeze-out exponent, a≈0.38 between
  > 1/3 and 2/5) has been blocked on since §5. Four open items, one method."
  > T2 is now answered (§5.1) and there is no exponent to supply a prefactor
  > for — `a` was a parameter of the wrong functional form. So the cluster is
  > down to **two** members and one of them (§12.1's ceiling factor) still has
  > not been checked across σ, which is the test §14 failed.

  **§14 was ejected from this cluster by its own kill test** → §14.1. The offset
  was measured at two further δ and falls ~40% (2.274 → 1.364 at n=16), so it is
  not a prefactor at all: `c ∝ δ²` is predicted, and the measured exponent is
  **2.27 at n=8 and 2.48 at n=16** against §2's 2.08 at n=2. The consolidation
  above was written down two commits before the test that broke it. (That left the
  table with two members; §12.1 then failed the same test, and the table has none.)

  **This is how the entry died, and it called its own shot.** The paragraph above
  used to end: *"§12's slope and §12.1's ceiling factor should be independent of the
  channel noise σ if they are prefactors. Neither has been checked across σ, and
  §14 is a warning that a 'constant' measured along one axis need not be constant
  along another."* §12.1's factor was then checked across σ and drifted 1.32×. The
  test was named correctly and the answer was still assumed for two commits before
  anyone ran it.

  **What survives.** §12's slope 0.74 is unexplained and is the only member never
  tested across σ — but it is now known to be the same *species* of defect (an
  exponent scale), so testing it as a prefactor would be testing the wrong thing.
  The live question is no longer "what amplitude is missing" but **"why is the
  saddle-point exponent off by 26% here and 7% there"** — see Q9a, which asks the
  same thing in δ.
- ~~**Q9. Why does the barrier exponent grow with alphabet size?**~~ **Measured at
  n = 32 and 64 → §14.1. It does not keep growing.** p = 2.08 / 2.27 / 2.48 / 2.53 /
  2.40 at n = 2 / 8 / 16 / 32 / 64; n=32 and n=64 differ by 0.45σ. Two things came
  out of it:

  **(a) p > 2 is solid** (3.9σ over the n ≥ 16 mean), so `c ∝ δ²` really does fail
  for large alphabets and §14's residual is a shape problem, not a prefactor.

  **(b) The *climb* was over-read.** Propagating the ±8% scatter gives ±0.21 per
  exponent, making the 2.08 → 2.47 separation only 1.9σ. §14.1 had presented the
  sequence as a clean trend before anyone propagated its uncertainty — the second
  time in this file that a sequence was read as a curve without error bars (see
  the depth-ceiling entry in §4).

  **Still open (Q9a):** p saturating is *probably* the same fact as c(n)
  saturating — §3's explanation (past n≈16 the champion's share has converged to δ
  and the field is fragmented, so the contest stops changing) predicts that nothing
  about the escape keeps evolving. That is a plausible unification, not a measured
  one. It predicts the *plateau value* ≈2.47 should be derivable from the same
  fragmented-field limit; nothing derives it yet.

  **Still open (Q9b), and it is the load-bearing one → §14.2.** Both the physical
  explanation for p > 2 and the suspicion that p is an artifact point at the *same*
  quantity: **molecules per rival**, which is `(1−δ)Ω/(n−1)` and falls to ~1 at the
  largest δ. Physics reading: van Kampen fails for few-molecule species, so p > 2
  is real. Artifact reading: `c` is an Ω-slope, each δ was measured in a different
  Ω band, and at n=64 those bands span 1–2 versus 4–9 molecules per rival, so the
  ratio between them is not an exponent.

  A linearity test (n=32, δ=0.18, Ω=150…540) found **no curvature** — F = 0.01 on
  1,5 dof — so `c` is well defined over **4.0–14.3** molecules per rival. The
  suspect band (1.0–2.3) is below that and untested, and is hard to test: reaching
  10 molecules per rival at δ=0.24, n=64 needs Ω≈800, an error rate ~1e-6, and ~1e7
  trials per point.

  **The cheap design I proposed here does not exist — retracted.** It was to hold
  molecules-per-rival fixed while varying δ, via Ω ∝ (n−1)/(1−δ). But `c` is
  *defined* as `−∂ln P/∂Ω` at fixed (n, δ), and `m` is proportional to Ω, so m
  varies along that derivative by construction. **You cannot differentiate in Ω
  while holding something proportional to Ω fixed.** The design is incoherent, not
  merely expensive.

  **And the claim is weaker than §14.1 states.** Refitting p from only the two δ
  whose bands sit in the verified-safe range gives **2.40 ± 0.32 — 1.2σ above
  quadratic, not significant.** The 3.9σ comes from the lever arm in ln δ, and the
  lever arm is the δ=0.24 point, i.e. the suspect one. **p > 2 is unproven.**

  **Parked deliberately.** Settling it needs a different observable than an
  Ω-slope, or ~1e7 trials per point. That is a lot of effort for a second-order
  correction to a barrier coefficient, when Q4 (asymmetric landscapes) is
  completely untouched and asks a first-order question. Anyone resuming this
  should start by finding an observable that is not a derivative in Ω.
- **Q8. Does the depth ceiling survive a better code?** `D_max ~ exp(δ*²/2σ²)` is
  for a bare repetition of one bit through one restoring stage per hop — the weakest
  possible code. Whether the ceiling belongs to the *chemistry* or the *encoding* is
  **still open**, and an attempt to answer it failed instructively:

  **Attempt 1 (rejected — it measured nothing).** The idea was R parallel vessels
  whose outputs are pooled, so channel noise averages down by `√R` and the ceiling
  should become `D_max^R`. Measured, it looked spectacular: at a fixed budget of 256
  molecules, one vessel reaches depth 9.14 while four vessels reach >3000.

  **Why it does not count.** No parallel vessels were ever modelled. "R vessels" was
  implemented as *dividing `noise_frac` by `√R` by hand*, so the result is the
  ceiling formula restated — `D_max` depends on `σ/δ*`, therefore reducing `σ` raises
  it — and not evidence about parallelism at all. **That, alone, is the reason it is
  rejected.**

  **And then I over-corrected.** The first version of this entry justified the
  rejection with a second claim: that "depth responds to noise and essentially not to
  molecule count," citing 355 → 489 (1.4×) for Ω 64→128. That is the *flattest step
  at the most saturated noise level*, chosen after I had already decided to reject.
  Measured properly across an 8× population sweep:

  | σ_ch/δ* | Ω = 16 → 128 (8× molecules) |
  |---|---|
  | 0.28 | 90.9 → 488.6 = **5.4×** |
  | 0.35 | 22.6 → 49.8 = **2.2×** |

  Molecule count moves the depth substantially; it saturates, and how quickly depends
  on where Ω sits relative to `Ω× = 1/(2κσ²)`. The rejection stands on the modelling
  grounds above, not on this — a cherry-picked step should not have been used to prop
  it up, and doubt applies to corrections as well as to claims.

  Worse, the setup quietly assumed a **free, perfect pooling operation** — and a
  pooler is itself a restoring element. That is exactly the error that killed Part C
  design 1, where a free `sign()` in the harness did all the restoring while the
  chemistry was decoration.

  **What a real test needs:** R vessel distributions propagated independently, each
  with its *own* channel draw; an explicitly modelled combining step with its own
  dissipation and its own noise; and a statement of whether the channel noise is
  independent per vessel or common-mode — because if it is common-mode, averaging
  buys nothing and the whole idea collapses.
- **Q2. Does σ's peak crossing γ_c mean anything?** §9.3: at Ω=30 the stationary
  dissipation rate peaks at γ=0.45; at Ω=60 and 120 it is still rising at γ=0.49.
  An Ω-dependent peak location that crosses the bifurcation point is either a real
  finite-size effect or a coincidence, and we cannot currently tell.
- ~~**Q3. n-winner reversible thermodynamics.**~~ **Answered** → §13. Not `n ln 2`
  and not the `3 ln n` I predicted before running: **γ_c(n) → n⁻³, so
  A_c(n) → 9 ln n**, i.e. 9 ln 2 ≈ 6.24 k_BT per bit of alphabet, exactly 9×
  Landauer. n=2 sits off that asymptote (ratio 3, not 9), so §9.1's famous case is
  the special one. **New (Q3a): why the cube?** The 9 factors as 3×3 — three
  reactions per cycle, and a γ_c suppressed by n³ — but nothing derives the cube.
  A derivation would turn a measured exponent (−3.02, still drifting) into a law.
- **Q4. Structured (asymmetric) landscapes.** Unequal rate constants, deformed
  basins — the honest bridge toward real chemistry. Note §9.2's closed-form EP
  identity **fails** there, which is exactly why `thermo.entropy_step` exists as
  the general primitive.
- **Q5. Does EIR's decode falloff share CRNL's mechanism?** §8 showed the radix
  penalty survives per-n hyperparameter tuning, so it is structural. Whether it is
  *basin crowding* (CRNL's mechanism) or partly a search effect needs a variant
  with early stopping disabled.
- ~~**Q6. Is the 0.35 noise fraction hiding a regime?**~~ **Yes** → §12.
- ~~**Q6a. Does the efficiency optimum appear on the wall side?**~~ **Answered, and
  my prediction was wrong** → §12.1. It does *not* appear at depth 30 even at
  σ_ch/δ* = 0.15: cost per bit stays monotone (28 → 887 k_BT). The reason I gave for
  §11.1's result — channel-limited saturation — was wrong. The real reason is that
  **information is bounded by one bit while cost is linear in Ω**, so an optimum
  needs a cascade deep enough that small systems fail. It duly appears at depth 300
  (Ω\*=10) and 1000 (Ω\*=12).

---

## 3. Retired — measured and moved

- ~~Does a landscape exist, and at what drive?~~ → §9.1, `γ_c = 1/2`, `A = 3 ln 2`.
- ~~What does a decision cost?~~ → §9.2.
- ~~What drive does remembering require?~~ → §9.3.
- ~~Can restoration be priced per bit without a comparator?~~ → §11.
- ~~Is the cascade channel-limited?~~ → §12, and with a parameter-free formula
  covering both regimes.

---

## 4. Disproven — kept deliberately

Confident, plausible, **wrong** results, kept with what killed them. Most were a
protocol artifact of one family — *something that scales with the landscape was held
fixed while the landscape shrank* — and the newest is a different family: a fitted
functional form that was never derived and never compared against a rival.

| Claim | What killed it |
|---|---|
| **"Restoration requires a minimum Ω as well as a minimum affinity."** | The passive control's dynamic range was absolute (±1) while its noise was scaled by δ*(γ). Rail the control to ±δ*(γ) and it becomes γ-independent (spread 0.003 vs 0.154); the crossover vanishes. Identical chemistry arm in both. §10.3. |
| **"Dissipation has a minimum near γ≈0.3"** — a clean U-shaped curve. | The decision threshold was held fixed while δ*(γ) shrank, so above γ≈0.42 the threshold sat *outside* the landscape (θ/δ* = 1.88 at γ=0.49) and "deciding" became "fluctuating past the attractor". Scaling θ with δ* makes the curve monotone. §9.2. |
| **"High-n AM collapses to all-blank."** | Simply ran it: single-winner probability ≈1.0 at every n from 2 to 48, under two conventions. The cost of radix is consensus *time*, not collapse. §3. |
| **"A direct SSA measurement at Ω=120 would take hundreds of hours."** | Measured throughput: 84,000 steps/s and ~0.4·Ω steps per unit time make one flip at γ=0.35 cost **5.5 minutes**. The hundreds-of-hours regime is γ ≤ 0.30 — exactly where the CME's own solve is rejected, so neither instrument reaches it. §9. |
| **"At σ/δ*=0.45 two populations die at *exactly* the same depth."** | True at Ω=64 and 128, and it is integer rounding across two adjacent points. At Ω=256 the integer crossing reads 10. Interpolated, the depth creeps 6.53→9.14 over Ω=16→256 with increments halving — convergent to ≈9.4, which is the real (and still strong) claim. Caught by a reader asking to run Ω=256. §12.1. |
| **"τ_SSA/τ_CME ≈ 0.8 from an arm-vs-attractor offset."** | Seed-averaged measurement gives 0.97 with no offset visible. The prediction conflated an MFPT-from-the-arm with the mean time between crossings of a long stationary trajectory, which is dominated by the dwell near the attractor. §10.1. |
| **"The cascade's decay length ξ grows like the restoration wall e^{cΩ}."** | Fit directly: ln ξ vs Ω has R²=0.69 and ξ saturates (24→76 for 15× Ω). Superseded by **T1**, which explains why — the cascade is channel-limited, so ξ cannot inherit the finite-count scaling. |
| **"Freeze-out is a genuine transition, `Hc ≈ 0.055`, `a ≈ 0.38`; in the Ω→∞ limit expand slower than Hc and you always decide."** | The expanding SSA is *exactly* ordinary SSA stopped at internal time `1/H` (bit-for-bit, 0/300 mismatches), so `H*` is one over the consensus time — which from a symmetric start diverges like `(3/2)lnΩ`. `Hc = 0`, and Kurtz's theorem makes `Hc > 0` **impossible**, not merely unsupported. Extended ×16384 in Ω: `dτ*/dlnΩ = 1.5005 ± 0.0023` where a positive Hc needs it heading to 0; `D` at H = 0.055 falls 0.988 → 0.268 where FSS demands a constant; §5's own fit extrapolates 26% wrong; the curvature is 21× too small; and a **zero-parameter** collapse beats the two-parameter FSS by 28×. A biased-start control makes `τ*` flat (slope −0.0022 ± 0.0003) with a real `Hc = 0.2102`, so the drift §5 read as criticality was the shrinking shot-noise seed. §5.1. |

**A second pattern, from the depth-ceiling correction.** A threshold observable
(the first *integer* depth below I=0.5) reported two different states as identical,
because 0.474796 and 0.495227 both round the same way — and the second was 0.005
from rounding the other way. Quantised observables hide trends and then present the
gap as an invariance. Quote the continuous quantity; round only at display.

**A third pattern: one spare parameter can eat a logarithm.** `Hc + C·Ω^{−a}` and
`1/((3/2)lnΩ + B)` are indistinguishable over ×32 in Ω — SSR 4.8e-7 versus 5.0e-7
on §5's own six crossings, with the *log* law using one fewer parameter — and they
even predict nearly the same transition width, so the obvious sanity check was no
check at all. Three lessons, all of which cost this project a headline result:

- **A data collapse tests interleaving, not functional form.** §5 read a good
  collapse as proof of a critical point. A collapse only says "some monotone
  reparametrisation lines these up", which a wrong form with a free offset can do.
- **A three-parameter fit that beats a two-parameter one by 3% is evidence
  *against* the extra parameter**, not for the model. §5 never compared.
- **The way out was deriving the form, not fitting harder.** Every ingredient
  needed (λ = ⅓, `D_δ = 1/(9Ω)`) had been sitting in `design.md` §9 since §2.

**A fourth pattern: the consolidation reflex.** Twice now, separate anomalies were
merged into "one underlying question" on the strength of a number looking constant —
and both times the constancy was an artifact of the axis that happened to be swept.
§14's 2.27 was constant in n and drifts in δ; §12.1's ≈3 looked flat at Ω=128 and
drifts in σ once the death depth is converged. Both times the real defect was a
multiplicative error in an *exponent* being forced into a prefactor, which is flat
over a short range in the exponent's argument and only then bends. **Constancy along
the axis you happened to sweep is not constancy.** Before merging two anomalies into
one cause, measure each along an axis you did not choose for it — and note that the
merge itself is what suppresses that test, because a shared explanation makes the
individual checks feel redundant.

**The pattern worth remembering.** In every case the error was invisible to the
guard that was supposed to catch it, because the guard watched the wrong quantity:
`θ/δ*` was correctly constant while the *input* relative to θ was not; the control's
rails were never checked against the chemistry's. And in every case the claim was
first written down as a number that had not been run.
