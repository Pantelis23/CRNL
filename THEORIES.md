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

    −ln p ≈ κ·Ω·δ*² / (1 + 2κΩσ²),    κ(γ) = (3/2)(1−2γ)/(1+γ) = λ(γ)/(2·D₀(γ))

Pooled collapse over 216 cells: **R² = 0.960** (0.933 with §12's κ, which scaled
the restoring gain with γ but left the diffusion at its γ=0 value — corrected in
§15). The population-limited side, which
nothing had explored, is dramatic — 11 orders of magnitude in p across Ω at
σ_ch/δ* = 0.10. §11's protocol sat on the channel-limited side, which is *why* it
found no efficiency optimum in Ω.

**What is still open here (T1a), now smaller.** The slope was 0.74 pooled,
drifting 0.80 / 0.63 / 0.42 / 0.50 with γ. **Part of that was a wrong exponent,
not a missing prefactor**: with the corrected `κ` of §15 the same 216 cells give
0.783 pooled and 0.81 / 0.68 / 0.51 / 0.68 per γ. So "the missing piece is the
prefactor and the Gaussian-tail correction" — what this entry used to say — was
wrong about the first ingredient it named.

What remains is a slope of 0.783 rather than 1, still non-monotone in γ, and **two
of its three candidate causes are now eliminated by measurement** (§22.1, §22.2):
the barrier's γ-dependence is right (exponent 2 to 0.1% asymptotically), and doing
§12's δ-minimisation exactly instead of by saddle point improves R² everywhere
(0.960 -> 0.970) while leaving the slope at 0.776.

**§22.3's answer is WITHDRAWN by §22.4 — the model it used is 10^3 out.**
Compared against the exact single-stage flip from `stage_kernel` (no fitting), the
convolution framework is off by 5-3688x with the exact barrier and 1-10x with the
quadratic one, because it assumes the chemistry RUNS TO COMPLETION where §12's
stage has t_stage = 16. A fitted slope from such a model says nothing about §12's
physics.

**What replaced it is better.** The exact barrier is shallower than `κδ²` away from
the saddle, so `κδ²` overestimates the barrier and suppresses predicted flipping --
partially cancelling the framework's overestimate from assuming completion. **§12's
formula fits as well as it does partly by error cancellation**, which explains why
correcting one ingredient made agreement worse.

**§22.5 measures it.** The framework error runs **5 to 3688x** and the barrier error
**0.001 to 0.78** -- two spans of about three decades cancelling to within one. `B`
approaches 1 as the barrier shrinks (0.78 / 0.74 / 0.70 at γ = 0.45) exactly as it
must, since §15 verified `κ` IS the true curvature **at** the saddle: the quadratic
is right there and wrong further out. The mechanism's functional form is
deliberately **left unclaimed** -- finite stage time predicts F should grow with the
barrier and it does, but the pooled log-log fit is only R² = 0.62 and the
within-γ behaviour differs (roughly exponential at γ = 0.30, sublinear at γ = 0.45).
Direction confirmed, law not; this thread has already produced one withdrawn
interpretation from over-reading a fit. §22.1 is untouched.

**The superseded reading:** the quadratic barrier WAS the main cause, and the
γ-dependence is not. Replacing `c(δ) = κδ²` with the exact ridge profile moves the
slopes from 0.51 and 0.68 to 0.90 and 1.25 — from 40% low to bracketing 1 — where
nothing else tried moved them at all. But R² *drops* (0.935 -> 0.810, 0.974 ->
0.954) and the γ-spread is untouched (ratio 1.35 -> 1.39). So the residual splits:
a large explained part from the quadratic assumption, and a smaller **still
unexplained γ-dependence** that survives every correction applied. A fix that
improves the slope while degrading the fit has the right first-order term and the
wrong second-order one. Of the three candidates named for that, the profile's own
Omega-error is **eliminated** -- R^2 is flat to the fourth decimal across a 3x
change in it -- leaving the omitted prefactor/transverse relaxation and the channel
model. The same check found the gamma = 0.45 slope is **not converged** (1.20 ->
1.34 and still climbing), so the overshoot there is larger than first reported and
the residual gamma-dependence is worse, not better.

**The original framing, now spent:** the one untested ingredient is `c(δ) = κδ²` —
the barrier quadratic in the displacement. §2's own table has `c/ε²` drifting 1.586 -> 1.809 across
ε ∈ [0.04, 0.20], and §14.1 measured `c ∝ δ^~2.5` at large n, so there is direct
evidence it is not quadratic over the range the channel samples. **Kill test in
§22.2:** put the exact ridge profile into the convolution and refit at γ = 0.30 and
0.45, which are both inside the quasipotential's window and are the two worst
slopes. A derivation that predicts 0.783 — or shows why the slope should
drift with γ at all — would turn a collapse into a law.
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

**§15 changes the coefficient but not this conclusion.** The corrected
`κ = (3/2)(1−2γ)/(1+γ)` still vanishes linearly in `(γ_c − γ)` — the `1/(1+γ)`
factor is finite and smooth at γ_c = 1/2, contributing only a constant `2/3`. So
the barrier still falls like `(γ_c−γ)²` and the population still diverges like
`1/(γ_c−γ)²`; only the prefactor moves.

**MEASURED -> §22. The exponent is 2 and T4's guess of 1 is dead.**
`dW = 3.09·(γ_c−γ)^1.9745` at R² = 0.99997 over seven γ, with the local slope
climbing monotonically to **2.0015** at the smallest gap — the approach an
asymptotic normal form must make. §15 passes the same test: an exponent of 1 would
have restored T4 and falsified the corrected κ, since that rests on λ vanishing
linearly at γ_c. They stood together.

**And it eliminates a suspect.** §12's slope 0.783 and §12.1's k = 1.0695 were both
plausibly a saddle point getting the barrier's γ-dependence wrong. It does not get
it wrong, so those residuals live in the other ingredients — §12's second saddle
over the flip location, §12.1's σ-dependence — and are not a shared failure. Unlike
Q7, this narrowing came from measuring the suspected common cause, not from
noticing that numbers looked alike.

**The original kill test, now spent:** sweep γ finely in [0.40, 0.499] and fit the
barrier's exponent in `(γ_c − γ)`. A measured exponent of 1 would restore T4's original guess and
falsify the κ correction; anything other than 2 means the gain's linear vanishing
is wrong, which would also break §15 — the two now stand or fall together, and
that is a cheap way to test both. **The exact quasipotential of §15 is the right
instrument for it, and its usable window (large γ, small barrier) happens to be
exactly this region.**

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
symmetry-breaking eigenvalue. (**That κ was later corrected to
`λ(γ)/(2·D₀(γ)) = (3/2)(1−2γ)/(1+γ)` in §15**: §12 scaled the gain with γ and left
the diffusion at its γ=0 value. It does not rescue T7 — the correction is a factor
in γ, and T7 died on the *n*-dependence.) §13 computes that eigenvalue for any n, so the hope
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


**T8: a coefficient is a gain over a noise, and BOTH depend on the drive.**
**CONFIRMED -> §15.** `kappa = lambda/(2 D_0)`, and every place this project
carried a coefficient to a new regime it scaled the gain and forgot the noise.
The correct `kappa(gamma) = (3/2)(1-2gamma)/(1+gamma)` was verified by two
instruments that share no machinery -- the exact quasipotential's ridge curvature
(1.0004 / 0.9987 / 1.0012 at gamma = 0.35 / 0.40 / 0.45) and the original
first-passage wall measurement (mean ratio 0.990 over 10 cells) -- and it was
sitting in the repo the whole time as `breaking_diffusion(2, gamma)`, written for
a different section.

**Where else to look:** anywhere a rate, gain, or eigenvalue was carried to a new
parameter and the diffusion was not recomputed. §13's `A_c(n) -> 9 ln n` and
§14's `lambda/D_0 -> 1` both *do* recompute it (that is what §14 is about), so
they are clean. §12 was the one that did not.

---


**T9: an asymmetry is decisive only when it beats the shot noise, and the
threshold falls as the system grows.** **CONFIRMED -> §18.** AM's `X + Y -> 2B` is
an annihilation, so a tilted AM in an expanding volume has all three of Sakharov's
ingredients. From an exactly symmetric start the surviving species is set by the
tilt when `g/lambda > sigma` and by chance below it, i.e. above

    beta * sqrt(Omega) = (sqrt3/2)(1-2 gamma)/(1-gamma)      [0.820 at gamma = 0.05]

and `P(X) = Phi(u)` with `u = (g/lambda)/sigma`, parameter-free, verified to
under 1% with a collapse holding to 0.0028 across a 4x population change. The
`sigma` uses §15's corrected `D_0(gamma) = (1+gamma)/9`, so this result rests on
that correction and would be visibly wrong with §12's.

**What is worth noticing:** the decisive asymmetry SHRINKS as `Omega^-1/2`. A
bigger system needs a *smaller* bias to have its outcome determined rather than
accidental -- the opposite of the usual "more molecules, more averaging" reading,
because what grows with Omega is the amplifier, not the noise.

**T9a, open: the analogy is not a mapping and the gap is the interesting part.**
Recruitment `B + X -> 2X` has no cosmological counterpart; there an asymmetry
survives linearly, here it is amplified. Nothing here says what a *passive*
annihilation network (no autocatalysis) would give, and that is the control this
entry needs before any of it is described as baryogenesis. **How to kill T9:**
run the same measurement on annihilation + pair production alone. If `Phi(u)`
survives with no restoring landscape, the landscape was never doing the work and
§18's framing is decoration.

---


**T10: the model had no temperature, and that is why expansion was only a clock.**
**CONFIRMED -> §19.** Every AM reaction is 2 -> 2, so dilution scales every
propensity identically and the landscape is invariant under expansion. §5.1's
exact reduction to "ordinary SSA stopped at internal time 1/H" is therefore a
theorem about **uniform-order kinetics with state-independent rate constants**,
not about restoration under expansion in general. Let the medium cool
(`gamma(s) = gamma0**((1-s)^-w)`, forward rates untouched) and the reduction fails
at once, because the drive profile is universal in `s = H*tau` while the number of
reactions inside it is not.

The payoff is an observable the fixed-drive model could not have: a relic
**abundance**. Conditioned on deciding, it rises 290x over a 4x range in H and
sits 10^5-10^8 above the equilibrium value at the drive it froze at, while the
fixed-drive arm is flat over the same range and simply equals its equilibrium.
**Abundance set by expansion versus abundance set by chemistry.**

**This does NOT overturn `Hc = 0`.** From an exactly symmetric start the
deterministic ODE stays symmetric under *any* `gamma(t)` while `beta = 0`, so
`D -> 0` for every H and §5.1's impossibility argument survives untouched. What
changed is the scope of its reduction, not its conclusion.

**T10a, open: does the cooling rate have an optimum?** The competition is real --
cooling deepens the landscape while dilution starves it -- but only H has been
swept, at one `w`. Sweeping `w` at fixed H asks whether there is a best cooling
schedule for getting a decision made, which is a different question from the
relic abundance and is unmeasured. **How to kill T10:** show the relic collapses
onto a function of `s_c` (the pitchfork crossing) alone, which would make the
cooling schedule a reparametrised deadline after all.

~~**T10b: the drive never depletes.**~~ **BUILT AND MEASURED -> §20.** The drive
is now a fuel species: `F` is consumed and `W` produced by every forward step, so
`gamma_eff = gamma_inf * w/f` RISES as the tank empties -- the mirror of §19's
cooling. Not bookkeeping: `n_F` is a genuinely independent coordinate, because a
full cycle `f1->f2->f3` returns (X,Y,B) exactly to its start while burning three
fuel, so **the fixed-gamma model is a projection that discards a coordinate which
must exist.** Anchored the way §19 was: with the tank frozen the drift equals
`am_reversible(gamma_inf w/f)` time-rescaled by `f`, to 1e-16.

**The second ceiling is real and has a different shape from the first.** The
fuel-limited memory lifetime is FLAT in Omega -- spread 1.08x and 1.16x over a 6x
population range at two fuel concentrations -- while the noise-limited lifetime on
the same clock is `exp(0.1215*Omega)` at R^2 = 0.984. They cross at Omega ~ 3-8,
and by Omega = 180 the noise ceiling is 1e9 times further away. **Above a
population of about ten, restoration is fuel-limited and more molecules buy
nothing** -- the exact mirror of §1's wall. §12.1's ceiling is untouched but no
longer alone: that one is a DEPTH set by channel noise and independent of the
drive, this one is a LIFETIME set by the budget and independent of the population.

**More fuel gives a SHORTER lifetime** (7.5 at concentration 10, 4.6 at 30). In
waste-fraction units the burn rate is exactly fuel-concentration-independent, so a
bigger tank buys no extra fractional runway; it only makes the chemistry fast
relative to the fuel clock, so the state tracks the collapsing landscape more
adiabatically and gives the bit up earlier in the tank's life.

**A prediction of mine failed in SIGN, and the miss is instructive.** I predicted
the bit is lost BEFORE the drive formally dies, since the barrier degrades
continuously on the way. It outlives it by 7-37%: once `gamma_eff` passes
`gamma_c` the landscape is gone but the state must still physically relax off the
old attractor, and that relaxation burns more fuel. I reasoned about the barrier
degrading and not about what happens after it vanishes.

**T10b-i, open: does the overshoot saturate?** `0.069*ln(Omega) + 1.018` fits at
R^2 = 0.969 and a saturating `c - a*Omega^-b` fits comparably at `c = 1.43 +- 0.03`.
Six times in Omega cannot separate them -- the same wall §17.2 hit, and the same
refusal to quote the flattering one. Needs a longer lever.

~~**T10b-ii: which ceiling binds for a CASCADE?**~~ **MEASURED -> §20.3, and the
premise of the question was wrong.** Both arms through one harness: at
sigma_ch/delta* = 0.40 they agree within 1.2x (noise binds both) and at 0.03 the
fueled cascade dies 36x earlier (exhaustion binds). So a cascade does have two ways
to die and the fuel concentration decides which.

**But they do not combine as a min().** The fuel-limited depth was predicted to be
independent of the channel noise and instead falls 2.2x-3.8x across the sigma
range: spent fuel raises `gamma_eff`, which shrinks `delta*` AND `kappa(gamma)`, so
a half-empty tank makes the same channel noise bite harder. **Exhaustion and noise
compound.** The right picture is not "restoration runs until the fuel is gone and
then stops" but "restoration degrades continuously as the fuel goes, and the
channel finishes it early."

~~**T10b-iii: how does D_fuel scale with the budget?**~~ **MEASURED -> §23, and the
kill test fired: `D_fuel` is NOT a budget property.** Five budgets x three channels x
three populations. The exponent is `0.6498 +- 0.019` at the quiet channel, so the
sub-linear two-point reading was right and §20's crossover argument (which assumed
linear) is wrong -- and the exponent **drifts with the channel**, 0.650 -> 0.495 over
sigma_ch/delta* = 0.03 -> 0.15 at 5.5 sigma, reproduced at 4.1-4.2 sigma at Omega 25
and 60. It drifts with Omega too (0.547 -> 0.650 -> 0.672).

**So the wording that has to go is "the" fuel-limited depth.** There is no
fuel-limited depth attached to a tank; there is one **for a given channel and
population**, with a budget exponent running 0.65 -> 0.26 across the range tested.
The compounding claim above is unaffected -- §23 is the same fact from the budget
side -- but any statement of the form "a tank of size Phi buys depth D_fuel(Phi)"
is now known to be incomplete.

Mechanism, and it is not the one predicted: fuel burned per stage is FLAT in the
budget (25-34 molecules over 16x, slope -0.01 to -0.03), exactly as the a-priori
argument requires (third-order burn rate goes like `f`, `stage_time` like `1/f`).
What collapses is the **usable burn fraction**, 0.198 -> 0.069. A bigger tank means
more stages spent while the drive is still strong, and every one of those stages
carries its own loss probability, so the bit goes at a shallower `gamma_eff`:
0.75 at Phi/Omega = 25 (past `gamma_c` -- the bit outlived the landscape) but 0.43
at Phi/Omega = 400 (landscape still alive). A louder channel raises the per-stage
loss probability, hence the sigma-drift.

~~**T10b-iii-a: what IS the exponent a function of?**~~ **TESTED -> §23.4, and the
mechanism above is mostly wrong.** The proposed kill test -- collapse the exponent
onto a curve in the healthy-tank per-stage loss probability `p` -- was badly posed
and §23.4 says why: at sigma_ch/delta* = 0.03 the healthy-tank `p` is at most
0.00125, which over the 44 stages the largest tank lasts accumulates to ~5%. It
cannot carry the effect, so there was nothing to collapse.

Tested properly instead, as an ABSOLUTE hazard integral with no free parameter
(measure `q(theta)` and `c(theta)` in single stages, integrate the survival product,
read off the median depth). Verdict, three parts:

- **The absolute depths come out right at the large budgets** -- 1.09x and 1.20x at
  Phi/Omega = 200-400, parameter-free. Real support for the picture.
- **Accumulated hazard explains the sigma-DRIFT and only that.** Forcing the hazard
  to zero changes the predicted depth by 0.0% at four of five budgets, but flattens
  the drift from -0.121 (measured -0.107) to -0.030.
- **It does NOT explain the sublinearity.** With hazard off the exponent is 0.9156,
  near the linear value; measured is 0.6474. The integral matches large tanks and
  under-predicts small ones by up to 1.75x, so the sublinearity is mostly **small
  tanks over-performing**, not big tanks under-performing.

Why they over-perform, measured at 200 trials/cell: the smallest tank holds the bit
for a median **3 of its 7 stages at `gamma_eff` > `gamma_c`**, i.e. in a MONOSTABLE
landscape with no rail, at a real separation (median delta = 0.19, ~8 molecules at
Omega=40); 97% of trials reach that regime against 9.5% for the largest tank.
Restoration is unavailable and the bit persists on kinetics -- the stage is two
relaxation times at `gamma_0` while `lambda(gamma_eff) -> 0` at `gamma_c`, so a tank
that drains in a few stages gives the state no time to follow the collapsing
landscape. NOT established: any general "the state lags the rail" statement. That
prediction's statistic was confounded and the matched-`gamma_eff` version scatters
0.70-1.65 with no ordering in Phi (§23.4).

~~**T10b-iii-b: can the kinetic excess be computed rather than named?**~~ **TESTED ->
§23.5, and the attribution is WITHDRAWN.** Two arms at the quiet channel:

- **Imposing quasi-staticity by hand** (re-seed X,Y,B onto the current rail every
  stage -- a DIAGNOSTIC with a free restoring element, rule 10, not a physical
  result) gives exponent `0.7077 +- 0.030`: **5.4 sigma from the hazard integral's
  0.8925 and only 1.6 sigma from the plain simulation's 0.6474.** Re-seeding on a
  median 52 of 52 stages closes essentially none of the gap. So the quasi-static
  state assumption is NOT what separates the integral from the measurement, and the
  integral's error is somewhere other than the state.
- **The adaptive-stage sweep proposed below as the kill test is a WEAK INSTRUMENT and
  its null carries nothing.** The exponent is flat across a 30x cap sweep
  (0.621/0.644/0.641/0.645) -- but the adaptive prescription stays within 2.4x of the
  fixed stage until `theta > 0.29` and only diverges in the last stage or two, and
  longer stages both help the state follow AND drain the tank faster, so a flat
  exponent could be two effects cancelling. Recorded as inconclusive. The `cap = 1`
  cell is a genuine control and passes: 0.6206 +- 0.012 against §23.1's 0.6474 +-
  0.022, 1.1 sigma.

**So: §23.4's phenomenon stands, its explanation does not.** Small tanks really do
hold real bits (delta ~ 0.19) past `gamma_c` for a median 3 of 7 stages, 97% of
trials against 9.5% for the largest tank. "Because `lambda(gamma_eff) -> 0` while the
stage stays at `gamma_0`'s relaxation time" is withdrawn.

~~**T10b-iii-c: why does the hazard integral over-predict the exponent at all?**~~
**TESTED -> §23.6. The hard stop was 39.6% of the error, and no more.** Seeding the
past-`gamma_c` cells at an imposed separation instead of the nonexistent rail moves
the exponent `0.8867 -> 0.7919` against a measured `0.6474`, leaving **+0.1445, 4.6
sigma** unexplained. Reported as the one-parameter model it is, with `delta_past`
swept rather than fitted -- and the sweep barely matters: 0.0103 of exponent across a
2.5x change, so the one imported number is nearly inert.

**An instrument fault fell out of this and is worth carrying forward.** The integrator
counted WHOLE stages, and at the smallest budget the predicted depth is ~5 stages, so
four visibly different models returned byte-identical exponents to four decimals --
which is the only reason it was caught. The measured past-`gamma_c` hazard really does
fall from `q = 0.429` to `0.302` across the sweep; the integer depth could not move.
Interpolating the survival crossing in `ln S` fixes it (rule 13: an approximation's own
numerical parameter is a second axis). **This changes no published number** -- the
continuous integrator puts §23.4's hard-stop exponent at 0.8867 +- 0.016 against the
printed 0.8925 +- 0.017, 0.25 sigma apart.

~~**T10b-iii-d: is the residue survivor selection?**~~ **MEASURED -> §23.7. NO, and I
had the sign backwards.** 3000 trials at each of three budgets, conditional hazard by
past-`gamma_c` survival index `k`, controlled by dividing out the unconditional
`q(theta)` the integral actually used. Predicted: ratio below 1 and FALLING with `k`,
carried separation RISING. Measured: the ratio climbs 0.80 -> ~1 and hovers, and the
separation **decays monotonically** at all three budgets (0.317 -> 0.117 at
Phi/Omega = 25). Survivors are not the ones that kept a big bit; the bit is being
ground down, which is what a monostable landscape does to it.

**The same data gives a better result than the suspect it killed.** The decay is a
function of `theta`, NOT of survival history: pooling 22 cells over three budgets and
eight survival indices, `delta ~ theta` alone gives weighted R^2 = 0.9177, and adding
`ln Phi` lifts it only to 0.9675 at coefficient -0.028 (~23% of the delta range, on
medians quantised at 1/Omega). At matched `theta ~ 0.35` every budget carries
`delta ~ 0.30` however many stages it took to get there. **So §23.6's closing note --
that any repair makes the model non-Markovian in `theta` -- is wrong in the useful
direction.** The constant `delta_past` was just the wrong constant.

**And the residue is two errors in two regions**, which is the reframing this thread
needed: at Phi/Omega = 25 the integral is -22% (past-`gamma_c`), at Phi/Omega = 400 it
is +10% (pre-`gamma_c`, where only 9.5% of trials ever reach `gamma_c`). No single
correction touches both.

**The Markovian REPAIR that §23.7 proposed off the back of that collapse is refuted
(§23.8), while the collapse itself stands.** Imposing the measured `delta_past(theta)`
moved the exponent 0.7919 -> 0.7880, closing 2.7%. It cancels: the curve raises the
imposed separation at low `theta` and lowers it at high `theta`, and the integral
traverses both. Already visible in §23.6 and under-weighted -- a 2.5x sweep of the
constant was worth 0.0103 of exponent. **The integral is not sensitive to the
separation it imposes past `gamma_c` at all**, so §23.6's 39.6% came from permitting
past-`gamma_c` stages at all, a structural correction, not from their value. Recorded
as a rule 17 violation committed in the same paragraph that reported rule 17 working.

**T10b-iii-e, open: the residue is STRUCTURAL, and one structural approximation is
left untested.** §23.9 reopened the thread on a different axis from the one closed
below. The integral makes exactly two structural approximations -- deterministic
`theta` and independent stages -- and the two surviving errors have their signs and
`Phi`-scalings. Testing the first WITHOUT CHANGING ANY HAZARD VALUE (same `q(theta)`,
same measured burn, `theta` propagated as an ensemble with per-stage burn drawn from
its measured spread) moves the exponent 0.7919 -> 0.7641, closing 19.2%.

Scoreboard against the 0.2393 gap from the hard-stop integral: **structural repairs
39.6% + 11.6% = 51%; parametric repairs 2.7% + ~0% = 3%.** That is the pattern, and
it points the remaining ~49% at the one structural approximation not yet tested:
**stages are treated as independent, but `delta` carries memory across them** -- a bit
knocked down stays down, so correlated trials die sooner than an independent product
predicts, which is an OVER-prediction and matches the +10% at the large tanks that no
past-`gamma_c` repair could touch.

**PARTLY ANSWERED -> §23.10: inter-stage memory is the largest single term, 56%.**
Arm A (the rail-reseeded diagnostic, which REMOVES memory) re-run at 400 trials and
paired against the plain arm from its own run gives `0.7217 +- 0.024` vs
`0.6325 +- 0.016`, a difference of **+0.0893 +- 0.0286 = 3.1 sigma**, localised at the
two largest budgets (+7%, +14%) and absent at the two smallest -- exactly where the
integral over-predicts. Decomposition, all paired: integral 0.7919, memoryless
simulation 0.7217, plain 0.6325, so **memory is 56% of the integral's error and the
integral's residual error against a MEMORYLESS simulation is 44%**, of which
theta-dispersion is 17%.

**§23.5's "1.6 sigma from doing nothing" is WITHDRAWN.** It compared Arm A against
§23.1's separate run instead of the plain arm in its own, and summarised at the
exponent level a structure that lives per-budget. Rule 9 one level up: the axis was
the budget, the statistic was a fit across it. Raised by an independent analysis of
the stored data, verified here before adopting. Still unexplained and flagged:
Phi/Omega = 50 moves the WRONG way under Arm A at both 80 and 400 trials.

~~**Also proposed by that analysis: `c(theta)` is a conditional-mean error.**~~
**Premise verified, consequence refuted (§23.10).** `hazard_at` really does pool
losers into `c`, and the stoichiometry really does make losers burn more -- `f1` fires
at a rate proportional to `n_X n_Y n_F`, maximal at `delta = 0`, and measured
`c_lost/c_surv` reaches **4.28**. But the contamination is `q*(c_lost - c_surv)`, and
`q` is tiny exactly where the ratio is large; peak contamination is 4.5% and
re-integrating with `c_surv` moves the exponent 0.7914 -> 0.7918. A real rule-12
defect inside the instrument with no measurable consequence, and another instance of
§22.4's cancellation motif.

**Still open:** the 44% the integral gets wrong even against a memoryless simulation,
of which theta-dispersion is 17%. **How to kill the rest:** measure the per-stage loss
probability CONDITIONED on the separation carried into the stage, `q(theta, delta)`,
and propagate `delta` as a second ensemble coordinate instead of collapsing it into
`q(theta)`. If the exponent falls to Arm A's 0.7217 the collapse onto `theta` alone
was the whole of the remaining instrument error. If it does not, the residue has
no surviving suspect, structural or parametric, and the integral should be abandoned
rather than repaired. **Note this is a genuinely different object** -- a two-variable
stage map, which is most of the way to just running the simulation -- so the result
worth having is the DECOMPOSITION (how much each approximation costs), not a working
model.

**Caveats to carry into it, both from §23.9.** (i) Part of the large-tank change under
the ensemble is trajectories crossing the `theta_max` cutoff, a model boundary rather
than physics. (ii) §23.9's own P4 fired: 0.0278 of exponent against a 0.032 combined
fit error is 0.87 sigma, so `theta`-dispersion is a real contributor of the right sign
and scaling but is NOT significant by the criterion set in advance. The paired
per-budget ratios are monotone across all five budgets, which is the argument for it
being real; both readings are on the record.

~~**T10b-iii-e as first framed (pre-`gamma_c` over-prediction), and any further
PARAMETRIC repair of the integral:**~~ **NOT OPENED, deliberately.** Three repairs have been tested -- imposed
quasi-static state (1.6 sigma from nothing), removing the hard stop (39.6%),
correcting the imposed separation (2.7%) -- and +0.14 survives as two errors in two
regions. Following §22.5's precedent the residue is left named and unexplained rather
than attacked a fourth time. The integral has done what an instrument should: it
reproduced the sigma-drift (-0.121 vs -0.107 measured), which is the only part of
§23.3's proposed mechanism that survived a test, and it eliminated two of three
candidate errors. §23's measured results (§23.1, §23.2, §23.4's phenomenon, §23.7's
collapse) require none of it. **If this is ever reopened, the kill test is the
large-tank +10%, which lives entirely below `gamma_c` where §23.5 already showed the
state assumption is innocent** -- so the next suspect there is the survival product's
independence, not anything about the dying landscape.

~~**Superseded framing of T10b-iii-d, left visible:**~~ The remaining error is a
SHAPE error across `Phi`, not a scale error -- predicted/measured runs 0.78, 0.72,
0.88, 1.00, 1.10 -- so no normalisation fixes both ends. *Suspect, stated as one:*
the survival product is mean-field and applies the unconditional `q(theta)` every
stage, but among trials that survive several past-`gamma_c` stages the separation is
selected upward, so survivors face a lower hazard than `q` says. The bias grows with
the number of past-`gamma_c` stages -- 3 of 7 at the smallest tank, 0 of 46 at the
largest -- i.e. exactly where the integral under-predicts. **How to kill:** measure
`q` conditioned on having already survived `k` past-`gamma_c` stages; if it falls with
`k`, selection is the residue. If it is flat in `k`, the independence assumption
itself is next. Either fix makes the model non-Markovian in `theta`, which is a
different object from the integral §23.4 set out to build -- so this may be the point
at which the integral has served its purpose rather than a defect to repair.

~~**Superseded kill test for T10b-iii-b, left visible:**~~ §23.4
attributes the 0.25 exponent gap (quasi-static 0.893 -> measured 0.647) to holding
past `gamma_c`, but nothing there predicts the number. The kinetic picture says the
excess should be controlled by the ratio of the stage time to the relaxation time at
the CURRENT `gamma_eff`, i.e. by `lambda(gamma_eff) * t_stage`, which vanishes at
`gamma_c` -- so it is a statement about a dimensionless group, not about Phi.
**How to kill:** re-run §23's five-budget sweep with the stage time set adaptively
from `lambda(gamma_eff)` instead of `lambda(gamma_0)`, so every stage is two CURRENT
relaxation times and the state cannot fall behind. If the exponent rises to the
quasi-static 0.89-0.92 the kinetic account is right; if it stays near 0.65 the excess
is not a lag at all and §23.4's attribution must be withdrawn. Note the adaptive
stage changes the fuel burned per stage, so `c(theta)` must be re-measured rather
than reused -- and the control has to move onto the same clock or this repeats
§10.3's mismatched-rail error.

---


**T11: coarse-graining restoration is a cliff, not a slope.** **MEASURED -> §21.**
Against an exact CME reference, every level that keeps ANY noise recovers the
restoration error exponent to 2-12% -- the chemical Langevin equation (real-valued
counts, Gaussian noise), tau-leaping (windowed Poisson firings) and the exact SSA
are all in one class. The ODE, which keeps none, reports p = 0 in all sixteen cells
where the truth spans 1.5e-3 to 1.6e-1, and has no refinement parameter that
improves it.

**So: the discreteness, the exact jump timing and the correct jump distribution are
all discardable for this observable; having noise at all is not.** The corollary is
about cost -- a cheap SDE gets the exponent, and the expensive exactness buys the
prefactor and the individual probabilities. A simulation that needs to know how
fast reliability grows with population can be cheap; one that needs the actual
failure rate cannot.

**Why this is not a numerics result.** Kurtz's theorem licenses the ODE limit on
finite time intervals and §5.1 uses it. It is true and it does not cover this
observable, because restoration lives in tails where the convergence is not
uniform. **A limit theorem cannot tell you what your simulation may throw away.**

**T11-REFINED -> §24: "having noise at all" is not the right axis, and §21 could not
have seen it.** Every level in §21's ladder keeps every species, so it retains the
noise AND the coordinates; the two were never separated. §24 separates them on one
measured kernel `K(theta, delta)`: keeping the coordinate `delta` while propagating
`theta` DETERMINISTICALLY buys `+0.1000` of exponent, while keeping `theta`-noise and
discarding `delta` buys `+0.0146` -- about 7x -- and once `delta` is kept, adding
`theta`-noise is worth `-0.0026`, the wrong sign. Paired ratios monotone across all
five budgets.

This RECONCILES with §21 rather than overturning it, under a sharper statement.
§21's ODE is deterministic in every coordinate INCLUDING `delta`; §24's working cells
are deterministic only in `theta` while `delta` still fluctuates. In every cell that
works the bit-carrying coordinate keeps its noise; in every cell that fails it has
been removed or collapsed. **Noise matters in the coordinate that carries the signal
and costs nothing in a bookkeeping coordinate** -- which also explains why §23.9's
theta-dispersion bought only 17%, a number that is awkward for "noise is what matters"
and natural here.

~~*Suspect, not law (rule 17).*~~ **KILL TEST RUN -> §24.1, and it GENERALISES --
more sharply in §21's own system than in the fuel network.** Same ladder, same
observable, same exact CME reference; `am_reversible` conserves n_X+n_Y+n_B so the
CLE noise splits cleanly into `delta = n_X - n_Y` (the signal) and `s = n_X + n_Y`
(the blank pool, bookkeeping). Projecting the noise while keeping the drift full:

| arm | delta-noise | s-noise | result |
|---|---|---|---|
| full CLE | 100% | 100% | correct to 0.2-5% |
| delta-only | 100% | 0% | correct to 2-18% |
| uniform 11% | 11% | 11% | wrong by 17-770x, or categorically 0 |
| s-only | 0% | 88% | **categorically 0, like the ODE, in all eight cells** |

**A model can keep seven-eighths of the noise and be as categorically wrong as one
that keeps none.** No graded middle -- the messy outcome named in advance did not
occur. And the `uniform 11%` arm, added expecting it to PASS and thereby show
amplitude was irrelevant, failed instead: keeping 11% of the SIGNAL's own noise is as
fatal as keeping none, because barrier crossing is exponential in the noise amplitude
along the crossing direction.

**So T11's headline is superseded: not "having noise at all" but "having the signal
coordinate's noise, at its own amplitude."** §21's measurements stand as printed.

Honest residue: dropping the blank-pool noise is not free -- delta-only runs +1.5% to
+4.5% high at eps = 0.20 and +10% to +18% at eps = 0.35, systematically over, and
worsening with the barrier. Bookkeeping noise is MOSTLY discardable, not exactly.
Also note §24's own full cell reached only 0.7026 against a measured 0.6325, so
`(theta, delta)` is not the whole state of the fuel system and 0.0701 of exponent lies
outside both coordinates -- that gap is untouched by this test.

**T13 -> §24.1b: §24's arms are NAIVE discarding, not correct reduction, and the
difference is a generalized force that resists closed form.** Prompted by the
coarse-graining literature (Zwanzig projection, arXiv:2512.03706): a rigorous
reduction returns a generalized force and a state-dependent diffusion, whereas §24
zeroes a noise component and leaves the drift alone. So §24-§25 measure what NAIVE
discarding costs -- what a practitioner does -- and not what a correct reduction
achieves. The scoping note is now in §24.

**Measured, not derived.** Running full and `delta-only` from an identical state and
differencing the mean delta-increment gives the missing force directly: **+0.0483 at
gamma = 0.05, -0.0015 at 0.30, -0.0041 at 0.45**, i.e. positive and large exactly
where §24.1a's residual is large (+13.2%) and ~0 where it vanishes (+0.4%). It also
BUILDS with the averaging window (+0.018 at 0.5, +0.048 at 2.0), which is a memory
term's signature. And it cannot be diffusion: `delta-only` preserves `a - b` exactly,
so delta's own diffusion is untouched. **The Zwanzig reading is confirmed in kind --
the cost of naive discarding is a DRIFT term.**

**Both closed forms fail, and rule 16 caught them.** The curvature term
`1/2 d2(b_delta)/ds2 Var(s)` is **identically zero**: `b_delta` is exactly linear in
`s` (second derivative 4.4e-16; a sweep across s+-3 is perfectly straight), because the
AM drift is bilinear in (delta, s). The `+force` arm came out bit-identical to
`delta-only`, which is how it was found. The cross-correlation term
`db_delta/ds * D_ds / lam_s` is **wrong in SIGN at gamma = 0.05** (-0.085 predicted vs
+0.048 measured) and **11x-140x too large** at 0.30 and 0.45. Fitting a coefficient
would have made it look successful at two of three gammas.

**T13-a, ATTEMPTED -> §24.1c, NOT closed.** The force §24.1b could not write down was
LEARNED instead (MLP + linear control, trained on the local force under common random
numbers, scored on exact P(error) which was never in the loss). Results:

- **A learned closure beats naive deletion**: at Omega = 80 the MLP takes `delta-only`'s
  +9.6% down to +4.0%; over nine cells, 14.1% -> 10.2%. The §24.1b direction holds.
- **The force is NONLINEAR** (P4 refuted). MLP train R^2 = 0.98-0.99 while linear
  swings 0.064-0.91, and averaged over cells linear (+22.2%) is WORSE than doing
  nothing (+14.1%).
- **Local accuracy is not tail accuracy** (P3 confirmed, and the point of the
  section): R^2 = 0.99 on the force buys about HALF the residual, not 99% of it. Rule
  16 restated for universal approximators -- which is why the scoring target was fixed
  in advance and kept out of the loss.
- **A' does not pass.** The fitted closure is window-dependent, improving monotonically
  toward tau = 4 in five of six model x Omega pairs (Omega=40 MLP: +18.1%, +13.2%,
  +1.6%). Either no Markovian closure exists, or tau = 4 is the first converged window
  -- **this run cannot separate them, and tau = 8 was not run for the closure.** That
  is the missing cell, and the +35.4% -> +2.1% jump between tau 2 and 4 is large for
  something supposedly converging.

Sampling floor: ~5% relative SE at p ~ 0.0094 and 40,000 trials, so single-cell
differences of a few percent are not individually significant and the claims rest on
the pattern across nine cells.

**T13-a-i, open: does the closure converge in the window?** Run the tau sweep out to
8 and 16 at fixed Omega. Saturation means a Markovian closure exists and the earlier
windows were undersampled; continued drift means it does not. This is the one cell
that decides between the two readings above.

~~**T13-a as first posed: is the correct reduction non-Markovian here?**~~ The measured force grows
with the correlation window, which no Markovian drift correction reproduces. **How to
kill:** measure the force as a function of window length out to several `1/lam_s` and
check whether it SATURATES. If it saturates, a Markovian generalized force exists and
neither candidate above is it; if it keeps growing, a correct reduction needs the full
memory kernel and this project's machinery cannot construct one. Note in advance that
none of this touches the CATEGORICAL failure -- no generalized force on a
deterministic coordinate produces barrier crossings.

**T11-REFINED-c -> §24.1a: the strongest competing explanation, tested and SPLIT.**
An independent review named an account I had not written down: `s` is a fast STABLE
variable (`lambda_sym = -(1+2g) = -1.60` against `lambda_antisym = +0.133` at
g = 0.30, ratio 12) started on its own nullcline, so §24.1 might be the textbook
large-deviation result that escape rates are set by diffusion along the UNSTABLE
direction -- a property of the saddle geometry, i.e. of the system, which is §24's own
named failure mode. The discriminating axis is `gamma`, since `3(1+2g)/(1-2g)` spans
3.7 to 57.

Measured: **`s-only` is exactly 0 at every gamma from 0.05 to 0.45**, including where
the separation is only 3.7 and the pool holds ~1.9 molecules and so cannot be
Gaussian-slaved at all. **The cliff is not slaving.** But `delta-only`'s residual cost
falls monotonically with the separation -- 13.2%, 8.1%, 3.1%, 0.4% -- which IS
subleading transverse diffusion. **So the categorical failure is about which subspace
carries the observable; the few-percent residual is about timescale separation.**
§24.1 left that residual unexplained; it now has a mechanism, and it is the
reviewer's, not mine.

**T12: PARTLY TESTED -> §25.1. The best candidate came back negative.** `Var(T)` at
first passage -- a PURE-noise observable, deterministic limit exactly zero, with two
credible and non-guessable mechanisms -- is recovered by `delta-only` alone (ratio
1.183, 8/8 inside a band fixed before looking) while `s-only` captures only 16% and
fails 8/8. The P4 control rules out bimodality: `Var(T | correct)` for `s-only` is
0.198 of full, so the shortfall is missing jitter, not missing slow error paths.
Exact reference added as `cme.first_passage_moments` and pinned against the SSA,
because a wrong factor in `Qtt m2 = -2T` still yields a plausible positive variance.

**So three observables have now been tested and all three live in `span(delta)`.**
Per the review's own framing that does NOT establish one stiff direction -- it
establishes that P(error), MFPT and Var(T) are the same question asked three ways.

**What survives as a positive result:** the pool's noise is not irrelevant to Var(T),
just insufficient. `delta-only` overshoots Var(T) by 18% where it overshoots P(error)
by 3.1% at the same gamma (§24.1a), so **the pool matters ~6x more to the timing
jitter than to the error probability**. Observable-dependence in the MAGNITUDE, where
the categorical requirement is identical. Weaker than T12 set out to show, and it is
what the data supports.

~~**T12 remaining: the two-target race.**~~ **RUN -> §25.2. THE REQUIREMENT REVERSES,
and T12's strong form is established.** Racing the decision target `|n_X-n_Y| >= thr`
against a pool target `n_B <= m`, with `m` tuned per cell to make it ~50/50 (at
Omega = 40, m = 5 against a starting pool of 9, so the pool must fluctuate DOWN BY 4):

| observable | needs | other subspace alone | cost of dropping the other |
|---|---|---|---|
| P(error) | `span(delta)` | `span(s)` -> exactly 0, 8/8 | +3.1% at g=0.30 |
| race | `span(s)` | `span(delta)` -> exactly 0, 4/4 | +9% to +13% |

**So the requirement varies in BOTH directions**, and §25's "downward only" is
superseded. The caveat, declared before the run and kept: each arm's categorical zero
is semi-definitional within its own observable -- a pool target is a pool-fluctuation
event, so an arm without pool noise cannot reach it. **The weight therefore rests on
the QUANTITATIVE halves, which no construction forces:** dropping the pool's noise
costs P(error) 3.1%, and dropping the signal's costs the race 9-13% (36 SE at 40,000
trials). `s-only` landed "partial" not "recovers" against a band fixed in advance,
and the band was not moved afterwards.

**This also vindicates §25.1's worry about the instrument.** P(error), MFPT and
Var(T) all live in `span(delta)` because all three are first-passage functionals OF
DELTA. The race is the first observable here whose absorbing set is defined on a
different coordinate -- which is what it took, and the limitation was the instrument,
not the chemistry.

**T12-RULE -> §25.3: the arc's one usable rule, stress-tested and standing.**

> **The required subspace is the one the observable is a functional of.**

P(error), MFPT and Var(T) are first-passage functionals of `delta` and need
`span(delta)`; the race's pool target is a functional of `n_B` and needs `span(s)`.
This is the practical payoff -- applicable by inspection, without running a
projection.

Attacked where it should break: the pool is `b* = g/(1+g)`, ~9 molecules at g = 0.30
but ~2 at g = 0.05, where its relative fluctuation is order one and it gates every
recruitment propensity. Predicted in advance that Var(T), a delta-functional, would
be forced to require `span(s)` there. **It is not.** `s-only`'s share of Var(T) is
0.140 at g = 0.05 against 0.160 at g = 0.30 -- unchanged -- and `delta-only` recovers
Var(T) at EVERY gamma from 0.05 to 0.45. The rule is not an artifact of the one gamma
it was found at.

**T12-RULE-a, open: does the rule survive an AMBIGUOUS observable?** Every observable
tested is transparently a functional of one coordinate. The rule has never been tried
where that attribution is unclear -- e.g. an absorbing set defined on a mixed
coordinate like `n_X` alone (which is `(s + delta)/2`, a genuine mixture), or on a
ratio. **How to kill:** race a target defined on `n_X` against one on `n_Y`; if the
required subspace is the mixture's own span the rule generalises to "the subspace the
functional depends on", and if instead it collapses to `span(delta)` or `span(s)` the
rule is really about the network's eigendirections and not about the observable at
all -- which would be a much weaker statement and would partly restore the slaving
reading §24.1a defeated.

**T12-SIGN, open and unpredicted: `s-only`'s error on Var(T) changes SIGN with gamma.**
It under-predicts by ~7x at g <= 0.30 and over-predicts by 2.4x at g = 0.45 (up to
4.1x), crossing between. The full-CLE control recovers in all cells at every gamma, so
this is not the instrument. *Suspect:* near `gamma_c` the landscape flattens, so with
delta deterministic the drift to threshold is weak and slow and pool fluctuations
modulating a weak drift produce enormous timing variance -- i.e. removing the signal's
noise INFLATES jitter near criticality rather than deflating it. **How to kill:**
measure the mean advance rate and its variance versus gamma under `s-only`; if the
inflation tracks the flattening drift it is this, otherwise something else changes
sign between 0.30 and 0.45.

~~**T12 remaining as first posed: the two-target race.**~~ Absorb on `|n_X - n_Y| >= thr` OR
`n_B <= m`, tuned to ~50/50; `cme.splitting_probability` already supports it.
**Caveat the review did not raise:** under `delta-only` the pool has no noise of its
own, so whether it reaches `m` is driven only by delta's fluctuations feeding the
drift -- not forced the way `Var(n_B)` is, but one step removed. Write the
pre-committed criterion before running. **And if this too comes back needing
`span(delta)`, the right conclusion is not "one stiff direction" but that every
observable reachable with this absorbing-set machinery is a question about delta,
which would be a limitation of the INSTRUMENT and should be recorded as one.**

~~**T12 as first posed: does any observable of this system require the POOL's
noise?**~~ §25
establishes that the requirement varies by observable but only DOWNWARD -- P(error)
needs `span(delta)`, MFPT needs nothing -- and an observable needing the EMPTY
subspace is weak evidence, since drift-dominated means are means. The strong form
needs an observable requiring `span(s)` or `span(delta)+span(s)`. **How to kill:** two
candidates, exact references, no definitional shortcut. (i) `Var(T)` at first passage:
the ODE gives 0, so it is a pure-noise quantity like P(error), and the jitter could
come from crossing-direction diffusion (-> delta) or from fluctuations in the RATE OF
ADVANCE, which `n_B` sets through the recruitment propensities (-> s). Neither is
guessable, and `E[T^2]` solves `Qtt m2 = -2T` beside the existing MFPT solve. (ii) A
two-target race absorbing on `|n_X - n_Y| >= thr` OR `n_B <= m`, tuned to ~50/50,
which `cme.splitting_probability` already supports with no new code. **If both come
back needing `span(delta)` that does NOT establish one stiff direction** -- it means
three saddle-dominated observables were chosen.

**A rule that falls out of this thread and is worth keeping:** *a projection arm is
informative only if its answer cannot be computed from the projection algebra alone.*
`Var(n_B)` under `s-only` is computable in advance, therefore dead as evidence --
though it is a good unit test of `project`. Stationary covariance is a 2x2 Lyapunov
equation, so solve it before running and only run what the algebra says discriminates.

~~**T11-REFINED-a: does the signal/bookkeeping split survive when the roles are not
obvious?**~~ **PARTLY ANSWERED -> §24.2. The SUBSPACE-level split holds at n = 3; the
within-subspace question is still open because the arm that tests it was invalid.**

At n = 3, with `gamma = 0.60 * gamma_c(3)` matching §21.4 and an exact CME reference:
`bookkeeping-only` (80% of the noise variance, all in the blank pool) is
**categorically 0 in all eight cells**, and `signal-only` (19-20%) carries the
exponent to 6.6-8.5%. So the split is not a two-species convenience -- it survives
where the signal subspace is two-dimensional.

**Wording consequence: "the signal COORDINATE" becomes "the signal SUBSPACE".** At
n = 2 the difference subspace is one-dimensional, so §24.1 could not distinguish the
two; n = 3 can, and it is the subspace that matters.

Honest cost: `signal-only` runs 13-28% high at eps = 0.25 and 88-157% high at
eps = 0.40, against n = 2's 2-18%. Dropping the pool's noise costs MORE with more
rivals, and worsens with the barrier in both systems.

~~**T11-REFINED-b: does the signal subspace decompose within itself?**~~ **MEASURED ->
§24.3. NO.** With the rival skew forced to a fixed value independent of Omega, the
parity artifact is gone and `decision-only` is smooth at every Omega -- and it
**never recovers**, under-estimating in all fifteen surviving cells across two skews,
averaging 0.666 of the exact answer at skew 2 and 0.716 at skew 4. Keeping only the
champion-vs-rivals direction loses about a third of the failure probability.

**So the coarse-graining picture is now bounded on both sides.** Removing the pool's
noise (80% of the variance) is categorical; removing one of the two signal directions
is NOT categorical but costs a third; only the full difference subspace reproduces the
answer. For a simulation: discard the bookkeeping noise entirely and keep the
exponent, but there is no further cheap truncation INSIDE the signal directions.

*Suspect, not established:* the order-statistic reading -- that the champion loses to
the best of two NOISY rivals, which sits higher than the best of two separated by
drift alone. Consistent with the sign and with the growth in the rival count, but
nothing here isolates it from other consequences of deleting d2.

**Limits on §24.3, carried forward rather than buried.** (i) The ratio drifts +0.063
between skew 2 and skew 4, so P7's DIRECTION is robust and its MAGNITUDE is not.
(ii) Exponents are not fittable from those runs: `setup_skewed` absorbs its
divisibility remainder into the MARGIN, so realised eps wobbles +-11% and the CME
reference goes non-monotone in Omega -- the very corruption
`approximation_hierarchy_nwinner.setup` was written to prevent, reintroduced while
fixing a different artifact. **The fix, for whoever runs this next:** absorb the
remainder into the SKEW instead, which §24.3 shows matters weakly, rather than the
margin, which §24.2 shows matters a lot. (iii) One cell is excluded by the stated gate
|full/CME - 1| < 0.25, because the full CLE control is itself 1.53x the CME there.

~~**T11-REFINED-b as first posed, left visible:**~~ The arm
that would answer it -- `decision-only`, retaining only the champion-vs-rivals
direction -- is INVALID as built. With the rival-vs-rival noise zeroed, `X2 - X3`
evolves deterministically from its initial value, so the whole arm is decided by an
integer-rounding parity in the start state: rivals tied -> p = 0 by construction,
rivals differing by 1 -> p = 0.054. Same class as §12's floor-division artifact and
§17's MI integer bias; more trials cannot fix it. **How to kill:** seed the rivals
with a deliberate fixed asymmetry independent of Omega parity and re-run. If
`decision-only` then recovers the exponent, an n-winner race needs noise in ONE
direction of n, which is a statement about simulation cost; if it fails while
`signal-only` works, the whole difference subspace is required and the split does not
decompose further.

**Also recorded, an error caught before it produced a result.** The first n = 3 run
copied `gamma = 0.30` from §21's AM defaults, but `gamma_c` falls with n and
`gamma_c(3) = 0.2023`, so it measured a network with NO LANDSCAPE -- width 0,
threshold at its floor, the state sitting symmetric. It would have read as clean
results. The experiment now takes `--gamma-frac` of `gamma_c` and refuses to run at
zero width. **A parameter copied across networks is a parameter that has not been
checked against the network it landed in.**

~~**T11a: is the cliff a property of restoration or of AM?**~~ **RUN AT n = 3 ->
§21.4. The cliff survives.** The ODE reports exactly 0 in all four cells where the
truth spans 2.1e-2 to 1.1e-1, and every noisy level lands within ~5% of exact
(kappa ratios 0.960 SSA / 0.947 CLE / 1.038 tau). §21.3 is not an AM artifact.

**T11a-i, still open: does the CLE's error GROW with alphabet size?** Not resolved.
The SSA is the exact chain, so its 4.0% deviation is the noise floor of that run,
and the CLE's 5.3% sits barely outside it -- where at n = 2 the floor was 0.2-0.5%
and the CLE's 2.5% was clean. Reading 5.3% against 2.5% as a trend would be reading
a difference smaller than the anchor's own scatter. **How to kill:** the same run at
n = 3 with the statistics of the n = 2 one (the cost is the CME at larger Omega,
C(Omega+3,3) states), or n = 4 where gamma_c = 0.068 makes the landscape shallower
and any breakdown should show sooner.

**T11b, WITHDRAWN as posed -> §21.3a.** At 40,000 trials the CLE's exponent matches
exact (ratios 0.980 and 1.001 at n = 2; per-cell agreement within 2 sigma at every
dt at n = 3), so the +2.5% and 5.3% deviations were sampling noise and there is no
sign to explain. Note the first draft of this withdrawal blamed dt-convergence,
which was itself a 4,000-trial artifact -- the CLE is converged AND correct.

**T11b-i, CLOSED AS UNRESOLVED -> §21.3a.** Measured over a 16x range in dt at
60,000 trials per cell. The excess is positive in 5 of 6 pooled cells (~+1.7% +-
0.5% at the well-resolved eps) but NON-MONOTONE in dt, scattering by more than its
binomial error -- so it is neither the monotone decay discretisation gives nor the
flat plateau a prefactor gives. The exponent is untouched at every dt, which is the
part §21.3 depends on. Closed rather than pursued: three rounds in, a ~2% effect
that will not hold still across dt needs a positivity-preserving integrator rather
than more samples, and that buys nothing for the cliff. What the entry asked:

**~~T11b-i, the residue worth keeping:~~** across all eight n = 2 cells the CLE sits
above the exact p (1.018 to 1.045, every z positive, ~3 sigma combined for a +3%
uniform excess). A uniform factor on p is a PREFACTOR effect and leaves the
exponent alone -- consistent with the CLE being the right diffusion limit with a
slightly wrong amplitude, which is what §15's saddle-curvature agreement would
predict. One step size and 3 sigma. **How to kill:** two more step sizes at 100k
trials; if the excess is dt-dependent it is discretization, if it is flat it is the
Gaussian-noise prefactor. What the original entry said:

**~~T11b, open: what sets the CLE's sign?~~** I predicted the CLE would OVERestimate
the failure probability from the 1-D birth-death comparison
(`ln r > 2(r-1)/(r+1)`), and it underestimates it by 2.5%. The scalar intuition
does not survive two dimensions and I cannot currently derive the sign. A
2-D Hamiltonian calculation would settle it, and would say whether 2.5% is
universal or accidental.

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

  **What survives, and it is now smaller.** §12's slope was 0.74 because part of
  its exponent was simply wrong: `kappa` scaled the restoring gain with gamma and
  left the diffusion at its gamma = 0 value, when `D_0(gamma) = (1+gamma)/9`
  (§15). Refitting §12's own 216 cells with `kappa = (3/2)(1-2g)/(1+g)` lifts the
  pooled collapse from R^2 = 0.933 to 0.960 and every slope toward 1 (pooled
  0.742 -> 0.783). **So one of the three "prefactors" turned out to be an
  ingredient left out of the exponent** -- not a Laplace amplitude, and not a
  mystery.

  What is left is smaller and still real: the corrected slope is 0.783, not 1, and
  still non-monotone in gamma (0.81 / 0.68 / 0.51 / 0.68). That residual most
  likely belongs to §12's *own* saddle point, which minimises a sum of two
  exponents and keeps only the minimum. **Do not merge it with Q9a** -- that
  reflex is what this section is a monument to.
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
- **Q4. Structured (asymmetric) landscapes.** **Opened -> §16.**
  `networks/am_asymmetric` tilts the two autocatalytic branches by `beta` while
  keeping every reverse at gamma x its own forward, so the cycle affinity stays
  `-3 ln gamma` and beta is a clean second axis costing no thermodynamic force.
  Three things are now measured:

  * a saddle-node at `beta_c(gamma)` past which the network answers X whatever it
    is shown; `beta_c` collapses as gamma rises (0.998 -> 0.050 over
    gamma = 0.05 -> 0.45), so near gamma_c a few percent of rate mismatch kills it;
  * the bias sits in the **saddle**, not the attractors, so reading a tilt off
    attractor positions reports no tilt at all;
  * mutual information falls monotonically in |beta| and **the penalty grows with
    Omega** -- the first place in this file where more molecules reliably hurt.

  ~~**Q4a: is the symmetric design still optimal for a BIASED source?**~~
  **ANSWERED IN FORM -> §17. No -- and the rule is the predicted one.** beta* > 0
  at every prior, and

      ln(e-/e+) at beta*  =  0.7625 * ln(p/(1-p)) + 0.018      R^2 = 0.999867

  i.e. **tilt until the log-ratio of the two error probabilities matches the
  prior log-odds**, the ratio holding to 2.7% over a 7.3x range in log-odds with
  a predicted-zero intercept measured at 0.018. That is the first statement in
  this project about how to *build* the chemistry rather than how it behaves.
  The optimal tilt is gentle: beta*/beta_c runs 0.016 -> 0.114, so the tilt that
  helps is 1-10% of the tilt that destroys the device.

  ~~**Q4a-i: is the coefficient 1?**~~ **ANSWERED: YES -> §17.3.** It was 0.76 at
  Omega = 200 and the 4x range available then could not pin the limit. The
  obstacle was cost, not principle: `beta*` was being found by maximising the
  exact mutual information, ~40 CME solves per cell. **Inverting the optimality
  condition** -- asking "given beta, which prior makes it optimal", a root find in
  p rather than an optimisation over beta -- yields the whole line from ONE error
  sweep, about 20x cheaper, and it reproduces the direct measurement to 0.0001 in
  slope. That buys the range: r rises 0.587 -> 0.911 over Omega = 100 -> 600 with
  the intercept collapsing 0.055 -> 0.0013 and R^2 reaching 0.999999.

  Over the 10x range the extrapolations still spread (0.978 / 1.029 / 1.131 /
  1.701), so two tests that need no chosen ansatz settle it: the **deficit**
  `1 - r = 27.9 * Omega^-0.902` is a clean power law at R^2 = 0.9978 (a limit
  below 1 would flatten it), and a **free-limit** fit gives `r -> 1.037 +- 0.030`,
  putting 1.0 at 1.23 sigma. So the rule is exactly

      ln(e-/e+) at beta*  =  ln( p/(1-p) )

  with a finite-population correction `1 - 28 Omega^-0.90`. The intercept,
  predicted zero, is 0.0002 at Omega = 1000.

  **The methodological point is worth keeping separately: the measurement was not
  impossible, it was badly posed.** An optimisation was being run where a root
  find would do, and the cost of that was mistaken for a limit of the method.

  **What was wrong in the prediction this entry used to carry.** It read: *"the
  optimum is set by matching the two barriers' difference to ln(p/(1-p))/Omega,
  which would make the optimal tilt shrink like 1/Omega."* The **matching
  condition is right**; the **1/Omega is not reached** (the measured exponent of
  `beta* ~ Omega^-x` is 0.51 / 0.60 / 0.68 / 0.75 across consecutive Omega, still
  climbing); and the barrier difference it names is the wrong K -- attractor-to-
  saddle barriers give `K = 0.874` at gamma = 0.35, which predicts a slope 1.9x
  off, because the errors here start from a *biased input*, not from the
  attractor. Those are different barriers and the entry conflated them.

  **A refinement that made things worse (recorded, because it is a pattern).**
  A more careful derivation keeping the log factors gives
  `p*h(e+) = (1-p)*h(e-)` -- equal deficit contributions from the two symbols.
  It is **refuted**: the ratio runs 1.21 / 1.48 / 1.90 / 2.80 / 4.05 across
  p = 0.60 -> 0.95, systematically worse the more extreme the prior, while the
  cruder form it was supposed to improve stays flat to 2.7%. See §5.1's
  three-parameter fit for the same lesson.

  **Still open (Q4b): what does the tilt cost in dissipation?** The affinity is
  beta-independent by construction, but the EP *rate* is force times flux and the
  flux moves. Unmeasured, and §16 is careful to claim only the force.

  Note §9.2's closed-form EP identity **fails** on these networks, which is
  exactly why `thermo.entropy_step` exists as the general primitive.
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

## 3a. State of the base

Twenty-six load-bearing closed forms are re-derived from the current code by
`experiments/verify_base.py` and checked against their published values; it runs as
part of the suite (`tests/test_verify_base.py`), so drift is a failure rather than a
later discovery. **All 26 agree** as of the §21 work.

This exists because the suite alone does not answer the question. The suite proves
the code is self-consistent *with itself*; the audit proves it still agrees with
what is **written down**, and the two come apart the moment a behavioural function
changes underneath sections already published. §15 changed
`information.wall_coefficient`, which feeds `predicted_exponent` and
`crossover_omega` and therefore §12 — those numbers are deliberately left as first
published, with the refit tabulated in §15.2 and `wall_coefficient_gain_only` kept
callable so they still reproduce.

Sampling-derived results are excluded from the audit on purpose: they cannot be
re-derived exactly, and their agreement is the suite's job, not this one's.

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
| **"The deterministic gain is 18.6× too small to explain the width growth"** (§94.1(c)). | Two independent faults. (1) The implied g² was solved out of σ₂/σ₁ **assuming σ_intr is stage-independent**, when §95–§96 show the operating point moves between stages and the intrinsic width follows it — so the "transmitted noise" being extracted was mostly stage 2's intrinsic width at a *different mean*. (2) The recursion itself is Thattai & van Oudenaarden (2002) Eq. 13 and carries a **time-averaging factor** τ = β_x/(β_s+β_x) that mine lacks; since §92–§93 measured the two timescales *equal*, τ = 1/2 and the deterministic side halves, making the gap 37× rather than rescuing it. §96's account reproduces the widths with **no transmitted-noise term at all**. §99(c). |
| **"The composition maximum sits at τ_up = τ_down, where the upstream correlation time meets the downstream response time"** (§93). | The knob was a **global speed scaling on stage 1**, which multiplies its relaxation time *and* its escape time by the same factor — so it could not say which was matched (rule 9), and the collapse onto s_up/s_dn was already an algebraic identity. Measured on Ω instead: at s_up = s_dn = 1 the upstream's correlation time is **0.6980** and the downstream's rail relaxation time **0.1511**, **4.62× apart**. The crossover sits at matched *speed scalings*, not matched *times*, because the two stages carry different boundary conditions. The crossover and the interior maximum (4.4419 at ratio 1, against 4.2102 slow and 1.5795 fast) are measured and stand; only the identification is withdrawn. §100.1. |
| **"The reflected upstream's correlation time is its rail relaxation time"** — implicit wherever §92–§98 called stage 1's clock a relaxation time. | The slowest mode of a stage reflected at its saddle is **diffusive relaxation across the whole box [saddle, cap]**, not the local curvature at the rail: measured gap 1.43–1.62 against |f'(r₃)| = 6.6195, flat in Ω and **not converging** to it. The wall does not merely delete the escape channel, it installs a timescale of its own — 4.6× slower — and that timescale is a property of where the boundary was put. §100 P4, predicted to converge to 6.6195 and refuted. |
| **"The reflected chain measures the fluctuation-transfer channel"** — implicit throughout §92–§98, and the basis for reading its penalty as a per-stage cost. | It **overstates that channel by 1.49× at D = 2 and 2.22× at D = 3** (Ω = 30, t₀ = 2.0; 1.63× and 2.64× at Ω = 14), because the wall bounces back probability that would have escaped and parks it against the saddle where it drives the next stage hardest. Separately it omits the failure channel entirely, which is 61–69% of the true error at D = 3. **The two errors have opposite signs and cancel to a 16–21% net**, which is why the construction looked sound: §94 reports 6.349 against a true total of 7.388. The inflation grows with depth and so does the omission; nothing makes the cancellation persist. §101. |
| **"P(last stage low) is comparable across depths at any window"** — assumed wherever §94 and §101 compare D = 2 with D = 3. | False below t₀ ≈ 2. Stage 1 seeds from its quasi-stationary law and is already spread; every later stage seeds as a **delta at its rail**, so until that delta propagates a deeper chain hands its last stage a *cleaner* input. At Ω = 30, t₀ = 0.5 the reflected D = 3 error is **ten times smaller** than D = 2's, and at Ω = 14, D = 3 stage 3 reads (2.7450, 0.6846) against stage 2's (2.6542, 0.7643) — the last stage cleaner than the one feeding it. §92.1(a) one level up: there the artifact reversed a trend in upstream speed, here it reverses the trend in **depth**. It removed §101's first headline (a 6.66× ratio at 91.8% contamination) before it was written down. §101.1. |
| **"The chain's operating point can be composed from the transfer map alone"** — §103's first version, withdrawn before publication. | It predicted 2.8501 → 3.0946 → 3.1173, each stage sitting *closer* to its rail than the one feeding it, against a measured 2.8042 → 2.6698 → 2.6113 that falls. **Not a 16% error but the wrong direction**, and the cause is structural: §91 built the coupling to be neutral at the rail, so `F(r₃) = 3.1827000000` against `r₃ = 3.1827` — **the rail is a fixed point of the transfer map**, and iterating `⟨F(x)⟩` can only converge toward it at any depth for any input. The degradation comes from each stage's own finite-Ω depression `d_intr = μ₁ − r₃`, which §96 includes and this omitted. Restoring it gives +3.45% and +4.32%. Same family as §96.1 — a term dropped from a composition — but caught by a **sign** rather than a magnitude, which no tolerance would have flagged. §103.1. |
| **"The early-start account is refuted — the offset needed grows in proportion to t, so it is a rescaling of the rate and not a shifted clock"** (§104). | **The refutation, not the account, was wrong.** §104 inverted `p(t) = form(k, t + Δ)`, modelling a head start as a **longer window**; a head start is the convolution `p(t) = 1 − e^{−kΔ}(1 − e^{−kt})/(kt)`. Under §104's form the required Δ spans **18.90×** across the five windows and does grow with t — a true reading of the wrong function. Under the correct one it spans **2.17×**, and a single Δ reconciles the curve. §105 then *derived* that Δ: rate-weighting the conditional occupation time below the saddle-node at x_up\* = 1.5795 (where the downstream descends 8.7× slower) gives **0.0686 against a fitted 0.0486**. The parameter-free worst residual falls 8.95% → 5.80% and the one-signed bias vanishes. **Rule 14 exactly — a withdrawal is a claim, and this one discarded a correct account.** §105. |
| **"P(stage 1 low at t) = 1 − e^{−k₁t}"** — §102's first ingredient, inherited by §103–§105. | A free stage has **no absorbing boundary**: its low-basin occupancy saturates at the stationary weight π_low, not at 1. π_low falls **0.9057 → 0.5247** over Ω = 14–70, and the one-way model correspondingly runs from 7% low to **59% high**, spanning 1.71× and crossing 1 near Ω = 35 — where §102 happened to validate it, accidentally accurate to 2.5%. Not a transient: at Ω = 55 the error *grows* with the window (1.2453 → 1.4267 over t₀ = 2 → 8). With `P = π_low(1 − e^{−λt})` the model is flat to **2.5% across four decades of probability**. Rule 9 on the axis the whole arc is parameterised by. §106.2. |
| **"At least two mechanisms are in play, one of them provably not an averaging effect"** (§108.3, withdrawn one commit later). | The Jensen argument was sound and the premise was an artifact. §108 measured stage 2 in a chain that seeds it as a **delta at its rail** and compared that against averages built from **stationary** laws. In a single stage with the upstream *pinned* — no fluctuating input, nothing to average — the delta seed carries a **25.2× transient** over a 16× window against the QSD seed's 1.12×, suppressing by 1.68× at t = 2. With matched seeding `true/k(⟨x⟩)` runs 0.9893 → 2.0257, and the one sub-unity cell is 1.1% below while the instrument itself under-reads a measured 5.3%; corrected, **1.0445 → 2.1387, none below**. The averaging family is not excluded. §109.1. |
| **"below 0.99 ⇒ §108's premise survives"** — §109's own P3 check. | A **fixed tolerance on a quantity whose own instrument carries a five-times-larger residual**. It printed *survives* off a single cell at 0.9893 when the QSD seed at t = 2 is 5.3% short of its converged value. Rule 20, in the test written to adjudicate a rule-14 withdrawal. Sixth broken criterion of the arc (§99.1, §101 P1, §102 P4, §104 P2, §104's refutation, this). §109.1. |
| **"stage i's escape rate is `escape_rate(Ω, mus[i])`"** — §102's `predict`, inherited by §103–§105. | `escape_rate(Ω, x_up)` is the rate of a stage whose **upstream** sits at x_up, so each stage was keyed to its own operating point rather than its input's. **Stage 1 has no upstream at all** — its rate is the free spectral gap, which §102's own P1 showed equals `escape_rate(Ω, r₃)`. The corrected indexing reproduces the measured gap **exactly at all seven Ω**; the as-coded one is **12.7% high at Ω = 14**. Correcting it improves the closure at every Ω without removing §106's divergence, so the published agreements were slightly lucky and the conclusions survive. §106.3. |

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

---

**T14: the sampling floor is the binding constraint, and it is now movable.**
Every tail claim here is capped by sampling rather than physics -- §24.1c's
conclusions rest on a nine-cell pattern precisely because per-cell error (~5% at
p ~ 0.0094) is comparable to the effects compared, and the founding claim concerns a
switch that errs at 1e-15 while every measured number sits between 1e-1 and 1e-2.
Exact CME reaches small probabilities but its state space grows as ~Omega^2/2;
sampling handles any Omega but floors at 1/N. Large Omega AND small probability is
reachable by neither.

**MEASURED -> §26.** An exact Gillespie SSA written in MLRift reproduces the exact
CME at the gate cell -- P(error) 0.232450 vs 0.233847 (1.5 sigma), MFPT 13.42224 vs
13.41524 (0.37 sigma) -- at **6,211 trajectories/sec on one core against Python's
351, an 18x single-core speedup**. That puts p ~ 1e-6 at 4.5 CPU-hours, or ~11
minutes threaded, against Python's 79 hours. **The GPU is not needed for the first
deep-tail result**, and that ordering is deliberate: the CPU implementation validated
against an exact reference is what a later GPU port must be diffed against.

**T14-a, open: does the barrier law survive five decades?** §22 measured the barrier
dying as `(gamma_c - gamma)^1.9745` and §12/§15 the collapse `P ~ exp(-Omega dW)`,
both over roughly ONE decade of probability. The transistor analogy that motivates
this project lives at 1e-15. **How to kill:** rerun the collapse at fixed gamma over
Omega wide enough to span 1e-2 to 1e-6, now affordable. If the exponent drifts with
decade the collapse is local and the extrapolation to real restoration is
unsupported; if it holds over five decades it is the strongest quantitative claim
this project has. **Gate:** every cell that HAS an exact CME reference must match it
before any cell that does not is reported.

**T14-b, open: the deep-tail coordinate result.** §24.1's categorical failure
(`s-only` reporting exactly 0) was measured where the truth is 1e-1 to 1e-2. Whether
"exactly zero" survives at 1e-6 -- or whether the pool's noise contributes a floor
that only becomes visible far into the tail -- is untested and is the sharper version
of the §24 claim.

~~**T14-a: does the barrier law survive five decades?**~~ **ANSWERED -> §27. It
survives 6.53 decades, and needed no sampling at all.** Exact CME from
`P = 9.833e-2` at Omega=40 to `2.927e-8` at Omega=620 (193,131 states, 13.6 s).
eps-controlled fit `ln P = -0.024904*Omega - 0.2145*Omega*(eps-dev) - 1.7280`,
R^2 = 0.999533; eps-corrected local slopes mean -0.025283 with 7.1% scatter and a
first-half/second-half difference of +7.3%, i.e. within noise. **The one-decade
extrapolation §12 and §15 made is supported at this gamma and eps.**

**Two corrections fall out and both are mine.** (i) §26 argued the deep tail needed a
faster sampler because the CME state space grows as ~Omega^2/2. It does, but Omega^2/2
is SMALL -- I assumed a limit instead of measuring one, and motivated a sampler on it.
(ii) The raw local slopes bounce by +-40% and that is the §24.3 integer-lattice
artifact again: realised eps swings 10.2% because the start margin `d0` is an integer.
Controlling for it takes R^2 from 0.998639 to 0.999533 and halves the max residual;
the threshold's rounding contributes nothing.

**Where the sampler DOES belong:** the exact state space is `~Omega^(n-1)/(n-1)!` in
species count -- fine at n=3, but `n_winner_reversible` at n=4 is ~Omega^4/24, i.e.
6.7e7 states at Omega=200 and out of reach. Multi-species networks and large Omega are
where sampling is the only instrument. §26's validated SSA was re-checked 265x deeper
than its gate (Omega=200, 8.7375e-4 vs exact 8.800e-4, 0.33 sigma).

**T14-c, open: does the collapse rate match `kappa(gamma) delta*^2` in ABSOLUTE
terms?** §27 measures the slope -0.025283 empirically. §15's closed forms give
`kappa = (3/2)(1-2g)/(1+g)` and `delta*`, so the barrier is computable with no free
parameter -- and rule 16 says a law that is only ever fitted is never tested. **How to
kill:** predict the slope from `kappa delta*^2` and compare the ratio; the sweep is
already cheap enough to repeat at several gamma, which turns one ratio into a curve.

~~**T14-c: does the collapse rate match `kappa delta*^2` in ABSOLUTE terms?**~~
**TESTED -> §28. At gamma = 0.30 yes, to 4.1% over 5.21 decades with NO free
parameter. At gamma = 0.15 no, off by 24-36%, and the drift is real.**

The prediction is the scale-function exponent `d(lnP)/dOmega = -2 V(x0)` with
`V = int mu/D` on the slaved 1-D manifold; its near-saddle limit is `-kappa x0^2`.
**Normalisation trap, stated because it would have read as physics:** `breaking_mode`
is a UNIT vector so `D_delta = 2 D_0`, making the near-saddle exponent `-kappa x0^2`
and NOT `-2 kappa x0^2`; the wrong factor turns 12% agreement into a 2.3x failure. A
numerical guard (exact integral / quadratic as x->0) returns 1.0000 at every gamma.

gamma=0.45 is EXCLUDED on a stated criterion: its collapse spans 0.40 decades and
three decades would need Omega ~ 3000 (~4.5M states). Under-determined.

**P4 fired and it is not finite-Omega.** Refitting on the upper Omega half moves the
ratios FURTHER from 1 (gamma=0.15: 1.240 -> 1.362; gamma=0.30: 1.041 -> 1.063), so
the small-Omega defence fails.

**T14-c-i, open: does the discrepancy collapse against the timescale separation?**
*Suspect:* the prediction's 1-D slaved reduction, whose error §24.1a measured as
shrinking monotonically with `3(1+2g)/(1-2g)` -- **5.6 at gamma=0.15 vs 12.0 at
gamma=0.30**, i.e. worse reduction exactly where the prediction is worse. **How to
kill:** measure at gamma in {0.20, 0.25, 0.35} and test whether discrepancy collapses
against separation. Two points fit anything; five make it a curve or kill it. **If it
does not collapse, the gamma-dependence lives in `kappa delta*^2` itself and §15's
closed form does not survive the absolute test** -- which would be the most consequential
withdrawal in the project, since §15 corrected §12 and everything after leans on it.

**T14-c-i -> §28.1: the excess does NOT collapse; it is scatter, not drift.** Five
gammas with criteria fixed in advance (>=2 decades, P >= 1e-12, eps-controlled):
ratios 1.155, 1.118, **0.952**, 1.061, 1.044. **§15's closed form predicts the
collapse slope to +-15% at every gamma across 5.5-9.9 decades each, with no free
parameter** -- the parameter-free test of the project's central law, and it passes.

**Two corrections.** (i) §28's gamma-drift was partly the INSTRUMENT: `p_cme` returns
`1 - split`, and §28's gamma=0.15 sweep reached ~1e-17, below double-precision
cancellation. With a 1e-12 floor its ratio is 1.155, not 1.240. §28 stands as printed
per rule 7 but its "24-36% off" overstates the failure. (ii) **My analysis script
printed `excess = 1.069 sep^-1.134, R^2 = 0.9966` and called it a clean power law. It
fitted 4 of 5 points** -- gamma=0.25's excess is negative so `log` dropped it, and the
dropped point is the one contradicting the pattern. No power law is claimed.

**The 1-D-reduction story is UNSUPPORTED.** The excess is not monotone in gamma and
changes SIGN at gamma=0.25; a reduction error tracking timescale separation cannot do
that. P1 and P3 both fail.

**T14-c-ii, open: is the +-10% scatter physics or unmatched grids?** Each gamma used a
DIFFERENT Omega list, so cells are unmatched and each carries its own realised-eps
pattern -- and that wobble alone moved §27's raw local slopes 40%. **How to kill:**
re-run with Omega chosen by ONE rule for every gamma (same decade span, same P range),
so the only difference between cells is gamma. If the scatter survives matched grids
it is physics and `kappa delta*^2` has a residual gamma-dependence; if it collapses,
§15's closed form is exact to the precision of the test and the whole T14-c thread
closes clean.

**T14-c-ii -> §28.2: MATCHED GRIDS REVIVE THE MECHANISM, and §28.1's central claim is
withdrawn.** One rule for every gamma -- same P window (1e-2 to 1e-6), ln P equally
spaced, 12 cells -- gives ratios 1.131, 1.137, 1.080, 1.055, 1.013. **All above 1**,
so §28.1's sign flip at gamma = 0.25 was a GRID ARTIFACT, not physics, and its
conclusion ("scatter, not drift; the reduction story is unsupported") is withdrawn.

What remains is a monotone DRIFT toward 1 as gamma rises, **11.7% against ratio
uncertainties of ~0.48%, i.e. 24x the measurement error**. With all five points and
nothing dropped: `excess = 5.987 * sep^(-2.031)`, R^2 = 0.887. §28's 1-D-reduction
story is back, and it is now consistent with §24.1a's independent measurement that the
reduction's error shrinks with the timescale separation.

**So §15 SURVIVES the absolute test.** The closed form predicts the collapse slope
with no free parameter to 1.3% at gamma = 0.35 and within 14% everywhere across ~4
decades per gamma, and the residual lives in the 1-D REDUCTION used to evaluate it
rather than in `kappa delta*^2` -- vanishing exactly where that reduction becomes
exact.

**Sequence worth remembering: §28 proposed the reduction story, §28.1 killed it on
unmatched grids, §28.2 revived it on matched ones.** The killing step was wrong for a
reason that had nothing to do with the hypothesis.

**T14-c-iii, open: what is at small gamma that a power law does not capture?** The
first two gammas differ by 0.006 in excess where `sep^-2` predicts ~37% between them
(7.00/5.57 = 1.26). The exponent -2.03 should be read as "roughly inverse-square", not
a measured constant, and R^2 = 0.887 over a factor of 3 in `sep` is not a strong test.
**How to kill:** extend to gamma = 0.05 and 0.10, where `sep` falls to 3.7 and 4.5 and
any real power law must predict a large excess. If the excess saturates instead of
growing, the reduction error has a floor and the inverse-square reading is wrong.

**T15 -> §29: §24.1's categorical zero is a THEOREM, and that bounds what §24 may
claim.** The projected CLE was ported to MLRift, gated against all four Python arms at
§24.1's cell, and run to 2,000,000 trajectories with `s-only` still exactly 0. The
reason is algebraic, not statistical:

    b_delta = delta * [ c_het*B - c_hom*(s-1) ]     -- no additive term

because `f1/r1` change X and Y equally (contributing 0), `f2/f3` give `c_het*B*delta`,
and `r2/r3` give `-c_hom*delta*(s-1)`. Verified numerically: `b_delta/delta` constant
to 4.6e-16 over a 25x range in delta. So under `s-only`, `d(delta)/dt = delta*g(t)`
gives `delta(t) = delta_0 exp(int g)`, and **sign(delta) is CONSERVED**. No number of
trajectories could ever have found a crossing.

**§24.1's "categorical" was right; §24's REASON for it was not.** §24 read the zero as
"noise in the wrong subspace is worthless for barrier crossing". The mechanism is that
AM's drift is exactly proportional to the signal coordinate. **A network whose
`b_delta` carries an additive term would give a small but NONZERO `s-only`**, so §24's
"seven-eighths of the noise and as categorically wrong as none" is a statement about
AM, not about coarse-graining in general. The QUANTITATIVE half of §24.1 is untouched:
`delta-only` at 11% of variance recovers to 2-18%, and `uniform 11%` is wrong by
17-770x while retaining reduced delta-noise.

~~**T15-a, open: does the identity hold at n >= 3?**~~ **CLOSED by §30 — yes, and in a
stronger PAIRWISE form; but the arm named below was the wrong one.** See the T15-a
block at the end of this file. The original question stands here as written:

**T15-a, open: does the identity hold at n >= 3?** §24.2's `rivals-only` also returned
categorical zeros, which the same structure would explain. **How to kill:** test
whether the drift along each signal direction of `n_winner_reversible` is exactly
proportional to that direction's coordinate, as it is at n = 2. If yes, §24.2 needs
the same qualification and the "subspace" result is likewise AM-family-specific. If
no, the n = 3 zeros have a different cause and are genuine evidence for the coordinate
reading.

**T14-c-iii -> §28.3: the eps axis CONFIRMS the attribution; the power law does not
survive.** Self-calibrating matched grids (every cell bisects Omega to span 1e-2 to
1e-6, 12 points, one rule everywhere) over gamma in [0.05,0.35] and eps in {0.35,0.50};
four cells exceeded the Omega<=900 cap and are reported, not dropped.

**P2, the discriminating test, passes decisively.** At fixed gamma the ratio moves only
2.6-3.0% between eps = 0.35 and 0.50, against a 16-point drift across gamma. The
residual is a property of the slaved MANIFOLD, not of where on it the trajectory
starts -- exactly what the 1-D-reduction attribution predicts and what a defect in
`kappa delta*^2` would not.

**P1 fails.** The excess keeps growing toward small gamma (0.1598 at gamma=0.10, so
saturation is dead too) but not as `sep^-2`: refitted over the wider range the exponent
is **-2.53, not -2.03**, and an exponent that moves with the window is not an exponent.
Decisively, **the excess reaches zero** (0.0044 at gamma=0.35) and no power law in
`sep` can. `excess = 0.2240 - 0.6276*gamma` fits at R^2 = 0.9905 crossing zero near
gamma ~ 0.357, but that is six points and two parameters -- a description, not a law.
**§28.2's exponent is withdrawn as a description; its ATTRIBUTION survives on P2.**

At gamma = 0.35, eps = 0.35 the parameter-free prediction is within **0.4%** -- the
closest agreement in the project.

**T14-c-iv, open: what vanishes near gamma ~ 0.357?** The excess crosses zero there and
one cell (gamma=0.35, eps=0.50) reads 0.975, i.e. slightly below 1 -- but with the
eps-spread at 3% the sign is NOT resolved. **How to kill:** measure at gamma in
{0.38, 0.41, 0.44} where the linear description predicts a clearly negative excess. If
the excess goes reliably negative, the reduction error changes sign and no current
account explains that; if it flattens at zero, the linear description is wrong and the
reduction error simply dies out. Note the Omega cap bites hardest here (shallow
collapses need large Omega), so the sampler of §26 may be the only instrument.

**T15-a -> §30/§30.1/§30.2: the identity is general, it is PAIRWISE, and it sorts
§24.2's arms -- but not the one T15-a guessed.** Derived by hand from the count-level
propensities, then verified against the network's own stoichiometry at n = 2..6 x four
gamma each (4,600 pairs; worst residual 4.4e-16 against the traffic that must cancel,
median 2.1e-16 strict, worst strict 2.9e-13 where |n_i - n_j| = 1 and the bracket is
itself near zero):

    d(n_i - n_j)/dt = (n_i - n_j) * (k/Omega) * [ n_B - sum_{l!=i,j} n_l
                                                  - gamma*(n_i + n_j - 1) ]

The conserved quantity is not "the signal coordinate" but **sign(n_i - n_j) for every
pair independently**, so a projection is covered exactly when it starves some pairwise
difference direction:

  * `bookkeeping-only` starves all C(n,2) -> the whole ordering freezes -> P = 0 is a
    THEOREM at every n. Demonstrated where it could not be barrier height: a champion
    leading by ONE count, exact CME 0.597, full noise failing 59.9%, and the arm still
    exactly 0 with zero flips in 40,000 trajectories. **§24.2's "subspace" result needs
    the same qualification §29 forced on §24.1 -- for this arm.**
  * `decision-only` starves d2 -> sign(n_2 - n_3) frozen. **This is §24.2's Omega-parity
    trap, which was filed under "floor-division artifacts" and blamed on integer
    rounding.** The rounding only chose the initial condition; the conservation law is
    what made it fatal. Confirmed in-run: 0 rival flips against 6331 champion flips.
  * `rivals-only` starves NONE of them, so the theorem does not cover it -- and its
    zero turns out not to be categorical at all (§30.2).

**Two mechanisms proposed and withdrawn in the same session, both about `rivals-only`
(rule 17).** First: it should be a mere sampling zero, since no sign is conserved. It
is not -- 440,000 trajectories, no champion flip, never closer than 0.51 of the initial
gap. Second: the exact erosion identity `du/dt = u*Gbar - (1-gamma)*delta_23^2/(4*Omega)`
(residual 1.3e-15, pinned in the suite) says the champion's noise-free mean margin has
a sink quadratic in rival spread, gated by sign(G_23) -- and the eps sweep ordered the
arm over four decades in exactly that sequence. **But eps sets the champion's margin,
which sets BOTH the barrier and G_23**, and I read the ordering as support before
noticing (rule 9). Holding the margin fixed and moving G_23 through n_B instead:
corr(ratio, G_23) = **-0.9957**, near-perfectly OPPOSITE, with the arm failing 18.2% at
G_23*Omega = -23.46 where the mechanism forbade failure. What breaks the confound is
that **the two sweeps confound it in opposite directions** and the paired ratio tracks
the barrier in both.

**The lesson is rule 16's, one level over.** The identity survived because it was
derived and checked in absolute terms to 4e-16; the mechanism died because it was
inferred from a monotone table. A table that orders perfectly is not evidence for the
quantity you happened to order it by.

**T15-b, open: is "the drift is multiplicative in every signal coordinate" what makes a
network restoring at all?** The identity says the deterministic flow can amplify or
contract a lead but can never REVERSE one, so in this family every error is a
fluctuation in a difference direction and never a drift. That is a candidate definition
of a restoring element and the closest thing this project has to a statement of the
digital abstraction: noise off the signal axis cannot flip a bit, whatever its
amplitude. **How to kill:** construct a restoring network whose pairwise drift carries
an additive term -- `am_asymmetric` and `am_fueled` are the two in the repo most likely
to, since neither is symmetric under the exchange that produces the cancellation. §29
predicts its `s-only` analogue would then be small but NONZERO. If every network that
restores has the multiplicative form, the identity is constitutive rather than
incidental; if one restores without it, the definition is too narrow and the general
claim in §24 stays dead.

**T15-b -> §31: CONFIRMED, and it is the first mechanism in this arc to survive its own
kill test.** §29 predicted that a network whose signal drift carries an ADDITIVE term
would give a small but nonzero `s-only`. `am_asymmetric` has one, exactly:

    d(delta)/dt = delta*(k/Om)*[n_B - gamma*(s-1)]
                + (k*beta/Om)*[n_B*s - gamma*((s^2+delta^2)/2 - s)]

verified against the network's own propensities at six beta, worst residual 1.8e-14.
The second term vanishes identically at beta = 0 and grows with beta, so beta is a knob
that turns off and on the very thing the zero is claimed to depend on.

At a barrier held to 7.1% (exact CME 0.323 -> 0.327, `full` flat at 0.324-0.333),
`s-only` goes from an EXACT zero at beta = 0 to 3.2% at beta = 0.30 -- while retaining
**0.891 of the variance at every beta, identical to three decimals**. The noise
amplitude is constant across the sweep and only the drift structure changes, so this is
not an amplitude effect. The claim is categorical, and no barrier change can turn an
exact zero into a nonzero.

**Two cautions kept in the record.** beta = 0.02 and 0.05 are BOUNDS (< 3.0e-5 at
100,000 trajectories), not zeros -- only beta = 0 is a theorem zero, and the transition
turns on below the instrument's floor rather than sharply. And P3 as written demanded
the effect under BOTH start rules; it failed, because the `fixed` rule's cells at
beta >= 0.05 are inadmissible (the saddle slides past the start and the noiseless ODE
already fails) and its one admissible beta > 0 cell sits under the resolution floor.
Confirmed under `matched`, UNTESTABLE under `fixed` -- not contradicted.

**The P4 control was load-bearing.** In the four inadmissible cells `s-only` reads
0.597, 1.000, 1.000, 1.000. Without the ODE check that reads as a spectacular
confirmation; it is `am_asymmetric`'s systematic tilt deciding the answer
deterministically, which is not a restoration failure at all.

**Why this one survived when §30's two did not:** it was derived algebraically, checked
in absolute terms to 1e-14, and then given a parameter that switches it off at fixed
barrier. §30.1's mechanism was inferred from a monotone table and died to the first
sweep aimed at it. Rule 16, one level up.

**T15-c/T15-d -> §42/§43: the identity is a THEOREM about exchange symmetry, and
conservation laws have nothing to do with it.** T15-b confirmed that the multiplicative
drift is what makes a network restoring, but neither §30 nor §31 could say which
structural feature *produces* the multiplicative form — every network in the AM family
has exactly one conservation law, so the two candidate explanations were perfectly
confounded (rule 9, at the level of the whole family rather than one sweep).

§42 broke the confound by construction. Cofactor pairs consumed and regenerated by the
signal-moving reactions add conservation laws without touching exchange symmetry: at
**two laws and at four**, the identity holds to 1.6e-13 and 4.5e-14. The load-bearing
row is the one that differs from the two-law arm **in a single rate constant** — same
species, same reactions, same orders, same laws — and there it fails by 1.9e1. So the
answer is symmetry, not conservation.

§43 then found that "symmetry" proves *more* than §30 measured, and proves it
algebraically. An antisymmetric polynomial in two variables is divisible by their
difference, so for **any** exchange-symmetric mass-action network,

    b_i - b_j = (n_i - n_j) * P(n),    P symmetric.

**delta = 0 is an invariant manifold of the deterministic flow at any reaction order and
any number of conservation laws — including none.** §30's constant *ratio* is the extra
condition that P not depend on delta, which needs the pair to enter at total degree <= 2;
a cubic pair term separates them (divisible, ratio affine in delta^2 with the coefficient
matching k/(4 Om^2) by hand to 1.00000008). On 200 randomly generated symmetrised
networks: 200/200 divisible, worst residual 3.0e-17, and **157 of them conserve nothing at
all**. Unsymmetrised: 2/200, so the probe has power.

**What stays a hypothesis.** The theorem gives no-reversal, NOT amplification — that
needs sign(P) > 0, which is a measurement and is negative below the separatrix. And it is
a statement about the deterministic drift only; the CME chain crosses delta = 0 freely,
which is exactly the "every reversal is a fluctuation" reading and nothing more. Both
limits are recorded in §43 rather than left for a later retraction.

~~**T15-e, open: is sign(P) > 0 characterisable, or is it just measured?**~~ **CLOSED by
§53: NO, and the two halves of restoration turn out to have different PREVALENCE, not just
different logical status.**

  * **P is the symmetry-breaking eigenvalue.** Differentiating b_i - b_j = delta*P at
    delta = 0 gives P(symmetric state) = d(b_i-b_j)/d(delta), which at a symmetric fixed
    point is T7/§14's lambda. Verified to 5.8e-10 over 52 states. §43's P and T7's lambda
    were the same object for twenty sections without anyone noticing.
  * **Closed form: P(symmetric fixed point) = (1 - 2*gamma)/3**, worst deviation 8.9e-12
    over 11 gammas from 0 to 0.9. Vanishes at gamma_c = 1/2 (the pitchfork), and at
    gamma = 0 equals 1/3 = 1/(2n-1) at n = 2 -- T7/§14's lambda(n) exactly. Sits in fixed
    ratio to §12's wall coefficient: P/kappa = (1+gamma)/4.5.
  * **The attractors lie on P's zero set** (|P| < 8e-16 at delta*), necessarily, since
    d(delta)/dt = delta*P. Not the separatrix -- §43 already made that delta = 0.
  * **sign(P) is NOT combinatorial.** 17 of 188 random topologies flip it under rate
    changes alone, spanning -38.26..+2.04 at one fixed state, and AM flips at gamma_c. One
    counterexample suffices. The 91% that do not flip say topology carries a lot of
    information without determining the answer.
  * **THE SHARPEST NUMBER: divisibility 200/200, amplification 21/200 = 10.5%.** Almost
    every exchange-symmetric network preserves the sign of a lead; almost none grows it.

**The pre-registered verdict rule for the kill test was badly designed and §53 says so:** it
demanded a MAJORITY of topologies flip before declaring non-combinatoriality, which is the
wrong logic for a universality claim -- one counterexample refutes it. The criterion, not the
data, is what changed.

**What this means for the founding claim.** Restoration = divisibility + sign(P) > 0. The
first is structural, universal, and needs no conservation law. The second is a rate-constant
condition that is RARE. **A restoring element is a tuned object, not a topological one**, and
"the transistor is a near-ideal restoring switch" is a statement about tuning rather than
about circuit topology.

**T15-f -> §54: WHERE the combinatorics runs out.** §53 said sign(P) is not combinatorial;
§54 says precisely where. Grouping reactions into mirror pairs and taking the X-heavy member,
with d_r = S_X(r) - S_Y(r),

    P = sum_pairs d_r * c_r * [ O_r (xy)^q sum_m x^m y^(p-q-1-m) ],   every bracket >= 0

verified to 1.9e-14 on AM, am_cubic and 120 random networks. **sign(P) is a rate-weighted
sum of INTEGERS with non-negative weights**, so:

    all d_r <= 0         -> P <= 0 everywhere, whatever the rates  113 nets, 0 amplify, 0 violations
    all d_r >= 0, some>0 -> P > 0 everywhere                         3 nets,            0 violations
    mixed                -> the rates decide, and ONLY here        183 nets, 33 amplify

**The cross-check against §53 is the sharp one.** §53 measured 17/188 topologies flipping
sign(P) under rate changes, for a different purpose. §54 predicts only MIXED ones can:
**0 of 168 unanimous topologies flipped, 24 of 132 mixed ones did (18%)** -- and 44% mixed x
18% = 8% reconciles with §53's independently measured 9%.

**d_r > 0 is NOT autocatalysis.** `B + X -> 2X` gives d = +1, but so does `2X + Y -> 2X + B`
where S_X = 0 and S_Y = -1 -- X catalysing Y's DESTRUCTION amplifies a lead just as well. The
governing notion is positive feedback on the DIFFERENCE.

**AM decomposed, and what gamma_c is.** Recruitment `B+X->2X` has d = +1 at weight k; its
reverse `2X->B+X` has d = -1 at weight gamma*k; the disagreement channel `X+Y->2B` is
SELF-MIRROR and contributes exactly zero. Evaluating the brackets gives **P = k(b - gamma*s)
-- §30's identity, recovered term by term.** So AM is a MIXED network, **gamma is literally
the weight on the contracting term**, and **gamma_c = 1/2 is where the mixed sum changes
sign.** §53's P = (1-2*gamma)/3 is that sum at the symmetric fixed point.

**The self-mirror zero is one fact wearing three hats:** §30's first cancellation, §51's
discovery that rho does not appear in mu, and §54's zero contribution. Three sections, one
structural cause.

**And it explains §53's "tuned, not topological" with a mechanism.** Of 11.3% amplifying,
**11.0 points come from MIXED networks and 0.3 from topology-guaranteed ones** -- the
guaranteed class is ~1%. Nearly all restoration in this family is rate-tuned rather than
topology-forced.

**T15-g -> §55: the decomposition predicts gamma_c(n) AND T7's lambda(n), with no fit.**
§54 worked only at n = 2; §30 proved the pairwise identity for every n, so the decomposition
follows. Sorting `n_winner_reversible` by the (i,j) swap gives exactly one d = +1
(recruitment, weight k, bracket b), one d = -1 (rev-recruitment, weight gamma*k), (n-2) more
d = -1 (rival disagreements X_i + X_k -> 2B), one SELF-MIRROR zero (X_i + X_j -> 2B), and
(n-2) p=q pairs that cancel. Summing at a symmetric state:

    P = k[ b - gamma*2x - (n-2)x ]

which is **§30's published n-winner bracket term by term**. Two absolute consequences:

  * **P(symmetric, gamma=0) = 1/(2n-1)** at n = 2..6, worst 8.9e-12 -- T7/§14's
    symmetry-breaking eigenvalue lambda(n), recovered from a stoichiometric decomposition
    that knows nothing about it.
  * **P vanishes exactly at gamma_critical(n)** (worst |P| 3.5e-12, three cells exactly 0),
    which `gamma_critical` finds by an independent bracketed root-find. gamma_c(3) =
    0.202226 against the published 0.2023.

**And it explains the radix penalty combinatorially.** The contributing counts are
**1 amplifying against n-1 contracting**, plus n-2 that vanish. gamma_c(n) falls with n
because the single d = +1 must outweigh one more contracting term at each step. **The radix
penalty (§3/§6.1), the symmetry-breaking eigenvalue (T7/§14) and the critical drive are the
same integer count evaluated in three places.**

**A reporting bug, mine, caught by a derivation written first.** P4's first pass read
1 / 2n-3 / 0 against a predicted 1 / n-1 / n-2. A p = q pair has an EMPTY bracket sum and
contributes zero, but `mirror_pairs` records its d as an arbitrary +-1 because "the X-heavy
member" is undefined when p = q. P1-P3 were unaffected -- they use the decomposition, where
the empty sum correctly gives zero. **The bug was in the reporting, not the algebra**, and
§55 records both counts.

**Two DIFFERENT routes to zero, and §54 saw only one:** self-mirroring (the reaction maps to
itself) and pairwise cancellation through p = q (distinct mirrors, opposite d, identical
propensity). Conflating them is what produced the bug.

**T15-h -> §56: THE RESTORATION TRICHOTOMY. Capability is combinatorial; realisation is
tuned.** §53 said "restoration is tuned, not topological". That is half right, and §54's
decomposition contained the other half.

> **Qualified by §62.** "Tuned" now has a criterion -- realisation is decided by
> `sum_r c_r d_r B_r(x*) > 0` at the symmetric steady state -- and the convexity below is
> a statement about **capability only**. The realising set of rate vectors is **not** convex:
> §62.2 gives two restoring c in AM's own classes whose sum does not restore. Nothing measured
> in §56 changes; what changes is that its cone must not be read as covering realisation.

P is LINEAR in the rate constants: with v_r(x) = d_r B_r(x), P(x) = <c, v(x)>. So a network
fails to restore exactly when <c, v(x)> <= 0 for every accessible x, and that is closed
under addition and positive scaling:

    **THEOREM. The non-restoring rate vectors form a CONVEX CONE** -- the polar of the cone
    generated by {v(x) : x accessible}.

Measured: 0 violations of "c1, c2 non-restoring => c1+c2 non-restoring"; linearity to
2.7e-16. Two corollaries complete the classification:

    all d_r <= 0          -> the WHOLE positive orthant is in the cone. Never restores.   0/157 capable
    all d_r >= 0, some >0 -> always restores.                                             2/2  capable
    mixed                 -> CAPABLE for some c; the failing c form a convex cone.        34/34 capable

**A network can restore for some rates <=> some d_r > 0. It restores for every rate <=> all
d_r >= 0.** Capability is decidable from stoichiometry; only realisation is tuned.

**Compactness does not bite.** The domination argument needs room to scale x, and a
conservation law confines the state to a compact simplex -- but 35/35 mixed networks are
capable there too. **Third time (with §42, §43) that conservation structure turns out not to
matter to this question.**

**§56 ALSO CORRECTS §54, and the correction sharpens it.** `classify` was counting p = q
pairs, whose bracket is empty and whose recorded d is an arbitrary +-1. §55 found exactly
this in the n-winner counter and **the fix was never propagated to `classify`**, which §54's
whole table rests on:

    class      §54 published    corrected
    all<=0          113             211
    all>=0        3 (1 amplifies)   7 (7 amplify)
    mixed       183 (33 amplify)   63 (27 amplify)

**148 of 232 misclassified.** And §54's one soft spot -- "two of three all>=0 networks read
P ~ 0 ... a threshold artifact" -- **was not an artifact, it was this bug**; corrected,
all>=0 amplifies 7/7 with no excuse attached. The flip cross-check strengthens: 0 of 263
unanimous flip, 24 of 37 mixed do. **Unchanged:** the decomposition (1.9e-14), the AM
term-by-term recovery of §30, and §55's gamma_c(n) and lambda(n) predictions -- all of which
use the brackets, where the empty sum correctly gives zero. The bug was in the CLASSIFIER,
not the algebra.

**Where this leaves the founding question.** 23% of random symmetric networks are capable of
restoration; **77% are forbidden by stoichiometry and no tuning can rescue them.** Of the
capable, ~49% realise it at randomly drawn rates. So "the transistor is a near-ideal
restoring switch" splits into a structural claim that is decidable and a tuning claim that is
not.

**~~T15-i, open: is the restoration cone's boundary computable in closed form?~~ -> §62:
YES, AND IT IS NOT THE CONE'S BOUNDARY.** The question as posed asked for the facets of the
polar of conv{v(x)} over the accessible set. That is the boundary of **capability**, and it
is the wrong object: the state that decides is not an arbitrary accessible x but the
**symmetric steady state x\***, the decision point the dynamics occupies. §55 had already
measured P there -- (1 - 2 gamma)/3 for AM, vanishing at gamma_c = 1/2 -- without noticing it
was a criterion. So:

> **RESTORES  <=>  sum_r c_r d_r B_r(x\*) > 0**, one linear inequality in the rate constants
> at fixed x\*, for any exchange-symmetric mass-action network.

Verified against the DYNAMICS rather than against another algebraic rule: 120 networks, 60
predicted restoring and 60 decaying, **0 disagreements**; reproduces AM's published
gamma_c = 1/2 to 1.1e-15. The trichotomy is now constructive -- pick d_r, then solve a linear
feasibility problem in c at the self-consistent x\*.

**And §56's convexity does NOT transport, which the original framing would have got wrong.**
The cone theorem needs P linear in c at FIXED x; x\*(c) moves with c. §62.2 exhibits two
restoring rate vectors **in AM's own two classes** whose sum does not restore
(P = +5.5e-2, +6.98e-2, sum -2.14e-1; confirmed by the ODE at x6236, x5659, x1.1e-12).
**Capability is a convex cone; realisation is not convex.**

**T15-k -> §63: THE RESTORATION THRESHOLD HAS A WIDTH, AND IT GOES AS Omega^(-1/2).**
The CME generator has an antisymmetric sector (exchange X<->Y commutes with Q whenever §43's
premise holds), so lambda_A -- its leading eigenvalue -- is the EXACT stochastic counterpart
of §53's Jacobian P, with no simulation and no first-passage definition. The chain is ergodic
at every finite Omega, so lambda_A < 0 always: **restoration at finite size is not a sign
change but a metastability, and the transition has a width.** Measured two ways, one anchored
to §62's exact deterministic rate and one mentioning no reference at all, the width goes as
**Omega^(-0.49 .. -0.51)**, stable across metastability levels and independent of rho.
Device reading: **to halve the blur on a chemical switch's threshold, use four times the
molecules.** The extrapolation beyond the measured decade is NOT resolved (three ansaetze
spread 90% at Omega = 1000) -- the exponent stands, the value past the data does not.

**T-THM -> §65: THE SYMMETRIC RESTORATION THEOREM, and the literature that was never
checked.** §43/§54/§55/§56/§62/§63 are one statement in six clauses: (1) invariance --
b_X - b_Y = (n_X - n_Y) P(n), so the sign of a lead is a deterministic invariant; (2)
decomposition -- P = sum_r c_r d_r B_r(x) with d_r = S_X - S_Y integral and every B_r >= 0;
(3) capability is combinatorial -- restores for SOME c iff some d_r > 0, and the non-restoring
c form a convex cone; (4) realisation is ONE linear inequality, sum_r c_r d_r B_r(x*) > 0, and
the realising set is a cone but NOT convex; (5) at finite Omega the generator block-
diagonalises and lambda_A is the exact stochastic P, negative always, with the transition
rounded to width ~ Omega^(-1/2).

**Prior art, checked at last (§65.1).** Clause 1 is the folk theorem on flow-invariance of
Fix(Z_2) in equivariant dynamics (Golubitsky, Stewart & Schaeffer 1988) -- §43's framing of it
as a discovery is withdrawn. **The nearest neighbour is the SPONTANEOUS MIRROR-SYMMETRY
BREAKING / homochirality literature**, not CRNT: an exchange-symmetric network whose two
species are enantiomers, asked when it leaves the racemic state, is that field's central
problem from Frank (1953) through Ribo & Hochberg to Montoya, Cruz & Agreda (Life 9:74, 2019)
and the Listanalchem tool. Clause 4 is, in substance, linear stability of the racemic state,
which they do routinely -- but their conditions come out SEMIALGEBRAIC AND HARD TO SAMPLE,
linearised where possible via Clarke's stoichiometric network analysis. Adjacent: Craciun-
Feinberg's species-reaction graph, Joshi & Shiu, and arXiv:1002.1054 on switching in mass
action networks by linear inequalities.

**So the claim is not that a criterion exists but that it has a FIXED SIGN STRUCTURE** --
clause 2 -- which separates topology from rates, reduces capability to the signs of a list of
integers, and attributes amplification to named reactions.

**T-THM-a -> §70: SETTLED for the two papers read, and it cost clause 4.** Montoya, Cruz &
Agreda (Life 9:74, 2019) Theorem 1, the "MM-condition": a racemic steady state is
symmetry-breaking iff the characteristic polynomial of A_Om - B_Om is unstable. For a
Z_2-symmetric Jacobian [[A,B],[B,A]] the spectrum splits into A+B and A-B, so **A_Om - B_Om IS
the antisymmetric block** and their Theorem 1 IS §62's clause 4 -- stated for k enantiomeric
pairs where CRNL has k = 1, so MORE general, not less. §65 called it "in substance" their
routine practice; it is more specific than that and §65 understated it.

**What is not in those papers:** no per-reaction signed decomposition (the criterion is
Jacobian-eigenvalue-based via Hurwitz-Routh on symbolic determinants), and no purely
combinatorial capability result -- their one structural theorem is negative and case-by-case.
So **the contribution narrows to clauses 2 and 3**, which is exactly what §65 hedged before
reading them.

Also worth recording: two independent routes reach the same object, and §63's lambda_A is the
exact CME counterpart of their deterministic A_Om - B_Om -- clause 5, where nothing comparable
was found.

**SCOPE: two papers, one of them poorly extracted. "Not found in two papers" is weaker than
"novel", and the homochirality field is seventy years old. §65 claims priority for nothing and
that still stands.**

**~~T-THM-a, open and the novelty claim rests on it: is clause 2's decomposition already in the
mirror-symmetry-breaking literature?** NOT SETTLED -- abstracts are not enough. **How to
kill:** read Montoya/Cruz/Agreda and the MATCH algebraic-analysis paper properly and check
whether the racemic-stability condition is ever written as a nonneg-weighted sum over
reactions carrying combinatorial signs. Until then §65 claims priority for nothing.

**T-DEPTH-a -> §74: WHAT DELTA COSTS -- closed elements pay affinity, open ones pay
material.** §73 collapsed the founding question to "what does it cost to place rails Delta
apart?", and §12's sigma = f*Delta convention had hidden it: under that convention
Delta/sigma = 1/f identically and every element gives the same depth. Holding sigma fixed in
PHYSICAL units instead:

  * **AM is conservative**, so delta* <= 1 and D_max SATURATES -- at sigma = 0.15 no drive
    pushes a conservative element past **D_max = 9.50e9**. Going from A = 8.99 to A = 27.6
    buys a factor of 8; going to infinity buys another 1.04.
  * **Schloegl is open**, and its affinity ln[e1 e2/e3] is EXACTLY invariant under r -> lambda r
    (measured drift 0.00e+00 over 2.5 decades) while Delta scales linearly. At matched affinity
    it passes 1e18 by scaling rails, paying in MATERIAL: rails at r1*Omega and r3*Omega.

**So closed and open elements pay different currencies for the same quantity, which is why §67
and §68 found no substrate-independent price -- they were pricing the wrong thing.** The cost
of depth is in Delta, and what Delta costs depends on whether the element is conserved.

**T-DEPTH-g -> §80: §75-§79 PRICED THE WRONG FAILURE MODE.** An element in a cascade can
fail two ways -- misreading its input (the Gaussian readout of the intrinsic noise,
eps ~ exp(-eta*Omega)) or spontaneously escaping its rail during the stage time
(eps ~ t*exp(-A*Omega)). Both are exponential in Omega, so the SMALLER exponent dominates
absolutely, and the arc never compared them. Measured on §75's own element:
**A = 0.026046 against eta = 0.059537**, so escape wins by e^432 at Omega = 12800 and readout
would bind only for stage times below 6e-188. **§75's own premise -- that the physical cascade
is chemically coupled -- selects the regime where §75-§79's eps is not the physical one.**

**What survives.** §76's D_max = c*/eps is mathematics about a binary channel and does not care
where eps comes from; only which eps to use changes. §77's eta is the right KIND of statement
(one exponent governs everything) with the wrong number for this regime; its measurement and its
transfer conclusion are unaffected, since A is just as substrate-specific. **§78's headline
survives with a different formula**: the escape action is an integral of the propensity
densities, A = -int ln(mu/lam) dx from saddle to rail, giving 0.026047 against the exact MFPT's
0.026046 -- still no master equation, but the quasipotential integral rather than the LNA
variance. SCOPE: that integral form is 1-D; a multi-dimensional element needs Hamilton-Jacobi,
which is harder than the LNA and not done here.

**§71/§72 are unaffected and this was checked.** They used an EXTERNAL channel 4.8x wider than
the intrinsic width, where ln eps_read = -6.15 against escape at -41.70, so readout binds by 35
nats. Engineered wiring and chemical coupling are physically different cascades.

**WHY IT WAS MISSED FOR FIVE SECTIONS, and the generalisable form:** the §73-§79 chain was
reductions each correct given the last, and the error entered at the ONE step that changed the
physical setup rather than simplifying the mathematics -- §75 replaced an external channel with
the intrinsic rail width, a change of MODEL, and that is exactly where a second failure mode
became available. **When a reduction replaces one physical mechanism with another, re-enumerate
the failure modes. Simplifying algebra cannot introduce new ones; changing what the model
represents can.**

~~**T-DEPTH-h, open: redo §77-§79 with A in place of eta.**~~ **-> §81. Done, and all three
redos give the same verdicts as their eta originals -- but §80 is not a relabelling.**

  * **§79 redone (out of sample) is the strongest version yet.** `A = -int ln(mu/lam) dx`
    predicts the exact mean first-passage time to **0.06% worst case on FOUR landscapes never
    used in §80**, including a different reaction order (the quartic 3X <-> 4X). §78's
    "deterministic side, no master equation" claim therefore holds for the exponent that
    actually governs depth, not just for the one that did not.
  * **§77 redone: A does not transfer.** 0.0101 .. 0.1358, a factor of 13, moving with
    landscape shape inside each substrate. Same conclusion as eta's, reached independently.
  * **A/eta is NOT constant: 0.2002 .. 0.7560, a factor of 3.78.** This is what makes §80 a
    correction rather than a change of name. eta is a local harmonic quantity at the rail; A is
    a global integral over the barrier; they respond differently to shape. Had the ratio been
    constant, every conclusion of §77-§79 would have carried over with a rescaled axis.

**T-DEPTH-h's own instrument failure, and it was inflating the headline.** AM at gamma = 0.20
first returned A = 5.7137, 3.8145, 2.8591 at Omega = 120, 180, 240 -- products with Omega of
685, 686, 686, i.e. exactly the ln(1e300) ~ 691 underflow floor of the stationary weights, so
"A" was 691/Omega and not a barrier at all. **That one cell was reporting A's span as a factor
of 282 instead of 13, and A/eta's as 7.79 instead of 3.78.** The tell is the shape, not the
size: a quasipotential barrier per molecule must CONVERGE in Omega, and anything decaying as
1/Omega is a saturated numerator. Floor cells are now excluded rather than fitted. Both verdicts
survived; both magnitudes did not, and the factor of 282 would have been the most quotable wrong
number in the record. This is rule 20's shape (never gate a converging quantity with a fixed
tolerance) turned around: **check that a quantity claimed to converge actually does, because a
saturated instrument produces a smooth, plausible, wrong series.**

**T-DEPTH-f -> §79: OUT OF SAMPLE. The ODE predicts systems the derivation never saw, and
§78's closing claim needs one qualification.** §78 was checked on the systems it was built from,
which is not a test (rule 16; §22 fitted a convolution for three subsections and was out by
3688x against an exactly-computable quantity). Predicting first from ODE + Lyapunov, then
checking against the exact CME:

  * AM gamma = 0.10 (unused value): 0.04% -> 0.02% -> **0.01%**
  * AM gamma = 0.35 (unused value): 13.96% -> 5.04% -> 1.98% -> **1.12%**
  * **Schloegl QUARTIC, 3X <-> 4X -- a DIFFERENT REACTION ORDER**, never touched by §73-§78:
    2.63% -> 0.58% -> **0.14%**
  * Schloegl with asymmetric rails: 0.70% -> 0.17% -> **0.04%**
  * AM gamma = 0.45, **chosen in advance to fail**: 10.76% -> 8.88% -> 23.72% -> 23.33%,
    NOT converging, 20.9x worse than the worst passing case.

All three out-of-sample systems converge; the case predicted to break, broke. That confirms
§78's P4 diagnosis on a system chosen in advance: the LNA fails on a shallow rail near gamma_c
and nowhere else.

**THE QUALIFICATION.** D_max ~ exp(Delta^2 Omega / 2V), so a relative error d in V becomes a
factor D_max^d in the depth. Measured: sigma errors of 0.01%, 0.04%, 0.14%, 1.12% become D_max
factors of 1.26x, 2.8x, 13.2x, 22.8x. **§78's "the exact CME is not needed to use it" holds for
eta and for ln D_max (predicted to 0.02-2%), NOT for D_max itself.** Anyone wanting a depth to
better than an order of magnitude still needs the master equation. That is a limit, not a
failure, and it matters because §72-§77 reported depths as numbers.

**PROCESS, recorded because it is about the rules themselves:** rule 20 was added two sections
before this one and I wrote exactly the gate it forbids -- a fixed 1% tolerance that called AM
gamma = 0.35 a MISS at 1.12% while its series ran 1.1396 -> 1.0504 -> 1.0198 -> 1.0112. Writing
the rule down did not stop the drafting error; it caught it on re-reading. And the quartic first
read non-monotone (0.966 -> 1.253 -> 1.047), which would have been evidence that the formula
breaks under a change of reaction order -- it was my landscape, with the rails only 2.4 sigma
apart so the basin second moment ate the other peak's tail. **Non-monotone convergence meant the
instrument, not the theory.**

**T-DEPTH-e -> §78: eta IS DERIVED, and §77's kill test named the wrong theory.** §77
proposed getting eta from WKB as the barrier action per molecule. WKB gives the ESCAPE
probability, but §75/§77's eps is not an escape -- it is the Gaussian readout of the rail's own
fluctuation, which is what a chemically-coupled cascade applies. The right object is the
LINEAR-NOISE variance: sigma^2 = V/Omega from the Lyapunov equation, so **eta = Delta^2/(2V)**.

Against §77's stored numbers with nothing fitted: 9.8797 vs 9.8813, 1.8329 vs 1.8346, 0.6781 vs
0.6813 for AM at gamma = 0.05, 0.20, 0.30, and 0.015577 vs 0.015617, 0.059537 vs 0.059574 for
Schloegl. **Worst 0.48%**, largest at gamma = 0.30 where the rail is shallowest and the LNA is
weakest -- as P4 predicted -- and eta(lambda=4)/eta(lambda=1) = 4.0000 exactly, as §75's
collapse requires.

**HOW MUCH IS CONTENT (stated because P2 is largely P1 restated):** eta_measured came from
eps = Phi(-Delta/sigma) with sigma from the exact CME, so the test is at bottom
V_LNA = V_exact plus algebra. The content is that **the linear-noise variance is the whole
story for this eps**, and that Phi's prefactor corrections do not disturb the slope at the 0.5%
level.

**THE FOUNDING QUESTION IS IN THE ODE.** Delta, the Jacobian and D = S diag(a) S^T are all
deterministic-plus-LNA. Chaining §73 -> §78: depth does not see the element's dynamics; the
element enters only through Delta; Delta/sigma is molecules at the rail; depth is c*/eps; eps
is fixed by eta; and eta = Delta^2/2V comes from the ODE and its linearisation. **The exact CME
was needed to ESTABLISH this and is not needed to USE it.** It also explains retrospectively
why every thermodynamic price failed to transfer: entropy production is not in the formula.
Affinity enters only by setting where the rails sit and how deep the wells are -- through Delta
and V -- and never in its own right.

**T-DEPTH-d -> §77: THE LAST FREE NUMBER IS eta, AND IT DOES NOT TRANSFER.** §76 left one
coefficient: eta = d ln(1/eps)/d Omega, nats of reliability per molecule. It EXISTS -- ln(1/eps)
becomes linear in Omega on both substrates, gated on convergence (rule 20): AM settles to 0.09%
and 0.01% at gamma = 0.20 and 0.05, Schloegl to 0.20% by Omega = 25600.

**P3's prediction held and it mattered:** Schloegl's ratios looked superlinear in §75's data
(2.288 against sqrt(4) = 2) and were predicted to be small-Omega contamination rather than
different physics, BECAUSE §75's exact lambda*Omega collapse cannot coexist with two exponents.
It straightens, and lambda = 4 gives 0.238185 / 0.059574 = 3.998 -- exactly the lambda-scaling
the collapse requires. The two sections agree where they must.

**eta spans a factor of 633** (AM gamma=0.05: 9.8813; Schloegl spread 0.6: 0.015617) and moves
with the landscape INSIDE each substrate -- gamma 0.30 -> 0.05 buys 14x, Schloegl's root spacing
0.6 -> 0.9 buys 3.8x. **It does not transfer, and it was never going to: eta is the one place
the chemistry survives.**

**CLOSING STATEMENT OF THE PROGRAMME.** Everything about a restoring element washes out of
composition -- reaction orders, symmetry, drive, dissipation, the shape of its commitment
function, even whether it is conserved -- **except eta, which is the element's own.** A
restoring element is a device for converting molecules into nats at an exchange rate set by its
landscape, and nothing else about it matters to what it can compute. The transistor's advantage
is not a special thermodynamic property; it is a good exchange rate.

**T-DEPTH-e, open: DERIVE eta.** Its landscape dependence is measured and not derived --
0.6813, 1.8346, 9.8813 at gamma = 0.30, 0.20, 0.05, plainly not linear. WKB should give it in
closed form as the barrier action per molecule. **How to kill:** compute the quasipotential
barrier for AM in closed form (§9.1 has the landscape; §63/§64 have the action machinery) and
check it against these three numbers ABSOLUTELY, no fit. And note the scope: eps here is the
Gaussian readout of the intrinsic noise, which §75 showed is what a chemically-coupled cascade
uses, but a different readout would move eta's absolute value even if the ratios survive.

**T-DEPTH-c -> §76: DEPTH *IS* THE ERROR RATE. D_max * eps = c* = 0.124266404564.** For a
binary symmetric channel the D-step bias is (1-2eps)^D and I falls through 1/2 at b* with
H((1+b*)/2) = 1/2, so D_max -> c*/eps with c* = -ln(b*)/2. Verified to 8.9e-9 over seven
decades of eps, and the SAME constant on every substrate measured -- Schloegl at two lambda and
two Omega, and §73's bare step function -- all within 0.8%, the residual being an asymmetry
correction rather than anything about the substrate. For asymmetric channels the ARITHMETIC
mean is the right eps (the decay rate 1 - e_hi - e_lo depends on the sum), exact in the
symmetric limit and drifting with the asymmetry ratio: 3x -> 1.050, 1000x -> 1.385.

**So §75's unification is an IDENTITY, not a coincidence: the depth ceiling carries no
information beyond the per-stage error rate**, and §12's whole cascade apparatus reduces to one
number the single element already determines.

**§12.1's reading of its own factor of 3 is WITHDRAWN.** It concluded "the factor ~3 is not a
prefactor -- it is a 7% error in the exponent". It is precisely a prefactor: exp(D^2/2s^2)/4 is
the Gaussian tail with its algebraic factor dropped, and the exact ceiling c*/Phi(-D/s) puts the
ratio at 4 c* sqrt(2pi) (D/s) -- RISING with D/s (3.204, 3.926, 4.758 at f = 0.45, 0.35, 0.28),
not a constant. §73's step-function measurements rise with it (2.71, 3.11, 3.80). The exponent
was never wrong.

**THE FOUNDING CLAIM, REDUCED.** Chaining §73 -> §76: depth does not see the element's dynamics,
the element enters only through Delta, Delta/sigma is fixed by molecules at the rail, and depth
is exactly c*/eps. **A restoring element is characterised, for every purpose this project set
out to measure, by ONE number: its per-stage error probability. Reliability is eps. Depth is
c*/eps. No second quantity, no trade-off, no thermodynamic price on either.** A deflationary
answer: the depth advantage of a restoring switch is not a property beyond its error rate but
the same property counted twice.

**T-DEPTH-b -> §75: ANSWERED, and it collapses §74 into something simpler. ONE CURRENCY:
MOLECULES AT THE RAIL.** In a chemically coupled cascade sigma is not a modelling choice --
stage 1's output species IS stage 2's input, so the inter-stage noise is stage 1's own rail
fluctuation, which the CME fixes. **My prediction that sigma ~ lambda was WRONG** (the
lambda^3-propensity argument with it): measured exponent 0.486. The rail fluctuation is
POISSONIAN in the count, sigma_n ~ sqrt(n_bar), so sigma_x ~ sqrt(lambda/Omega) and
**Delta/sigma ~ sqrt(lambda*Omega)**.

**The collapse is exact: Delta/sigma depends on the product lambda*Omega ALONE, to 0.000%
spread over 32x in each.** Scaling an element's rails and increasing its molecule count are the
same act. AM behaves identically -- rail width ~ Omega^(-0.50), so delta*/sigma ~ sqrt(Omega)
without bound.

**So §74 is deflated in both halves** (true at fixed sigma, which §74 flagged): AM's "maximum
composition depth" was a fact about the convention and not about conservation, and Schloegl's
"material" currency is not separate from Omega. What replaces it:

> **Depth is bought with molecules at the rail and with nothing else: D_max ~ exp(c*N_rail).**

That is the same currency and the same functional form as reliability (§1/§35/§38's exp(-c*Omega)),
so **the founding question's two halves collapse into one statement: a restoring element buys
both its error rate and its composition depth with molecules, exponentially, and buys neither
with free energy.** Affinity is a GATE, not a price -- required for rails to exist (§9.1, §68)
and setting where they sit, but once they exist more depth is a matter of counting.

**~~T-DEPTH-b, open and it is load-bearing: which channel-noise convention is physical?** §74's
entire content hangs on sigma being fixed in physical units rather than scaling as f*Delta. The
case for fixed sigma is that the inter-stage channel belongs to the wiring, not the gate; but a
real chemical cascade may well have noise scaling with the signal, in which case §12's
convention is right and no element beats any other. **How to kill:** build a two-stage cascade
where the coupling is itself chemical -- the output species of stage 1 IS the input of stage 2,
with no abstract channel -- and measure how the effective sigma scales with the rail
separation. That is a modelling question this project has never posed, and it decides whether
§74 is a result or an artifact of a convention.

**T-CASC-b -> §73: §72's TRANSFER IS DEFLATED -- the ceiling does not see the element.**
§72's own reduction has the chemistry entering only through p_cross, a monotone commitment
function near 1 at the low rail and 0 at the high one. Holding rails, channel and criterion
fixed and swapping ONLY that function: Schloegl exact gives 2.37/3.04/3.69, an exact Langevin
double well 2.71/3.04/3.79, a slope-matched sigmoid with NO dynamics 2.37/3.04/3.66, and a
**STEP FUNCTION with no element at all** 2.71/3.11/3.80. The sigmoid matches the chemistry to
three significant figures. **So the Schloegl-AM agreement is real and is not evidence about
substrates: both agree with a function that has no substrate.**

**Corrected claim, sharper than the one it replaces: composition depth is fixed by the READOUT
GEOMETRY -- rail separation against channel noise -- and the restoring element's entire
contribution is where it puts its rails.** Reaction order, symmetry, drive, chemostatting and
the shape of the commitment function do not survive into D_max once Delta is fixed. For the
purposes of depth a restoring element is a rail-placer.

**Consequence for the founding question:** composition depth is not a fifth currency, because
nothing is being bought. Three costs were tested across substrates and none transferred (§67,
§68); the fourth was not a cost. **No cost of restoration measured here is
substrate-independent, and the one substrate-independent quantity is not a cost.** This does
NOT say chemistry is irrelevant to restoration -- the element still sets Delta, whether rails
exist at all (§68), how fast it decides (§58) and how sharply it thresholds at finite Omega
(§63). All of that is upstream of the geometry.

**Instrument note (rule 10, and it nearly went the other way):** the Langevin column first read
0.68/0.14/0.01, comfortably OUTSIDE, which would have SAVED §72. It was wrong twice -- the
scale density is exp(U/D) and I wrote exp(2U/D), and a plain cumsum over an integrand spanning
e^9450 gave 0/0 that np.clip turned into a plausible column. **The bug pointed toward the
conclusion just published.**

**T-CASC-a -> §72: CLOSED. THE DEPTH CEILING TRANSFERS -- the first quantity in this
project that does.** Both routes through the matrix exponential wall out (§71.2), including the
second-eigenvalue route, which still needs exp(Qt)v. The reduction that works avoids the
exponential entirely and is exactly the regime §12's formula describes: in the fully-restoring
limit the chemistry contributes one number per rail, the probability of committing to the wrong
one, which is the EXACT birth-death splitting probability §61/§69 already compute in closed
form. O(cap) per cell. The element is asymmetric so eps_hi != eps_lo and the cascade is a binary
ASYMMETRIC channel, handled by a two-state chain.

All three channels saturate in Omega up to 28800, and the ratio of measured D_max to
exp(Delta^2/2 sigma^2)/4 lands on AM's published ratios: **2.71 / 3.04 / 3.75 against AM's
3.00 / 3.38 / 3.33** at sigma/Delta = 0.45 / 0.35 / 0.28. §12's ceiling, which AM overshoots by
~3x, is overshot by ~3x on an element with one species, no symmetry, chemostats and a different
reaction order.

**CAVEAT that keeps this at "tens of percent", not 10%:** this is the t_stage -> infinity
idealisation and AM's numbers are finite-t. The gap is measured -- at Omega = 400, f = 0.45 the
finite-time cascade plateaus at 5.10 against the saturated 6-7 -- so a like-for-like correction
moves the disagreements to roughly 25/25/10%. Still a different category from a quantity with
no counterpart (§67) or one varying 4x inside a family (§68).

**Standing statement: the THERMODYNAMIC accounting of restoration is substrate-specific; the
INFORMATION-THEORETIC accounting is not.** What makes a restoring element good at composing is
how far apart its rails are relative to the noise, and that transfers across chemistry with
nothing in common.

**T-CASC-a's kill test is REPLACED by §71.2, because the obvious route makes things worse.**
Going sparse -- exp(Qt)v by expm_multiply plus a Gaussian convolution instead of a dense
kernel -- was built and VALIDATED against the exact instrument to 3e-6 in I with D_max
identical to four decimals, and it is **50-90x SLOWER** (95s against 2.1s). The reason is
structural: propensities are extensive, so ||Qt|| ~ Omega (8.5e4 at Omega = 200, 1.5e6 at
3600), and every Krylov/Taylor exponential costs O(||Qt||) while dense scaling-and-squaring
absorbs the norm in log time but costs O(n^3). **Both wall out near 1e4 states, for opposite
reasons.** So T-CASC-a needs a different instrument, not a bigger machine: **D_max is set by
the second eigenvalue of the stage map C*K**, and a sparse partial eigensolve is O(nnz) per
iteration with no dependence on ||Qt|| at all. That is the reachable form of the measurement.

**T-CASC-a, open -> §71: the depth ceiling is the most promising transfer candidate found,
and it is NOT DECIDED.** The founding claim is about composition, and §12.1's ceiling
D_max ~ exp(delta*^2/2 sigma^2)/4 is INFORMATION-THEORETIC, not thermodynamic -- its formula
mentions only the rails against the channel noise. So unlike §67's cost per e-fold and §68's
affinity floor it has a reason to transfer, and the prediction was written that way on purpose.

Run on §67/§69's asymmetric Schloegl element with both gates holding (kernel stochastic to
7e-13 and holding its rail; with the channel OFF, I stays at 0.9637 through 600 stages, so the
decay measured is the channel's). **But D_max keeps rising with Omega** -- 13.6%, 37.6%, 84.8%
spread, still climbing at Omega = 1300 -- and §12's ceiling is an Omega-SATURATED quantity. Not
the same object yet, so the ratio is withheld rather than reported.

**The qualitative contrast is worth keeping and is NOT a result (rule 17):** Schloegl's
measured/predicted runs 0.97-1.70 against AM's 3.0-3.4 -- the same quantity within a factor of
~2, on a formula with no free parameter -- where §67's dissipation comparison missed by 20x
with no intensive definition available and §68's floor varies 4x inside one family. **How to
kill / what it costs:** f = 0.45 is decelerating toward ~5.7 and is reachable at Omega of a few
thousand; cap = 2*r3*Omega states with a dense expm puts the practical wall near 1e4 states,
and the narrower channels need much more.

**T-THM-b -> §67: PARTLY ANSWERED, and the answer is that the COST QUESTION does not
port.** Schloegl's model -- one dynamic species, chemostatted, no exchange symmetry to break,
and 1-D so the chain is exact -- was run through this project's rulers.

**What transfers: the affinity floor exists in both, and the values are close.** As the
landscape dies, Schloegl's cycle affinity tends to **ln 9 = 2 ln 3 = 2.1972**, independent of
where the degenerate point sits (derived, then confirmed by the engine's `cycle_affinity` to
8.9e-16). AM's is **3 ln 2 = 2.0794** (§9.1). **Not universal -- 5.66% apart -- but both are
ln(small integer) from chemistry with nothing in common**, and both fit `(pairs) x ln(pairs+1)`.

**What does NOT transfer, and not as a number: §38's cost per e-fold has no counterpart.** AM
is closed and conservative, so all dissipation belongs to the decision; Schloegl is
chemostatted and burns entropy sitting still (Sigma = 14293 against a cycle affinity of 2.26).
Two repairs were tried and both refuted -- a housekeeping subtraction came out NEGATIVE in
every cell, and a dimensionless cycles-per-molecule measure grows with Omega. The reason is
structural: a 1-D birth-death chain has zero stationary probability current, so the
non-adiabatic part is the system term, which is negative along the decision path. **The price
of restoration is defined relative to closed conservative bookkeeping and does not survive the
move to a driven device.**

**~~T-THM-c, open: is the affinity floor (pairs) x ln(pairs+1)?~~ -> §68: REFUTED, one
commit after it was opened, and by derivation rather than by fitting.** Generalising Schloegl's
autocatalysis to order p at fixed pair count and imposing a triple root gives
**A_c(p) = 2 ln[(p+1)/(p-1)]** -- confirmed against the engine's `cycle_affinity` to 4.4e-16
and x0-independent to 8.9e-16. At p = 3 that is 2 ln 2 = 1.3863, killing both the
(pairs)x ln(pairs+1) reading (predicts 2.1972) and the (pairs) x ln(max order) reading
(predicts 2.7726) with a single counterexample. **The floor runs 2.1972 down to 0.5026 within
ONE family at FIXED pair count and tends to 0 as p grows**, so §67's 5.66% agreement between
3 ln 2 and 2 ln 3 is a coincidence of two points. Worse, the "pairs" reading was malformed:
at three or more reversible pairs the cycle space is no longer one-dimensional and a single
affinity is not defined at all.

**Standing statement:** every restoring element examined has an affinity floor -- a nonzero
drive below which no bistable landscape exists -- but its value is set by the element's own
stoichiometry and varies by at least a factor of four. **There is no universal price of
admission.** §68 also closed a gap §67 left: A_c is a genuine MINIMUM (A rises monotonically
away from the degenerate point), which §67 had assumed rather than checked.

**T-THM-c (original entry, kept per rule 3), open: is the affinity floor (pairs) x ln(pairs+1)?** Two substrates fit it: AM with
3 reversible pairs gives 3 ln 2, Schloegl with 2 gives 2 ln 3. Two points fit almost any
two-parameter form, so this is a PATTERN, not a law (rule 17). **How to kill:** find a
restoring element with 4 reversible pairs and check for 4 ln 5 = 6.44, or one with 2 pairs and
a different topology and check it still gives 2 ln 3. A single counterexample settles it, so
do not write a tolerance.

**~~T-THM-b, open: does the theorem survive loss of exchange symmetry?~~ (original entry
below, kept per rule 3 -- note it asked about the THEOREM and §67 answered about the COST;
whether clauses 1-4 have an asymmetric analogue at all is still open.)** Every quantitative
result in this project -- G ~ 2 k_B per molecule per e-fold, the frontier exponent, gamma_c --
is measured on ONE symmetric network, while the founding object (an inverter driven toward
one rail) is ASYMMETRIC. Clause 1 fails outright without the symmetry (§42: residuals 1.9e1,
2.9e2). **How to kill:** run the rulers on a driven Schloegl switch (A + 2X <-> 3X, B <-> X,
chemostatted), which restores with no exchange symmetry at all and is 1-D so everything is
exact. Per this project's three failures of exponent transfer (§39.2, §46, §59), the prior is
that G != 2 and the frontier exponent differs; if instead they transfer, that is a
substrate-independent price of restoration and a far bigger result.

**~~T15-l -> §63.2: the exponent is nearly, but not exactly, §9.1's pitchfork.~~ ->
§64: WITHDRAWN, and the pitchfork is not excluded.** §63.2's estimator was gated on the DRIFT
between consecutive local slopes while its REMAINING error was 4x that drift, admitting a
point with ~6% error at gamma = 0.45 against ~0.1% at gamma = 0.38 -- a bias systematic in
gamma and signed exactly to lower nu. Extrapolated, nu moves UP by 0.15 as §64 predicted
before running. But three routes then disagree on where it lands (width of the transition,
which needs no extrapolation: 1.95-2.03; stationary distribution: 1.99; extrapolated action:
2.10-2.19), so **nu ~ 2 +- 0.1 and nothing here excludes 2.** §63.1's width exponent is
untouched -- it never used the action. The original entry, for the record:

**T15-l -> §63.2 (as first written):** w ~ Omega^(-1/nu)
where A(gamma) ~ (gamma_c - gamma)^nu, so -1/2 says nu = 2 -- which is what §9.1's pitchfork
predicts independently (barrier a^2/4b with a ~ gamma_c - gamma). Fitted with gamma_c held
FIXED at the value §62 proves exactly, **nu = 1.9496 +- 0.0026, flat and non-monotone across
nested windows, with 2 lying 9.7x the scatter away.** Note the agreement 1/nu = 0.5129 against
the width's 0.5133 is NOT an independent confirmation -- both are readings of the same surface
lambda_A(gamma, Omega) related by the scaling ansatz, so it confirms the ANSATZ, not nu.

~~**T15-n, open: WHY DO THE THREE ROUTES TO nu DISAGREE?**~~ **-> §84. CLOSED. nu = 2 exactly,
and the three routes were not disagreeing about a value -- each was measuring how far its own
window sat from gamma_c.**

The original kill test proposed running the extrapolation machinery on §61's closed-form 1-D
chain to measure the bias. **That was not needed, and the reason it was not needed is rule 16:
all three routes were FITS, and nobody had computed what the pitchfork normal form predicts.**

Eliminating the fast variable b (symmetric eigenvalue -(1+2g)) from the exactly-factoring lead
drift du/dt = k u [b(1+g) - g], with b0 = 1/3 exactly for every gamma:

    eps = (1-2g)/(3(1+g)),  b2 = -(1-g)/(2(1+2g)),  D_u(0) = 2k(1+g)/9
    **A = 9 eps^2/(4|b2|) = (1+2g)(1-2g)^2 / (2(1-g)(1+g)^2)**

k, D_u and the prefactor all cancel. **The noise had to come from the propensities, not the
drift** -- a normal form built on the drift alone would be exactly §83's error.

  * **P1**: k(1+g)eps reproduces `lambda_antisym` to 5.6e-17 at eight gammas; the normal form's
    delta* converges to the exact one, residual 0.1835 -> 0.0001.
  * **P2, ABSOLUTE**: against the exact first-passage action the ratio runs 0.8630 -> 0.9776
    over gamma = 0.30..0.46, monotone from below, residual/(g_c-g) = 0.685 -> 0.559. Zero cells
    excluded, zero fitted parameters.
  * **P3/P4**: the closed form's effective exponent is 1.8300 on §63.2's window [0.20, 0.45] and
    1.9898 on [0.48, 0.499]. The MEASURED A gives 1.942 -- squarely on §64's width route
    (1.95-2.03) and stationary route (1.99).

> ⚠ **DEFLATED BY §85.2, and heavily.** nu = 2 at a supercritical pitchfork is a STANDARD
> universal exponent (Dykman, Maloney, Smelyanskiy & Silverstein, PRE 57:5202, 1998: xi = 2 at a
> pitchfork, 3/2 at a saddle-node), and it has been MEASURED at 2.0 +- 0.1 and 2.00 +- 0.03
> (Chan & Stambaugh, PRL 99:060601, 2007). **It also does not discriminate** -- transcritical gives
> 2 as well. And the "every finite window reads below 2" account is the named, standard
> effective-exponent / corrections-to-scaling phenomenon (Wegner 1972; Riedel & Wegner 1974;
> Kouvel & Fisher 1964), documented for escape near bifurcations specifically by Dykman & Ryvkine,
> PRL 92:080602 (2004). **Worse, it is algebra inside §84's own formula**: with eps = 1/2 - gamma,
> A = (32/9) eps^2 (1 - (5/3) eps + ...), so nu_eff = 2 - (5/3) eps, which reproduces both window
> fits on inversion. §84's P3 presented as a finding the leading correction term of its own closed
> form. What survives is the closed form for reversible AM (not found elsewhere -- which per §70 is
> not novelty) and the absolute validation.

**So §63.2's "2 is excluded over [0.20, 0.45] with no drift toward it" was reading the window
bias.** §64 had already withdrawn the exclusion; §84 says what it was. And §64's present
statement -- "nu ~ 2 +- 0.1, not determined more precisely than that by any instrument here" --
is superseded: **nu = 2 exactly, by derivation, validated in absolute terms.** §63.1's width
exponent of exactly 1/2 follows, as §64 said it would if nu climbed to 2.

**T15-n's own two broken criteria (§84.1), in a section about broken criteria.** P3 first demanded
the closed form land within 0.05 of §64's width route and printed FAILS off 1.8300 -- but those
routes fit a width, a stationary distribution and an extrapolated action over their own grids, and
**a threshold between effective exponents of different observables fitted by different protocols
cannot be satisfied only by the thing it claims to test** (rule 19). P4 then gated the
measured-vs-predicted exponent gap at 0.10 and passed at 0.091 -- a fixed tolerance on a
converging quantity, rule 20 verbatim, two sections after rule 20 was last invoked.

~~**T15-n.1, open: why is the closed form's window bias LARGER than the true one?**~~ **-> §85.
CLOSED. Three suspects, two killed, and the survivor confirmed by a test that isolates it.**

  * **S1, truncation of the slow manifold: REFUTED, and backwards.** Solving ds/dt = 0 exactly at
    every u makes agreement WORSE -- residual (gc-g)^0.750 against the truncated form's
    (gc-g)^1.121. §84's closer numbers are a CANCELLATION: A_trunc/A_exact is non-monotone
    (1.0169 .. 1.0401 .. 1.0246), and a systematically better approximation would be one-signed.
  * **S2, the Fokker-Planck vs jump-process action: REFUTED, and it was the motivating suspect.**
    It has the right SIGN -- ln r - 2(r-1)/(r+1) ~ (r-1)^3/12 > 0, so the diffusion form must
    underestimate -- but on the same manifold the two agree to 0.30% against a 2-10% deficit.
    lam_u/mu_u stays near 1 across this barrier.
  * **S3, the adiabatic elimination: CONFIRMED.** The pair X+Y <-> 2B has Delta u = 0, so scaling
    it by M sharpens the timescale separation, leaves the cycle affinity exactly unchanged, and
    makes the elimination asymptotically exact. The residual falls monotonically to zero:
    0.0977 -> 0.0130 at gamma = 0.40 and 0.0651 -> 0.0100 at gamma = 0.44, over M = 1..16.

**And the sign forbids the obvious repair.** The true action is LARGER than the reduced one, while
fast-variable noise feeding the lead would ADD to its diffusion and so LOWER the barrier. Whatever
the elimination loses, it is not captured by widening D_u -- so "add the fast variable's noise" is
refuted before it is tried.

~~**T15-n.2, open: what does the elimination lose?**~~ **-> §86. ANSWERED: it integrates the
right formula along the WRONG CURVE.**

In this regime the barrier is only 6-20 nats, so the exact stationary distribution does not
underflow and the 2-D quasipotential Phi(u,b) = -ln pi(u,b)/Omega is readable outright. The 1-D
reductions use the DETERMINISTIC slow manifold b_det(u); the curve realising the marginal
quasipotential is the RIDGE b_ridge(u) = argmax_b pi(u,b), and they differ by a resolved amount.

  * **Resolved, not quantisation.** At gamma = 0.40, u = 0.38 the displacement is 0.00512,
    0.00474, 0.00447, 0.00436 over Omega = 200..600 -- constant in PHYSICAL units while growing
    from 1.02 to 2.62 LATTICE units. An artifact would be a fixed number of lattice units.
  * **Integrating along the ridge closes the deficit**: 0.9027 -> 1.0088 (gamma = 0.40) and
    0.9348 -> 1.0107 (gamma = 0.44). So a 1-D chain of this form CAN reproduce the action given
    the right b(u) -- the failure is the curve, not the projection. This is a DIAGNOSIS: the
    ridge is read from the exact answer and predicts nothing in advance.
  * **The displacement tracks the deficit toward gamma_c** (0.00254 -> 0.00104 against
    0.0977 -> 0.0651), so it is not a bystander.

~~**T15-n.3, open: is there a PREDICTIVE curve?**~~ **-> §87. CLOSED, affirmatively -- and it
corrects §86's mechanism while leaving §86's numbers standing.**

The tilted-generator route works from the rate functions alone: write the 2-D WKB Hamiltonian in
(u, s), send the FAST HAMILTONIAN FLOW to its fixed point at frozen (u, pu) via dH/dps = 0 and
dH/ds = 0, let H = 0 fix pu(u), and integrate A = -int pu du.

  * **Predicted action 1.0257 and 1.0151** against the exact first passage, at gamma = 0.40 and
    0.44, versus the deterministic manifold's 0.9027 and 0.9348. Grid-converged
    (1.0246 -> 1.0257 over n = 30..240). **No master equation, no stationary solve, no lattice,
    no Omega.**
  * The residual overshoots and shrinks toward gamma_c, consistent with finite timescale
    separation. **No mechanism is attached to it** (rule 17), and note that §85.3's convexity
    argument is about a DIFFERENT reduction and does not license a sign claim here.

**§86's MECHANISM IS WITHDRAWN; its measurement is not.** The predicted instanton curve is
displaced DOWNWARD (-0.00061, -0.00024) where §86's measured ridge is displaced UPWARD (+0.00254,
+0.00104): opposite sign, ~4x smaller. **The ridge of the stationary distribution is not the
instanton's curve.** §86's ridge integration does close the deficit of the 1-D functional it was
fed to, and that number stands; the account -- "the reduction integrates along the wrong curve, and
the ridge is the right one" -- was wrong in its second half. Two curves displaced in OPPOSITE
directions both land near the measured action because they are evaluated in two different
functionals. **Fifth mechanism in this project to die to the first test aimed at it.**

~~**T15-n.4, open: what is the remaining 1.5-2.6% overshoot?**~~ **-> §88. It IS finite timescale
separation -- but the 1/M law predicted for it is REFUTED.**

  * **P1, validation on all eight gammas §84 measured** (not §87's two): overshoot runs +5.61% at
    gamma = 0.30 down to +1.02% at 0.46, monotone, shrinking toward gamma_c as critical slowing
    down requires. Had it grown toward gamma_c the finite-separation diagnosis would be dead.
  * **P2, the M axis §86.1 could not use** (the stationary instrument failed at M > 1; §87 needs
    no stationary solve): the overshoot falls 2.56 -> 0.79% and 1.52 -> 0.59% over M = 1..16. So
    the reduction is exact in the timescale-separation limit.
  * **But |overshoot| x M grows by 5.0x and 6.2x** over that range, and the series is NON-MONOTONE
    at the first step (2.56 -> 2.99 at gamma = 0.40). **The decay is slower than 1/M and its law
    is not identified.** Reported rather than fitted (rule 15).

~~**T15-n.5, open: what IS the decay law of the reduction error?**~~ **-> §89. CLOSED, and no law
was fitted: the error IS the term the reduction drops.**

The slow-fast reduction discards the fast pair's action contribution C = int p_s ds
(Assaf-Roberts-Luthey-Schulten, PRL 106:248102, 2011, name it; §85.3 recorded that the literature
agrees on its magnitude and claims no sign). §87's solve already returns p_s(u) and s(u), so C is
integrable along the curve and comparable in ABSOLUTE terms (rule 16).

  * **C accounts for 91-98% of the overshoot at every gamma**, with the right sign, improving the
    agreement in 8/8 cells: a +1.02..+5.61% error becomes **+0.07..+0.18%**.
  * **It reproduces the non-1/M shape §88 could not fit**, including the NON-MONOTONE first step
    where the raw overshoot rises: corrected ratio 1.0006-1.0012 at every M = 1..16.

**THE ARC'S RESULT, stated plainly.** The escape action of a two-species restoring element is
computable to ~0.1% from its rate functions alone -- solve the fast Hamiltonian fixed point,
integrate -int p_u du, subtract int p_s ds -- with no master equation, no stationary solve, no
lattice and no Omega, validated against exact CME first passage over eight gammas and ten (gamma,
M) cells. §84's closed form was the leading term of this; §85-§88 were the diagnosis.

~~**T15-n.6, open: the residual is FLAT, and that is a suspect**~~ **-> §90. RESOLVES AS
UNRESOLVABLE HERE, and §89's quoted precision is withdrawn.**

The suspect was right about a mechanism and wrong about its size. §84's two-point local slope IS
biased -- the slopes are linear in 1/Omega (correlation -0.99 at all three gammas), and Richardson
extrapolation puts the bias at **+0.234%, +0.297%, +0.415%, LOW**. §89 needed +0.07%. So correcting
it OVERSHOOTS and flips the residual negative (1.0007 -> 0.9984, 1.0006 -> 0.9977, 1.0007 ->
0.9965).

**Then rule 15.** Varying the extrapolation ansatz -- last slope, linear in 1/Omega over all
points, over the last three, over the last four, quadratic -- gives a spread of up to **0.407%** in
the measured A, and §89's corrected ratio against each **STRADDLES 1** at every gamma
(1.0006, 0.9965, 0.9980, 0.9996, 0.9975 at gamma = 0.44). **The measured A is not determined to
better than the residual it was being used to judge.**

  * ~~§89: "the corrected ratio is 1.0007 +- 0.0004 across every cell"~~ -- **the +- is withdrawn.**
  * **§89.1's flatness argument is weaker than it looked, and this is the transferable part.**
    Every cell used the SAME two-point estimator, so its bias is **COMMON MODE**, and a flat
    residual is precisely what that produces. **Flatness across cells sharing an estimator is not
    evidence about physics.** §89.1 read it as evidence pointing at the instrument; it was not
    evidence either way.
  * What survives is better stated than what was claimed: §89's correction brings theory and
    measurement into agreement **to within the measurement's own resolution (~0.4%)**. §88's
    +1.02..+5.61% overshoot is far larger than that and is untouched.

**Correction to §84/§88's published numbers (rule 7):** their measured A carries a two-point bias
of +0.23..+0.42%, low. The originals stand as printed; no re-fit is offered, because the ansatz
spread is as large as the bias.

**T15-n.7, open: can the measured A be pinned below 0.1% at all?** Not with this instrument. **How
to kill:** the bias is O(1/Omega) with a coefficient set by the WKB prefactor b, and T14-e already
localises b (start-side factor, exactly theta-independent, a function of (gamma, epsilon)). If b
can be computed rather than fitted, the two-point bias becomes a KNOWN subtraction rather than an
extrapolation, and the ansatz spread collapses. That makes T14-e -- open since §61 -- the blocker
for a 0.1% test, which is a use nobody had for it.

**§88's PAYOFF, and it is the reason the route was worth building.** §81/§82 could not measure
AM's action at deep barriers -- stationary underflow (§81.1) and negative first-passage times with
ln T saturating at ~35 (§82's scope note) -- so §82's AM arm was confined to the shallow regime.
§87's route never builds a lattice, and reaches gamma = 0.05 (A = 0.8719), 0.10 (0.6151) and
0.20 (0.3099). **These are EXTRAPOLATIONS and the bias runs the WRONG WAY**: it grows away from
gamma_c, so they carry more than the validated +1.02%..+5.61% band, not less. Direction known,
size not.

**And §80's conclusion extends to barriers it was never measured at.** Against §77's measured eta
(9.8813 at gamma = 0.05, 1.8346 at 0.20, 0.6813 at 0.30), A/eta = 0.0882, 0.1689, 0.1944 -- escape
binds throughout. **This does not depend on the extrapolation being accurate**: A would have to be
wrong by 415% to flip it, and P1's bias is a few percent and of the wrong sign to help.

**§86.1: the M-axis cross-check is UNRESOLVED, and was not gated until it agreed.** Reading a
ridge off the stationary distribution at M > 1 produced three distinct contaminations -- a whole
slice at the 1e-300 floor (which drove the mean displacement to -0.078 and flipped the verdict);
partial underflow leaving a spurious peak at the edge of the survivors; and np.interp
extrapolating flat below the first traced slice, **rule 19's own named trap from §59, hit in the
same session the rule was quoted**. After all three were fixed, M = 1, 2 and 8 give 1.0147, 1.0095
and 0.9977 -- P3 reproduced on a second axis -- while **M = 4 gives 0.7082 and is unexplained**;
its coverage is no worse than M = 8's, and a gate written to exclude it excluded M = 8 instead.
**Writing a third gate after seeing which cell disagreed would have been fishing.** P3 and P4 do
not depend on this axis.



**T15-j, open (replacing T15-i): which fixed point decides when the symmetric steady state is
not unique?** 7.4% of drawn networks had more than one symmetric steady state and are outside
§62's claim. **How to kill:** for those networks, evaluate the criterion at every symmetric
fixed point and ask which one's sign matches the ODE from a small kick at each. The natural
guess is that the one that decides is the one reached from the relevant basin, which would
make the criterion basin-dependent rather than network-dependent -- a strictly weaker
statement than §62's, and worth knowing before §62 is cited beyond its stated scope.

**T-OPT -> §57: an exhaustive search REDISCOVERS AM, and the optimum is what §54
predicted.** The founding question, with §56 making it tractable. On {X,Y,B} with X+Y+B
conserved and bimolecular reactions there are 30 conservative reactions in 16
exchange-symmetric classes; networks from 1-3 classes number 696, and AM is one of them.
Every network parameterised as AM is: forward 1, reverse gamma. Figure of merit is §40's
Q = (Var(T)/<T>^2)*<Sigma>/2, floor Q >= 1.

**GATE: AM reads Q = 5.4750 at gamma = 0.05 against §40's published 5.39** -- 1.6% apart
through separate code.

**And the gate caught a screening bug that had discarded the answer.** `delta_star_of`
brackets sign changes between fixed grid points, but `slaved` returns None as delta -> 1, so
a large delta* falls in a dead zone: at gamma = 0.03 and 0.05 the closed form gives 0.971 and
0.952 and it returns None. AM's Q minimum sits at gamma ~ 0.05, so **the first pass threw
away exactly the best cells** and put AM's optimum at gamma = 0.10, Q = 5.79. Caught only
because §40's published number disagreed. Fixed locally, gated to 2.1e-13; `delta_star_of`
untouched since §36 and §39.2 rest on it. Scope limit stated: below gamma ~ 0.05 `slaved`
itself dies, so the grid starts there -- which is also §40's grid.

**THE ANSWER: nothing beats AM by more than 9%, and every one of the top ten CONTAINS AM.**
Best Q = 5.0045 against AM's 5.4750. The search finds no alternative motif, only AM with a
decoration. 0 of 39 cells read Q < 1, so §40's instrument warning never fired.

**The winner is what §54 predicted.** It is AM plus `2B -> X+Y` as a separate forward class,
which makes the disagreement channel **detailed-balanced** (both directions rate 1) while
recruitment stays driven at gamma. §54: `X+Y -> 2B` is SELF-MIRROR, so d = 0 and it
contributes identically zero to P. **A channel that carries no signal but carries a drive is
spending entropy for nothing**; removing the drive buys the 9% without touching the drift.
This converges with §44 from a different direction -- §44 found speeding that same channel a
free lever worth 43-50%, for the same structural reason, and an exhaustive search told
nothing about §44 arrives at the same channel.

**Three real limits.** 39 valid cells from 696 networks -- most capable networks never form a
landscape on this slice, which §56 EXPLAINS: capability is combinatorial, but realisation
needs rates inside the complement of the non-restoring cone, and a single gamma traces a
1-D curve through a high-dimensional rate space. The family is 1-3 classes, bimolecular,
3 species. Omega = 200, where Q is converged to ~1.6% on AM but only ~25% at gamma = 0.30.

**T-OPT-a, open: does the FULL rate space change the answer?** §56 characterises the
restoring rates as the complement of a convex cone, so per-network rate optimisation is
well-posed and this search did not do it -- it sampled one curve, which is the same class of
mistake as the screening bug above. **How to kill:** optimise Q over independent per-class
rates for the top networks. If AM's margin survives, the claim strengthens from "on a slice"
to "in the family"; if something overtakes it, the slice was hiding the answer.

**T-OPT-a -> §58: AM IS ON THE SPEED-OPTIMALITY FRONTIER, and §40's ruler had a factor of
1e57 in it.**

**(1) The one-sided setting is inapplicable, quantified.** §40 flagged in advance that its
absorbing set is TWO-SIDED where the TUR's standard statement is one-sided, and nothing had
measured the cost. Measured: one-sided first passage requires waiting out excursions into the
WRONG basin, which need a barrier crossing back, so it diverges exponentially --
**ln(T_one/T_two) = 0.630*Omega + 5.06**, a factor of **1e57 at Omega = 200**. The two-sided
set is not a sloppy substitute; the one-sided quantity is exponentially large and a committed
bistable element never waits for it. §40's Q = 5.39 is measured against a bound whose
derivation assumes a setting this system exponentially cannot occupy.

**(2) My own optimiser failed, and a direct scan caught it.** Nelder-Mead over free log-rates
reported AM's optimum at rho = 0.997, Q = 4.84 -- i.e. that §44's rho lever does not help Q.
A direct scan: Q falls 4.856 -> 2.919 over rho = 1 -> 20. **A 40% improvement along a single
axis, missed.** Third instrument failure in two sections (§55's counter, §57's screen, this),
and like both, caught only by an independent check.

**(3) THE FRONTIER, which is the answer.** Grid search over (rho, gamma_dis, gamma_rec), 274
cells, drives Q to **1.25-1.33** against the floor of 1 and AM's 5.475. But the optimum rides
the bifurcation (gamma_rec = 0.49, gamma_c = 0.50) with mean decision time **783 against AM's
4.09**. Q is dimensionless and therefore TIME-BLIND. As (time, Q):

    time   4.04  5.91  9.88  18.8  45.5  113   241   783
    Q      5.40  3.38  2.43  2.04  1.89  1.84  1.37  1.25

**AM sits ON this frontier at its fast end** -- AM is (4.09, 5.475), the frontier's fastest
point is (4.04, 5.400). **AM is Pareto-optimal at its own operating speed.** Buying 4.3x
closer to the bound costs **194x in time**.

**This reframes §40 and corrects §57.** §40's "5x from the bound" is not slack -- it is the
price of deciding in 4 time units instead of 800. §57's "nothing beats AM by 9%" was true
only on the slice; with free rates plenty beats AM ON Q, all of it by being slower, and none
of it at AM's speed. §57 is qualified in place, and its conclusion survives in a STRONGER
form: Pareto-optimal, not merely near-optimal.

**And the methodological lesson repeats §38's one level up: Q is not a design objective,
because a system can improve it arbitrarily by slowing down.** §38 found "cost per nat of
reliability" is not a quantity since reliability is bought with free input margin. **A figure
of merit improvable by doing nothing faster measures a trade-off, not a quality.** The
frontier is the quantity; Q alone is a coordinate on it.

**T-OPT-b, open: is the frontier's shape universal?** The points are close to
Q - 1 ~ t^(-0.5). **How to kill:** fit per network family and check whether the exponent
transfers. §39.2 and §46 both found coefficients that did not transfer between axes, so the
prior is that it will not, and a shared exponent would be the surprise.

**T-OPT-b -> §59: NO OTHER TOPOLOGY BEATS AM. The winners are AM re-parameterised.**
§54 buys back the search dimensions -- classes with d_r = 0 contribute identically zero to
the drift, so the search collapses to (rho_ns, gamma_ns, gamma_s) whatever m is. 14 of the 16
conservative classes are signal-carrying.

**My P2 test had an extrapolation flaw that inflated the result:** it scored networks against
AM by np.interp, which returns the ENDPOINT outside the traced range, so anything slower than
AM's slowest point automatically "beat" it. Restricted to AM's traced range, three of five
candidates fall away.

**The two survivors are AM.** cls1 = {2X->B+X, 2Y->B+Y} is AM's recruitment pair with the
forward direction RELABELLED -- its reverse is B+X->2X, the recruitment itself -- and
{cls1, dis} generates a reaction set IDENTICAL to {dis, rec} (checked, pinned by a test). So
they are AM with the recruitment rates decoupled from the shared gamma, which is rate freedom
§58 did not have. **Every candidate with genuinely different chemistry -- AM+revdis,
AM+cls0, AM+cls2 -- is worse at every overlapping time, by 2x to 90x.**

**§58's headline strengthens: from "AM is on its family's frontier" to "no other topology in
the enumeration reaches it."** And §58's far end is corrected -- pooling every rate
assignment of AM's reaction set reaches **Q = 1.115 at t = 3747**, not 1.253 at 783. Buying
5.5x closer to the bound costs 950x in time.

**THE FRONTIER EXPONENT DOES NOT TRANSFER -- A THIRD FAILURE.** Fitting Q - 1 = a*t^(-b):
AM 0.584, AM+revdis 0.078, AM+cls0 0.300, {cls1,dis} 0.402, AM+cls2 0.493 -- **0.078..0.584,
136% spread.** There is no universal time-cost law for approaching the bound; the approach
rate is a property of the network, not of the bound. **P4 predicted this in advance so the
confirming outcome would not be the flattering one.** §39.2's coefficient did not transfer
between axes, §46's SCALING did not, and now §59's exponent does not. **Three attempts,
three failures -- the prior is now explicit: in this system, exponents fitted on one axis do
not carry to another, and any future claim of one needs a cross-axis test BEFORE it is
written down.**

**T-TUR-d RESCOPED BY §66: the closure below is established for EXCHANGE-SYMMETRIC elements
only.** §60's mechanism -- the two absorbing boundaries are exchange images with equal
stationary weight, so the system term cancels outcome by outcome -- is a fact about beta = 0,
and every real restoring device is tilted. Repeated on `am_asymmetric`, |r_e - 1| rises 40x
from beta = 0 to beta = 0.40; whether that is the tilt or the solve is NOT settled, because the
instrument's precision fails on the same variable (the error boundary carries e^-38.9 of the
correct one's stationary weight, below double precision, so Phi_all itself returns 0.81). What
does not depend on settling it: **§60 claimed a general closure from the one case whose
symmetry guarantees it.** For tilted elements the factorisation is untested, not established,
and the founding question's sharpest form is OPEN where real devices live. See §66.

**T-TUR-f -> §69: CLOSED. The factorisation is DEEPER than symmetry, and §60's mechanism
is withdrawn while its result is strengthened.** Asked on §67's Schloegl element -- one dynamic
species, so no exchange symmetry exists to appeal to, and 1-D so the tilted generator collapses
to the chain with lambda and mu swapped and shifted, tridiagonal and free of §66's enormous
factors. **Phi_o = p_o to 2e-13 across skews from 1:1 to 4:1, with the two boundaries'
stationary weights differing by up to e^32.** §60 attributed the cancellation to the boundaries
being exchange images of equal weight; that account is unavailable here and the factorisation
is exact anyway, so **it is wrong, and what does produce the identity has no candidate
mechanism (rule 17: it gets none until something independent supplies one).**

**This also resolves §66 against its own prediction.** §66 saw |r-1| rise 40x with the tilt and
could not separate it from the solve because the nuisance grew with the cause. §69 breaks that
confound by construction -- ln w spans 0 to 32.3 with |r-1| at 1e-13 -- so **§66's rise was the
conditioning**, and its within-band correlation of +0.579 was reading solve error. §66's
refusal to call it either way was correct.

**Standing statement on the founding question's sharpest form:** the fluctuation theorem
factorises over outcomes on a symmetric two-species element at every gamma and on an asymmetric
one-species element at every skew tested, so the error rate is not obtainable from the entropy
production by this route, and the obstruction is NOT a symmetry accident. Combined with §68
removing the affinity floor's universality, the surviving cross-substrate statement is that
reliability is bought with MOLECULES, not free energy. Scope: single elements, two outcomes,
one separatrix; cascades and n > 2 outcomes are untouched.

**~~T-TUR-f, open: does the outcome-wise factorisation survive a tilt?~~ (original entry
below, kept per rule 3)**

**T-TUR-f, open: does the outcome-wise factorisation survive a tilt?** §66 could not reach
past beta/beta_c ~ 0.6 because the pi(n)/pi(n0) boundary convention exhausts double precision
exactly where the predicted effect is largest. **How to kill:** redo the tilted generator in
LOGS -- the boundary weight is a log-sum-exp, as §61's 1-D splitting probability already is --
so the e^-38.9 ratio costs nothing. That is a rewrite of the boundary assembly in
`outcome_split*.py`, not new physics, and it would open the whole beta < beta_c range. Until
then no claim is made for tilted elements in either direction.

**T-TUR-d -> §60: THE FLUCTUATION THEOREM FACTORISES OVER OUTCOMES, so it cannot bound the
error.** §41 verified <e^(-S_tot)> = 1 at absorption to 5.5e-14 and never split it by
outcome. Splitting it -- Phi_o = p_o * <e^(-S)|o>, computed with §41's tilted generator under
§35's outcome-selective boundary -- was the founding question's sharpest form: if error paths
carried exponentially little entropy, the identity would pin p_e to the dissipation.

**The hypothesis Phi_e = p_c is REFUTED by three orders of magnitude** (Phi_e/p_c spans
3.75e-4..0.205, drifting with every axis, so P4's constant-ratio fallback fails too).

**What is true instead: <e^(-S_tot)|o> = 1 for EACH outcome separately** -- 2.8e-3 on the
correct branch, 5.0e-2 on the error branch, with RANDOM SIGN, the signature of numerical
error. So **Phi_o = p_o identically and the aggregate identity is p_c + p_e = 1 in disguise.
The IFT carries no information about which outcome occurred and cannot bound the error
probability.** Reliability is not dissipation in this exact sense -- dead, not unconfirmed.

**The trap worth naming: the aggregate's exactness does NOT validate the split.**
Phi_c + Phi_e = 1 with p_c + p_e = 1 forces p_c(dev_c) + p_e(dev_e) = 0 BY CONSTRUCTION, so
§41's 5.5e-14 says nothing about the outcome-wise identity, which is established here only to
~5% -- eleven orders worse.

**Gate note:** the first grid failed with Phi_c ~ 1e279 at gamma = 0.10, Omega = 90 -- the
conditioning failure §41 itself documented at exactly gamma = 0.10. Restricted to gamma >=
0.20, Omega <= 60 the gate holds at 5.7e-9.

**AND A PATTERN THAT NOW NEEDS ITS OWN RULE.** P6 declared "suspect a sign error in
sigma_local" on a premise that was simply wrong -- error paths do NOT carry strongly negative
total entropy production, because the two absorbing boundaries are exchange images and carry
equal stationary weight, so the system term cancels. **The criterion was wrong, not the
instrument.** That is the FOURTH badly-designed verdict rule this session: §53 demanded a
MAJORITY of topologies flip before declaring non-combinatoriality (one counterexample
suffices); §55 counted p=q pairs whose recorded d is arbitrary; §59 used np.interp, which
extrapolates FLAT beyond the traced range so anything slower than AM automatically "beat" it;
and this. In every case the measurement was sound and the SUMMARY RULE was not.
**Pre-registering a prediction does not protect against pre-registering the wrong test of
it**, and the rules in CLAUDE.md are all about the measurement -- none of them covers the
verdict criterion.

**T15-e, open: is sign(P) > 0 characterisable, or is it just measured?** Divisibility is
settled; restoration is divisibility plus a positive P over the relevant region, and
nothing yet says which symmetric networks have one. **How to kill:** take the cofactor
and cubic families of §42/§43, compute P's zero set, and check whether it coincides with
the separatrix each network's own ODE produces. If P's sign region is readable off the
stoichiometry the way divisibility is, restoration becomes a structural property end to
end; if it depends on the rate constants in a way no combinatorial statement captures,
then the theorem covers the *cannot-reverse* half only and the *does-amplify* half stays
a per-network measurement.

**T-COST-e -> §44: temperature has TWO jobs, and the second one is a free lever.**
`cooling.py` maps the drive to a temperature but puts T only in the reverse rates, and
says so. §44 gave every rate an Arrhenius form. A *uniform* activation energy turns out
to be provably inert -- Q -> lambda Q and sigma -> lambda sigma together, so Sigma is
exactly invariant (measured: 4.6e-15 over a 1000x rescale) -- so the upgrade has content
only through the RATIO of the two forward channels, rho = k_dis/k_rec. **Every result in
FINDINGS before §44 is at rho = 1, and nothing had ever checked it.**

Under the delta* >= 0.40 admissibility floor, T* splits by sign(dEa): -8.3%/-8.5% at
dEa = +0.6 against +0.7%/+6.1% at -0.6. **Temperature does not act through gamma alone**,
so cooling.py's minimal model is not sufficient to support an optimal-temperature claim.
The pre-registered form of that test FAILED and its verdict is kept in §44.3; it failed
because its argmin sat at delta* = 0.227 with the landscape 4% from death, which is §9.2's
withdrawal pathology in a new guise.

**The bigger finding is the one nobody was looking for.** rho is a lever that dominates
gamma and costs nothing: from AM's rho = 1 to the asymptote, cost falls 43-50%, reliability
DOUBLES (L: 17.2 -> 34.4 nats at gamma = 0.16), and mean time halves -- with delta* frozen
to 0.16%, so it is not geometry. **gamma trades cost against reliability, which is what
makes §38's optimum an optimum; rho does not trade.** And physically rho is what a
CATALYST sets: a catalyst lowers an activation barrier without touching dG, so it is
thermodynamically free, and it accelerates both directions equally, which is why gamma and
the affinity 3ln(1/gamma) are untouched.

**§44.1 is a refutation worth keeping.** The prediction was rho* < 1, on the reasoning that
the disagreement reaction moves delta by exactly zero (§30's first cancellation) and so
produces entropy without signal. The sign was backwards: a_f/a_r for that pair is
independent of rho, so driving it fast pushes it onto its own local equilibrium where
ln(a_f/a_r) -> 0 and it carries unbounded flux at bounded dissipation. **Flux is not
dissipation.** The replacement account predicted G must asymptote, and it does.

~~**T-COST-f, open: does rho act through the timescale separation?**~~ **PARTLY CLOSED by
§45 -- it survives in its weak form and is REFUTED in its quantitative form.** The original
statement is kept below.

**T-COST-f, open: does rho act through the timescale separation?** rho raises reliability
while delta* is frozen, which points at slaving -- a faster disagreement channel is a
faster pool, which is §36's on-manifold condition and §39.2's 1/sep law. **How to kill:**
compute sep(gamma, rho) and check whether §39.2's law predicts the L improvement
QUANTITATIVELY, not just in sign. Absolute test, not a fit (rule 16). If it does not, the
mechanism is wrong and only §44.2's measurement stands.

**T-COST-f -> §45: sep matters, sep does not govern, and §39's 6% closes anyway.** The
design scan supplied the control before any prediction was written: sep is NON-MONOTONE in
rho with a minimum near rho ~ 1.5, while §44.2's cost is monotone, so "tracks sep" and
"tracks rho" make opposite predictions.

  * **The residual peaks at rho = 1.5, exactly where sep bottoms out.** Non-monotone in rho
    while cost is monotone -- so it tracks sep, not rho. A monotone residual would have
    looked like a confirmation and meant nothing.
  * **It closes.** |residual| 0.083 at sep < 8 vs 0.0138 at sep > 20; at rho = 32 §39's
    residual is 0.7%, down from 5.3%. **First thing that has ever moved T-COST-c**, and it
    confirms §39.2's "exact in the slaved limit" by reaching that limit with a new knob.
  * **But the sharp tests fail.** The 8x matched-sep pair disagrees 3x (0.022 vs 0.067 at
    sep 11.7 vs 10.6); the two knobs give different curves in 1/sep (slopes 0.441 vs 1.147,
    opposite intercepts); resid*sep drifts 4.5x. **sep is not the governing variable.**

So "rho works by deepening slaving" is INCOMPLETE rather than established, and §44.2's
lever keeps its measurement while its account stays a suspect. §45 also declines to claim
the sign change the raw numbers suggest along gamma: it rests on one cell whose two Omega
disagree by 8%.

~~**T-COST-f2, open: is the governing variable the PATH separation rather than the point
separation?**~~ **CLOSED NEGATIVE by §46 — and §45's framing was corrected in the process.**
All three harmonic path conventions (uniform, time-weighted, cost-weighted) leave the 8x
matched pair matched to within 1.11x against a 3.11x residual gap, and leave the two knobs
split by 2.79-2.89x — marginally WORSE than the point convention. The structural reason was
visible beforehand: over rho = 0.5 -> 4 both eigenvalues scale together (3.4x and 3.7x), so
no ratio, at a point or along a path, can manufacture a 3x.

**§45's framing was wrong in two ways and §46 says so.** T-COST-c was closed by §39.1, so
there were no "remaining candidates" to hand the residual back to; and §39.2 had ALREADY
recorded, under rule 9, that the 1/sep coefficient does not transfer between axes — so
"the two knobs give different slopes" was that published finding on a third axis, a
confirmation read as a refutation. §45 also tested the COST ratio where §39.1/§39.2 test
the TIME ratio.

**T-COST-f2 -> §46: the sharper result is that the 1/sep SCALING is axis-dependent, not
just its coefficient.** On `(T_det/MFPT - 1)*sep`, with §39.2's own resolution cut:

    T axis (§39.2)   0.6465   spread 12.3%   holds
    gamma axis       0.5963   spread 22.2%   holds
    rho axis         0.7465   spread 72.1%   FAILS

and it fails by changing SHAPE, not drifting: gap*sep is U-shaped in rho, bottoming near
0.545 at rho = 2-3 and rising to 0.946 and 1.082 at the ends. **§39.2's headline survives
as a statement about the slaved LIMIT; the RATE of approach is axis-dependent and on rho it
is not 1/sep at all.** The gate reproduced §39.2's published cell (+0.0914 at gamma = 0.20,
sep = 7.00) with +0.0805 / +0.0938 bracketing it.

**T-COST-h -> §47: YES, and the parameter-free form caught a broken instrument.** Singular
perturbation gives the lag with no fitted constant at all:

    eps(delta) = (dmu/ds)(ds*/ddelta)/(dnu/ds),   T_det/MFPT - 1 ~ <eps>_time

Measured against the exact gap: pred/gap = 0.921 over 11 rho cells and 0.931 over 5 gamma
cells -- **a parameter-free prediction of a quantity three sections had fitted, right to 8%
on average.** And <eps>*sep on the gamma axis is 0.627 against the 0.5963 §46 fitted and the
0.6465 §39.2 fitted on the T axis: the constant that would not transfer between axes is now
COMPUTED, and it lands between the two fitted values.

**Then it disagreed with §46 and §46 lost.** The prediction is flat in rho (7% spread) where
§46 measured 62% scatter and declared the 1/sep scaling axis-dependent. §46 computed that at
Omega = 200 without checking Omega-convergence of the gap, against a §39.2 result
established over Omega = 400-800 with convergence checked. Re-run: spread 46.7% -> 44.4% ->
22.9% over Omega = 300/500/700, mean settling at 0.686 against the predicted 0.654.
**§46.1's headline is withdrawn.**

**The withdrawal is not "the scaling holds" (rule 14).** Several cells still move at
Omega = 700 -- rho = 6 shifts 26% from Omega = 500 -- and rho = 32's gap has fallen to 0.0095,
near §39.2's 0.008 floor. The rho axis is UNRESOLVED. §46's gamma-axis 22.2% is under the
same suspicion, same Omega.

**Why the parameter-free form mattered.** §46 applied rule 13 to the model (three averaging
conventions) and never to the measurement. A FITTED model would have absorbed the scatter
into its constant and reported agreement; only a prediction with no free parameter could
say "your data is wrong". That is rule 16's actual content.

~~**T-COST-i, open: what is the 8% systematic shortfall?**~~ **CLOSED by §48: it is
finite-Omega absorption bias and it extrapolates away. The lag model is EXACT as
Omega -> infinity.** Five cells, endpoints matched to the lattice, all converge
monotonically upward in Omega and extrapolate to pred/gap = 1.00 (0.984..1.004, mean 0.996).
That identifies §39.1's candidate (iv) -- absorption selecting the leading edge of the
packet -- which §39.2 left live for the time and which nothing had measured.

**Two things in §48 are worth more than the closure.**

  * **The endpoint fix bought MONOTONICITY, not accuracy.** Matching `T_det`'s limits to the
    lattice endpoints the CME actually uses (rule 11) moves T_det by under 1% and does NOT
    reduce the shortfall -- the pre-registered per-cell prediction failed outright. But the
    nominal series BOUNCE with Omega (0.936, 1.052, 1.000, 0.991) while the matched ones are
    monotone in all five cells. Without that, no convergence could be read off at all, and
    P3/P4 were untestable.
  * **Rule 18 caught in the act.** Averaged over cells, the shortfall looks like a clean law:
    (1 - mean ratio)*sqrt(Omega) = 3.842, 3.876, 3.866, 3.708 -- constant to 4% over 3.3x in
    Omega, and it would have been reported as a discovered 1/sqrt(Omega) absorption law. Per
    cell the same coefficient spans 0.7 to 6.5. **The mean is an averaging artifact and the
    law is not claimed.**

**Rule 15 applies to the closure itself: the DECAY LAW is unresolved.** Free-exponent fits
give intercepts 0.984..1.004; a fixed 1/sqrt(Omega) gives 0.843..1.098 and sends gamma = 0.35
to 0.843, not 1. Exponents span 0.20..0.82 from four points with two parameters. **The
intercept is robust, the exponent is not**, and convergence to exactly 1 is established only
under a free exponent.

**T-COST-j -> §49: the exponent is DERIVED as 1, not 1/2 and not fitted.** The Laplace
expansion of the first-passage integral, with D0 = (up+dn)/2 read straight off `updown`,
gives `T_det/MFPT - 1 = <eps>_time + K/Omega`. Confirmed where testable: (gap - pred)*Omega
is constant to 6% (rho=32) and 17% (gamma=0.20). **§48's fitted 0.20..0.82 was averaging
asymptotic cells together with pre-asymptotic ones** -- the pre-asymptotic ones being mostly,
though not entirely, the longest-traversal cells (that ordering FAILED at gamma = 0.07 and
§49 says so).

**T-COST-k -> §50: a COMPLETE 1-D account supplies only ~21%. The rest is 2-D.** The 1-D
slaved birth-death chain contains bulk, boundary layer and jump discreteness at once and is
one tridiagonal solve. In the two cells converged by Omega = 2000 it gives bd/cme = 0.213
and 0.217, **agreeing to 2%**. So ~79% of the absorption coefficient is not one-dimensional
at all.

  * **§49's K is withdrawn by §50.** At the same cell the bulk term alone gives 0.293 while
    the complete account gives 0.217 -- a whole smaller than one of its parts. K is an
    asymptotic expansion resting on a near-cancelling integral (g' changes sign mid-path,
    and the cancellation flips K NEGATIVE at gamma = 0.07); bd_coeff is exact. §49's
    surviving contribution is the exponent, which does not depend on K.
  * **Rule 11 fired twice in three sections, in the same experiment family.** §50's first
    pass integrated T_det over unrounded limits while the chain ran between rounded lattice
    points -- an O(1/Omega) mismatch in delta which, times Omega, is exactly the size of the
    effect. It gave sign-flipping coefficients spreading 3029%. Matching the endpoints
    removed every oscillation.
  * **§50's cell selection is POST-HOC and it says so** -- the convergence criterion was
    chosen after seeing which cells converged, unlike §47's pre-registered delta* >= 0.40.
    Provisional pending larger Omega.

~~**T-COST-l, open: is the missing ~79% the pool-fluctuation Jensen term?**~~ **REFUTED
EXACTLY by §51, and by §30's identity.** mu = k*delta*(1 - (1+gamma)s), verified against the
network's own fluxes to 9.8e-16 over 81 states x 3 gamma x 3 rho -- and that bracket IS §30's
identity at n = 2, k(b - gamma*s), in concentrations. **The drift is exactly LINEAR in the
pool coordinate**, so d2mu/ds2 = 0 identically and the Jensen term vanishes: J = 0 to
machine precision in all five cells, (bd+J)/cme unchanged at 0.213 and 0.217. §39.1's
candidate (iii) is now dead for the TIME as well as the cost. The theorem that opened this
session is what closes the question.

  * rho does not appear in mu at all -- the disagreement channel moves delta by exactly
    zero, §30's first cancellation -- which is why the whole rho family shares one drift law.
  * **§51's P1 gate "failed" as a division-by-zero artifact** and is reported as one: it
    measured the RELATIVE convergence of an exact zero. P2's sign test never fired either;
    d2mu/ds2 is neither positive nor negative.
  * The NOISE is curved where the drift is not -- d2(up+dn)/ds2 = k(gamma-2) != 0 -- but a
    Jensen term in the diffusion reaches the MFPT only through the correction itself, so it
    is O(1/Omega^2) and cannot be the missing 79% either.

**§52: the 21% is now solid, the Omega-asymptotics is not.** Before proposing a third
mechanism for the missing 79%, §52 verified the measurement it would be fitted against --
§50's figure rested on two POST-HOC selected cells, one with a top point §49 had excused.
Adding Omega = 1400 and three new cells: **gamma = 0.28, rho = 8 converges independently and
gives bd/cme = 0.202** against rho=32's 0.214, mean 0.208, spread 6.1%, and **seven of eight
cells give 0.194-0.239 regardless of convergence status.** So "~79% is two-dimensional" is on
much firmer ground.

**But only 2 of 8 cells are Omega-converged at Omega = 1400, and they are the two FASTEST
POOLS.** Every slow-pool cell is still climbing. The coefficient is asymptotic in both Omega
AND sep, and at moderate sep, Omega ~ 1000 is not asymptotic. Consequences, all recorded in
place:

  * **§49's "numerical floor" diagnosis of gamma = 0.07 is WITHDRAWN.** The series is 2.608,
    2.494, 2.601, 3.653, 4.870 -- it genuinely rises, it is not a floor.
  * **§49's exponent-1 confirmation holds at rho = 32** (0.1% last step, 6.7% total drift
    over 4.7x in Omega) **and was premature at gamma = 0.20** (14.1% last step).
  * **§48's extrapolations to pred/gap = 1 were fitted on non-asymptotic data** in six of
    eight cells. The Omega -> infinity claim is established at rho = 32 and merely consistent
    elsewhere.

**Running this before T-COST-m was the right order:** a third mechanism fitted against six
non-asymptotic cells would have been fitted against a moving target. T-COST-m stays open and
should be tested at rho = 32 and gamma = 0.28/rho = 8 -- the two cells where the target holds
still -- rather than on the axis-spanning grid.

**T-COST-m, open: is it noise-induced drift from the delta-s cross-correlation?** Adiabatic
elimination produces an effective slow drift with a term set by <delta-fluct * s-fluct>,
distinct from both the deterministic lag and the Jensen curvature. It survives here
precisely because **dmu/ds = -k(1+gamma)*delta != 0 even though d2mu/ds2 = 0** -- the drift
is linear in s, not independent of it. **How to kill:** solve the stationary 2x2 Lyapunov
equation about the manifold, take the delta-s entry, and test (dmu/ds)*<delta s>/mu in
ABSOLUTE terms against the missing 79% (rule 16). The sign is again forced positive.

**Budget after §51:** of the absorption coefficient C, ~21% is exact 1-D discreteness (§50),
**0% is Jensen (§51)**, ~79% unaccounted. Two named candidates eliminated rather than merely
unfavoured.

**T-COST-l, open: is the missing ~79% the pool-fluctuation Jensen term?** The pool b
fluctuates about its slaved value, so delta feels E[mu(delta,b)] rather than mu(delta,b*) --
an O(1/Omega) term no 1-D reduction can carry. **This is §39.1's candidate (iii)**, withdrawn
there as an explanation of the COST while explicitly left live for the TIME: "Either may
still explain the time gap, which is now a separate and cleaner question." Nine sections
later it is the leading explanation of exactly that gap. **How to kill:** compute
E[mu(delta,b)] - mu(delta,b*) with the pool variance from the linear-noise approximation
about the manifold -- one 2x2 Lyapunov solve per delta -- and test the coefficient in
ABSOLUTE terms against the missing 79% (rule 16, not a fit). The SIGN is forced: absorption
shortens the MFPT, so the term must be positive.

**T-COST-j, open: what sets the absorption exponent?** The shortfall decays as Omega^-p with
p measured at 0.20..0.82 across five cells -- not the 0.5 a naive packet-width argument
gives, and not constant. **How to kill:** a packet-width account predicts the coefficient
scales with the local drift steepness at the threshold and the diffusion there, both of which
are computable exactly on the slaved manifold; compute them and test the coefficient in
ABSOLUTE terms (rule 16), not by fitting p. gamma = 0.35, the slowest converger at p = 0.20,
is the discriminating cell.

**T-COST-i, open: what is the 8% systematic shortfall?** pred/gap ~ 0.92 on both axes. It
was predicted in advance as next-order quasi-steady-state, which is O(gap) in the same small
parameter -- **but it does not behave like one.** A higher-order term must shrink where the
gap is small; measured, pred/gap is 0.936 at gamma = 0.07 (gap 0.150) and 0.780 at
gamma = 0.35 (gap 0.051), which is the wrong direction. **How to kill:** gamma = 0.35 is also
the cell nearest the resolution floor, so re-measure shortfall-vs-gap at Omega = 700+ before
treating the sign as real. If it survives, higher-order QSS is eliminated and the shortfall
is something else.

**T-COST-h, open: is the lag set by manifold velocity over relaxation rate, not by the
eigenvalue ratio?** The three axes move the two eigenvalues in OPPOSITE proportions — rho
takes fast x26 and slow x4, gamma takes fast x1.4 and slow /3.3 — so a law in their ratio
alone cannot capture both. The slow-manifold lag should physically go as
`|db*/ddelta * mu| / |lambda_fast|`, which is dimensionally a different object from
`|lambda_slow|/|lambda_fast|` and coincides with it only when the manifold's motion is set
by the slow eigenvalue. **How to kill:** compute it on all three axes (T, gamma, rho) and
check constancy against §39.2's 12.3% bar. A partial check already discourages the simplest
version -- `gap*|lambda_fast|` alone drifts 5x across the rho axis -- so the
`db*/ddelta * mu` factor must carry the whole difference or the idea is wrong.

**T-COST-f2, open: is the governing variable the PATH separation rather than the
point separation?** `sep_of` takes the eigenvalue ratio at the symmetric point (x = y,
b = b*), but the traversal runs over delta in [eps*delta*, theta*delta*], away from it. A
separation measured at one point need not represent the separation along the path, and the
rho and gamma knobs deform the manifold differently -- which would explain BOTH the
matched-pair failure and the two-curve split without abandoning slaving. **How to kill:**
compute a path-averaged separation over the actual traversal and re-run §45's P3 and P4
against it. If the 8x matched pair then agrees and the knobs collapse onto one curve, the
governing variable was the path separation. If they still split, slaving is not the
mechanism at all.

**T-COST-g, open: is selective catalysis the real substrate criterion?** The question that
prompted §44 was whether basing logic on atoms rather than molecules would be cheaper. §44
says the dominant free lever is selective acceleration of the non-signal channel -- which
requires a substrate with REACTION-SPECIFIC catalysts. Enzymes have that; electronic or
nuclear transitions do not, since you cannot catalyse one atomic transition among many the
way an enzyme picks one reaction among thousands. **This is a suspect, not a result**, and
it is not tested by anything in §44. **How to kill:** it predicts that the achievable rho
range, not the energy scale dE, is what separates substrates -- so a substrate comparison
holding dE/kT matched and varying only the achievable rho should reproduce the whole
effect. Until that is run, §44 supports only the claim that rho matters more than gamma,
not any claim about atoms.

**T16, open: does CONCATENATED AM show a threshold, or does the ceiling survive?**
FINDINGS 1's wall and 12.1's depth ceiling say a single AM stage has a fidelity ceiling
-- more molecules buy exponentially less, and past a point nothing. The threshold
theorem of fault-tolerant computing says the opposite about a restoring code: below a
threshold physical error rate, CONCATENATION buys arbitrary fidelity at polylog
overhead. Both are statements about restoration, they disagree, and the reason may
simply be that nothing in this project has ever concatenated -- every result here is
level-0.

**The construction must be chemistry-only, and this is where the experiment can be
faked without noticing (rule 10).** The obvious version -- run k tanks and take the
majority of their answers in numpy -- inserts a free, noiseless, perfectly reliable
majority gate, the exact class of error that has already cost three withdrawn results.
The honest version is a POOL MERGE: run k independent AM tanks of size Omega, then
physically combine their contents into one tank of size k*Omega and run AM there. The
merged tank's initial margin is the sum of the k committed margins -- positive iff a
majority answered correctly -- and the combining stage carries its own noise because it
IS an AM tank. No sign(), no free comparison, nothing the chemistry could not do.

**How to kill.** Measure p_L at concatenation level L and fit the level-to-level map.
  * THRESHOLD (QEC): with k = 3 the merge fails only if two of three units fail, so
    p_{L+1} = C*p_L^2 and p_L falls doubly exponentially below p* = 1/C. **The sharp
    test is the EXPONENT of a log-log fit of p_{L+1} against p_L: 2, not 1.** A trend
    alone proves nothing (§30.2).
  * CEILING: p_{L+1} ~ p_L, or p_L -> p_inf > 0, the merged tank's own error floor
    cancelling whatever the vote bought.

**Two traps already visible.** (i) The merged tank has k*Omega molecules, so comparing
it to a level-0 tank of Omega compares different systems; the control is a SINGLE tank
of size k*Omega at the same margin (rule 18), without which a "threshold" is just the
exp(-Omega V) scaling FINDINGS 1 already published. (ii) The k units must fail
independently for the vote to help. They share nothing by construction, but the merge
correlates them afterwards, so independence must be MEASURED, not assumed.

**Why it is worth doing.** If the ceiling survives concatenation there is a structural
difference between chemical and quantum restoration worth naming. If a threshold
appears, this project has an independent classical instance of the deepest theorem in
fault-tolerant computing, in a model where the thermodynamic cost of every stage is
already measured (9.1's affinity floor, 20's fuel lifetime) -- which the quantum
version does not have.

**T16 -> §32: CLOSED. The ceiling survives concatenation, and the reason is an exponent
count that passes an absolute test.** The pool-merge construction was used, so no free
`sign()` does the voting: k tanks of size Omega commit, their contents are physically
combined into one k*Omega tank, and that tank runs AM itself. Everything exact, no Monte
Carlo. The merged stage is reliable -- p_merge = 0 / 0.0030 / 0.9970 / 1.0000 for j = 0,
1, 2, 3 wrong -- so the vote genuinely works and has no floor of its own, and
p_1 ~ p_0^1.85 follows near-definitionally from that.

**But voting loses to POOLING the same molecules, in 9 of 10 cells, by a factor that
grows exponentially.** The reason is countable in advance: p_0 ~ exp(-Omega*c) gives
3*p_0^2 ~ exp(-2*Omega*c) for the vote against exp(-3*Omega*c) for one k*Omega tank --
**voting squares the error, pooling cubes the exponent** -- so ln(p_1/p_pool) must be
linear in Omega with slope c, the collapse rate measured independently in §27/§28. The
eps-controlled slope is 0.0231 (attractor readout) and 0.0208 (threshold), against
-2*V_exact = 0.0211 from §15's closed forms with no free parameter: **within 1.4% and
9.1%**, and the two conventions bracket the prediction. Rule 16's absolute test, passed
on a quantity that has nothing obviously to do with the closed form it was predicted
from.

**The structural difference from quantum error correction is now nameable, and it is not
that one restores and the other does not.** In QEC the physical error rate is FIXED --
there is no operation that lowers it by using more of the same qubit -- so concatenation
is the only lever. Here error falls exponentially in Omega, so a bigger tank is a lever
and a better one. **Chemistry does not need the code because it has a cheaper knob.**
That answers the question T16 was opened to ask, and it answers it against the QEC
analogy rather than for it.

**T16-a, open: does PERIODIC RE-MERGING beat the single-tank hold?** §32 is one-shot
voting -- each tank runs once to commitment -- and that is not what QEC does. Real fault
tolerance is time-extended: fresh ancillas repeatedly remove errors that accumulate
during storage, and the threshold theorem is a statement about that repetition, not
about a single vote. **§12.1's depth ceiling is the right target**, since it concerns a
bit HELD over time, and the fair comparison is a k*Omega tank holding a bit against
periodic re-merging of k smaller tanks at interval tau. **How to kill:** measure the bit
lifetime under both at matched molecule count and matched total dissipation, sweeping
tau. If re-merging wins at some tau, the exponent count in §32 is a one-shot artifact
and concatenation does buy something over time. If the single tank still wins at every
tau, the conclusion generalises and the cheaper-knob statement is about restoration in
this family rather than about one protocol. **The dissipation match is what makes it
honest** -- re-merging burns drive every cycle, and §9.1's affinity floor plus §20's
fuel lifetime already price that, so this is a comparison the project can make and the
quantum version cannot.

**T16-a -> §33: CLOSED, and it CORRECTS §32.** Re-merging beats the single-tank hold on
lifetime in 10 of 13 cells -- 5.0x longer at (k=3, Omega=8) while burning 29% LESS, and
3.7x longer at (k=5, Omega=14). §32's one-shot conclusion did not generalise to the
time-extended protocol, which is exactly what §32's own scope note flagged as untested.
**§32's closing sentence, "chemistry does not need the code because it has a cheaper
knob", was too strong; it is left standing in FINDINGS with the correction beside it.**

**The model was validated before any lifetime was quoted.** The two-state renewal
reduction is legal only if first passage is near-exponential, and `first_passage_moments`
settles it exactly: std/mean = 0.9890 -> 1.0000 as the barrier deepens, which is Kramers.
L_remerge ~ 1/tau holds to 0.6% at small tau. **P1 FAILED** -- ep_rate/N is not
size-independent (10.8% over 4x in N) -- so dissipation was accounted explicitly rather
than assumed, and the burn-rate ratio spans only 0.71-1.10 against lifetime ratios of 5x
to 20x. Dissipation cannot explain a difference an order of magnitude larger than itself.

**THE WIN IS CONDITIONAL ON CYCLE SPEED and that is the dominant sensitivity.** Since
L_remerge ~ 1/tau exactly, the winning region tracks it: at tau = t_relax re-merge wins
up to Omega = 18 (k=3), at 5 t_relax only to Omega = 8, and at tau >= 10 t_relax it
never wins at any Omega tested. Quoting the win without the cycle time is quoting half
the result.

**P5 as it came out.** Predicted slope ratio k - ceil((k+1)/2) = 1 and 2; measured 1.0755
(MATCH) and 2.3139 (**MISMATCH**, twice the tolerance fixed in advance). Both run high
and in proportion, consistent with the Kramers power-law prefactor the pure-exponential
ansatz omits. Comparing the two measured slopes to each other cancels that shared bias:
2.1514 against a predicted 2.0000, 7.6%. So the integer structure across k is real and
the absolute values are not clean enough to claim it exactly. Recorded as failed-as-
stated with the diagnosis, not restated to fit.

**THE CONTRAST WITH QEC IS SHARPER THAN §32 HAD IT, and this is the durable statement.**
In QEC, below threshold, concatenation's advantage GROWS WITHOUT BOUND with level,
because the physical error rate is fixed. Here re-merging's advantage occupies a BOUNDED
WINDOW and then REVERSES, because the physical error rate itself falls exponentially in
Omega -- growing the tank eventually outruns the code. **Chemistry has a knob QEC lacks,
and the code wins only until that knob is turned far enough.** That is a statement about
why the digital abstraction looks different in the two substrates, and it came out of a
runnable experiment rather than an analogy.

**T16-b, open: where does the crossover sit in the (Omega, tau, k) volume, and is there a
closed form for it?** §33 has the crossover at Omega ~ 18 for k=3 at tau = t_relax and
~14 for k=5, but only over the range the exact MFPT solve reaches -- N <= 72 before
`first_passage_moments` stops being trustworthy, which is why k=7 could not be fitted at
all. The exponent count predicts the crossover satisfies (k - m)*c*Omega = ln(T(Omega)^?
/ tau) + const, so it should be extractable in closed form from c and t_relax alone.
**How to kill:** derive the crossover from the exponent count, then test it in ABSOLUTE
terms against measured crossovers at a gamma NOT used to fit it (rule 16 -- and §30.1 is
the standing reminder that a formula matching one sweep proves nothing). The MFPT ceiling
is the binding constraint, so §26's MLRift sampler is the likely instrument for the
larger tanks, where an exact solve is unavailable but the lifetime is still measurable.

**T16-b -> §34: CLOSED. §33's crossover is a consequence of the exponent count, not an
empirical boundary.** With m = ceil((k+1)/2) and ln T(N) = c*N + a,

    Omega_x = [ (m-1)(a - ln tau) - ln C(k,m) ] / [ c (k-m) ]

and because m-1 = k-m = (k-1)/2 for EVERY odd k, this is (a - ln tau)/c minus a small
combinatorial term -- **the leading term contains no k**. `c` and `a` come from a
straight-line fit to the HOLD protocol's own MFPTs; no crossover measurement enters the
prediction anywhere, which is what makes it an absolute test in §28's and §32's sense
rather than a fit to the thing it explains (§30.1 being the standing reminder).

**P1, the structurally surprising prediction, CONFIRMED:** the crossover is nearly
k-independent -- spread 3.00%/3.61%/4.39% across k = 3,5,7 -- while the win MARGIN at
fixed Omega differs by more than 2x between k=3 and k=5. Margins differ, crossings
coincide. No fit to crossover data would have suggested that.

**P2 at tau = t_relax: within 3.3% across three gamma**, over a 1.9x range in Omega_x
and a 2.6x range in c. **But P2 as stated FAILED over the full grid** (0.9191 +- 0.0976,
range 0.685-1.033, well outside the +-10% predicted), and the failure tracks the
derivation's own tau << T assumption: correlation with tau/T is -0.696. Solving the
crossover condition EXACTLY -- same c and a, still no crossover data -- moves the mean
from 0.9191 to 0.9864 and doubles the within-2% count, so the linearization was the bias
and the exponent count under it is sound. **Scatter falls only 1.14x and residual
disagreement persists at large tau and small Omega_x, where the crossover lands near
Omega ~ 7 and both renewal assumptions are marginal. That residual is NOT explained and
no cause is named for it** (rule 17).

**P3: d(Omega_x)/d(ln tau) = -1/c exactly**, no k, no combinatorics -- measured 0.78 /
0.87 / 0.996 of that across gamma = 0.25/0.30/0.35, the agreement improving as the
reachable range widens, and exact to 0.4% at gamma = 0.35.

**P5's DIRECTION WAS BACKWARDS and the reason is instructive.** I predicted the high-
gamma end would run out of reach first, since shallower barriers give smaller c and
larger Omega_x -- true (15.5 -> 28.2 from gamma 0.25 to 0.35) but the wrong effect. The
MFPT validity ceiling moves FASTER than Omega_x does, because a shallow barrier makes
T(N) grow slowly and keeps exact solves trustworthy to N = 126 at gamma = 0.35 against
N = 72 at 0.30. So gamma = 0.35 gave 10 reachable cells and gamma = 0.15 gave 2. **When
the instrument's reach and the physics move along the same axis, predict the
instrument's reach explicitly** -- it decided which cells existed, and I had it exactly
inverted.

**T16-c -> §34.1: CLOSED, and instrumentally -- it is the Kramers prefactor.** Refitting
the hold data as ln T = c*N + b*ln N + a (still hold-only, no crossover data) and
re-running the identical test moves the mean from 0.9864 to **0.9988** and cuts the
scatter **3.86x**, from sd 0.0854 to 0.0221, with 22 of 31 cells inside 2%. Candidates
(i) and (ii) are not needed and are withdrawn rather than left standing. The fitted b is
gamma-dependent -- -0.653, -0.459, -0.285, -0.069, +0.074 across gamma 0.15..0.35,
crossing zero near gamma ~ 0.32 -- and **no mechanism is proposed for that** (rule 17);
it is the obvious next thing to explain. Killing the instrumental candidate first cost
one refit and settled the question outright.

**Superseded, kept for the record:** **T16-c, open: what is the residual scatter at large tau?** Removing the linearization
fixed the mean but left the spread (0.773-1.196), concentrated where Omega_x ~ 7 and
tau/T > 0.05. Three candidates, none tested and none preferred: (i) the committed state
is not meaningful at ~7 molecules, so T(Omega) itself is the wrong object there;
(ii) tau is no longer long against t_relax, so merged portions do not re-amplify before
the next merge and the two-state reduction loses its justification even though
std/mean still reads ~1; (iii) ln T(N) curvature -- the Kramers prefactor -- which a
straight-line fit absorbs into `a` and which would bias small-Omega extrapolations
most. **How to kill:** (iii) is separable without new machinery by fitting
ln T = c*N + b*ln N + a and re-running the same absolute test; if the spread collapses,
the prefactor is the cause and (i)-(ii) are not needed. Do that BEFORE proposing a
physical mechanism -- the cheap instrumental explanation has to be eliminated first,
which is the §28.1/§28.2 lesson.

**T14 -> §35: THE PROBABILITY FLOOR WAS AN IMPLEMENTATION ARTIFACT. It is gone.**
T14 has said since §21 that "large Omega AND small probability is reachable by neither
instrument". The probability half was never physics: `p_cme` returns `1 - split`, an
error probability computed as a difference of two numbers near 1, so it dies to
cancellation near 1e-12 -- which is what cost §28 its gamma = 0.15 cells. Naming the
WRONG outcome as the favoured set in `splitting_probability` solves for the small number
DIRECTLY, with no subtraction anywhere: **P = 6.35e-33 at Omega = 2000 in 115 s**, 25
orders below anything previously measured and through the founding claim's own regime.

Validated three ways, because a result this convenient is exactly what this project has
been burned by: identical to the subtractive route to 7-8 digits across the whole
overlap; componentwise relative correction of **1.0e-13 at h = 6.35e-33** (exact
per-row summation plus one refinement step -- a NORM residual would not have noticed a
garbage small component); and it is the M-matrix property, so the LU solve carries no
subtractive cancellation and relative accuracy survives to arbitrarily small values.
**The remaining limit is the Omega^2/2 state space and double-precision underflow near
1e-308, not the probability.**

**T14-a/T14-c-iv -> ANSWERED as a side effect, and the answer costs a headline.** Every
collapse slope this project has published is a finite-Omega EFFECTIVE slope. The local
slope drifts monotonically (-0.051603 -> -0.049895 over 29.21 decades at gamma = 0.20,
eps-controlled), so P ~ A(Omega) exp(-c Omega) with an algebraic prefactor, and a
straight-line fit returns a contaminated c. Three ansaetze all reported (rule 15); the
prefactor form beats the pure exponential by **10-45x in rms** for one extra parameter,
which is §5.1's test run in the opposite direction and read that way deliberately. The
prefactor exponent sits at -0.45..-0.40, drifting toward WKB's -1/2 as the lever arm
lengthens, but **fixing it at -1/2 costs 2-6x in rms so it is not exactly that, and no
mechanism is asserted** (rule 17); c shifts under 1.3% either way, so the rate is robust
to the prefactor's form.

**Against §15's closed form the asymptotic disagreement is 7.5-15.5%, not §28.3's
0.4-11%**, and much flatter in gamma: the asymptotic excess falls 0.155 -> 0.075 across
gamma = 0.20-0.35 where §28.3 measured 0.110 -> 0.004. **§28.3's zero crossing near
gamma ~ 0.357 is withdrawn as an asymptotic statement** -- an 18.7x correction at
gamma = 0.35, where its "closest agreement in the project, 0.4%" was almost entirely
finite-Omega contamination and was flattering precisely where the window was shallowest.
§28.3's numbers stand as measured effective slopes and its P2 eps-independence is
untouched, being a fixed-window comparison. **T14-c-iv needed gamma = 0.38-0.44 and no
longer does: the excess does not cross zero because it was never as small as the shallow
window made it look.** gamma = 0.35 is UNRESOLVED by the 3% ansatz-spread criterion
fixed in advance (4.99%, the shortest lever at 8.79 decades) and is reported so.

**T14-d -> §35.3: ILL-POSED, and it can be PROVED rather than asserted.** Read from the
LOCAL SLOPE instead of ln P -- where s(Omega) = -c + b/Omega makes b a two-parameter
line with the constant differentiated away -- the curvature is real (adding 1/Omega^2
lifts R^2 from 0.94 to 0.99) but **b does not converge with model order**: -0.356 ->
-0.643 -> -0.340 as linear -> quadratic -> cubic, with q flipping +44 -> -78. The reason
is structural: over Omega = 150..1950 the reciprocal spans only 13x, and on that range
corr(1/Om, 1/Om^2) = **+0.961**, corr(1/Om^2, 1/Om^3) = **+0.986**, with design
condition numbers 6.7e2 -> 3.6e5 -> 2.8e8.

**The FUNCTION is determined while the DECOMPOSITION is not** -- the three orders
extrapolate to 1/Omega -> 0 within 0.69/0.85/1.46% and to Omega = 4000 within ~1%, while
their coefficients differ by 2x. **The decades of P are irrelevant to this question**;
only the range in 1/Omega matters, and no depth in P widens it. That retroactively
explains the whole arc -- §35's b values, §35.2's failure to separate gamma from lever
arm, §35.2's P4 scatter, the 40-139% sliding-window spreads -- as ONE ill-conditioned
projection reported four ways. §35.2 concluded "the basis was wrong"; truer is that
EVERY basis is wrong for this question over a bounded range.

**§35.1's b values are withdrawn as measurements of anything, and nothing in this
project currently constrains b, including whether it is -1/2.** What that buys is the
opposite of a loss: §35's rate c is now verified across four bases and three model
orders to **0.03-0.62%**, so the 7.5-15.5% asymptotic disagreement and the withdrawal of
§28.3's zero crossing stand more firmly than when published -- the one quantity they
depend on is the one this arc proved robust while everything around it was not.

**T14-e, open: DERIVE b instead of measuring it.** The splitting probability is a ratio
of scale-function integrals S(x) = int exp(2 Omega int mu/D), and Laplace asymptotics on
that RATIO gives the algebraic prefactor in closed form -- the Gaussian widths at the
dominant endpoints partially cancel, which is precisely why the exponent need not be
-1/2 and why assuming it was, was never safe. **How to kill:** derive it, then test in
absolute terms against the FUNCTION s(Omega), which §35.3 shows is pinned to ~1% even
though its coefficients are not. Predict the curve, not the coefficients. Measuring b
directly would need 1/Omega decorrelated over ~100x -- Omega ~ 15,000, ~1e8 states --
out of reach for the exact solver and not worth reaching for when the analysis exists.

**T14-e, NARROWED by §61: the prefactor discrepancy is a START effect, and the exponent
is exact.** §35.3's proof that fitting cannot separate prefactor from exponent is about
the 2-D problem. §50's exact 1-D slaved chain has a CLOSED-FORM splitting probability
(scale-function sums evaluated in logs -- no solve, no fit, no cancellation), and §44's
rho is a lever that drives sep -> infinity, which §39.2 says makes the reduction exact.
Together they separate the two: an exponent error is a SLOPE in ln(P_1D/P_2D) vs Omega,
a prefactor is an INTERCEPT. At rho = 1024 the slope vanishes (**exponent error 0.006%**,
span 0.0138 nats against ln P ~ 70) and a constant intercept survives. That constant is
**not universal** -- it moves with gamma (1.26-1.60) and strongly with eps (1.37->1.83) --
**but it is exactly theta-independent**, identical to four decimals while the threshold
moves 14% (thr 98 vs 86, 328 vs 287). **So the reduction gets the quasi-stationary escape
mode right and its excitation wrong**: the error is a function of (gamma, eps) alone.
**How to kill / what would settle it:** an Assaf-Meerson dissipative-WKB prefactor for
the 2-D problem must reproduce a ratio-to-the-1-D-chain that is a function of the START
and independent of the absorbing threshold. That is now a falsifiable target with known
arguments, which it was not before §61. **Still open:** the prefactor itself.

**Superseded, kept for the record:** **T14-d -> §35.2: NOT CLOSED. Both dependencies survive, and that is the answer for
now.** Matched-decade sweeps do not collapse the gamma-spread (at the one genuinely
matched target, 21 decades with a 7.5% span, b still runs -0.4519/-0.4352/-0.4172
monotone in gamma), AND b moves 6.9-12.0% with window length at fixed gamma. Neither is
eliminated. **P4 failed informatively:** adding a `d/Omega` term was meant to stabilise b
if the ansatz were merely incomplete, and instead scattered it from +0.0624 to -0.7362,
because 1/Omega and ln Omega are collinear over these windows. By P4's own terms that
makes the missing-term reading unsafe too. My own decade-matching tolerance (0.35 x
target) was too loose and left three of four targets uninterpretable -- disclosed rather
than dropped. **Present statement: b ~ -0.45 +- 0.05, consistent with WKB's -1/2, with
neither dependence eliminated.**

**THE POINT OF RUNNING IT, AND IT SUCCEEDED: §35 does not depend on b.** The rate c is
stable to **0.12-0.19%** across every window and ansatz while b moves up to 12%, so
§35's closed-form disagreement moves only in the fourth decimal (1.1540-1.1561 against
the quoted 1.1553 at gamma = 0.20). **§35's 7.5-15.5% asymptotic disagreement and the
withdrawal of §28.3's zero crossing stand independently of the prefactor being
unresolved** -- rule 14, a withdrawal verified as carefully as an assertion.

**Superseded, kept for the record:** **T14-d, open: is the prefactor exponent -1/2, and why is it gamma-dependent?**
Measured -0.4484 / -0.4394 / -0.4089 / -0.3964 at gamma = 0.20/0.25/0.30/0.35, i.e.
drifting with gamma AND with lever-arm length, which are confounded here because deeper
gamma gives more decades at the same Omega cap. **How to kill:** break that confound the
§30.2 way -- hold the decade count fixed across gamma by capping the fit window rather
than the Omega grid, so lever arm is constant while gamma varies, then sweep the
opposite way (fixed gamma, varying window length). If the exponent tracks gamma at fixed
lever arm it is physics; if it tracks lever arm at fixed gamma it is the fit absorbing
higher-order curvature and -1/2 stands. **Do not propose a mechanism before that
sweep** -- §30.1 died for exactly this.

**§36: THE 7.5-15.5% DISCREPANCY WAS AN INITIAL CONDITION. Fourteen sections
misattributed it.** `V = int mu/D` integrates along the SLAVED MANIFOLD -- it describes
delta evolving with the pool on its own nullcline. Every exact run since §12 has been
started by `_setup`, which puts the pool at the ATTRACTOR's `gamma/(1+gamma)`. At
gamma = 0.25 those differ by 36.6%; at gamma = 0.20 by 46.3%.

Same network, same gamma, eps, threshold, Omega grid, solver and parameter-free
prediction -- only the pool's start moves:

    OFF the manifold (what §22-§35 compared): pred/meas = 1.0805..1.1678, mean 1.1201
    ON  the manifold (like-for-like)        : pred/meas = 0.9822..1.0047, mean 0.9942

a factor of **14.8**. And the mechanism is quantitative: the pool gap orders exactly
with the excess across gamma (46.3% -> 16.8%, 20.0% -> 8.1%). **That is the
gamma-dependence** §28.3 could only fit with a straight line and §35 re-measured
asymptotically without explaining -- it was never a property of the closed form.

**The honest framing is neither "§15 was wrong" nor "§15 was right".** Both starts are
legitimate preparations; `_setup`'s is a modelling convention, not a law. The error was
comparing a prediction about ONE against a measurement of the OTHER, for fourteen
sections, while attributing the gap to successively more refined approximations -- the
1-D reduction (§28.3), the finite-Omega window (§35), the Gaussian truncation (§35.4).

**Every number in §22.4, §28-§28.3 and §35 stands. What they were measurements OF
changes.** Untouched: §35's instrument unlock (which is how this was found), §35.1's
drift, §35.3's ill-posedness proof, §35.4's elimination of the Gaussian truncation, and
the §29-§31 identity arc.

**How it was found, because the route is the lesson.** Not by looking for it. The
`slaving_axis` experiment was built to test §28.3's attribution on a separation axis
independent of gamma, and its **P0 gate FAILED** -- T = 1 disagreed with §35's published
number. The gate existed to catch a broken instrument; it caught a fourteen-section
misattribution instead. The sweep it was gating had R^2 = 0.29 and would have been the
weakest result of the session. **Gate every new instrument against the established one
at a cell where they must agree** -- that is the entire content of this finding.

**T14-f, open: is the remaining 0.5-1.8% the 2-D minimum-action correction?** Right size,
only candidate left -- but its sign FLIPS across gamma (0.9822, 0.9907, 0.9993, 1.0047)
where a path-minimisation correction must be one-signed, since a minimum over paths
cannot exceed the value along the slaved one. **How to kill:** compute the 2-D geometric
minimum action with the full WKB Hamiltonian and test it absolutely against the
on-manifold rates. One-signed and ~1% closes it; a residual that keeps flipping sign is
numerical, and the closed form is exact to the precision of this test.


**T-COST -> §37: an optimal drive EXISTS and is protocol-robust; the cost per nat is NOT
a constant.** R = Sigma/L, the k_B spent per nat of reliability, with both sides exact
linear solves on the same generator (Q_tt Sigma = -sigma_local is the MFPT system with
the local entropy rate as its source). R is Omega-independent to 0.25%, so it is a
property of the chemistry. It diverges at both ends of the drive with an explicit
mechanism -- s/ln(1/gamma) is constant at 0.77-0.80 while c saturates at ~0.19 and then
collapses at gamma_c -- so **R ~ 0.79 ln(1/gamma)/c(gamma)** and the minimum is forced,
not fitted.

**The deciding test SPLIT.** gamma* moves only **1.6% across theta** -- which is exactly
the failure mode that killed §9.2's "dissipation minimum near gamma ~ 0.3" -- and 26%
non-monotonically across eps, a factor of 1.32 over the whole 3x3 grid. **But R* varies
203%**, collapsing with eps, because eps sets how hard the decision is and a wide-margin
decision is cheap per nat. **P4 is withdrawn: there is no universal cost per nat, and
"transistor-grade reliability costs ~564 k_BT" is eps-specific (123 at eps = 0.50, 1548
at eps = 0.25) and must not be quoted as a constant.**

The optimum is BROAD: within 5% of minimum over gamma in [0.03, 0.08], a factor 2.7 in
drive. Quoting gamma* to three digits would overstate what a minimum this flat locates.
§9.2's withdrawal STANDS -- its optimum was at gamma ~ 0.3 and was a threshold artifact.

**T-COST-a, open: is there a cost that IS margin-independent?** R divides by the
reliability bought from a particular start, which is why eps drags it. The candidate
invariant is Sigma for the FULL traverse, saddle to attractor, which has no free start
point. **How to kill:** compute it and check whether Sigma/L is eps-free by construction.
If it is, that is the founding question's number and §37's R is a projection of it. If it
is not, the cost of reliability is preparation-dependent as a matter of physics rather
than of protocol, which is itself worth stating clearly.


**T-COST-a -> §38: R WAS NOT A QUANTITY. Restoration is priced per E-FOLD OF GAIN, and
that corrects §37's optimum.** Measured directly: as the input margin eps rises,
Sigma FALLS (1094 -> 205 k_B) while L RISES (15.0 -> 101.7 nats), so R = Sigma/L collapses
by a factor of **36** for a trivial reason -- a start nearer the threshold needs fewer
reactions AND is more reliable. **Reliability is bought with input margin, which is free.
Dissipation buys GAIN.**

    G = Sigma / (Omega * ln(theta/eps))   [k_B per molecule per e-fold]

preparation-free by construction. Naive spread across the eps x Omega grid falls from R's
**3600% to 7-27%**. Estimated as the SLOPE of Sigma against Omega*ln(theta/eps) fitted
jointly, since dividing leaves Sigma's offset drifting through ln(gain).

**An interior minimum at gamma = 0.20 that is theta-INVARIANT** -- gamma* = 0.20 at
theta = 0.70, 0.80 and 0.90, spread 0.0%, with G* spanning only 4.8% (2.0358 / 1.9927 /
1.9395) against R*'s 203%. **G* ~ 1.94-2.04 k_B per molecule per e-fold of gain.**

**§37 IS CORRECTED: the optimal drive is gamma ~ 0.20 (A ~ 4.83), not 0.07.** §37's
theta-robustness stands as measured; what it robustly located was the optimum of a
construction that conflates gain with margin.

**Kept honest: the gate MARGINALLY FAILS.** Joint-fit R^2 = 0.9940-0.9995 against the
0.995 threshold fixed in advance, failing at gamma = 0.20 and 0.30, and the fitted
intercepts are NEGATIVE (-13 to -28), unphysical as Omega*ln(gain) -> 0 and marking where
the linear form gives out. **G is a very good description of how restoration is priced; it
is not a universal constant.** G* ~ 1.99 sitting near 2 is left as an observation and NOT
claimed -- the theta-trend is monotone through 2 rather than converging on it, and
§28.2's power law and §35.1's -1/2 were both structure read into fitted quantities.

**T-COST-b, open: what is the residual 7-27%?** The negative intercepts point at a
finite-size offset the linear form omits. **How to kill:** fit
Sigma = G*Omega*ln(gain) + A*Omega + B*ln(gain) + C and check whether G stabilises and the
intercept turns physical. **§35.3 is the standing warning** -- if the extra terms are
collinear over the available range, G will swing without converging and the decomposition
will be ill-posed rather than incomplete, exactly as the prefactor exponent was.


**T-COST-b -> §39: the cost has a CLOSED FORM good to ~6%, and the OPTIMUM is predictable
analytically.** Physics forbade §38's proposed four-term fit -- Sigma -> 0 as eps -> theta
kills the constant and bare-Omega terms -- and exposed ln(theta/eps) as the pure-
exponential approximation to the traversal. The unfitted prediction is
Sigma_pred = Omega * int sigma(d)/mu(d) dd from the network's own fluxes.

**P1 and P2 both FAIL**: pred/exact = 0.9050..1.1674, mean 1.0583 against a 5% gate, and
the eps-spread is 7.4-10.0%, no better than §38's. **P3 fires and is the interesting
part: the residual is FLAT in Omega** (0.0784/0.0773/0.0771 over a 3x range), so it is not
finite-size -- and a gate against the CME's own sigma_local shows my sigma's error DOES
converge (8.0% at Omega=150 -> 1.1% at 1200), so it is not sigma either. **A ~6% gap
between the deterministic path cost and the exact stochastic cost survives Omega -> infinity
and is unexplained.**

**Process failure worth recording: that gate was run AFTER the comparison, not before.**
§36 was found by exactly such a gate failing, and its stated lesson -- gate every new
instrument against the established one where they must agree -- was available and not
applied. It happened to exonerate sigma; it might not have.

**P4 HOLDS and is what survives.** Minimising the closed form over gamma with NO CME solve
puts the optimum at gamma = 0.240 against §38's measured gamma* ~ 0.20. Two computationally
independent routes agree on the optimal drive, which upgrades §38's design principle from a
measurement to a prediction.

**T-COST-c, open: what is the flat ~6%?** Not finite-size, not sigma. Untested candidates,
none preferred (rule 17): (i) the MFPT differs systematically from the deterministic
traversal because absorption selects early-fluctuating trajectories; (ii) the exact path
leaves the slaved manifold under noise, making the manifold the wrong contour;
(iii) a Jensen gap E[sigma(state)] != sigma(E[state]) that does not close. **How to kill
(i):** compare the exact MFPT to threshold against int d(delta)/mu directly -- one linear
solve, and it separates the time from the entropy.


**T-COST-c -> §39.1: CLOSED. The cost residual is entirely a TIME residual.** Candidate
(i) was named as the one to kill first because it separates time from entropy in one
solve, and it is (i). Comparing T_det = int d(delta)/mu against the exact MFPT:
T_det/MFPT = 1.0209..1.1764 (mean 1.1069) against Sigma_pred/Sigma_exact =
0.9367..1.1683 (mean 1.0796), **correlation +0.9513**, with the difference shrinking
monotonically with Omega in EVERY cell (0.0897 -> 0.0207, 0.0843 -> 0.0236, ...) to
0.7-2.4% at Omega = 700.

**The entropy RATE along the path is right; the CLOCK is wrong.** int sigma/mu
overestimates the cost by exactly the factor by which int d(delta)/mu overestimates the
first-passage time. **Candidates (ii) (the path leaving the slaved manifold) and (iii)
(a Jensen gap in sigma) are withdrawn as explanations of the COST**, though either may
still explain the time.

**Sigma = Omega * sigma_bar * T with sigma_bar right and T wrong** means the cost and the
first-passage time are now ONE problem, not two -- any improvement to the MFPT carries
straight through. It also explains why §39's optimum survived while its magnitude did
not: a smooth multiplicative factor varying 1.16 -> 1.08 across gamma moves a minimum's
location far less than its value.

**T-COST-d, open: why is the MFPT persistently BELOW the deterministic traversal?** The
gap converges to a nonzero limit (~1.16 at gamma=0.07, ~1.08 at 0.30) and grows with
landscape depth, so it is not a finite-count effect. Live candidates: (ii), (iii), and
new (iv) -- absorption at a threshold selects the leading edge of the packet, putting
first passage below mean arrival by an amount set by packet width against drift, which
need not vanish if the threshold sits on the steep part. **How to kill (iv):** compare
MFPT against int d(delta)/mu for thresholds at different theta. If the gap tracks the
local drift steepness AT the threshold rather than the path as a whole, it is an
absorption effect and the deterministic traversal is exact away from the boundary.


**T-COST-d -> §39.2: CLOSED. The gap is the SLOW-MANIFOLD LAG and the closed form is
EXACT in the slaved limit.** The manifold is defined by ds/dt = 0 AT FIXED delta, but as
delta evolves the manifold moves and the pool lags -- an O(1/sep) error, not O(1/Omega),
which is precisely why the gap survived Omega -> infinity. Tested on §36's independent
separation axis: the gap runs +0.0914 -> +0.0002 as sep goes 7 -> 621, with intercept
**1.00089 at 1/sep -> 0, R^2 = 0.9977**, and mean |ratio(400) - ratio(800)| = 0.0048 so
it is not finite-count.

**Sigma = Omega * int sigma/mu d(delta) is exact where the reduction is**, with the
entropy rate exact (§39.1) and the traversal exact as sep -> infinity.

**The coefficient does NOT transfer between axes, and I checked (rule 9).** Along the T
axis (gap)*sep = 0.6465 +- 0.0285, constant to 12.3% over 9x in sep. Carried to the gamma
axis it gives 16.2% at gamma=0.07 against 16% measured, but **5.4% at gamma=0.30 against
8%, 33% off** -- because scaling the pool pair deforms the network as well as the
separation. **The 1/sep SCALING is established on both axes; the coefficient 0.6465 is a
T-axis value and is NOT universal.** Quoting "1 + 0.65/sep" as a general law would repeat
§28.2 exactly.

**BOTH HALVES OF THE FOUNDING QUESTION ARE NOW CLOSED-FORM:**

    reliability   -ln P = 2 Omega V(x0),  V = int mu/D     0.5-1.8% on-manifold (§36)
    cost          Sigma = Omega int sigma/mu d(delta)      exact as sep -> inf (§39.2)

Neither held at the start of this session: the first was believed 7.5-15.5% wrong, and
the second did not exist.


**T14-f -> §36.1: REFUTED, and its premise was my own reasoning error.** Two candidates
eliminated first: the realised pool wobble (0.00-0.53%, adding a b regressor moves rates
under 0.03%) and the fitting window (six Omega windows give a 1.30-2.29% spread, LARGER
than the 0.5-1.8% residual -- at short windows the residual is entirely instrument).

Redone at §35 grade (Omega to 1800, 11 points, on-manifold, eps AND b controlled, with a
half-split precision check): ratios 0.9963 / 0.9975 / 0.9921 / 0.9820, **0 of 4 above 1**,
and the residual GROWS with separation (-0.37% at sep 7 -> -1.80% at sep 17) where a
deviation from the slaved manifold must SHRINK.

**The premise was a category error and it was mine.** I argued the 1-D slaved result
bounds the 2-D action from above because a minimum over paths cannot exceed the value
along the slaved one. **The 1-D reduction is not the action along a 2-D path** -- it is
the exact WKB action of a DIFFERENT process, a birth-death chain from rates projected
onto the manifold, with momentum conjugate to delta alone rather than to both coordinates.
No variational inequality relates them, so a residual below 1 was never forbidden and the
sign test the experiment was built around had no content.

**Precision against the claim:** half-split 0.06-1.39%, so gamma = 0.25 and 0.30
(-0.25%, -0.79%) are NOT resolved. Only gamma = 0.35 is clearly resolved (-1.80% against
0.06%) and it has the fewest decades (11.6 against 41.5), exactly where a systematic is
most suspect.

**§36 STANDS and is strengthened:** on-manifold 0.9820-0.9975 mean 0.9920 against
off-manifold 1.0805-1.1678. §15's closed form is right to within 2% like-for-like. Only
the EXPLANATION of the last percent is withdrawn.

**T14-g, open: what is the ~1% one-signed deficit?** It grows with sep and gamma, i.e. it
is largest where the landscape is shallowest and the fit shortest, so the leading
candidate is now INSTRUMENTAL: 11.6 decades at gamma = 0.35 against 41.5 at 0.20. **How
to kill:** hold the decade count fixed across gamma using §35.2's matched-decade
construction (already written) and re-measure. Flattens -> lever arm. Survives -> physics,
and it needs a mechanism that is not the 2-D action.


**§40: THE FIRST EXTERNAL STANDARD. AM sits ~5x from the thermodynamic bound.** Every
cost number in §37-§39 was a measurement with no external reference. The first-passage
TUR supplies one -- Var(T)/<T>^2 >= 2/<Sigma>, so Q = (Var(T)/<T>^2)*<Sigma>/2 >= 1 -- and
both sides were already exact and already built (`first_passage_moments` for Var(T),
§37's Q_tt Sigma = -sigma_local for <Sigma>, same generator and absorbing set).

**The bound HOLDS at all 32 cells**, Q in [5.39, 180.9]. That is genuine external
validation of the §37-§39 apparatus: an independent inequality the exact entropy solve
and the exact first-passage moments had to satisfy together, and did.

**AM is not near it.** Closest approach Q = 5.39 -- roughly 5.4x more dissipation than the
thermodynamic minimum for the timing precision achieved. Not a thermodynamically optimal
decision element; its ubiquity is presumably about robustness, simplicity or speed.

**The TUR optimum and the COST optimum are DIFFERENT.** Q falls monotonically toward small
gamma with its minimum at the grid EDGE (0.05, unbounded below), where §38's cost per
e-fold has an INTERIOR minimum at 0.20. Two objectives, two optima: harder drive makes the
decision more deterministic (timing variance falls faster than Sigma rises) which the TUR
rewards without limit, while gain-per-dissipation pays a ln(1/gamma) penalty that
eventually dominates. **§38's optimum is real and is NOT the thermodynamic one.**

**Kept honest:** Q's Omega-spread is 1.7-4.7% for gamma <= 0.20 but 16.6-35.6% for
gamma >= 0.25, so the shallow-landscape cells are NOT converged and are indicative only.
The small-gamma end, where the conclusion lives, is converged.

**T-TUR-a, open: is the two-sided bound tighter?** The 5.4x gap rests on a ONE-SIDED
inequality applied to our TWO-SIDED absorbing set (|delta| >= thr). The bound HOLDING is
robust to that -- applicability could only matter if it were violated -- but reading 5.4
as the distance from optimal assumes the form is tight here. **How to kill:** re-run with
a one-sided absorbing condition (delta <= -thr only) so the standard TUR applies verbatim,
and compare Q. If Q drops toward 1 the gap was the boundary convention, and AM is far
closer to thermodynamically optimal than §40 reports.


---

**T-COST-n -> §82: CLOSED. The escape action has no thermodynamic price, and the closure is
two-way.** §80 showed the exponent governing reliability and depth is the escape action A;
§3's four earlier closures all priced quantities that do not govern it, so the founding question
was reopened by §80 whether or not anyone said so. §82 asks it of A directly, using the fact that
both networks have a ONE-DIMENSIONAL cycle space, so the entire non-equilibrium drive is a single
number that can be pinned by construction:

  * **At the affinity held fixed to 4.9e-14 nats, A spans a factor of 926** (Schloegl, exact
    quadrature). At AM's affinity held fixed to 4.4e-16 nats, a factor of 5.0.
  * **At A held fixed to 0.0000%, the affinity spans 0.92 nats.** So the drive does not set the
    action AND the action does not report the drive. Earlier closures established only the
    first direction.

**How to reopen:** find any network whose cycle space is one-dimensional and whose A is constant
on the affinity level set. The measurement is exact and cheap in 1-D, so this is falsifiable by
one counterexample rather than by a fit.

**T-COST-n's own rule-9 failure, in the section whose method IS rule 9.** P5 asked what does set
A if the affinity does not. Swept over r1 in [0.05, 0.35] -- a 1.4x move in the saddle's relative
position -- A spanned 3% and the verdict printed "Delta alone reproduces A". Widened to the full
admissible range (32.7x) at three separate Delta, A spans up to 1.754x. **The narrow sweep was
not wrong about its own data; it was wrong that its range was the range that mattered, and
nothing in the output said so.** Both numbers are kept and both are in the test.

The corrected reading is sharper than either prediction: **Delta and the affinity both fail to
determine A, by margins three orders of magnitude apart -- 1.75x against 926x.** The drive has
essentially no grip; the landscape's headline number has most of one.

~~**T-COST-o, open: is A a functional of the field alone?**~~ **-> §83. CLOSED, and by a
sharper test than the one it was opened with, which settles more than it asked.**

The kill test as opened was to match Delta and the saddle placement and vary the curvature. The
sharper version: **A depends on lam and mu separately; the ODE depends only on their difference.**
So a reaction pair adding the SAME function to both leaves the deterministic dynamics exactly
unchanged while moving A -- and mass action supplies such pairs directly, since a birth and a
death with identical REACTANT COMPLEXES have identical propensities. Schloegl: X -> 2X and
X -> 0 at equal rate, both c*n. AM: X+Y -> 2X and X+Y -> 2Y, both c*n_X*n_Y/Omega.

  * ODE invariance gated on BOTH substrates: 4.9e-15 relative for Schloegl over 41 points x 6
    rates, 1.1e-16 for AM over 200 random simplex points x 3 rates.
  * Quadrature checked against the EXACT first passage of the modified chain at every rate,
    worst 0.215% (rule 16 -- otherwise a quadrature drifting with c fabricates the result).
  * **A falls 5.6x and eta falls 4.4x, monotonically and in the same direction**, by two
    unrelated routes. AM falls 3.52x at gamma = 0.35 and reproduces §82's A to four figures at
    c = 0, pinning the two sections' instruments together.

**The consequence is larger than T-COST-o.** This project's method is the gap between the ODE and
the CME. Two networks sharing an ODE exactly and differing in A means **the ODE does not determine
reliability, and does not determine the maximum composition depth** -- at Omega = 6400 the ceiling
spans e^137 between elements no ODE measurement can distinguish. §78/§81's operational claim
survives (A is still quadrature, no master equation); **its gloss did not.** "A property of the
deterministic field" was too loose: the Lyapunov route needs D = S diag(a) S^T, which is the
propensities and not the drift, and so does A. A is a functional of the PROPENSITY PAIR.

**And it sharpens §82 rather than qualifying it.** Both neutral pairs are irreversible as written,
so their dissipation is formally infinite: the element pays an UNBOUNDED thermodynamic cost and
gets WORSE reliability for it. Drive and reliability move in OPPOSITE directions inside a single
network, which is stronger than §82's independence.

**T-COST-p, open: how much of the depth ceiling is propensity and how much is drift?** §83 shows
the drift does not determine A but says nothing about how the two divide. **How to kill:** hold
lam + mu fixed while varying lam - mu, the orthogonal sweep to §83's, and see whether A is
constant on THAT level set. If it is, A is a functional of (lam+mu) alone and the pair collapses
to one scalar function; if not, both halves matter and the ceiling needs the full pair.

**T-CASC-c -> §91: COMPOSITION IS NOT A PROPERTY OF THE ELEMENT. It is a property of the
transfer function between stages, and the figure of merit is the NOISE MARGIN in units of the
upstream rail width.**

The founding claim -- error does not accumulate as you compose restoring gates -- was the half of
this project that had never been tested. Every cascade built here (§12, §71, §72) coupled stages
through a Gaussian readout channel with sigma imposed by hand; §75 argued the physical cascade has
no readout at all, and §76-§79 then computed epsilon for ONE element and used D_max = c*/epsilon.
**The chain was never built.** §91 builds it.

Three couplings, ALL exactly neutral at the rail (each reproduces the isolated element to 0.0e+00
when the upstream is correct), behave completely differently:

  * **source** (drives the influx): does NOT transmit. The downstream keeps its high rail when
    the upstream is at r1, because the autocatalysis carries the landscape and the source term is
    a correction.
  * **catalytic** (drives the autocatalysis): transmits, margin **0.88 sigma** -- less than one
    upstream fluctuation -- and the per-stage error is **17.7-63.5x** the isolated element's.
  * **hill** (saturating, half-max at the saddle): transmits, margin **3.39 sigma**, per-stage
    error **1.71-2.29x**.

**The controlling variable is the margin, tested by collapse (P6, rule 17).** Moving it two
independent ways -- Hill exponent and half-max -- over OVERLAPPING ranges, 14 points:
log(per-stage penalty) = **-0.952 x margin/sigma**, RMS residual 0.165 in log, exponent knob 0.127
and half-max knob 0.196. Margin 1.81 sigma -> 16.0x; 4.70 sigma -> 1.06x. **Each added sigma of
noise margin costs a factor of e.**

All three chains use the SAME element, rails and Omega. **This is what a transistor's saturating
gain curve supplies and what the external channel imposed by hand.**

**Consequence for §76-§81:** their depths use an isolated element's epsilon and are right only in
the large-margin limit -- overestimating ~2x per stage at 3.39 sigma and ~46x at 0.88 sigma.

**§91.1: A NEUTRALITY GATE IS HALF A GATE, and the missing half produced a beautiful false
result.** The first version gated only neutrality, measured source coupling, and found accumulation
strongly SUBLINEAR -- an added stage costing 0.095 of an isolated element's error, total 0.53-0.61x
the union bound. That reads as restoration filtering upstream errors: the founding claim confirmed.
Decomposing the joint law exactly gives **P(stage 2 low | stage 1 low) = 0.0017 .. 0.136**: the
flip never propagates. **"Error does not accumulate" was true because nothing was connected.** No
accumulation number means anything without the transmission probability printed beside it.

**§91.3: two more broken criteria.** P6 first demanded STRICT monotonicity of the penalty in the
margin and failed on a tied pair (margins 2.59 and 2.61 giving 4.845 and 4.860). And both half-max
points fell OUTSIDE the exponent sweep's range, so the collapse -- the whole content of P6 -- was
never tested; rule 19's extrapolation trap in a new costume.

~~**T-CASC-d, open: is the margin law derivable?**~~ **-> §92. PARTLY: it is a FROZEN-upstream
statement, and there is a second control variable the margin does not contain.**

The downstream barrier depends on the upstream concentration, A = A(x_up), vanishing at the
collapse. Two averages are possible and they are NOT the same:

    frozen upstream   <exp(-[A(x_up)-A(r3)] Omega)>   = 4.845
    fast upstream     exp(-[<A(x_up)>-A(r3)] Omega)   = 1.139

  * **P1**: A(r3) equals §80's isolated quadrature to 1.1e-16 for all three couplings.
  * **P2**: A ~ (x_up - x_crit)^p with p_eff(window) = 1.452 (catalytic), 1.214 (hill),
    bracketing the saddle-node 3/2. Window-dependent by construction (rule 21).
  * **P6**: scaling stage 1's propensities leaves its landscape, barrier, rail width and
    stationary law IDENTICAL and changes only its clock. The penalty runs
    **4.244, 4.385, 4.442, 4.155, 3.302, 2.356, 1.645** over speed 1/8 .. 32 -- a factor of 2.7
    at FIXED margin. The slow plateau (4.24-4.44) matches the frozen formula (4.845) to 10% with
    no fitted parameter; the fast tail falls toward the mean-landscape rate (1.139).

**The frozen formula is the SLOW limit, not the fast one -- the two labels were backwards when
this was opened -- and the fall as the upstream speeds up is MOTIONAL NARROWING**: fluctuations
faster than the downstream's response average out before it can act. The crossover sits where the
upstream correlation time meets the downstream response time; here they are the same element, so
speed ~ 1 sits at the top of the curve.

**So §91's margin law holds at fixed timescale ratio.** All of §91's variants ran at speed 1,
inside the plateau, which is why one variable sufficed there. **A cascade whose upstream stages are
FAST relative to the downstream is protected by motional narrowing** -- a design principle
contained in neither the element nor the margin.

**§92.1: three confounds, each of which changed the answer, two of which had already printed a
verdict.** (a) Seeding stage 1 as a delta at its rail left the slow cells un-SPREAD, and the trend
read 1.171 -> 2.269 RISING; seeded from the quasi-stationary law it reads 3.010 -> 2.524, FALLING.
**The equilibration artifact reversed the sign.** (b) P(s2 low, s1 high) conditions on stage 1
surviving to the end, excluding exactly the trajectories that dipped, with the exclusion growing
with speed -- worst where the signal is. (c) Pre-equilibrating the joint chain AGES stage 2: at
speed 1/8 it had 14.6% error before the window opened. All three produced smooth, plausible,
monotone numbers.

~~**T-CASC-f, open: where is the crossover, and can it be placed?**~~ **-> §93. CLOSED. It sits at
tau_up = tau_down, and composition is a function of TWO single-element numbers.**

  * **P1, and the first version failed it.** Held at fixed WALL time while scaling the downstream
    clock, the diagonal s_up = s_dn = c ran 26.295, 4.442, 3.340, 2.941 -- because scaling the
    downstream changes how much of ITS OWN clock has elapsed. With the window set in the
    downstream's PROPER clock (t = t0/s_dn) both stages see identical proper time and the diagonal
    is 4.442 at every c, to machine precision. Nothing was read from the grid until it was.
  * **P2 is ALGEBRA and is labelled as such.** Q = s_up Q1 + s_dn Q2 with the coupling inside Q2,
    so Q t = t0[(s_up/s_dn) Q1 + Q2] depends on the speeds only through their ratio, exactly. The
    measured collapse (every ratio group at spread 1.000) is a WIRING CHECK -- it would catch the
    coupling in the wrong block -- **not a discovery**. §84 published its own Taylor term as a
    phenomenon; this is the same trap one level up.
  * **P2b/P3, the physics, which the algebra does not fix.** Penalty against tau_up/tau_down:
    4.210, 4.223, 4.293, **4.442**, 3.302, 1.835, **1.580** over 1/64 .. 64. **Plateau for a slow
    upstream (flat to 5.5%, matching §92's frozen formula 4.845 to 10% unfitted), a turn at ratio
    ~ 1, a fall above it.** The crossover is not adjustable: it sits where the upstream correlation
    time meets the downstream response time, so §92's curve was centred at speed 1 as a consequence
    of its two elements being identical, not by coincidence.

**THE COMPOSITION RESULT, as it now stands.** The penalty is a function of the **noise margin**
(§91) and the **timescale ratio** (§92, §93), and both are properties of SINGLE elements. The exact
joint CME is needed to validate that, not to compute it -- which is what made §91's D = 2 ceiling
look binding. **A cascade is safest when its upstream stages are fast and its transfer function
saturates**: motional narrowing and noise margin, neither of them a property of the restoring
element itself.

~~**T-CASC-g, open: does the two-number description survive depth?**~~ **-> §94. The LAW composes;
the WIDTH SEQUENCE it needs does not, and three readings of that sequence are all wrong.**

  * **P1**: stage 1's (mean, sd) is identical at D = 2 and D = 3 to 1e-9 -- the coupling is
    one-way as built.
  * **P3, and it HOLDS.** §91's law, fed stage 2's measured width, predicts stage 3's penalty as
    5.6514 against a measured 6.3494 -- **11% low, inside §91's own 18% scatter**. So the
    two-number description composes stage to stage.
  * **The penalty GROWS with depth**: 4.4419 -> 6.3494, because the margin is measured in units of
    the upstream rail width and each output is wider than its input.

**§94.1: three readings of the width sequence, all refuted.** (a) Geometric erosion --
sigma_2/sigma_1 = 1.0855 extrapolated per stage, giving a hard depth cap -- refuted BEFORE running,
since it contradicts what a restoring element is (LNA: sigma_out^2 = sigma_intr^2 + g^2 sigma_in^2
with g < 1). (b) Convergence to the LNA fixed point -- refuted by measurement: unconditioned widths
are **0.49922, 0.54181, 0.70814**, and sigma_3 blows past both the recursion (0.54904) and the
fixed point (0.55071). (c) That the DETERMINISTIC gain explains the growth -- it does not:
g^2 = 0.0109 and 0.0135 at the measured operating points against the 0.1783 the width ratio
implies, **thirteen times too small**.

**Robustness**: the overshoot of the fixed point grows with Omega -- 1.015, 1.034, 1.086, **1.286**
at Omega = 14, 16, 20, 30 -- so it is not a small-Omega artifact. (The ACCELERATION is an Omega = 30
statement only: the ratio-of-ratios runs 0.941, 0.954, 0.999, 1.204.)

~~**T-CASC-h, open: a stochastic mean-shift compounds down the chain**~~ **-> §95. PARTLY
CONFIRMED, and the problem is now localised to the MEAN.**

§94's sharpest anomaly was that the deterministic gain is 13x too small to explain sigma_2/sigma_1
at the FIRST step. The missing term is CURVATURE. The stage-to-stage map F (the downstream high
rail given a pinned upstream) is concave near the rail because the Hill coupling is near
saturation, so by Jensen <F(x_up)> falls short of F(<x_up>) by (1/2) F'' sigma^2, lowering the
operating point onto a flatter part of the landscape. **Second order in sigma, which is why a
first-order gain missed it.**

  * F(r3) = r3 to 4.4e-16; F' = 0.1042 and **F'' = -0.17265**, converged in the finite-difference
    step to five decimals over h = 1e-2..5e-4 (rule 13). Concave, so Jensen pushes the right way.
  * **The mean, nothing fitted:** mu_2 predicted 3.00750 without Jensen, **2.98599 with**, against
    a measured 2.95635. The residual goes -0.05115 -> -0.02964: **curvature closes 42%**, right
    sign. The other 58% is not explained.
  * **THE WIDTH IS NOT A SEPARATE PROBLEM.** Placed at its exactly-measured operating point the LNA
    gives sigma_2 = 0.53385 against 0.54181 -- **1.5%** -- and the LNA's own error on stage 1,
    whose input is a noiseless chemostat, is 2.7% the other way. **The width follows from the
    operating point within the LNA's own accuracy, so the width sequence IS the mean sequence.**

Both fractions are reported and not tuned: §22's convolution was fitted for three subsections and
was out by 3688x against a quantity computable exactly.

~~**T-CASC-i, open: what supplies the remaining 58% of the mean shift?**~~ **-> §96. 89% of it was
the SECOND-ORDER TRUNCATION, not physics.**

§95 used F(mu) + (1/2) F'' sigma^2 -- the frozen-upstream limit, truncated. Computing the frozen
average EXACTLY over stage 1's own law:

    §95's truncation      mu_2 = 2.98461
    <F(x_up)> exact       mu_2 = 2.94799      measured 2.95165
    G(<h>) fast limit     mu_2 = 2.96416

**The residual falls from -0.03296 to +0.00366 -- 0.12% of the mean -- and changes sign.** The two
limits are genuinely different objects (the coupling enters only through the Hill factor h, so the
rail is G(h); frozen averages the OUTPUT <G(h)>, fast averages the INPUT G(<h>), and Jensen applies
to a different function in each), and **all eight measured means lie between them** as the upstream
clock is swept 256x. So the finite correlation time is what the remaining 0.12% is made of.

**Consequence:** since §95 showed the width follows from the operating point within the LNA's own
accuracy, and §96 predicts the operating point to 0.12%, **the width sequence is predicted end to
end** -- which is exactly what §94.1(b) wrongly denied.

**§96.1: TWO INSTRUMENT BUGS in §94/§95, one of which reversed a headline.**

  * **(a) THE SEED.** `seed()` indexed the last stage by its raw molecule count -- right for a FREE
    stage (range from 0), wrong when every stage is reflected (indexed by position in the reflected
    range). The all-reflected run seeded stage 3 at a mean of **1.1333 against a rail at 3.1827**,
    essentially at its saddle, and it relaxed UPWARD through the window. **§94's sigma_3 = 0.70814
    and the Omega-scaling table built on it are artifacts.** Corrected: widths 0.50120, 0.54965,
    0.58375, and the LNA recursion predicts sigma_3 to 4.3%. ~~§94.1(b): "the widths do not
    converge and BOTH simple readings are refuted; the width sequence is not predicted by anything
    here"~~ -- **WITHDRAWN.** §94's own P4, which had been marked superseded, was right.
  * **(b) THE STATISTIC.** `stage_stats` applied a high-side filter to every stage. Right for a
    free stage (condition on not having escaped), wrong for a reflected one, where it merely drops
    the boundary lattice site: stage 1 read 3.02222 against an exact stationary 3.02117. **Caught
    by a P1 gate in §96 that FAILED** -- two independent builds of the same two-stage system
    disagreed, including on stage 1, which is provably invariant.

**What survives unchanged:** the one-way coupling; §94 P3 (now 8.5% low, was 11%); §94.1(a) and
(c) -- the deterministic gain is **18.6x** too small, was 13x; and §95's structure, with curvature
closing 40% of the mean shift instead of 42%.

~~**T-CASC-j, open: does the end-to-end prediction hold at a second element?**~~ **-> §97. THE
DERIVATIONS TRANSFER; THE FITTED SLOPE DOES NOT.**

Run forward on roots 0.20/1.0/4.3505 (A = 0.439977 vs 0.190241) with Hill n = 6, K = 1.3 (vs
n = 4, K = 1), nothing measured on the cascade, against one exact joint solve:

    mean      predicted 4.26779   measured 4.26781    -0.00%
    width     predicted 0.51622   measured 0.50236    +2.76%
    penalty   predicted 1.15325   measured 4.47424   -74.2%

**The mean is right to five decimals on an element it was never calibrated on and the width is
inside the LNA's own accuracy. §91's slope -0.952 is 74% wrong, so it is a FIT and not a law.**

**But the derivation behind the fit does transfer.** §92's frozen average carries Omega and the
barrier explicitly where §91's slope carries only the margin, and this element's barrier is
A*Omega = 13.2 against the calibration's 5.7:

    frozen 8.92   |   measured 4.47   |   fast 1.01

The two limits bracket it, exactly as in §92. **So §91's margin law is a ONE-ELEMENT
PARAMETERISATION of §92's average, and the margin alone is not the controlling variable -- the
barrier depth A*Omega enters too.**

**§91's headline needs that qualification, and it is rule 9 at the top level.** The 14-point
collapse onto log(penalty) = -0.952 x margin/sigma was real, but every one of those points had the
SAME ELEMENT and therefore the same A*Omega. **Constancy along the axes swept is not constancy, and
the axis not swept was the barrier depth.** It is also rule 19's extrapolation trap: the new
element's margin is 5.38 sigma against a calibrated range of 1.81-4.70.

**WHERE THE COMPOSITION ARC STANDS.** Mean and width are predicted end to end from single-element
quantities, out of sample. The penalty is BRACKETED by two computable limits whose separation is
the upstream's correlation time. **There is no closed form for the penalty**: §91's slope is not
it, and §92's average gives two numbers rather than one.

~~**T-CASC-k, open: is there a closed form for the penalty?**~~ **-> §98. CLOSED. The position
between the two limits is a function of the BARRIER DEPTH, and it transfers.**

Omega is the clean axis: it moves A*Omega without touching the landscape, the transfer map, the
margin in concentration units or the coupling. So the curve was traced on ONE element and the
OTHER checked against it, nothing fitted.

  * **P2**: position (log space, 0 = fast, 1 = frozen) runs **1.3089, 1.1369, 0.9400, 0.8221,
    0.7096, 0.6352, 0.5805** over A*Omega = 2.66..16.17 -- monotone falling. Deeper barriers sit
    further from the frozen limit, as an excursion-duration argument gives: a deeper barrier needs
    a deeper, rarer and SHORTER upstream excursion, so fewer last long enough for the escape to
    complete and the rate-average over-counts.
  * **At A*Omega = 2.66 the position EXCEEDS 1** -- the measurement is outside the bracket on the
    shallow side. The two limits are asymptotic and the bracket is not guaranteed at 2.7 nats.
    **-> §99(a) supplies the reason, and it is not asymptotics.** There are THREE regimes, not two:
    the mean rate `<exp(-dA*Omega)>` is the INTERMEDIATE one, and the truly quenched limit is
    `<1/k>`, dominated by the highest barrier the upstream presents. A measurement above the
    mean-rate "limit" is allowed and expected. So "position between the limits" over-reads the
    variable: it is a coordinate that runs through [0,1] over most, not all, of the range. The
    finding -- that it is a reproducible function of A*Omega which transfers out of sample --
    is unaffected.
  * **P3, out of sample**: §97's element (different landscape, different coupling, barrier 2.3x
    deeper) sits at A*Omega = 13.2 with position 0.6831; the first element's curve predicts
    **0.6383, off by 6.6%**.

**So the composition penalty CLOSES: frozen, fast and the position between them are all computable
from single-element quantities**, and §97's bracket becomes a number. And §91's slope is now fully
accounted for -- it was a one-element parameterisation whose missing variable is A*Omega, invisible
to it because all 14 of its points shared one element.

**§98.1: np.interp printed FAILS off the right data, for the third time in this project.** P3's
first verdict was FAILS with a "predicted position" of 1.3089 -- **exactly the first row of the
sweep**. np.interp requires increasing x; `aoms` already increases with Omega and it was passed
reversed, so it returned fp[0] and the comparison was against the SHALLOWEST cell. The tell was
that the prediction was a row of the table rather than a value between rows. Same family as
§86.1(3) and rule 19's §59. **np.interp is the single most frequent source of broken comparisons
here and has now produced three wrong verdicts.**

**SCOPE OF THE COMPOSITION RESULT, as it now stands:** depth is computable **given the widths**;
the width sequence is open.

**T-CASC-e, open: does the margin law survive depth?** §91 is D = 2. The exponential-in-margin
penalty compounds if each stage contributes independently, which would make D_max = c*/(penalty x
eps). **How to kill:** D = 3 is 1.77M joint states and reachable (a run was started and stopped
for time, not feasibility). Check whether the per-stage penalty at stage 3 equals stage 2's.

## 5. Where this sits in the literature

The project has run without citations. That is a defect, not a style: several results
here have known counterparts, and three open questions have published frameworks that
are better routes than anything attempted so far. Recorded with what each bears on.

### The network itself

**Approximate Majority is not ours and is well studied.** Angluin, Aspnes & Eisenstat
(2008) introduced it as a population protocol and proved O(log n) convergence with
correctness for an initial gap omega(sqrt(n) log n); Condon et al. (UBC) give a
tri-molecular analysis improving this to Omega(sqrt(n) log n). Cardelli &
Csikasz-Nagy (Sci. Rep. **2**, 656, 2012) showed the eukaryotic cell-cycle G2/M switch
computes AM, and Dodd et al. (2007) show the histone M/U/A epigenetic switch is
structurally identical to it.

**What this means for CRNL's scope.** That literature analyses CONVERGENCE TIME and
CORRECTNESS PROBABILITY. **It does not price the switch thermodynamically** -- there is
no affinity floor, no cost per e-fold, no dissipation-reliability trade in it. §9.1,
§20, §37-§40 are therefore complementary rather than duplicative. **But note the regime
mismatch:** their guarantees are for a gap ~ sqrt(n) log n, i.e. eps ~ log(Omega)/sqrt(Omega)
-> 0, whereas every CRNL result uses FIXED eps. Claims should not be transported between
the two without saying so.

### T14-e (the prefactor) -- the framework exists and numerics provably cannot substitute

§35.3 proved the algebraic prefactor is unextractable numerically: the asymptotic basis
functions are 90-99% correlated over any bounded Omega range. **Assaf & Meerson**
(J. Phys. A **50**, 263001, 2017 -- review; Phys. Rev. E **81**, 021116, 2010) develop
exactly the dissipative-WKB machinery for master equations and explicitly "yields both
entropic barriers to extinction AND PRE-EXPONENTIAL FACTORS, and holds for a general set
of multistep processes WHEN DETAILED BALANCE IS BROKEN" -- our case precisely (gamma < 1).
**This is the route for T14-e**, and it is analysis, not fitting.

### T14-g / §39.2 (the slow-manifold gap) -- it has a name, and it is a known phenomenon

**Borner, Deeley, Romer, Grafke, Lucarini & Feudel (arXiv:2311.10231), "Saddle avoidance
of noise-induced transitions in multiscale systems":** noise-induced transitions in
NON-GRADIENT systems with FAST AND SLOW degrees of freedom "may deviate significantly
from the FW instanton even for noise so weak that transitions become extremely rare",
and they "highlight the link between saddle avoidance and timescale separation".

**That is our situation exactly** -- gamma < 1 breaks detailed balance (non-gradient),
the pool is fast against the signal (multiscale), and §39.2 measured a deviation that
persists as Omega -> infinity and closes as 1/sep. **Our O(1/sep) law is very likely an
instance of saddle avoidance**, which would explain why it is not a finite-count effect.
**Berglund & Gentz** (J. Diff. Eq. **191**, 1-54, 2003) give the rigorous slow-manifold
concentration results that would make the correction derivable rather than fitted.

### T-TUR-a (the two-sided bound) -- partially CLOSED, and a sharper test exists

Two results bear on it:

  * **Pal, Reuveni & Rahav (arXiv:2103.16578)** derive an FPT TUR valid for arbitrary
    initial conditions AND absorbing states: CV^2 >= 1/(Sigma/2 + 1). §40's `2/Sigma` is
    its large-Sigma limit; measured, the refinement moves Q by 0.13-0.22% (§40.1).
    **So §40's bound choice is justified and its conclusion robust.**
  * **Neri (SciPost Phys. 12, 139, 2022)** derives bounds for first passage of a current
    with TWO thresholds -- our geometry -- expressing uncertainty via the SPLITTING
    PROBABILITY rather than the variance. Crucially the bounds become **EQUALITIES when
    the thresholded current is the stochastic entropy production**, via the martingale
    property of exp(-S), giving **p_- = exp(-l_-)** independent of everything else.

**T-TUR-b, open: does the martingale equality hold here?** CRNL computes splitting
probabilities exactly and to 1e-33 (§35), so p_- = exp(-l_-) is directly checkable -- an
EXACT parameter-free identity of a kind this project has never tested. **The obstacle is
construction, not precision:** the threshold must be on the accumulated ENTROPY
PRODUCTION, not on delta, which requires an augmented chain (state x accumulated S). That
is a real build. **How to kill:** construct it at small Omega where the augmented state
space is affordable and check p_- against exp(-l_-). Agreement would validate the entropy
solve against an exact external identity; disagreement at small l_- is expected since the
result is asymptotic in the thresholds.

### §37-§39 (the cost) -- the logarithmic divergence is a known general result

**Ouldridge, Govern & ten Wolde (Phys. Rev. X 7, 021004, 2017)** show a canonical
biochemical readout network faces an accuracy-dissipation trade-off strictly worse than
the thermodynamic bound, with "the thermodynamic cost diverging LOGARITHMICALLY as
accuracy approaches 100%". **§37's R diverging and §38's ln(gain) structure are instances
of this**, which is reassuring for the measurements and deflating for their novelty --
what is new here is the OPTIMAL DRIVE (§38, gamma ~ 0.20-0.24), not the divergence.
Rao & Peliti (JSTAT P06001, 2015) and Sartori & Pigolotti (PRL 110, 188101, 2013) give
the kinetic-vs-energetic discrimination framework that §38's gain/margin split resembles.

### §80-§83 -- checked EARLY this time, and the check came back badly (§85.1)

**§83's construction is published.** Plesa, Zygalakis, Anderson & Erban, *Noise control for
molecular computing*, J. R. Soc. Interface 15:20180199 (2018), arXiv:1705.09392, name these
**"zero-drift networks"**; their Eq. (5) `s + s_bar -> 2s` / `s + s_bar -> 2s_bar` at equal rate is
§83's AM pair verbatim, and they state the deterministic-invariance/stochastic-difference property
in the same words.

**The framing has a settled vocabulary §83 did not use.** Identical mass-action ODEs =
**dynamical equivalence** (Horn & Jackson 1972, *macro-equivalence*) or **confoundability**
(Craciun & Pantea, J. Math. Chem. 44:244, 2008 -- canonical). That such networks differ
stochastically: Enciso, Erban & Kim (arXiv:2006.02272, EJAM 2021) prove the CTMC rates determine
the network uniquely, so the stochastic model is strictly finer; Faul, Hoessly & Xia
(arXiv:2505.07638) give the sharp criterion -- the ODE sees y'-y, the generator also sees
(y'-y)(y'-y)^T -- and **their Example 5.12 is §83's one-species construction**.

**§80/§81's action formula is textbook.** `A = -int ln(mu/lam) dx` is Assaf & Meerson, *WKB theory
of large deviations in stochastic populations*, J. Phys. A 50:263001 (2017), Eqs. (43)-(44). It was
derived here independently and should have been cited from §80. The same review states §85's S2 in
general: the Fokker-Planck reduction is not a controlled approximation to the action.

**A published near-instance of §83's phenomenon.** Assaf & Meerson §IV's bursty Verhulst model has
"the deterministic rate equation coincides with" the non-bursty one, with an exponentially
different mean time to extinction. **§83 is not the first instance.**

**A scope statement §83 should have carried.** Anderson, Craciun & Kurtz (arXiv:0803.3042):
complex-balanced networks have product-form Poisson stationary distributions AT the deterministic
equilibrium, so for that class the ODE DOES determine the stationary distribution. §83's
construction escapes only because X -> 2X, X -> 0 is neither weakly reversible nor complex-balanced.

**What is left, and it is narrow.** No source was found computing an escape action or switching
rate for a dynamically equivalent PAIR, nor drawing a composition-depth consequence. Per §70's own
rule that is **not novelty** -- it is a one-line corollary of two literatures that have not been
joined. The residual claim is the quantification and the CRNT framing, not the construction and
not the mechanism.

### §84-§85 (nu = 2 and the escape action near the pitchfork) -- mostly known, and one part is our own algebra

**nu = 2 is a standard universal exponent.** Dykman, Maloney, Smelyanskiy & Silverstein, PRE
57:5202 (1998): the activation barrier goes as |mu - mu_b|^xi with xi = 2 at a pitchfork and 3/2 at
a saddle-node, after Dykman & Krivoglaz (1979). Measured at 2.0 +- 0.1 and 2.00 +- 0.03 by Chan &
Stambaugh, PRL 99:060601 (2007). **It does not discriminate**: transcritical also gives 2.

**The Langevin-to-CME transfer is citable, not a gap closed here.** Assaf & Meerson (2017) §II: the
van Kampen/Fokker-Planck Hamiltonian, wrong in general for large deviations, COINCIDES with the
exact WKB Hamiltonian near a bifurcation -- which is why §85's S2 found FP and WKB agreeing to
0.30%. Their §VII A 1 works the transcritical case by §84's own method.

**The effective-exponent bias is named and standard.** Kouvel & Fisher (1964); Riedel & Wegner,
PRB 9:294 (1974); Wegner, PRB 5:4529 (1972). For escape near a bifurcation specifically: Dykman &
Ryvkine, PRL 92:080602 (2004). Chan & Stambaugh's own fits give 2.0 near the pitchfork and
1.43/1.53 further out -- the same story as §84's 1.99 vs 1.83.

**What was not found:** any escape action or quasipotential for Approximate Majority, reversible or
not -- checked across the molecular-programming literature (Cardelli & Csikasz-Nagy 2012; Cardelli
et al. 2016; Angluin-Aspnes-Eisenstat 2008; Condon et al. 2017/2020) and the voter-model analogues
(Castello, Eguiluz & San Miguel 2006; the noisy/partisan voter models). Per §70 that is not novelty.

### The methodological point

**Three of the open questions had published frameworks and one had a published
phenomenon, and the project found none of them because it never looked.** The measurements
stand on their own -- that is what the discipline bought -- but §35.3 spent a session
proving numerically that a prefactor is unextractable when Assaf & Meerson had published
how to derive it, and §39.2's 1/sep law was measured without knowing it likely has a name.
**Search the literature when a question is NAMED, not when it is finished.**


**T-TUR-b -> §41: CLOSED. An exact identity validates the entropy machinery.** §37-§40
all rest on sigma_local, the reverse pairing and the sign convention, and nothing had
checked that object externally -- §40 tested an INEQUALITY, which an entropy wrong by a
constant factor would still satisfy.

The integral fluctuation theorem <exp(-S_tot)> = 1 holds at ANY stopping time, including
the absorption already used, and needs no augmented chain: tilting each transition by its
entropy weight collapses to a_j exp(-Delta s_j) = a_rev(n'), so the tilted generator is
built from REVERSE propensities at the same cost as every other solve.

**Verified to 5.5e-14 at the best-conditioned cell, median 1.33e-9 over 36 cells**, with
solve residuals ~1e-15. The medium term alone spans 3.7e-11 to 0.92, nowhere near 1 (P1).
**sigma_local, the pairing and the sign convention are validated together against an exact
identity** -- the first equality, as opposed to inequality, the entropy machinery has met.

**Fails at gamma = 0.10 as CONDITIONING, not physics:** deviation by gamma runs 1.3e-12 ->
1.5e-9 -> 4.2e-3 as gamma falls 0.30 -> 0.10, with one cell overflowing to 5.9e+280,
tracking the dynamic range of E[exp(-S_med)] (0.92 -> 3.7e-11) exactly. The tilted solve
carries the full range of pi, so a deep landscape destroys it.

**The test caught a real convention error -- mine.** The first pass used the boundary
pi(n_0)/pi(n) instead of pi(n)/pi(n_0) and returned 0.0008 instead of 1, a 1000x
deviation, with a 1e-15 solve residual proving the algebra exact and the CONVENTION wrong.
That is exactly the failure mode §37-§40 were exposed to and never tested for.

**T-TUR-c, open: does p_- = exp(-l_-) hold?** Still needs thresholding on accumulated
entropy production rather than delta -- an augmented chain (state x accumulated S). Better
motivated now: §41 shows the entropy bookkeeping is correct, so the remaining question is
about first-passage geometry rather than conventions. **How to kill:** build the augmented
chain at small Omega and check p_- against exp(-l_-), expecting agreement only
asymptotically in l_-.


**T14-g -> §36.2: CLOSED. The 1-D reduction is ~1% shallow and it is REAL.** Three sweeps,
because each held a different thing fixed. §36.1's "instrumental" reading is WITHDRAWN:
the deficit SURVIVES at matched decades (spread 1.64 vs §36.1's 1.55), and the two sweeps
confound oppositely -- one fixes the Omega window and varies decades, the other fixes
decades and varies the Omega-ratio -- with the deficit tracking gamma in both.

**But both shared the window's lower edge at Omega = 150**, where -lnP is 8.6 at
gamma = 0.20 and only 2.3 at gamma = 0.35 (P ~ 0.10, not a tail). Every large-gamma fit
reached into a region where the Laplace/WKB asymptotic has not taken hold. Starting every
window at P ~ 1e-4 instead HALVES the spread, 1.64 -> 0.78 points. **Half the apparent
gamma-dependence was asymptotic validity and it took a third sweep to see it.**

**Result: §15's 1-D closed form runs 0.9911 +- 0.0039 against the exact 2-D collapse --
one-signed, ~1% shallow**, resolved above half-split precision at 3 of 4 cells. A negative
deficit is PERMITTED for the reason §36.1 established: the projected-rate chain is a
different process, not a constrained path, so no variational inequality forces >= 1.

**§36 is strengthened:** on-manifold 0.9911 +- 0.0039 on properly matched windows against
1.0805-1.1678 off-manifold. Like-for-like agreement is ~1% where fourteen sections read
7.5-15.5%.

**Methodological, and it cost three runs:** "lever arm" was TWO hidden variables -- span
and starting depth -- and matching the first left the second free. Rule 9 says measure
along an axis you did not choose; this says count how many axes you did not choose,
because a matched-looking comparison can still share an unexamined edge.

### Quantum mechanics -- where this project sits, and the one thing a quantum substrate buys

Checked against the quantum-thermodynamics literature. Three findings matter, one of them
placing the project's machinery precisely and one answering a question asked earlier in the
session about basing logic on atoms.

**1. CRNL is the INCOHERENT LIMIT of quantum stochastic thermodynamics, and that is exactly
where the TUR §40 used is exact.** Nishiyama & Hasegawa, *Quantum thermodynamic uncertainty
relations without quantum corrections: a coherent-incoherent correspondence approach*
(arXiv:2505.09973, U. Tokyo). Quantum TURs generally need "quantum correction terms" because
coherence can beat the classical bound. Their CIC framework maps a Lindblad system S with
Hamiltonian H onto S-null with H set to zero, and shows entropy production and dynamical
activity are INVARIANT under the map -- so in S-null the classical-form relations

    Var[J]/E[J]^2 >= 2/Sigma        and        Var[J]/E[J]^2 >= 1/A

hold with no correction at all. **A classical Markov jump process is realisable as exactly
that H = 0 sector** (jump operators L_ij = sqrt(w_ij)|i><j|, coherences decaying
independently, the diagonal obeying the CME). So the chemical master equation is the
incoherent sector of a Lindblad equation, and §37-§41's entropy apparatus is the H = 0 case
of quantum stochastic thermodynamics.

**This retroactively justifies §40's choice of bound**, which §40/§40.1 flagged as the least
secure part of that section (T-TUR-a, "partially CLOSED"): the uncorrected classical form is
the right one here *because* there is no coherence, not by assumption. It is a placement,
not a new result, and the embedding is not unique -- so it should not be read as "CRNL is
secretly quantum".

**2. Coherence can go BELOW the classical bound -- and that is the one thing a quantum
substrate buys that chemistry provably cannot.** arXiv:2510.20873 (*Quantum Coherence as a
Thermodynamic Resource Beyond the Classical Uncertainty Bound*), arXiv:2501.00627
(*Violation of the TUR in quantum collisional models*), arXiv:2604.05747 (kinetic UR in
collective dissipative many-body systems), arXiv:2607.12264 (*Heisenberg Scaling in
Many-Body Kinetic Uncertainty Relation via Quantum Feedback*, U. Tokyo).

**§40 measured Q_min = 5.39 -- AM sits about 5x above the classical floor Q >= 1.** The
literature says a coherent element can go under Q = 1 entirely. So the honest answer to
"would atoms be better" gains a third term beyond §44's: not the energy scale (§44 says
that is not the lever), not selective catalysis (§44's suspect, which chemistry WINS on),
but **coherence, which removes the floor AM is being measured against.** That is a real
asymmetry and it is the first thing in this session that a chemical substrate cannot match
even in principle. **It is also entirely outside this rig** -- nothing here can test it, and
saying so is the point.

**3. The closest external analogue to §12/§35 is a cat qubit.** arXiv:2507.18714,
*Non-perturbative switching rates in bistable open quantum systems: from driven Kerr
oscillators to dissipative cat qubits* (Alice & Bob et al.) computes switching rates in a
bistable open quantum system by path integral -- **the same WKB-with-a-prefactor problem as
T14-e**, in a system engineered for exactly the property this project measures: a
restoring element whose bit-flip rate is exponentially suppressed in a size parameter.
Worth reading against §35.3, which proved the prefactor decomposition ill-posed
numerically; if their path-integral route gives the prefactor analytically, it is the
method §35.3 said numerics could not substitute for.

Also noted, not yet used: Zurek, *Decoherence, einselection, and the quantum origins of the
classical* (quant-ph/0105127) -- einselection is a restoration mechanism in this project's
sense, the environment destroying off-basis components and leaving a discrete pointer set,
which is "why discrete" one level below chemistry. And arXiv:2512.03770, *Quantum
Simulations of Opinion Dynamics*, builds exactly-solvable quantum consensus models -- the
nearest thing to a quantum AM.

**Nothing here changes a measurement in FINDINGS.** It places §37-§41 in a named framework,
retroactively secures §40's bound choice, and identifies the single quantum advantage that
is real and that this rig cannot reach.

### The composition arc (§91–§98) — recorded §99

Four searches were run against §91–§98 after the arc closed. **The mechanisms are known,
most for decades. The object is not found.** Full account in FINDINGS §99; what each
bears on:

**Resonant activation — §97's bracket is Pechukas & Hänggi (1994).** The two-limit
heuristic (average the rate vs rate at the average, with the Jensen argument) is stated
in those terms in PRL **73**, 2772 (1994); reviewed by Reimann & Hänggi (Springer LNP
**484**, 1997); limits named by Boguñá et al. (PRE **57**, 3990, 1998); the slow side is
Zwanzig's dynamical disorder (1990). **Already done for master equations:** Assaf,
Roberts, Luthey-Schulten & Goldenfeld, PRL **111**, 058102 (2013) — bistable gene, CME,
WKB action, OU noise at finite correlation time, optimal correlation time — is the closest
prior art in this project, and their adiabatic limit is my mean-rate limit. Also Spalding,
Doering & Flierl (PRE **96**, 042411, 2017); Flomenbom & Klafter (2004).

  * **What it corrects:** three regimes, not two (see T-CASC-k above).
  * **What it leaves open:** in every reference the modulator is **exogenous** — an OU or
    dichotomous process whose statistics are chosen. Here the modulator is a
    chemically-coupled upstream bistable stage whose statistics come from the same
    chemistry as the downstream's, which is *why* the two timescales cannot be tuned
    independently (§93). **That is the part of §91–§98 with no prior art attached.**

**Noise propagation — §94's recursion is Thattai & van Oudenaarden (2002).**
`σ²_out = σ²_intr + g²σ²_in`, its fixed point and the g < 1 condition are Biophys. J.
**82**, 2943 (2002), Eq. 13, with g the "differential amplification factor" and
attenuation length n₀ = 1/ln(1/g²). §94 rederived a 24-year-old result. **Their form
carries a time-averaging factor τ = β_x/(β_s+β_x) that mine lacks**, valid only under
separated timescales — which §92–§93 measured to be *equal* here. This killed §94.1(c)
(§4 catalogue).

**The Jensen shift is the piece that is not standard *in this combination*.**
`⟨F(x)⟩ < F(⟨x⟩)` is textbook away from cascades — Paulsson's stochastic focusing (PNAS
2000), Ochab-Marcinek & Tabaka (PNAS 2010), Hernández-García & Velázquez-Castro
(arXiv:2307.03057, whose Eq. 17 is §95's second-order truncation) — but it does not appear
combined with the variance recursion in the cascade literature. **§95 + §96 — operating
point from the exact `⟨F(x_up)⟩`, then LNA about the shifted point — is what makes §96
work where §94 did not**, and that combination was not found.

**Noise margin.** The term is used in synthetic biology (Weiss & Basu; Yaman et al. 2012)
but always statically and deterministically. Closest quantitative work is Privman & Katz
(arXiv:0803.4197, 2008): σ propagated through gate surfaces, depth ≈ 11 gates — but
*monostable* gates against an arbitrary tolerance band, not a bistable element's own
escape. **A σ-normalised margin-to-failure for a bistable cascade out of an exact CME was
not found.**

**The object.** Nothing found combining all four of: *directed* cascade (not symmetric
mutual coupling), genuinely *bistable* elements, *exact* CME, and a *depth* question.
Nearest misses: Ungarelli et al. (1999), QCA wires; Maity et al. (2024), concatenated
bistable networks, deterministic; a 2026 CMOS-bits-in-series analysis whose coupling is
symmetric and whose bits form one collective degree of freedom — the opposite sign of the
question.

**T-CASC-l, open: is §93's crossover at a relaxation time or an escape time?** §93 placed
the penalty maximum where τ_up ≈ τ_down and read it as *matched relaxation times*. Reimann
& Hänggi's review distinguishes **type-I** potentials, where the modulation moves the
barrier *top* and the resonant-activation optimum matches an **escape** time, from
**type-II** ("breathing") potentials, where the modulation moves the well *bottom* and the
optimum matches a **relaxation** time. §93 never asked which one the Hill coupling
produces. **The measurement stands; the account is a suspect (rule 17). How to kill:**
compute where the coupling enters the landscape — does raising the upstream move the
downstream's barrier position, its well depth, or only the region between? That is a
derivative of the WKB landscape already computable with `stage_kernel`, and it decides the
reading without a new run. If type-I, §93's "matched relaxation times" is withdrawn and
the crossover has to be re-placed against an escape time.

**T-CASC-m, open: is §98's A·Ω the real variable, or is it §91's margin/σ in disguise?**
Near a saddle-node the two are not independent — margin/σ ∝ μ^(3/4)/√D and the action
∝ (margin/σ)² — and **§98's Ω-sweep cannot separate them at all**, since within one
element margin/σ ∝ √Ω. §99.1 discriminated them post-hoc on the second element, where the
ratio (margin/σ)²/(A·Ω) differs by 22% (1.7931 vs 2.1945): A·Ω predicts to −6.6% and
(margin/σ)² to −15.0%. **That is one out-of-sample point with no error bar and it is not
decisive.** **How to kill:** a THIRD element chosen so the ratio differs from *both*
existing values, which pushes the two candidate predictions far enough apart that one
point falsifies one of them. The element must move the landscape and the coupling
together — Ω will not do it, by the argument above.

~~**T-CASC-l, open: is §93's crossover at a relaxation time or an escape time?**~~ **-> §100.
SETTLED, and not by the route the question named.** The reflecting wall at stage 1's saddle removes
its escape channel entirely, so **only one of the two candidate timescales exists in this chain**.
Measured on the axis §93 lacked — Ω, which moves an escape time like exp(A·Ω) and leaves a
relaxation time alone — the reflected stage-1 gap varies **1.13×** over Ω = 14→55 while the free
gap varies **1708×**, tracking exp(−A·Ω)'s 2440×. Three orders of separation. The type-I/type-II
classification does not bite: not because the coupling is of one type, but because the construction
admits one timescale.

**But the same measurement killed the identification §93 attached to the crossover.** The reflected
gap is **1.43–1.62**, flat in Ω — *not* the macroscopic rail rate |f'(r₃)| = 6.6195, and not
converging to it. **The wall installs a timescale of its own**: the slowest mode of a reflected
stage is diffusive relaxation across the whole box [saddle, cap], a domain-scale mode. The box is
Ω-independent in concentration (width 2.978), which is why the gap is flat and why the
discriminator still works. So at s_up = s_dn = 1 the upstream's correlation time is **0.6980** and
the downstream's rail relaxation time is **0.1511** — **4.62× apart**. See the §4 catalogue.

**And §100 priced the wall itself.** Freed, the chain has 12.2% more stage-2 error *at t₀ = 2.0* —
but that net hides two effects of opposite sign, each 33–45% of `P_refl`: the wall *removes* the
failure channel (−7.29e−3) and *inflates* the surviving branch by **+49%** (+5.32e−3), because
probability that would have leaked out is pushed back against the saddle where it drives stage 2
hardest.

**§100.2 then showed 12.2% is one cell of a trend.** Swept over the window, free/reflected runs
**1.9615 → 1.0673** for t₀ = 0.5 → 8: **at short windows the wall hides half the chain's stage-2
error**, because the failure channel is then almost the only route to a low stage 2, while at long
windows stage 2's own escape swamps it. The optimism in `D_max` is therefore a window-dependent
factor between 1.07 and 1.96 over the range measured — **not a constant that can be absorbed**.
Rule 18 again: sound measurement, unsound summary across cells, caught only by sweeping an axis the
prediction did not name.

> **SCOPE, restated for §92–§98.** What they measure is the transfer of the upstream's *rail
> fluctuations*, on a clock set partly by the wall. §91, unreflected, measured the other thing —
> accumulation of upstream *errors*. **`D_max = c*/(penalty × ε)` has been built from the reflected
> penalty alone and is therefore conditional on every upstream stage surviving.**

**T-CASC-n, open: what is the depth ceiling when upstream stages are allowed to fail?** §100's P3
was refuted in the direction that makes this urgent: `P(stage 2 low | stage 1 low)` is **0.7254 at
the shortest window measured and rises monotonically to 0.9830**, never falling. *(It is an endpoint
co-occurrence, not a transmission probability, and it must go to 1 as t₀ → ∞ carrying no information
— §100.2. The claim that survives is that there is no window in which the coupling buys any
protection.)* **The Hill coupling's saturation protects against fluctuations and not at all against
failures** — a
transfer function that flattens a 1σ wobble has no flat region once the input has moved to the other
rail, because transmitting that is what it is for. So a chain with free upstream stages has a second
error channel that transmits at ~94% per stage, and the noise-margin law never priced it. **How to
kill:** build D = 3 with no stage reflected, seeded from the joint stationary law, and compare the
per-stage error against the reflected chain's at the same Ω and t. If the free channel dominates,
the ceiling is set by upstream escape (a quantity §80–§90 can already compute exactly) rather than
by margin/σ, and §91–§98's law is the *sub*-dominant term. The instrument exists — §100's
`build_free` generalises to D = 3 at 121³ = 1.77M states, the size T-CASC-e already identified as
reachable. **The kill test must sweep t₀ rather than evaluate at one** (§100.2): the free channel's
weight relative to everything §92–§98 measured varies by 1.8× across the windows already tried, and
a single-window verdict would be the same error §100.2 corrects.

**T-CASC-l', open: which downstream timescale must the upstream beat?** §100 established the
crossover is *not* at the downstream rail relaxation time (4.62× off). It did not establish what it
is at. **How to kill:** the wall's box width is a free parameter — move the reflecting boundary
below the saddle and the upstream correlation time moves with it, at fixed landscape, fixed
coupling and fixed downstream. If the crossover tracks the *upstream* correlation time, the ratio
at the peak is a constant of the downstream alone and that constant names the timescale. This is
the same confound-breaking move as §100 itself: a knob that moves the nuisance against the
suspected cause.

~~**T-CASC-n, open: what is the depth ceiling when upstream stages are allowed to fail?**~~
**-> §101. ANSWERED for D = 3, and the answer moves the arc's centre of gravity.** With every stage
free, the failure channel is the **majority in every admitted cell at D = 3 — 55.7%, 61.3%, 68.8%** —
against a **minority in every one at D = 2 — 32.6%, 40.1%, 47.3%** (two barrier depths, two
windows). **One added stage moves the majority from the channel §91–§98 priced to the channel it
never saw, and it does so everywhere it was measured.** The exact share is window-dependent and
falls as the window lengthens (61.3% → 55.7% at D = 3 over t₀ = 2 → 8); the crossing does not move.

**And the reflecting construction is wrong in two directions at once**, which is why the damage was
invisible. At Ω = 30, D = 3 it reports 6.349 where the true total is 7.388 — 16% low, a minor-looking
correction — while **overstating the fluctuation-transfer channel it claims to isolate by 2.22×**
(2.64× at Ω = 14). The inflation grows with depth, 1.49× → 2.22×, and so does the omission; the
small net is a cancellation of two ~120–160% errors and there is no reason for it to persist further
along the chain.

> **SCOPE, restated again.** `D_max = c*/(penalty × ε)` is not merely conditional on upstream
> survival (§100). **The term it omits overtakes the term it keeps within one added stage.** What
> limits depth here is upstream *escape* — the quantity §80–§90 already computes exactly from the
> rate functions alone, with no master equation.

**T-CASC-o, open: is the depth ceiling just the upstream escape rate?** The suspect (rule 17, NOT
confirmed): at depth the cheapest route to a downstream error is one upstream escape followed by
near-certain propagation (§100 measured transmission at 0.73–0.98), rather than the last stage
escaping its own barrier. It predicts the contaminated share grows with the number of upstream
stages, which it does. **It is not confirmed and the naive version predicts the wrong Ω-dependence**,
since both routes cost one escape. **How to kill:** the suspect gives
`contaminated/pure ≈ (D−1) × [upstream escape rate] / [last-stage escape rate at its degraded
operating point]`. Both rates follow from §80–§90's escape action with no CME at all, and §96
already supplies the degraded operating point. A mismatch kills the account while leaving §101's
measurement standing.

**T-CASC-p, open: does the depth trend have a shape?** §101 has two depths. D = 4 at Ω = 14 is
10,556,001 joint states; it was started and abandoned after ~50 minutes without completing a single
window. **No `D_max` is extrapolated from two points** (§90 and rule 15 are why), so the *shape* of
the growth — and therefore any ceiling — is open. **How to kill:** the free chain is a product-form
candidate if the suspect above holds, in which case D = 4 is unnecessary and the ceiling follows in
closed form; if it does not hold, D = 4 is required and needs a solver that does not build the joint
generator explicitly.

**§101.1 is a methodological result and constrains §94 itself.** A cell is only admissible for a
depth comparison if P(last low) *increases* with D — a necessary condition, not a tolerance. At
t₀ = 0.5 it fails at both Ω: the reflected chain's D = 3 error is **ten times smaller** than its
D = 2 error, because stage 1 seeds from its quasi-stationary law while every later stage seeds as a
delta at its rail, so until that delta spreads a deeper chain hands its last stage a *cleaner* input.
Measured at Ω = 14, D = 3, t₀ = 0.5: stage 3 is (2.7450, 0.6846) against stage 2's (2.6542, 0.7643)
— the last stage cleaner than the one feeding it. **This is §92.1(a) one level up: there the
equilibration artifact reversed a trend in upstream speed, here it reverses the trend in depth,
which is the entire subject.** It removed §101's first headline (Ω = 30, D = 3, t₀ = 0.5: a
6.66× ratio at 91.8% contamination) before it was written down.

> **The usable window band is squeezed between equilibration below and saturation above**, and cells
> with P_free > 0.5 are excluded from both the ratio and the contaminated-share readings, since both
> are forced toward 1 there. At Ω = 14 only t₀ = 2.0 survives both bounds. **§94's own D = 3 result
> was taken at t₀ = 2.0 — admissible, but one step from the lower edge of the band**, and nothing in
> §94 records that the band exists.

~~**T-CASC-o, open: is the depth ceiling just the upstream escape rate?**~~ **-> §102. THE SUSPECT
SURVIVES ITS KILL TEST.** Reading each stage's escape rate off a **1-D** pinned generator at its
measured operating point, and predicting `contaminated ≈ P(any upstream failed) × p_transmit` and
`pure ≈ P(none failed) × (1 − e^{−k_last·t})`, the contaminated/pure ratio comes out at
**1.159 / 1.040 / 1.359** times the measured value across (Ω, D) = (14, 2), (14, 3), (30, 2) —
**all inside the pre-registered factor-of-two gate, and closest at D = 3**, the cell the question is
about. Wiring was exact: `escape_rate(Ω, r₃)` equals §100's free stage-1 gap to **0.00e+00**, because
`h(r₃) = 1` makes a pinned downstream stage literally the same element.

> **So the channel §92–§98 never measured is accounted for, to tens of percent, by escape rates —
> with no joint master equation on the prediction side.** What limits depth in this cascade is the
> escape action, which §80–§90 computes from the rate functions alone. §101 predicted this
> reconnection and §102 supplies it.

**§102's P4 was a broken criterion — the third this session.** It asked whether the measured
contaminated/pure tracks `k₁/k_last` across Ω, but the model's own formula sums over the upstream
stages, so the quantity depends on D while `k₁/k_last` does not. Measured, `k₁/k_last` is
0.9374 / 0.9066 / 0.9632 while contaminated/pure runs 0.8987 / 2.2001 / 0.6703 — **read literally,
P4 kills a suspect that P3 confirms.** Rule 19's own kill test would have caught it: the data that
makes P4 print "tracks" is a table at fixed D, and the table written spans depths. (§101's P1
demanded bitwise agreement between two formulations of the same arithmetic; §99.1 needed a post-hoc
label. The measurements in this session all survived; three of the criteria attached to them did
not.)

**T-CASC-q, open: is the frozen/fast bracket the whole of §102's residual?** P3's three misses are
all one-signed — measured `pure` always exceeds a rate-at-the-mean model — and `escape_rate(x_up)` is
steeply convex, so the rate at the average and the average of the rate differ. **Measured, the
effective rate is bracketed by those two limits in every cell**, at position **0.2200 / 0.0866 /
0.1783** in log space (0 = fast, 1 = frozen) — near the *fast* end. The limits separate by 2.44× at
Ω = 14 and **4.76× at Ω = 30**, which is why the residual is worst there. **That the value is
bracketed is measured; that the position means "the upstream is fast against the downstream's escape"
is a suspect (rule 17)** — and §100 showed the upstream's correlation time here is set partly by
where a boundary was put. **How to kill:** slow the upstream with §93's speed knob at fixed
landscape; the position must move toward 1. If it does not move, the bracket is a coincidence of two
numbers rather than an averaging limit.

**T-CASC-r, open: can the chain be predicted with no joint solve at all?** §102 takes the operating
points from §101's joint solve, so it shows escape rates *at* the operating points account for the
split — not that the operating points need no CME. §96 already predicts them from single-element
quantities to 0.12%. **How to kill:** feed §96's predicted operating points into §102's rate curve
and compare the resulting contaminated/pure against §101's measured values. If it survives the same
factor-of-two gate, **the depth ceiling of a chemically-coupled bistable cascade is computable from
one element's rate functions**, which is what the whole composition arc has been reaching for.

~~**T-CASC-r, open: can the chain be predicted with no joint solve at all?**~~ **-> §103. YES, and
the composition arc closes on single-element quantities.** Stage 1's law from a 1-D restricted
generator, each later mean from `⟨F(x)⟩ + d_intr`, each width from the LNA, every escape rate from a
1-D pinned generator — **no joint generator is built on the prediction side.** The contaminated/pure
split comes out at **1.184 / 1.059 / 1.370 / 1.169** times measured across (Ω, D) = (14,2), (14,3),
(30,2), (30,3), all inside §102's own pre-registered factor-of-two gate.

**And discarding the joint master equation costs about two percent.** §102, which read its operating
points off the joint CME, got 1.159 / 1.040 / 1.359 on the same cells; predicting them instead gives
1.184 / 1.059 / 1.370.

> **So the depth ceiling of a chemically-coupled bistable cascade is computable from one element's
> rate functions.** §101 showed the contaminated channel is the majority of the error by D = 3; §102
> showed that channel is set by escape rates; §80–§90 computes the escape action from the rate
> functions alone; §103 supplies the operating points those rates are evaluated at. The whole
> prediction side runs in under a second against a joint solve that costs ~25 minutes at Ω = 30,
> D = 3 and 1,771,561 states.

**One empirical input, and it cannot carry the result.** `p_transmit` is §100's measured 0.9376.
Swept over its entire measured range (0.7254–0.9830 across a 16× window, §100.2), **every cell stays
inside the gate** — 0.82 to 1.44. The model is not parameter-free and does not claim to be.

**T-CASC-s, open: where does the Ω-expansion break, and does the one-signed miss have a name?**
§103's P1 and P2 predicted the operating points within 1% and 2%; both **held at Ω = 30 (0.13%,
0.50%) and failed at Ω = 14 (1.64%, 3.45%, 4.32%)**, and every miss is in the same direction —
predicting a *less* degraded chain than the one measured. The QSD, the LNA width and the
static-transfer average are all asymptotic in Ω, and Ω = 14 is the A·Ω = 2.66 cell §98 already flagged
for reading outside its own bracket. **How to kill:** sweep Ω at fixed D and fit the residual's decay.
If it falls like 1/Ω the miss is the leading LNA correction and is nameable; if it saturates, one of
the three ingredients is wrong rather than merely truncated, and the conditioning choice (dropping
inputs where F has no high rail) is the first suspect since it biases the mean upward by exactly this
sign.

**T-CASC-t, open: can `p_transmit` be predicted too?** It is the last empirical number in the chain.
§100 measured it as an *endpoint co-occurrence* that rises 0.7254 → 0.9830 with the window and must
reach 1 as t → ∞, so it is not a constant of the coupling and a predicted version has to be a
function of the window. **How to kill:** it should follow from the same 1-D machinery — the
probability that a stage whose input has collapsed to the low rail has itself crossed its saddle by
time t is `1 − exp(−k(x_low)·t)` with `k` read off §102's curve at the *low* rail rather than the
high one. If that reproduces the measured window dependence, the model becomes parameter-free.

~~**T-CASC-t, open: can `p_transmit` be predicted too?**~~ **-> §104. PARTLY, and the model is now
parameter-free.** The window dependence follows from a one-parameter closed form —
`p_transmit(t) = 1 − (1 − e^{−k_low·t})/(k_low·t)`, stage 1 falling at a near-uniform time in [0, t]
and stage 2 then descending at rate `k_low` — with `k_low` computed from the same 1-D pinned
generator §102 used, evaluated at the **low** rail. Predicted against §100.2's five measured
windows, nothing fitted: **−8.95%, −7.51%, −3.01%, −1.22%, −0.58%**, one-signed and worst at the
shortest window exactly as predicted.

**The payoff (§104.1).** Feeding that predicted `p_transmit` into §103 in place of §100's measured
number leaves all four cells inside §102's gate: **1.1488 / 1.0266 / 1.3290 / 1.1337**.

> **So the depth ceiling of a chemically-coupled bistable cascade is computable from one element's
> rate functions, with no free parameter and no joint master equation anywhere.**

**Do not read the improvement as improvement.** Every cell is closer to measured with the predicted
`p_transmit` than with the measured one, because §104's value is ~3% low while §103's model runs
15–37% high — **two errors of opposite sign, partly cancelling**, which is the third time that
structure has appeared in this arc (§100's net, §101's net, now this).

**T-CASC-u, open: why is the descent rate ~30% too slow?** §104's P2 was **refuted**: the computed
`k_low` is 5.5169 against 7.07–8.77 implied by the windows, and against the measured curve the
disagreement is 0.6–9.0%. Two candidates were tested and neither suffices — an off-by-one in the
absorbing boundary, fixed in the *unflattering* direction (5.7660 → 5.5169), and starting the descent
from §103's degraded operating point rather than the pristine rail (**+3%**). **The remaining
suspect** is P1's second finding: the downstream is already monostable when the upstream merely
reaches its *saddle* (one root at 0.1351), so stage 2 begins descending **before** stage 1 formally
falls and the model starts its clock late. Right sign — **but a fixed delay is refuted too**: the
offset needed grows as 0.13, 0.78, 2.2 over t₀ = 0.5, 2, 8, i.e. roughly proportional to t, which is
a rescaling of the rate and not a shifted clock. **How to kill:** replace the indicator "upstream
below its saddle" with the continuous Hill factor and compute the descent under a *moving* landscape;
if that closes the 30% the account holds, and if it does not, the near-uniform-fall approximation
itself is the suspect rather than the rate.

**A criterion note, the fourth in this arc.** §104's P2 compared the computed rate against a rate
obtained by **inverting the model's own closed form on the data**. That tests form × rate jointly
and reports the verdict as a rate test; if the form carries a systematic error the implied interval
absorbs it. P3, against the measured values themselves, is the independent test — and it puts the
disagreement at 0.6–9.0% rather than 22–37%. Both describe the same shortfall; only one is a
measurement. (§99.1, §101's P1, §102's P4, now this.)

~~**T-CASC-u, open: why is the descent rate ~30% too slow?**~~ **-> §105. IT IS NOT. `k_low` was
never wrong; the clock was, and §104's refutation of that account was itself an algebra error
(rule 14).**

§104 modelled a head start as `form(k, t + Δ)` — a **longer window** — and concluded from the
resulting 18.90× spread in required Δ that "the offset grows in proportion to t, so it is a
rescaling of the rate and not a shifted clock". A head start is a different convolution:
`p(t) = 1 − e^{−kΔ}(1 − e^{−kt})/(kt)`. Under it the required Δ spans **2.17×** across the same 16×
window, and the early-start account was never refuted.

**And the head start turns out to be derivable.** The naive candidate fails twice — the conditional
traversal below the bistability edge is **τ = 0.4464**, 9.19× the fitted Δ and *longer than stage 2's
entire descent time* 1/k = 0.1813, which would put p_transmit at 0.9711 where the measurement says
0.7254. But x_up\* = 1.5795 is a **saddle-node**, so just below it the downstream descends 8.7×
slower than at the low rail, and time spent there cannot count at full rate. Rate-weighting the
conditional occupation time gives **Δ_eff = 0.0686 against a fitted 0.0486 — ratio 1.414**, inside
the pre-registered factor-of-two gate, with nothing taken from the measured p_transmit.

> **With that derived Δ the curve stays parameter-free and improves**: worst residual **8.95% →
> 5.80%**, the one-signed bias gone, and three of five windows inside 0.3%. §104's published Δ = 0
> numbers stand as printed (rule 7); §105 supplies the missing term without adding a parameter.

**§103's closure is untouched** — all four (Ω, D) cells stay inside §102's gate under §100's measured
p_transmit, under §104's Δ = 0, and under §105's derived Δ.

**T-CASC-v, open: is the remaining 5.8% at the shortest window the near-uniform-fall approximation?**
That is the one large residual left, and §104's own P3 predicted the shortest window would be worst
for exactly this reason: there `k₁t` is smallest but stage 2's own independent escape and stage 1's
chance of returning are least negligible *relative to the transmitted signal*. **How to kill:**
replace the uniform fall-time density with the exact conditional density `k₁e^{−k₁s}/(1 − e^{−k₁t})`
and add stage 2's independent-escape channel, both already computable from §102's curve. If the
5.80% closes, the account is complete; if it does not, the head start is absorbing an error that
belongs elsewhere.

~~**T-CASC-s, open: where does the Ω-expansion break, and does the one-signed miss have a name?**~~
**-> §106. IT IS NOT AN Ω-EXPANSION PROBLEM AT ALL, and the scope statement is now precise.**

Swept over Ω = 14–70 at D = 2, the operating-point residuals converge **super-algebraically** —
local exponents −2.53, −3.50, −4.47, −6.05, −8.49, **−11.43**, reaching 0.0001% at Ω = 70. §103's
transfer model is far better than the 1/Ω decay predicted. **But the closure ratio built from those
points runs AWAY from 1**, 1.1845 → 1.8158.

**Three candidates cleared, one convicted.** Not the operating points (exact), not the rate ratio
k₁/k₂ (→ 1), not `p_transmit` (measured flat to 1% across the sweep — and §105's derived-Δ model
tracks it to 2% *including its non-monotone shape*, an unplanned confirmation of a section built at
Ω = 30 alone). The convicted term is the crudest in the model: `P(stage 1 low) = 1 − e^{−k₁t}`, one
free stage, no coupling. **It spans 1.71× and crosses 1 near Ω = 35 — and §102 validated it at
Ω = 30, where it is accidentally accurate to 2.5%.** Rule 9, on the axis the whole arc is
parameterised by. Not a transient either: at Ω = 55 the error grows with the window (1.2453 → 1.4267
over t₀ = 2 → 8).

**The fix, and why it does not help.** A free stage has no absorbing boundary, so its low-occupancy
saturates at the stationary weight below the saddle, not at 1: `P(low at t) = π_low(1 − e^{−λt})`.
π_low falls 0.9057 → 0.5247 over the sweep, which **is** the missing 1.71×; with it the model is flat
to **2.5% across four decades of probability**. But π_low(stage 1) ≈ π_low(stage 2), so it
**cancels in contaminated/pure** — a correct fix to a real error that the reported quantity is blind
to.

**§106.3, an indexing error in §102/§103.** `escape_rate(Ω, x_up)` is the rate of a stage whose
*upstream* sits at x_up, but `predict` keys each stage to **its own** operating point. Stage 1 has no
upstream: its rate is the free gap, i.e. `escape_rate(Ω, r₃)` — which §102's own P1 established.
**The corrected indexing reproduces the measured gap exactly at all seven Ω**; the as-coded version
is 12.7% high at Ω = 14. Correcting it improves the closure everywhere (1.0416 → 0.9714 at Ω = 14)
without removing the divergence, so §102's and §103's published agreements were slightly lucky and
their conclusions survive.

> **SCOPE, now precise.** The transfer model is superb, the escape model was wrong in a way that
> cancels, and what actually bounds §103's closure is the **averaging of the escape rate over a
> fluctuating input** — §102.1's rate-at-the-mean approximation, whose two limits separate 2.44× at
> Ω = 14 and 8.80× at Ω = 55. With the two-state form in place the whole residual sits in `pure`,
> under-predicted ~2× at Ω = 70. **T-CASC-q is no longer a residual; it is the binding constraint.**

**T-CASC-w, open: does §103's closure survive once the escape rate is properly averaged?** §102.1
showed the effective rate is bracketed by `k(⟨x⟩)` and `⟨k(x)⟩` at position 0.09–0.22, and §106 shows
the bracket widens with Ω. **How to kill:** replace `k(⟨x⟩)` in `split_from` with the position-
interpolated rate and re-run §106's sweep. If the closure ratio then flattens toward 1, the arc's
model is complete and its only approximation is a named, measured interpolation; if it does not, the
position itself is Ω-dependent in a way §102.1's three cells could not see — which is the same
mistake §106 just caught, one level down.

~~**T-CASC-w, open: does §103's closure survive once the escape rate is properly averaged?**~~
**-> §107. THE LEVEL IS FIXED; THE DRIFT IS NOT.**

Four corrections applied — §106.3's indexing, §106.2's two-state occupancy, §102.1's averaging, and
§105's derived head start. **The correct fast-limit average is the GEOMETRIC mean of the rate**:
with `A(x) = −ln k(x)/Ω`, `exp(−Ω⟨A⟩) = exp(⟨ln k⟩)`. §102.1 used `k(⟨x⟩)`, which is a different and
smaller number; the geometric mean reproduces §102.1's measured effective rate to **±9% at Ω = 14,
30 and 55** where `k(⟨x⟩)` is ~21% low and the frozen limit is 2–7× off.

**The ablation lands where §106 said it would.** At Ω = 30 the ratio runs 1.3702 → 1.2879
(indexing) → 1.2575 (two-state) → **1.0138** (geometric) → 1.0142 (derived p_t). The geometric
correction carries almost all of it — it is the only one that touches `pure`, where §106 localised
the entire residual — and two-state is nearly invisible, as predicted, because π_low cancels.

**But the drift survives every correction.** The model is now *centred*, crossing 1 near Ω = 30
instead of sitting above 1 everywhere, and 6/7 cells land in [0.7, 1.4] against 4/7. **It still
trends monotonically and the span gets slightly worse: 1.533× → 1.880× over Ω = 14–70.**

> **WHAT SURVIVES, precisely.** Every cell across Ω = 14–70 is inside §102's pre-registered
> factor-of-two gate under the corrected model (0.78 to 1.47), so §103's closure holds across the
> whole swept range at the criterion it was registered against. **What does not hold is any claim of
> asymptotic correctness**: the model is calibrated near Ω ≈ 30 and drifts 1.88× across a 5× change
> in Ω.

**T-CASC-x, open: what drives the residual Ω-drift?** Every identified approximation has now been
corrected or shown to cancel, and the drift is unchanged. **The last remaining candidate has the
wrong sign**: §102.1's position moves *toward* the fast end as Ω grows (0.2200 → 0.1072 over
Ω = 14 → 55), so a fast-limit average should get *more* accurate at large Ω, while the measured ratio
gets worse (1.01 at Ω = 30, 1.47 at Ω = 70). **How to kill:** the drift is monotone and smooth over
1.88× — fit it against the candidate scalings that are actually available (A·Ω, √Ω, the rail width
σ = √(V/Ω), π_low itself) and see whether any reproduces it. If none does, the missing term is not a
finite-size correction to anything already in the model, and the honest scope statement is the
calibrated one. **Report every candidate, not the flattering one (rule 15).**

**T-CASC-x, open -> §108: NOT CLOSED, but reduced from "no candidate" to "two mechanisms, one of
them provably outside the family of fixes tried."** The goal was to explain and remove the drift or
prove it irreducible; neither was achieved, and this records how far it got.

**The drift is entirely in one channel.** With §107's corrections applied, `contam` measured/model is
**flat to 0.4%** across Ω = 14–70 (1.1988 … 1.2151) while `pure` runs **0.9474 → 1.7604**. The one
thing stage 2 has that stage 1 does not is a fluctuating input.

**Six candidates cleared**, each by measurement: the operating points (exact to 0.0001%), the
predicted input law and the averaging of λ (both 0.4% against the measured law), splitting π_low from
λ (1.3% → 0.4%), `p_transmit` (flat), and **return trips** — pre-registered, and refuted. `pure`
conditions on stage 1 being high *at the end*, so dip-and-return trajectories count as pure while
having dragged stage 2 down, and π_low falling 0.906 → 0.525 makes returns go from rare to common:
right sign, right Ω-dependence. Making stage 1 **absorbing** below its saddle removes them outright.
**The drift survives, 0.5926 → 1.2142, span 2.05×.**

**The result: one mechanism cannot be enough, and this is an inequality rather than a fit.**
Inverting stage 2's effective rate from the absorbing run gives true/k(⟨x⟩) = **0.7880, 0.9599,
1.0784, 1.1380, 1.2459, 1.4138** over Ω = 14–70. And `ln k` is **convex in the input at every Ω**
(d²/dx² ∈ [+0.408, +2.230] at Ω = 14, [+2.779, +14.493] at Ω = 70).

> **By Jensen, convexity of `ln k` puts both `exp⟨ln k⟩` and `⟨k⟩` at or above `k(⟨x⟩)`. No average
> of the rate over the input can land below `k(⟨x⟩)`. The Ω = 14 and Ω = 20 cells are below it.**
> Those cells are outside the reach of *any* averaging prescription — not the ones tried, the whole
> family, by an inequality.

**So there are at least two mechanisms**: one suppressing, dominant at small Ω and provably not an
averaging effect, and the averaging enhancement, dominant at large Ω. **That is why §107's four
corrections recentred without flattening** — every one addresses the enhancing mechanism and none can
reach below `k(⟨x⟩)`. The whole class of fixes §102.1–§107 were drawn from is retired.

**T-CASC-y, open: what suppresses stage 2's escape at small Ω?** It must be something that does *not*
act on stage 1, since `contam` is flat. Candidates worth a kill test, in order: (i) the delta seed —
stage 2 starts at its rail and has to spread, which suppresses early flux, though the deficit is a
fixed absolute time and so should be Ω-independent as a fraction of the window; (ii) correlation
between stage 2's position and its input, which every averaging treatment here assumes away by
evaluating a rate against a *static* input law; (iii) the finite response time of stage 2 to input
excursions — the input's deep excursions are exactly the ones that are brief, which is §92's motional
narrowing acting on the escape rather than on the penalty. **(ii) and (iii) are the same suspect
stated two ways and neither is an average over a static law, which is what §108 proves is needed.**

**T-CASC-z, open: the flat 21.5%.** `contam`'s measured/model sits at 1.215 at every Ω — the model
under-predicts stage 1's escape by that much with an *exact* gap, an *exact* stationary law and no
coupling at all. It does not drift, so §108 set it aside, but it is unexplained, and a transient from
the delta-seeded start would have the wrong sign.

~~**T-CASC-z, open: the flat 21.5%.**~~ **-> §109.2. CLOSED. It was a seed mismatch.** Stage 1 is
seeded from `stage1_stationary`, the **reflected** stage's stationary law, while the model assumes
the **free** stage's QSD. The reflected law is not depleted near the saddle, so it escapes
**1.2460 / 1.2770 / 1.2768 / 1.2493 / 1.2402** times faster across Ω = 14–70 — flat, which is why
the residual was flat. **With matched seeding the two-state model predicts a free stage's escape to
3–5% at every Ω, spanning 1.0188×.**

**T-CASC-y(i) CONFIRMED, and §108's dismissal of it was reasoned from the wrong seed.** §108 argued
a delta-seed transient "would have the wrong sign". True of stage 1's reflected law; **false of stage
2's delta**, which is the other mismatch and runs the other way. In a *single* stage with the
upstream **pinned** — no fluctuating input, nothing to average — the delta seed carries a **25.2×**
transient over a 16× window against the QSD seed's **1.12×**, and at t = 2 it suppresses by 1.68×.

> **§108.3's central claim is WITHDRAWN (see §4).** Its Jensen argument was sound; its premise —
> that the measured rate lies below `k(⟨x⟩)` — holds only for a delta-seeded chain. With matched
> seeding the sequence is 0.9893, 1.2825, 1.5190, 1.8161, 2.0257, and the one sub-unity cell is 1.1%
> below while the QSD instrument itself under-reads by a measured **5.3%** at t = 2. Transient
> corrected: **1.0445, 1.3540, 1.6038, 1.9175, 2.1387.** No cell is materially below `k(⟨x⟩)`; the
> averaging family is **not** excluded.

**But the drift survives both seed fixes** — span **2.048×** with matched seeding against the
delta-seeded 1.794×. §106–§108's *localisation* stands; only §108's explanation falls.

**T-CASC-x, restated more sharply than §108 left it.** With matched seeding the effective rate runs
from **at the fast limit** at Ω = 14 (0.9893 against a bracket top of 2.3937) to **past the frozen
limit** at Ω = 70 (2.0257 against 1.8346). The system traverses the entire fast-to-frozen bracket as
Ω grows, and exits it. **The naive timescale argument predicts the opposite**: the escape time grows
exponentially in Ω while the upstream's correlation time is macroscopic, so the upstream should
become *relatively faster* and the position should move *toward* the fast limit. **How to kill:** the
position is now a measured curve over five Ω; compute the upstream correlation time and the
downstream instanton traversal time independently (both 1-D, both available from §100 and §105's
h-transform) and check whether their ratio moves the way the position does. If it moves the other
way, the position is not set by that ratio and the resonant-activation framing inherited from
§92/§93 does not apply to the escape rate at all.

**A correction to this arc's own running commentary (§109.4).** §102–§109 each carried a banner
counting "broken criteria" and listing §99.1, §101's P1, §102's P4, §104's P2, §104's refutation and
§109's P3. **The count is right; the category is not.** Three are rule-19 gates that cannot be
satisfied by what they test (§101's P1, §102's P4, §104's P2). **One** is a rule-20 fixed tolerance
against a larger instrument residual (§109's P3 — and it is 12× the gate, not the 5× first quoted,
because that 5% was measured at Ω = 30 and extrapolated onto Ω = 14). **One is plain arithmetic**
(§104's refutation: `form(k, t+Δ)` for a head start instead of the convolution — verified against
quadrature, §105's form agrees to 0.00e+00 and §104's is off 2.56e−02). **One was never a broken
gate at all** (§99.1, a rule-2 disclosure that was correctly disclosed). Grouping them made one
recurring fault out of four different ones, and pointed at the wrong remedy for three of them.

**§109.5 -- the four documented-but-unfixed defects are now closed in code.** §106.3's indexing,
§106.2's one-way occupancy and §102.1's rate-at-the-mean were written up and demonstrated in separate
functions while `predict` and `split_from` kept carrying all three as their default. They now
**default to the corrected model**, with `legacy=True` for the original; every caller reproducing a
published table (§102's, §103's and §106's `main()`, and the tests pinning their values) passes it
**explicitly, naming the section it reproduces**, so rule 7 is honoured *visibly* rather than by a
silent default. The seed mismatch is handled the other way round on purpose: the delta-at-the-rail
seed is §101's intended initial condition and both its arms share it, so the default stands, but
`solve()` now carries a SEEDING WARNING with the measured cost and `matched_seed=True` supplies the
QSD seed any stationary-model comparison needs. No published number changes.
