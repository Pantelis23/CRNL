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
