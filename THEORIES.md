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
correcting one ingredient made agreement worse. §22.1 is untouched.

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

**T10b-iii, open: how does D_fuel scale with the budget?** Measured sub-linear
(11 -> 23 for a 3x budget, ~Phi^0.67) where the crossover argument assumed linear.
Two fuel concentrations cannot pin an exponent, and with the ceilings compounding
there is no reason to expect a clean power law at all. **How to kill:** four or more
budgets at a fixed quiet channel; if the exponent drifts with sigma, D_fuel is not a
budget property and the whole two-ceiling framing needs rewording.

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
