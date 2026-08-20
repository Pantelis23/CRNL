# What restoration costs

*A synthesis of CRNL §1–§105. Every claim here points at the section that measured it and
carries the scope that section stated. Where a claim was withdrawn, the withdrawal is here
too — the retraction record is part of the result, not an appendix to it.*

---

## The question

Binary computation won. The usual explanation is that the transistor is a near-ideal
**restoring** switch: it pushes a marginal input back toward 0 or 1, so error does not
accumulate as you compose gates. That is a story about physics told by engineers, and it has
never been priced.

**Is signal restoration a measurable physical quantity, and does it have a thermodynamic
cost?**

CRNL asks this in chemistry, where the accounting is exact. A reaction network is run two
ways — the deterministic mass-action ODE and the exact chemical master equation — and the gap
between them *is* the noise the restorer must fight. Entropy production is computed per jump,
exactly, with no coarse-graining.

The workhorse is **Approximate Majority** (`X+Y→2B`, `B+X→2X`, `B+Y→2Y`, each reversible at
γ× the forward rate): a two-symbol restoring element that stops restoring at γ_c = 1/2. Later
sections add **Schlögl's model** (`2X⇌3X`, `∅⇌X`, chemostatted), which restores with no
symmetry at all and is exactly solvable.

---

## The answer in two sentences

> **A restoring switch is a tuned object, not a topological accident.** Whether a network
> *can* restore is decidable from its stoichiometry alone; whether it *does* is a condition on
> its rate constants.
>
> **Restoration has a thermodynamic cost. Reliability does not have a thermodynamic price.**
> Error rates are bought with molecules, not with free energy — and this is not a limitation of
> the measurement but a theorem-shaped obstruction.

The second sentence is the one that took sixty sections and three withdrawn attempts to reach,
and it is a better answer than the constant everyone expects. The transistor's 10⁻¹⁵ error rate
is bought with electron count, exactly as AM's is bought with Ω.

---

## 1. What a restoring element *is*

For any mass-action network symmetric under exchanging two species X and Y, the drift of the
lead `n_X − n_Y` factors exactly:

    b_X − b_Y = (n_X − n_Y) · P(n)

so the diagonal is invariant: **the sign of a lead is a deterministic invariant, and every
reversal is a fluctuation** (§43). P is the amplification. It decomposes as

    P(x) = Σ_r  c_r · d_r · B_r(x),    d_r = S_X(r) − S_Y(r) ∈ ℤ,    B_r(x) ≥ 0

with every bracket nonnegative and `d_r` fixed by stoichiometry alone (§54). Two consequences:

* **Capability is combinatorial** — the network restores for *some* rates iff some `d_r > 0`;
  the non-restoring rate vectors form a convex cone (§56).
* **Realisation is one linear inequality** — it restores iff `Σ_r c_r d_r B_r(x*) > 0` at the
  symmetric steady state, reproducing AM's γ_c = 1/2 to 1.1e−15 and agreeing with the ODE on
  120/120 networks (§62). The realising set is a cone but is **not convex**: two restoring rate
  vectors can sum to a non-restoring one, inside AM's own family (§62.2).

At finite molecule count the deterministic sign change becomes a **metastability transition**
with a width, going as **Ω^(−1/2)** (§63) — *four times the molecules for half the blur on the
threshold.*

**Prior art, checked late and honestly (§65.1, §70).** The invariance is a folk theorem of
equivariant dynamics. The realisation criterion is Theorem 1 of Montoya, Cruz & Ágreda (2019)
in the homochirality literature, where it is stated more generally than here. What was not
found there is the **sign decomposition** and the combinatorial capability criterion — and
"not found in two papers" is weaker than "novel". This project claims priority for nothing.

---

## 2. The four currencies, and which ones survived

Restoration is not bought with one thing. Four distinct costs were isolated, and the last two
sections destroyed the universality of two of them.

| you buy | you pay | status |
|---|---|---|
| **existence** of a landscape | a nonzero drive — an affinity floor | **real everywhere, universal nowhere** |
| **gain** (e-folds of amplification) | ~2 k_B per molecule per e-fold | **AM-specific bookkeeping** |
| **reliability** (error rate) | **molecules**, as exp(−cΩ) | **survives; the only cross-substrate claim** |
| **speed** at fixed distance from the TUR bound | time | untested off AM |
| **composition depth** | **molecules at the rail** (§75) | same currency as reliability |

**Existence.** Every restoring element examined has an affinity floor: below some drive, no
bistable landscape exists at all. AM's is exactly `3 ln 2 = 2.0794` (§9.1); Schlögl's is
exactly `2 ln 3 = 2.1972` (§67). They sit 5.66% apart, both `ln(small integer)`, from chemistry
with nothing in common — and that agreement is a **coincidence of two points**. Generalising
Schlögl's autocatalysis to order p at fixed pair count gives `A_c(p) = 2 ln[(p+1)/(p−1)]`,
which runs from 2.1972 down to 0.5026 and tends to zero (§68). *There is no universal price of
admission.* The floor is a property of the element.

**Gain.** §38 measured ~2 k_B per molecule per e-fold for AM, in closed form. On Schlögl that
quantity **does not exist**: the element is chemostatted, so it dissipates while sitting still,
and Σ-to-absorption is dominated by housekeeping proportional to the decision *time*. Two
repairs were tried and both refuted — a housekeeping subtraction goes negative in every cell,
and a dimensionless cycles-per-molecule measure grows with Ω (§67.2). *The price of gain is
defined relative to closed conservative bookkeeping, and does not survive the move to a driven
device.*

**Reliability.** The error probability falls as exp(−cΩ) (§1, §35, §38). No free-energy cost
attaches to it, and §3 below says why that is structural.

**Composition depth — and it is not a currency at all.** §12's ceiling
`D_max ≈ exp(δ*²/2σ²)/4` is set by the inter-stage channel against the rails. §72 measured it
on Schlögl, found the ratios landing on AM's, and read that as the one quantity that transfers.
**§73 deflated it**: a *step function* — commit to whichever rail is nearer, no chemistry, no
dynamics — reproduces the same ratios, and a slope-matched sigmoid matches Schlögl to three
significant figures.

> The corrected statement is sharper. **Composition depth is fixed by the readout geometry, and
> the restoring element's entire contribution is where it puts its rails.** Reaction order,
> symmetry, drive, chemostatting — none of it survives into D_max once Δ is fixed. So a
> restoring element, for the purposes of depth, is a *rail-placer*.
>
> And the search for a substrate-independent price is over: three costs were tested across
> substrates and none transferred (§67, §68); the fourth was not a cost. **No cost of
> restoration measured here is substrate-independent, and the one substrate-independent
> quantity is not a cost.**

**§74 says why, and it is the closest thing here to an answer.** If depth is bought with Δ,
then the cost of depth is the cost of Δ — and that depends on a structural property nobody had
isolated:

> **A conservative element has a maximum composition depth.** AM's δ\* ≤ 1 because
> X + Y + B = Ω, so at fixed channel noise no amount of drive pushes it past D_max = 9.5×10⁹.
> Going from 9 to 28 k_BT of affinity per cycle buys a factor of 8; going to infinity buys
> another 1.04.
>
> **An open element has none.** Schlögl's affinity `ln[e₁e₂/e₃]` is *exactly* invariant under
> rescaling its rails (measured drift 0.00e+00 over 2.5 decades), so at the same affinity it
> passes 10¹⁸ — paying in **material**, not free energy.

**Closed and open elements pay different currencies for the same quantity.** That is why no
substrate-independent price was ever found: the earlier sections were pricing dissipation, gain
and affinity floors, all of which are upstream of Δ and invisible to depth.

**§75 settled the convention, and §74 did not survive it.** In a chemically coupled cascade σ
is not a choice: stage 1's output species *is* stage 2's input, so the inter-stage noise is
stage 1's own rail fluctuation. It is **Poissonian in the molecule count**, σ_n ~ √n̄, so
Δ/σ ~ √(λΩ) — and measured, **Δ/σ depends on the product λΩ alone, to 0.000% spread over 32× in
each.** Scaling an element's rails and increasing its molecule count are the *same act*. AM
behaves identically (σ ~ Ω^−0.50), so its "maximum composition depth" was a fact about the
convention, not about conservation.

> **So there is one currency, not four.** Depth: `D_max ~ exp(c·N_rail)`. Reliability:
> `p_err ~ exp(−cΩ)`. **Same currency, same functional form.** A restoring element buys both
> its error rate and its composition depth with molecules, exponentially, and buys neither with
> free energy.
>
> Affinity is a **gate, not a price**: drive is required for rails to exist at all (§9.1, §68)
> and it sets where they sit — but once they exist, more depth is a matter of counting. §74's
> own numbers showed it without my seeing it: A = 9 → 28 k_BT bought a factor of 8 in depth,
> and A → ∞ bought another 1.04.

**§76 closes it: the two are not merely the same currency, they are the same number.** For a
binary channel with per-stage error ε, the depth at which the bit dies is exactly

> **D_max = c\*/ε,  c\* = 0.124266404564**

verified to 8.9e−9 over seven decades and identical across every substrate measured, to under
1%. So the cascade apparatus — kernels, mutual-information decay, the ceiling formula — carries
no information beyond ε.

> **A restoring element is characterised, for every purpose this project set out to measure, by
> one number: its per-stage error probability.** Reliability *is* ε. Depth *is* c\*/ε. There is
> no second quantity, no trade-off between them, and no thermodynamic price on either — ε is
> bought with molecules, exponentially.
>
> This is a **deflationary** answer to the founding question, and it is worth stating as such:
> the depth advantage of a restoring switch is not an additional property beyond its error
> rate. It is the same property counted twice.

**§77 names the one number that is left.** If reliability and depth are both ε, the whole
question is how fast ln(1/ε) grows with molecule count — `η = d ln(1/ε)/dΩ`, **nats of
reliability per molecule.** It exists: ln(1/ε) becomes linear in Ω on both substrates, settling
to 0.01–0.20%. And it does **not** transfer — spanning a factor of **633**, and moving with the
landscape inside each substrate (γ from 0.30 to 0.05 buys 14×).

> **Everything about a restoring element washes out of composition — reaction orders, symmetry,
> drive, dissipation, the shape of its commitment function, even whether it is conserved —
> except η, and η is the element's own.**
>
> A restoring element is a device for converting **molecules into nats** at an exchange rate set
> by its landscape, and nothing else about it matters to what it can compute. The transistor's
> advantage is not a special thermodynamic property. **It is a good exchange rate.**

**§78 derives η, and it lands somewhere unexpected.** Not WKB — that gives the *escape*
probability, a different exponent. The ε here is the Gaussian readout of the rail's own
fluctuation, so the relevant object is the **linear-noise variance**: σ² = V/Ω from the Lyapunov
equation, giving

> **η = Δ²/(2V)**

Against §77's measured values with nothing fitted: 9.8797 vs 9.8813, 1.8329 vs 1.8346, 0.6781 vs
0.6813 (AM), 0.015577 vs 0.015617 and 0.059537 vs 0.059574 (Schlögl). **Worst 0.48%**, largest
exactly where the rail is shallowest and the LNA weakest.

> **Δ, the Jacobian and D = S·diag(a)·Sᵀ are all deterministic-plus-LNA quantities.** No master
> equation, no stationary solve, no simulation, no entropy production.
>
> **So the whole founding question — how deeply can you compose a noisy restoring element, and
> what does it cost — is answered by the ODE and its linear-noise correction.** The exact CME
> was needed to *establish* that; it is not needed to *use* it.

That also explains, retrospectively, why every thermodynamic price failed to transfer: **entropy
production is not in the formula.** Affinity enters only by setting where the rails sit and how
deep the wells are — through Δ and V — and never in its own right.

**§79 tested that out of sample, which is the only test that counts.** Predicting from the ODE
alone and then checking against the exact CME, on systems the derivation never saw: AM at two
unused γ (0.01%, 1.12%), Schlögl with **quartic autocatalysis, 3X⇌4X — a different reaction
order** (0.14%), and Schlögl with asymmetric rails (0.04%). All converge. And a case chosen *in
advance to fail* — AM at γ = 0.45, where the rail is shallowest — duly failed, not converging at
all and landing 20.9× worse.

> **One qualification, and it is real.** `D_max ~ exp(Δ²Ω/2V)`, so a relative error δ in V
> becomes a factor `D_max^δ`. Sub-percent σ errors become depth factors of 1.3–23×. **The ODE
> route predicts the exponent — η and ln D_max — to well under a percent; it does not predict
> D_max itself to better than an order of magnitude.** For that, the master equation is still
> required.

### The correction that matters most

**§80 found that §75–§79 priced the wrong failure mode.** An element in a cascade can misread
its input *or* spontaneously escape its rail during the stage time. Both are exponential in Ω,
so the smaller exponent wins outright — and on §75's own element the escape exponent is
**A = 0.026046 against the readout's η = 0.059537**. Escape is more likely by **e⁴³²** at
Ω = 12800; readout would bind only for stage times below 10⁻¹⁸⁸.

> **§75's own premise — that a physical cascade is chemically coupled — selects the regime in
> which §75–§79's ε is not the physical one.**

What survives: §76's `D_max = c*/ε` is mathematics and doesn't care where ε comes from. §77's
"one exponent governs everything" is the right shape with the wrong number. **And §78's headline
survives with a different formula** — the escape action is also a deterministic-side quantity,
`A = −∫ln(μ/λ)dx` from the rate functions alone, matching the exact first-passage time to 1e−4.
So the founding question is still answered without a master equation; it is the quasipotential
integral rather than the linear-noise variance.

### The ODE does not determine reliability

**§83 is the sharpest single result of this arc, and it constrains the project's own method.**
A depends on λ and μ separately; the mass-action ODE depends only on their difference. So a
reaction pair adding the same function to both — a birth and a death with identical *reactant
complexes*, like `X→2X` with `X→∅`, or AM's `X+Y→2X` with `X+Y→2Y` — leaves the deterministic
dynamics **exactly** unchanged and still moves A.

> Gated on both substrates: the drift is invariant to **4.9e−15** (Schlögl) and **1.1e−16** (AM).
> Across the sweep A falls 5.6× and η falls 4.4×, monotonically, by two unrelated routes. **At
> Ω = 6400 the depth ceiling spans e¹³⁷ between elements no ODE measurement can tell apart.**

§78/§81's operational claim survives — A is still quadrature, no master equation needed. Its
gloss did not: *"a property of the deterministic field"* was too loose. **A is a functional of the
propensity pair, and the ODE is strictly coarser than the thing that sets reliability.**

> ⚠ **And almost none of that is new — checked early, and the check came back badly (§85.1).**
> The construction is published as **"zero-drift networks"** (Plesa, Zygalakis, Anderson & Erban
> 2018), with §83's `X+Y→2X` / `X+Y→2Y` pair as their own worked example. Networks with identical
> mass-action ODEs are **dynamically equivalent** / **confoundable** (Horn–Jackson 1972;
> Craciun–Pantea 2008), and that they differ stochastically is established (Enciso–Erban–Kim 2021;
> Faul–Hoessly–Xia, whose Example 5.12 is §83's one-species case). The action `−∫ln(μ/λ)dx` used
> since §80 is Assaf–Meerson (2017) Eqs. (43)–(44) and should have been cited there.
> **What is left is the quantification, not the phenomenon** — and per §70, "not found" is not
> novelty.

Both neutral pairs are irreversible, so this sharpens §82 rather than qualifying it: the element
pays an *unbounded* thermodynamic cost and gets **worse** reliability for it. Drive and
reliability move in opposite directions inside one network.

§71/§72's external-channel cascade is unaffected — checked, not assumed: there the readout term
binds by 35 nats. **Engineered wiring and chemical coupling are physically different cascades.**

**§81 redid all three sections against A, and the answers hold in their new form.** The integral
predicts the exact first-passage time to **0.06% on four landscapes it never saw**, one of them a
different reaction order — so §78's claim (the coefficient is available from the deterministic
side alone) is now established for the exponent that actually governs depth. A does not transfer
either: a factor of 13 across shapes and substrates. And **A/η spans a factor of 3.78**, which is
what makes §80 a correction rather than a relabelling — a constant ratio would have carried
§77–§79's conclusions over intact.

---

## 3. Why reliability has no thermodynamic price

The natural hope is an exact relation between error rate and dissipation. It was pursued five
ways and closed each time — and the fifth is the one that counts, because the first four priced
quantities §80 later showed govern nothing.

**§82 asks it of the escape action, the exponent that does govern reliability and depth.** Both
networks have a *one-dimensional* cycle space, so their entire non-equilibrium drive is a single
number that can be pinned by construction while the kinetics move:

> With the affinity held fixed to **4.9e−14 nats**, A spans a factor of **926**. With A held
> fixed to **0.0000%**, the affinity spans **0.92 nats**.

Neither implication holds. The drive does not set the action, and the action does not report the
drive — and unlike the earlier closures, this one runs in both directions. What *does* constrain
A is the landscape, though not tightly: at fixed Δ, across a 32.7× move in the saddle's position,
A still varies by 1.75×. **Both fail, by margins three orders of magnitude apart.**

The four earlier closures, each correct about its own observable:

The integral fluctuation theorem `⟨e^(−S_tot)⟩ = 1` holds at the decision time — best cell
5.5e−14, **median 1.33e−9** over 36 cells (§41; quoting only the best cell would be the
flattering half of a spread this project's rule 15 exists to prevent). Split it by outcome, and the hope is that the odds of being right appear as the
exponentiated entropy of the error paths. **They do not: the identity factorises.**
`⟨e^(−S_tot)|o⟩ = 1` for *each* outcome separately, so the split carries no information about
the error rate (§60).

That closure was first proved only where symmetry guarantees it (§66 caught this), and §60's
stated mechanism — the two boundaries being exchange images of equal stationary weight — is
available only at zero tilt. So it was re-asked on an element with **one species, where no
symmetry exists to appeal to**:

> `Φ_o = p_o` to **2e−13**, across skews from 1:1 to 4:1, with the two boundaries' stationary
> weights differing by up to **e³²** (§69).

**The result generalises; §60's mechanism for it is wrong.** What produces the identity has no
candidate explanation, and gets none until something independent supplies one.

> **The obstruction is not a symmetry accident.** The error rate is not obtainable from the
> entropy production by this route, on either substrate family.

Scope: single elements, two outcomes, one separatrix. Cascades and n > 2 outcomes are
untouched.

---

## 4. Why this should matter outside chemistry

* **Fault-tolerant computing.** §7's cascade and §12.1's depth ceiling
  `D_max ≈ exp(δ*²/2σ²)/4` are von Neumann's (1956) and Evans–Schulman's noisy-circuit depth
  bounds — *with the gate thermodynamically priced*, which that literature has never had. (The
  measured ceiling ran ~3× the prediction, later read as a 7% error in the exponent rather
  than a missing prefactor, §12.1.)
* **Quantum error correction.** §32–§34: chemistry has a knob QEC lacks, and the bounded-window
  result says where the analogy stops.
* **Statistical mechanics.** P is a mean-field order parameter, γ_c a critical point, §63 its
  finite-size rounding. Spin systems have no analogue of *which reaction* carries the
  amplification — `d_r` is precisely that.
* **Thermodynamic uncertainty relations.** §40's first-passage Q has floor 1; AM sits at
  **5.475 at a mean decision time of 4.09**, and is **Pareto-optimal at its own speed** (§58).
  Networks reaching Q = 1.25–1.33 — a 4.3× improvement — do it *by being slower*, at a mean
  decision time of **783**. Roughly 190× the time to get 4.3× closer to the bound.

---

## 5. What this project is actually for

Results have been withdrawn after publication repeatedly, each with the mechanism that killed
it left visible (the catalogue is `THEORIES.md` §4; `FINDINGS.md` carries ~46 withdrawal
markers). Several were killed by the *next* section — §64 withdrew §63.2's exponent, §68
refuted §67's affinity-floor pattern one commit after opening it, and §69 overturned §66's
reading of its own data. The failure modes were stable enough to become rules:

* a **mechanism** proposed in the same breath as the measurement it explains is a hypothesis,
  and three in two sessions died to the first test aimed at them;
* a **misread statistic** is more dangerous than a wrong mechanism, because it is a plausible
  sentence with correct numbers attached;
* a **broken verdict criterion** is worse still — it prints a confident verdict *from correct
  numbers*, and nothing in the output looks wrong. Seven of these were caught, four in a single
  session, which is why every verdict rule is now unit-tested on data engineered to trigger
  each branch **before** the experiment runs. That convention has caught a defect in four
  consecutive sections, at a cost of about ten minutes each.

The measurements in this project survive. The explanations attached to them mostly have not.

---

## 6. What is open

* **Clause 2's status** — not found in two papers is not novelty (T-THM-a, §70).
* **Whether the theorem has an asymmetric analogue at all.** Clause 1 fails outright without
  exchange symmetry (§42: residuals 1.9e1, 2.9e2), and the founding object is asymmetric.
* **Why the fluctuation theorem factorises over outcomes.** Measured twice, explained never.
* **The WKB prefactor** (T14-e): localised to a start-side factor, exactly θ-independent, a
  function of (γ, ε) — §61 gave it a falsifiable target that analysis must now reproduce.
* ~~**The barrier exponent ν**, where three extraction routes disagree between 1.95 and 2.19~~
  — **closed by §84: ν = 2 exactly.** Eliminating the fast variable gives the escape action in
  closed form with no fitted parameter, `A = (1+2γ)(1−2γ)²/(2(1−γ)(1+γ)²)`, which matches the
  exact first-passage action to **0.9776 and rising** at γ = 0.46. The three routes were never
  disagreeing about a value: ν is a *limit*, every finite window reads below it, and the measured
  effective exponent (1.942) is exactly §64's width route. §63.1's width exponent of exactly 1/2
  follows. **T15-n.1 closed by §85**: the residual is the adiabatic elimination — truncation and
  the diffusion approximation are both refuted, and sharpening the timescale separation drives the
  residual 0.0977 → 0.0130.
  > ⚠ **But almost none of §84 is new (§85.2).** ν = 2 at a supercritical pitchfork is a standard
  > universal exponent (Dykman et al. 1998), measured at 2.00 ± 0.03 (Chan & Stambaugh 2007), and
  > it **does not discriminate** — transcritical gives 2 too. The "every window reads below 2"
  > account is the named effective-exponent/corrections-to-scaling phenomenon (Wegner 1972;
  > Dykman & Ryvkine 2004) — and it is algebra inside our own formula: `A = (32/9)ε²(1 − (5/3)ε)`,
  > so `ν_eff = 2 − (5/3)ε`. What survives is the closed form for reversible AM and the absolute
  > validation against exact CME first passage.
  **§86–§87 then located and fixed the reduction's error.** The deficit is the adiabatic
  elimination (§85), and the tilted-generator route predicts the escape curve from the rate
  functions alone — **1.0257 and 1.0151 against the exact first passage, with no master equation,
  no stationary solve and no Ω** (§87), where the deterministic slow manifold gives 0.9027/0.9348.
  §86's ridge measurement stands; §86's *mechanism* was withdrawn by §87, since the instanton's
  curve is displaced opposite to the ridge.
  **§89 closes the arc**: the remaining overshoot is the fast-pair term `∫p_s ds` the reduction
  drops — computed, not fitted, it accounts for 91–98% of it at every γ and reproduces the
  non-1/M shape §88 could not. **The escape action of a two-species restoring element is
  computable to ~0.1% from its rate functions alone**, validated absolutely against exact CME
  first passage over eight γ and ten (γ, M) cells. What began as "ν ≈ 2 ± 0.1, not determined
  more precisely by any instrument here" (§64) ends as a quadrature.
  > ⚠ **§90 withdraws §89's quoted precision.** The measured A is itself only determined to ~0.4%
  > once its extrapolation ansatz is varied, and the candidates straddle the theory value. §89's
  > 0.1% residual is *below the resolution of the instrument that would test it*. The agreement
  > stands — to within the measurement's own resolution — but "1.0007 ± 0.0004" does not, and the
  > flatness that made it look physical was a common-mode estimator bias.
* ~~**Cascades.** Everything above is one element.~~ **§91 built the chemically-coupled chain
  the founding claim is actually about, and the answer is that composition is not a property of the
  element at all.** Three couplings, all *exactly* neutral at the rail, using the same element and
  rails: one does not transmit, one transmits with a noise margin of 0.88σ and amplifies per-stage
  error 18–64×, one saturates with a 3.39σ margin and costs only 1.7–2.3×.
  > **The figure of merit is the noise margin in units of the upstream rail width, and the penalty
  > is exponential in it — `log(penalty) = −0.95 × margin/σ`, with two independent knobs collapsing
  > onto one curve.** That is what a transistor's saturating gain curve supplies, and what §12/§71/
  > §72's external channel imposed by hand as σ = f·Δ instead of deriving.
  §76–§81's depths, computed from an isolated element's ε, are right only in the large-margin
  limit.
  **§92 then found a second control the margin does not contain.** Scaling the upstream stage's
  rates leaves its landscape, barrier, rail width and stationary law identical and changes only its
  clock — and the penalty moves by 2.7×, from a slow plateau of 4.24–4.44 (matching the
  frozen-upstream average `⟨exp(−ΔA·Ω)⟩ = 4.845` to 10%, unfitted) down toward the mean-landscape
  rate 1.139. **That fall is motional narrowing**: upstream fluctuations faster than the
  downstream's response average out before it can act on them.
  > **So composition is governed by two numbers — the noise margin and the timescale ratio — and a
  > cascade whose upstream stages run fast is protected by narrowing.** §91's law is the
  > frozen-upstream case; its variants all sat at speed 1, at the top of the curve.
  **§93 placed the crossover: it sits at τ_up = τ_down**, where the upstream correlation time meets
  the downstream response time — a plateau below it (matching the frozen formula to 10%, unfitted)
  and a fall to 1.58 above. *(The collapse onto the ratio is an algebraic identity,
  `Q·t = t₀[(τ ratio)·Q₁ + Q₂]`, and is labelled a wiring check rather than a result.)*
  > **So the composition penalty is a function of two numbers — the noise margin and the timescale
  > ratio — and both are properties of single elements.** The joint master equation validates that;
  > it is not needed to compute it. **A cascade is safest when its upstream stages are fast and its
  > transfer function saturates** — motional narrowing and noise margin, neither of them a property
  > of the restoring element.
  > ⚠ **§99(a)–(b): the bracket and the crossover both have prior art, and one reading is a
  > suspect.** The frozen/fast bracket is the standard heuristic account of **resonant activation**
  > (Pechukas & Hänggi, PRL 73, 2772, 1994), already done for master equations with a WKB action
  > and finite-correlation-time noise by Assaf et al. (PRL 111, 058102, 2013). And there are
  > **three** regimes, not two: `⟨exp(−ΔA·Ω)⟩` is the *intermediate* one, not the quenched limit.
  > Separately, "τ_up = τ_down" is a suspect (rule 17) — the literature's optimum matches an
  > *escape* time for barrier-top modulation and a *relaxation* time only for "breathing"
  > potentials, and §93 never asked which the Hill coupling produces. **The crossover is measured
  > and stands; its account is open (T-CASC-l).** What has no prior art is that the modulator here
  > is a chemically-coupled upstream stage rather than an imposed process — which is exactly why
  > the two timescales cannot be tuned independently.
  > ⚠ **§100 settles T-CASC-l and withdraws §93's identification.** The reflecting wall at stage 1's
  > saddle removes its escape channel, so **only one of the two timescales exists in this chain** —
  > confirmed on the axis §93 lacked: over Ω = 14→55 the reflected gap varies **1.13×** while the
  > free gap varies **1708×**, tracking exp(−A·Ω). But the reflected gap is **1.43–1.62, not** the
  > rail rate |f′(r₃)| = 6.6195, and is not converging to it: **the wall installs a box-scale
  > timescale of its own, 4.6× slower.** So at s_up = s_dn = 1 the upstream's correlation time is
  > 0.6980 and the downstream's rail relaxation time 0.1511 — **4.62× apart. The crossover sits at
  > matched *speeds*, not matched *times*.** The crossover and the interior maximum stand; the
  > identification is withdrawn, and *which* downstream timescale must be beaten is open again.
  **§94 took it to D = 3.** The law *composes*: fed stage 2's measured width it predicts stage 3's
  penalty as 5.65 against a measured 6.35 — 11% low, inside §91's own scatter. But the width
  sequence it needs is not predicted by anything: unconditioned widths run 0.499, 0.542, **0.708**,
  blowing past the LNA fixed point (0.551), and the deterministic gain at the operating points is
  **13× too small** to explain even the first step.
  > ⚠ **§99(c) withdraws that last clause and dates the recursion.** `σ²_out = σ²_intr + g²σ²_in`,
  > its fixed point and the g < 1 condition are **Thattai & van Oudenaarden (2002), Eq. 13** — 24
  > years old. Their form carries a **time-averaging factor** τ = β_x/(β_s+β_x) that mine lacks,
  > valid only under separated timescales, which §92–§93 measured to be *equal* here. And the
  > "too small by 13×" gap came from solving g² out of σ₂/σ₁ **assuming σ_intr is
  > stage-independent** — when §95–§96 show the operating point moves and the intrinsic width
  > follows it. §96's account reproduces the widths with **no transmitted-noise term at all**, and
  > is what survives.
  **§95–§96 then closed the width question.** The mean shift is the *static-transfer average of a
  concave map*: computed exactly rather than to second order it predicts the operating point to
  **0.12%**, with the residue bracketed by the frozen and fast limits as the upstream clock is
  swept 256×. And the width follows from the operating point within the LNA's own accuracy.
  > **So the chain is predicted end to end** — mean, then width, then penalty — from single-element
  > quantities.
  **§97 ran it forward on an element it was never built from**, and split the result cleanly: the
  **mean is right to −0.00%** and the width to 2.76%, while **§91's fitted slope is 74% wrong**.
  The derivation behind that slope — §92's frozen/fast average — still brackets the measurement
  (8.92 / 4.47 / 1.01), so **the margin law was a one-element parameterisation, and the barrier
  depth A·Ω enters too**. §91's 14-point collapse was real but every point shared one element:
  rule 9 at the top level.
  **§98 then closed it.** Ω moves the barrier depth without touching the landscape, so the curve
  could be traced on one element and the other checked against it: the position between the two
  limits falls monotonically with A·Ω, and §97's element — different landscape, different coupling,
  barrier 2.3× deeper — lands on that curve to **6.6%**, nothing fitted.
  > **So the composition penalty is computable from single-element quantities after all**, with
  > A·Ω as the variable §91's margin law could not see. What began as "the chain was never built"
  > ends with mean, width and penalty all predicted out of sample. *(Two instrument bugs were found and corrected on the way, one of which had
  > reversed §94's headline; the originals stand in FINDINGS §94 with the correction beside them.)*
  > ⚠ **§99.1 tests the obvious objection and does not fully settle it.** Near a saddle-node
  > A·Ω ∝ (margin/σ)², and §98's Ω-sweep **cannot separate them at all** — within one element
  > margin/σ ∝ √Ω. Only the second element breaks the tie, where the ratio differs by 22%: A·Ω
  > transfers to −6.6% and (margin/σ)² to −15.0%. That is **one out-of-sample point with no error
  > bar, run post-hoc**, and it is recorded as a discrimination rather than a prediction. A third
  > element settles it (T-CASC-m).

**§100 then priced the reflecting boundary that §92–§98 all rest on**, and found the arc's scope
narrower than its sentences. Freed, the chain carries 12.2% more stage-2 error at t₀ = 2.0 — but that
net hides two effects of opposite sign, each 33–45% of the reflected error: the wall *removes* the
upstream failure channel and *inflates* the surviving branch by **+49%**, pushing back probability
that would have leaked out. **And the 12.2% is one cell of a trend**: swept over the window,
free/reflected runs **1.9615 → 1.0673** for t₀ = 0.5 → 8, so at the short windows where a cascade
actually has to hold its value **the wall hides half the error**.
  > **So §92–§98 measure the transfer of the upstream's *rail fluctuations*, not the accumulation of
  > its *errors*, and `D_max = c*/(penalty × ε)` is conditional on every upstream stage surviving.**
  And the prediction that mattered was refuted: **`P(stage 2 low | stage 1 low)` is 0.7254 at the
  shortest window and rises monotonically to 0.9830, never falling.** *(An endpoint co-occurrence,
  not a transmission probability — it must reach 1 eventually. What survives is that no window buys
  the downstream any protection.)* **Saturation protects
  against fluctuations and not at all against failures** — a transfer function that flattens a 1σ
  wobble has no flat region once the input has reached the other rail, because transmitting that is
  what it is built to do. The noise-margin law prices the first channel and never saw the second,
  which is now the arc's largest open question (T-CASC-n).

**§101 took the free-upstream chain to depth 3, and the arc's centre of gravity moved.** With every
stage able to fail, the failure channel is the **majority in every admitted cell at D = 3 — 55.7%,
61.3%, 68.8% across two barrier depths and two windows — against a minority in every one at D = 2**
(32.6%, 40.1%, 47.3%). **One added stage moves the majority from the channel
§91–§98 priced to the channel it never saw.** And the reflecting construction turns out to be wrong
in two directions at once, which is why the damage stayed invisible: it *understates the total* by
16–21% while **overstating the very channel it claims to isolate by 2.2–2.6×**. §94's 6.349 against
a true 7.388 reads as a minor correction; the same number is more than double the fluctuation
transfer it is supposed to measure. Both errors grow with depth and nothing makes the cancellation
persist.
  > **So `D_max = c*/(penalty × ε)` is not merely conditional on upstream survival (§100) — the term
  > it omits overtakes the term it keeps within one added stage.** What limits depth here is upstream
  > *escape*, the quantity §80–§90 already computes exactly from the rate functions alone. That is a
  > reconnection, not a demolition: the arc's own earlier machinery is what the ceiling needs.

  **§101.1 is the methodological half, and it deleted the headline before it was written.** The first
  reading was a 6.66× ratio at 91.8% contamination (Ω = 30, D = 3, t₀ = 0.5). The tell was that the
  *reflected* chain's D = 3 error was **ten times smaller** than its D = 2 error — a longer chain
  cannot be more reliable in steady state. Stage 1 seeds from its quasi-stationary law and is already
  spread; every later stage seeds as a delta at its rail, so until that delta propagates a deeper
  chain hands its last stage a cleaner input. **§92.1(a) one level up: there the equilibration
  artifact reversed a trend in upstream speed, here it reverses the trend in depth, which is the
  entire subject.** The gate that follows — admit a cell only if P(last low) increases with D — is a
  necessary condition rather than a tolerance, and it applies to §94's cells too, which nothing in
  §94 records.

**§102 then closed the loop §101 opened.** §101's suspect — that the failure channel is upstream
escape followed by near-certain propagation — was given its kill test, and survived it. Reading each
stage's escape rate off a **1-D pinned generator** at its operating point predicts the
contaminated/pure split to **1.16 / 1.04 / 1.36** of the measured value, inside a pre-registered
factor-of-two gate and closest at D = 3. Wiring was exact: pinned at the rail, a downstream stage
*is* stage 1, agreeing to 0.00e+00.
  > **So what limits depth in this cascade is the escape action — which §80–§90 already computes from
  > the rate functions alone, with no master equation.** The composition arc's ceiling turns out to be
  > set by the arc's own earlier machinery rather than by the noise-margin law that replaced it.

  And the residual is one-signed with a named cause: the escape rate is steeply convex in the
  operating point, so the rate at the average and the average of the rate differ — **§92's frozen/fast
  pair, arriving at a rate instead of a penalty.** The measurement is bracketed by both limits in
  every cell, at position 0.09–0.22, near the fast end; the limits separate 2.44× at Ω = 14 and 4.76×
  at Ω = 30, which is exactly where the residual is worst. *(That it is bracketed is measured; what
  the position means is a suspect with a kill test — T-CASC-q.)*

**§103 closed the arc.** §102 still read its operating points off the joint CME, so it had not shown
the chain is computable without one. §103 predicts them instead — stage 1's law from a 1-D restricted
generator, each later mean from the exact static-transfer average plus the stage's own finite-Ω
depression, each width from the LNA, every escape rate from a 1-D pinned generator — and **builds no
joint generator at all.** The contaminated/pure split lands at **1.184 / 1.059 / 1.370 / 1.169** times
measured across four (Ω, D) cells, all inside §102's own gate. **Discarding the joint master equation
costs about two percent.**
  > **So the depth ceiling of a chemically-coupled bistable cascade is computable from one element's
  > rate functions.** §101: the contaminated channel is the majority of the error by D = 3. §102: that
  > channel is set by escape rates. §80–§90: the escape action follows from the rate functions alone.
  > §103: so do the operating points those rates are evaluated at. The prediction side runs in under a
  > second against a joint solve costing ~25 minutes and 1.77M states.

  Two honest limits. `p_transmit` remains empirical — but swept over its whole measured range
  (0.7254–0.9830) every cell stays inside the gate, so it does not carry the result. And the scheme is
  an **expansion in Ω**: P1 and P2 held at Ω = 30 (0.13%, 0.50%) and failed at Ω = 14 (up to 4.32%),
  every miss one-signed toward a less degraded chain. §103 found where it frays.

  **§103.1 is the one that would have been embarrassing.** The first version composed the operating
  point from the transfer map alone and predicted the chain getting *better* with depth —
  2.8501 → 3.0946 → 3.1173 against a measured sequence that falls. The cause was structural and
  visible: §91 built the coupling neutral at the rail, so **`F(r₃) = r₃` exactly and the rail is a
  fixed point of the transfer map** — iterating it can only climb toward the rail, never degrade. What
  degrades a chain is each stage's own finite-Ω depression below its rail. **Caught by a sign, not a
  magnitude**; no tolerance would have flagged it.

**§104 removed the last free parameter.** `p_transmit` was §103's one empirical input. Its window
dependence follows from a one-parameter closed form whose rate is the *descent* rate of a stage whose
input has collapsed — the same 1-D pinned generator, read at the low rail instead of the high one.
Predicted against §100.2's five windows with nothing fitted: −8.95% to −0.58%, one-signed and worst
at the shortest window as predicted. Fed back into §103, all four cells stay inside the gate at
**1.1488 / 1.0266 / 1.3290 / 1.1337**.
  > **So the depth ceiling of a chemically-coupled bistable cascade is computable from one element's
  > rate functions — no free parameter, no joint master equation.**

  Two things not to over-read. The cells got *closer* with the predicted `p_transmit` than with the
  measured one, but that is two errors of opposite sign partly cancelling — the third time that
  structure has appeared in this arc. And the rate behind the form is **~30% too slow** (5.5169
  against 7.07–8.77 implied): an off-by-one was found and fixed in the unflattering direction, the
  degraded-start correction is worth 3%, and the natural suspect — that the downstream is already
  monostable when the upstream merely reaches its saddle, so the clock starts late — is refuted in
  its simplest form, because the offset needed grows in proportion to t rather than staying fixed.
  T-CASC-u carries it.

**§105 corrected §104 — and the thing that needed correcting was a *withdrawal*.** §104 named the
right mechanism for its 30% shortfall (stage 2 starts descending before stage 1 formally crosses)
and then refuted it, because the required offset appeared to grow in proportion to the window.
**It grew because §104 modelled a head start as a longer window.** A head start is a different
convolution, `p(t) = 1 − e^{−kΔ}(1 − e^{−kt})/(kt)`; under it the required Δ spans 2.17× rather than
18.90×, and the account was never refuted.
  > **`k_low` was never wrong. The clock was.**

  And Δ is derivable. The naive candidate fails twice — the conditional traversal below the
  bistability edge is 0.4464, nine times too big and *longer than stage 2's whole descent*, which
  would put p_transmit at 0.9711 where the measurement says 0.7254. But x_up\* = 1.5795 is a
  **saddle-node**, so just below it the downstream descends 8.7× slower and that time cannot count
  at full rate. Rate-weighting the conditional occupation time gives **0.0686 against a fitted
  0.0486**, with nothing taken from the measured p_transmit. The curve stays parameter-free and
  improves: **worst residual 8.95% → 5.80%, the one-signed bias gone, three of five windows inside
  0.3%.**

  **One thing this session's record makes plain.** Every measurement in §98–§102 survived. Three of
  the *criteria* attached to them did not: §99.1 needed a post-hoc label, §101's P1 demanded bitwise
  agreement between two formulations of the same arithmetic, and §102's P4 was ill-posed across
  depths and would have killed a suspect its own section confirms. That is rule 19's class, three
  times in five sections, and in each case the numbers were right and the sentence about them was not.

**What §99 leaves standing.** The composition arc's *measurements* survive intact; its
*mechanisms* took the damage, which is now the project's most reliable pattern (rule 17, and
§99(a)–(c) are the fourth, fifth and sixth instances). The arc's genuinely unclaimed ground is
narrower than it looked and sharper: **a directed cascade of bistable elements solved by exact CME,
where the modulator is itself chemistry rather than an imposed noise process, and where the
operating point is taken from the exact static-transfer average before the LNA is applied.** No
reference was found combining a directed cascade, genuine bistability, an exact CME and a depth
question.
