# What restoration costs

*A synthesis of CRNL §1–§70. Every claim here points at the section that measured it and
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

What that leaves open is η itself: its landscape dependence is measured, not derived. WKB should
give it in closed form as the barrier action per molecule (T-DEPTH-e).

---

## 3. Why reliability has no thermodynamic price

The natural hope is an exact relation between error rate and dissipation. It was pursued four
ways and closed each time — the last closure being the sharpest.

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
* **The barrier exponent ν**, where three extraction routes disagree between 1.95 and 2.19
  (§64, T15-n).
* **Cascades.** Everything above is one element. The founding claim is about *composition*.
