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

> **Read with §5.1.** The construction has an exact reduction the module notices
> and does not use: under the internal clock `dτ = e^{−Ht}dt` this is *ordinary*
> AM stopped at internal time `τ = 1/H`. The "finite integrated propensity" is a
> finite time budget. Everything below is true at each finite Ω, but the critical
> rate is `1/(consensus time)` and therefore drifts to **0** as Ω grows.

**There is a critical expansion rate at any given Ω.** Below it consensus completes
and the system
reaches a clean rail; above it the reaction freezes mid-decision, locking in a relic
minority abundance that grows toward 0.5. This is the chemical analogue of
cosmological freeze-out (the Γ-vs-H competition that set the relic dark-matter
abundance and primordial helium) — shared mathematics, not a claim about cosmology.

At the transition the frozen composition is maximally broad (critical fluctuations);
below it the distribution is sharply bimodal at the rails, above it a single spike
at the initial 50/50.

---

## 5. Freeze-out sharpens with Ω — `freezeout_scaling.py`

> **CORRECTED by §5.1, and the heading used to read "is a genuine transition".**
> The collapse below is a real fit; its reading is wrong. `Hc = 0`, there is no
> critical point, and `a ≈ 0.38` is a parameter of the wrong functional form.
> Read §5 as the measurement and §5.1 as what it means.

Six system sizes spanning ×32 collapse onto a single master curve under

    D(H, Ω) = F((H − Hc)·Ω^a)

with **Hc ≈ 0.055** and **a ≈ 0.38** (`results/freezeout_fss.json`). Finite-size
crossings drift toward Hc as predicted:

| Ω | 40 | 80 | 160 | 320 | 640 | 1280 |
|---|----|----|-----|-----|-----|------|
| H*(D=0.5) | 0.1658 | 0.1409 | 0.1219 | 0.1078 | 0.0972 | 0.0879 |

A crossover would not collapse; this does. **The conclusion drawn from that here —
"finite-size scaling about a critical point, and in the Ω→∞ limit expand slower than
Hc and you always decide" — is wrong, and §5.1 replaces it.** A collapse tests
whether one coordinate makes the curves interleave; it does not test whether the
coordinate has the right functional form, and here it does not.

**The exponent was flagged as not-to-be-over-read, and that was the right instinct
for the wrong reason.** a ≈ 0.38 comes from a two-parameter grid search with no error
bars, and sits between 1/3 and 2/5. The universality class is not identified because
there is no universality class. "Connecting it to the quasipotential of §2" was
indeed the next step — done in §5.1, and what the quasipotential says is that the
whole power law is spurious.

### 5.1 There is no critical rate: Hc = 0 — `freezeout_law.py`

**`Hc = 0`.** From a symmetric start there is no critical expansion rate; in the
Ω→∞ limit *every* positive H freezes the decision, which is the opposite of what §5
concluded. Four routes agree and the first is exact rather than fitted.

**Route 1 — the expanding SSA is ordinary SSA with a finite clock.** `expanding.py`
notices that a purely bimolecular network's total propensity decays as
`a₀(t) = a₀(n)·e^{−Ht}`, and stops there. Put in the internal clock
`τ(t) = (1−e^{−Ht})/H`, so `dτ = e^{−Ht}dt`. The next-event condition becomes

    −ln u = ∫₀^Δt a₀(n)·e^{−H(t+s)}ds = a₀(n)·[e^{−Ht} − e^{−H(t+Δt)}]/H = a₀(n)·Δτ

which is *exactly* the ordinary Gillespie increment `Δτ = −ln(u)/a₀(n)`; reaction
choice is untouched by the overall 1/Ω scaling. And `τ(∞) = 1/H`. So

> **The expanding SSA at rate H is ordinary AM stopped at internal time τ = 1/H.**

Freeze-out is not a separate dynamics — it is the ordinary chain failing to finish
inside a finite budget, and "the remaining integrated propensity never reaches the
target" is that budget running out. Verified **bit-for-bit**: identical seeds give
identical state sequences and identical frozen/absorbed verdicts, 0 mismatches over
5 × 60 seeds at H ∈ {0.05, 0.12, 0.3, 1, 3}, and it holds for n-winner too (the
argument needs only uniform reaction order). Two consequences: `H*(Ω) = 1/τ*(Ω)`
with `τ*` the plain consensus time, and **one SSA pass measures every H at once** —
which is what makes the range below affordable.

**Route 2 — the consensus time diverges logarithmically.** `design.md` §9 already
holds every ingredient. `dδ/dt = δ·b` holds **exactly**, not to linear order, so
with `b → 1/3` the decision axis grows at `λ = 1/3`; and the effective spread at
the saddle is `σ² = D_δ/λ = 1/(3Ω)`. From an *exactly* symmetric start that
`Ω^{−1/2}` shot noise is the only seed, so

    τ*(Ω) = (1/λ)·ln(1/σ) + O(1) = (3/2)·ln Ω + O(1)

— slope **3/2, no fitted parameter** — and `H*(Ω) = 1/((3/2)lnΩ + B) → 0`.

The Ω→∞ statement is stronger than a fit. At any *fixed* internal time the density
process converges (Kurtz) to the mass-action ODE, which from an exactly symmetric
start stays exactly symmetric forever, so `D(H,Ω) → 0` for **every** H > 0. A
positive Hc is not unsupported; it is impossible.

**Measured, 15 sizes over ×16384 in Ω** (16 replicates × 1250 trials each;
`results/freezeout_law.json`):

| Ω | 40 | 160 | 640 | 2560 | 10240 | 40960 | 163840 | 655360 |
|---|---|---|---|---|---|---|---|---|
| τ* = 1/H* | 5.999 | 8.208 | 10.363 | 12.469 | 14.532 | 16.635 | 18.717 | **20.772** |
| H* | 0.1667 | 0.1218 | 0.0965 | 0.0802 | 0.0688 | 0.0601 | 0.0534 | **0.0481** |

`H*` passes straight through §5's `Hc ≈ 0.055` and keeps going. Fits of
`τ* = A·lnΩ + B`:

| range | A (D = 0.5) | vs 3/2 | χ²/dof |
|---|---|---|---|
| all 15, Ω ≥ 40 | 1.5110 ± 0.0016 | +6.7σ | 5.44 |
| Ω ≥ 640 (11 pts) | **1.5005 ± 0.0023** | **+0.2σ** | **1.10** |
| Ω ≥ 2560 (9 pts) | 1.4999 ± 0.0035 | −0.0σ | 1.34 |

The full-range fit is 0.7% high because the local slope drifts *down* to 3/2 from
above; drop the small-Ω transient and the parameter-free prediction is confirmed to
**0.03%**. The 14 local slopes `dτ*/dlnΩ` run 1.621, 1.566, 1.551, 1.558, 1.522,
1.517, then scatter about 3/2 with no trend — mean of the top 8 = **1.497 ± 0.020**,
where a positive Hc requires them to be heading for **0**. All four crossing levels
agree (A = 1.497 / 1.511 / 1.515 / 1.519 at D = 0.25 / 0.5 / 0.75 / 0.9), as they
must if the whole distribution is translating.

**Route 3 — the exact CME, no sampling at all.** Integrating `dp/dτ = Qᵀp` on the
conserved simplex gives τ* with zero sampling error. Against the SSA:

| Ω | 40 | 80 | 160 | 320 | 640 | 1280 |
|---|---|---|---|---|---|---|
| exact τ* | 5.9934 | 7.1185 | 8.2144 | 9.2887 | 10.3488 | 11.4002 |
| SSA τ* | 5.999 | 7.123 | 8.208 | 9.283 | 10.363 | 11.417 |
| diff | −0.006 | −0.004 | +0.007 | +0.006 | −0.014 | −0.017 |

Agreement to **0.15%** between sparse linear algebra and an SSA vectorised across
trials, sharing no propensity or selection code. (Ω ≤ 1280: the exact route costs
~Ω³ and takes 235 s at 1280 against 0.15 s at 40.)

**Route 4 — the ODE, no stochastics at all.** Displace the mass-action ODE from
`(1/2, 1/2, 0)` by `δ₀ = Ω^{−1/2}` and read off when it decides. Local slopes per
factor 4 in Ω converge to **1.5000**, and the loser-clearing time to **2.5001**
(see below). This route cannot predict the seed *amplitude* — only its `Ω^{−1/2}`
scaling — so it pins the limits, not the size of the finite-Ω excess.

**Five things that kill Hc > 0 without fitting a functional form.**

1. **D at a fixed H, as Ω grows.** Under FSS, `D(Hc, Ω) = F(0)` must be
   Ω-*independent*. At H = 0.055 it is not:

   | Ω | 40 | 640 | 2560 | 10240 | 40960 | 163840 | 655360 |
   |---|---|---|---|---|---|---|---|
   | D at H = 0.055 | 0.988 | 0.946 | 0.893 | 0.801 | 0.644 | 0.449 | **0.268** |
   | D at H = 0.04 | 0.999 | 0.995 | 0.989 | 0.978 | 0.955 | 0.908 | 0.832 |

   A 3.7× fall where a critical point demands a constant, and still accelerating.
   (H = 0.04 is falling too but has not crossed, so it is not by itself a
   refutation — the log law only forces D → 0 once `(3/2)lnΩ > 1/H`, which at
   H = 0.04 needs Ω ~ 10⁷. Stated because quoting the H = 0.04 column as
   evidence would be overreach.)

2. **§5's own fit, extrapolated.** Refitting `Hc + C·Ω^{−a}` on Ω ≤ 1280 reproduces
   §5 (`Hc = 0.0580, a = 0.371`) and then predicts `τ*(655360) = 16.42`. Measured
   **20.77** — 26% out. On the full range the same form drifts to `Hc = 0.0364`,
   `a = 0.233` with χ²/dof = **85.7** against the log law's **5.4**, using one more
   parameter. Pin `Hc = 0.055` and χ²/dof = **4535**.

3. **The curvature is 21× too small.** A positive Hc forces τ* to bend over toward
   the ceiling `1/Hc`. The measured quadratic term in lnΩ is
   **−0.00405 ± 0.00065** — a bend of 0.09 against a rise of 14.8. §5's own
   power-law form requires **−0.0845**, i.e. 20.9× more.

4. **The transition width in `1/H` is constant.** The log law translates the curve
   without reshaping it, so the width is fixed; FSS with `Hc > 0` cannot manage
   that, because its width in H falls like `Ω^{−a}` while `H*` saturates at `Hc`:

   | Ω | 40 | 320 | 1280 | 5120 | 20480 | 81920 | 655360 |
   |---|---|---|---|---|---|---|---|
   | measured τ*(.75) − τ*(.25) | 5.28 | 5.66 | 5.73 | 5.70 | 5.69 | 5.69 | **5.70** |
   | required by §5's own params | 5.28 | 5.77 | 5.20 | 4.16 | 3.03 | 2.06 | **1.06** |

   Flat to 0.7% from Ω = 320 onward, where FSS demands a **5.4× collapse**.
   **This is the check that was worthless over §5's range** — there the two forms
   agree to ~1% (see below) — and it is decisive here only because the range is 64×
   longer.

5. **The "exponent" is not constant, and the log law predicts its drift with no
   free parameter.** For *any* pure power law `−dlnH*/dlnΩ` is a constant; the log
   law says it equals `A·H*`, hence drifts to 0:

   | Ω pair | 40→80 | 320→640 | 2560→5120 | 20480→40960 | 327680→655360 |
   |---|---|---|---|---|---|
   | measured −dlnH*/dlnΩ | 0.248 | 0.159 | 0.113 | 0.098 | **0.074** |
   | log law, `(3/2)·H*` | 0.229 | 0.153 | 0.116 | 0.093 | **0.074** |

   A 3.4× drift, tracked to a few percent at all 14 pairs by a parameter-free
   expression. Since the exponent is not constant, there is no exponent.

**A zero-parameter collapse beats the two-parameter one.** Same
Bhattacharjee–Seno residual, same points, `D` sampled at `τ = 1/H` on a common grid:

| collapse coordinate | free params | residual |
|---|---|---|
| `τ − (3/2)·lnΩ` | **0** | **2.09e-05** |
| `τ − A·lnΩ`, A fitted | 1 | 1.21e-05 (A = 1.517) |
| `(H − Hc)·Ω^a` (§5's) | 2 | 5.94e-04 |

28× worse with two more parameters — and, refit on the wider range, §5's own form
**chooses `Hc = 0`** (it returns Hc = −0.0000, a = 0.141).

### The control that settles it: give the system a bias

Same instrument, same network, same observable — start with an Ω-*independent*
pairwise margin `δ₀ = 0.10` instead of an exactly symmetric one. Then the seed does
not shrink, so `τ*` should be Ω-independent and `H*` should tend to a positive
constant: a real critical rate, with no finite-size scaling at all.

| Ω | 320 | 1280 | 5120 | 20480 | 81920 | slope dτ*/dlnΩ |
|---|---|---|---|---|---|---|
| symmetric start, τ* | 9.283 | 11.417 | 13.485 | 15.545 | 17.621 | **+1.5005 ± 0.0023** |
| δ₀ = 0.10, τ* | 4.829 | 4.767 | 4.762 | 4.758 | **4.758** | **−0.0022 ± 0.0003** |

Flat to 0.15% over ×256 in Ω, giving `H* = 0.2102` — and the power-law fit, handed
this data, returns `Hc = 0.2101` with no drift. `D` at H = 0.055 reads **1.0000 at
every Ω**, which is precisely the Ω-independence FSS predicts and the symmetric
start refuses. The collapse ranking inverts too, as a real discriminator should.

**So the object §5 mistook for a critical point is the shrinking initial
condition.** `H*(Ω)` drifts because the shot-noise seed does.

### A second, independent number: absorption costs (5/2)·lnΩ, and where it misses

Full absorption needs the symmetry to break *and* the last molecules to be cleared
off the rail. Near the rail `dY/dt = (Y/Ω)(B − X) ≈ −Y`: unit rate, so clearing
from `O(Ω)` to one molecule costs its own `lnΩ`. Prediction: **5/2 = 3/2 + 1**.

Measured (median absorption time, symmetric start): global fit 2.636 ± 0.002,
falling to **2.589 ± 0.005** over Ω ≥ 2560 and 2.550 at the top pair. The ODE route
gives exactly 2.5001. So 5/2 is approached but **not** reached — 3.5% high and
13σ from 5/2 on the restricted fit.

The biased control locates the discrepancy exactly. With the symmetry-breaking term
removed by a fixed δ₀, the absorption slope measures the *clearing* term alone:

    clearing coefficient  = 1.0895 ± 0.0176        (prediction 1)
    breaking coefficient  = 1.5005 ± 0.0023        (prediction 3/2)
    sum                   = 2.590                  measured together: 2.589 ± 0.033

**The decomposition closes to 0.05%.** The 3/2 is right to 0.03%; the entire
residual sits in the clearing term, which is **9% above 1** at 5σ. That is a
last-molecule stochastic effect the deterministic `y → 1/Ω` proxy cannot see, and
it is not derived. Small, specific, and open.

### Why the original collapse looked so good

`Ω^{−a}` and `1/lnΩ` are hard to separate over ×32, and the obvious cross-check
does not separate them either. The width in `1/H` is `ΔH/H*²` with `ΔH ∝ Ω^{−a}`
and `H*` crossing over from `∝ Ω^{−a}` to `Hc` — so under FSS it rises, peaks and
then falls, and it happens to be **flat to ~1% right across §5's range** (its own
fitted parameters give 5.28 → 5.79 → 5.20 over Ω = 40…1280). Measured there:
5.28 → 5.73. Consistent with both. It only becomes a discriminator ×512 further
out, where FSS demands 0.20× and the measurement gives 1.08×. Recorded because I
nearly
quoted the flat width as evidence *for* the log law when it was evidence for
nothing. What separated the two was deriving the functional form instead of
choosing it, plus a ×16384 range.

**What survives §5.** The sharpening is real, the collapse quality is real, and
"expand fast and the decision freezes half-made" is real at every finite Ω. What
does not survive is `Hc > 0`, the power law, and therefore the hunt for a
universality class: **open question 1 and THEORIES T2 are void, not open.**

**Caveats.** All of this is the *symmetric* start, which is what §4–§6 use; with a
bias there is a genuine positive critical rate (the control). Ω ≤ 655360 by cost,
so `Hc = 0` rests on the Kurtz argument plus a slope that has not budged over
×16384 — not on reaching the limit. One order-parameter convention (§5's, in which
all-blank counts as decided), and the headline fit is one level (D = 0.5), though
all four are reported and agree. The clearing coefficient's 9% excess is measured,
not explained.

---

## 6. Bigger alphabets freeze more easily — `expansion_radix.py`

Running n-winner AM under expansion ties §3 and §4 together. The critical rate falls
monotonically with alphabet size (Ω=160, D=0.5 crossing):

| n | 2 | 3 | 4 | 6 | 8 | 12 | 16 |
|---|---|---|---|---|---|----|----|
| H*(n) | 0.121 | 0.110 | 0.101 | 0.089 | 0.079 | 0.074 | 0.071 |

H*(n) appears to saturate toward a floor, mirroring c(n) — **but that is a
small-Ω artifact; see §6.1.** The frozen relic is richer: fast expansion leaves
≈n coexisting species, versus 2 for AM.

**So a larger alphabet is penalised twice** — a lower restoration margin *and* a
lower tolerance for expansion. The restoration-margin penalty does bottom out
(§3, §14); the expansion penalty does **not** (§6.1).

### 6.1 The radix freeze-out penalty is unbounded, and 1/λ(n) sets it

§5.1's reduction applies unchanged (the argument needs only uniform reaction order,
which n-winner has), so `H*(n, Ω) = 1/τ*(n, Ω)`. **Check: the table above is
reproduced to 1–3% by an *ordinary, non-expanding* SSA measuring nothing but
consensus time** (`--sec6-check`, Ω=160, same D=0.5 level, 4×2000 trials):

| n | 2 | 3 | 4 | 6 | 8 | 12 | 16 |
|---|---|---|---|---|---|----|----|
| §6's H*(n) | 0.121 | 0.110 | 0.101 | 0.089 | 0.079 | 0.074 | 0.071 |
| 1/τ*, no expansion | 0.1221 | 0.1097 | 0.0993 | 0.0863 | 0.0801 | 0.0729 | 0.0692 |
| ratio τ*·H*_§6 | 0.991 | 1.003 | 1.017 | 1.031 | 0.986 | 1.015 | 1.025 |

So §4–§6 are all one measurement: the consensus-time distribution.

**Prediction, before running.** At γ=0 the n-winner symmetric point has
`x_i* = 1/(2n−1)` and symmetry-breaking eigenvalue `λ(n) = 1/(2n−1)`
(THEORIES T7, §14). Shot noise starts the *relative* asymmetry at
`√((2n−1)/Ω)`, which must grow to O(1) at rate λ, so

    dτ*/dlnΩ = 1/(2λ(n)) = (2n−1)/2      → 1.5, 2.5, 3.5, 5.5 at n = 2, 3, 4, 6

— unbounded in n, because λ(n) → 0 like 1/(2n). Measured (Ω = 360…46080, 8
replicates × 1250 trials, `results/freezeout_law_n*.json`):

| n | 2 | 3 | 4 | 6 |
|---|---|---|---|---|
| predicted (2n−1)/2 | 1.5 | 2.5 | 3.5 | 5.5 |
| measured, top-4 pairs | 1.492 ± 0.034 | 2.471 ± 0.040 | 3.373 ± 0.041 | 4.921 ± 0.118 |
| ratio | 0.995 | 0.988 | 0.964 | 0.895 |

and the local slopes **rise toward the prediction from below**, more slowly the
larger n: n=4 runs 2.97 → 3.45 across the sweep, n=6 runs 4.21 → 5.01. So the law
holds, approached from below, and the approach is not finished at n=6.

**Why from below, measured directly.** The collective route needs all n species
alive while the slow mode grows, and a species that fluctuates to zero is gone for
good (`B + X_i → 2X_i` needs `X_i > 0`). Counting survivors at the crossing, n=6:

| Ω | 60 | 120 | 240 | 480 | 960 | 3840 | 46080 |
|---|---|---|---|---|---|---|---|
| Ω/(2n−1) per species | 5.5 | 10.9 | 21.8 | 43.6 | 87.3 | 349 | 4189 |
| species alive at D=0.5 | 4.15 | 4.75 | 5.34 | 5.73 | 5.90 | 5.98 | **6.00** |

At a handful of molecules per species a third of the alphabet is already extinct
when the contest is decided — a **faster** route than the collective mode, which is
why the slope falls short. By ~90 molecules per species it is gone.

**That is exactly the regime §6 measured in.** At Ω=160, n=16 there are
`160/31 = 5.2` molecules per species. So §6's apparent saturation of `H*(n)` is the
extinction route taking over as n grows at fixed Ω, not a floor. At large Ω the
penalty is much stronger and still growing:

| n | 2 | 3 | 4 | 6 |
|---|---|---|---|---|
| H*(2)/H*(n) at Ω = 160 (§6) | 1.000 | 1.100 | 1.198 | 1.360 |
| H*(2)/H*(n) at Ω ≈ 46080 | 1.000 | **1.361** | **1.689** | **2.231** |

**And Hc = 0 at every n**: D at H = 0.055 falls monotonically to 0.207 (n=3), 0.094
(n=4), 0.045 (n=6) across the sweep, and refitting `Hc + CΩ^{−a}` on Ω ≤ 1280
extrapolates 18% low at n=3 exactly as it did for AM.

**The absorption decomposition holds at every n, with one shared constant.** §5.1
found absorption = breaking + clearing, clearing measured at 1.0895 ± 0.0176:

| n | 2 | 3 | 4 | 6 |
|---|---|---|---|---|
| measured absorption slope | 2.565 ± 0.029 | 3.573 ± 0.075 | 4.495 ± 0.021 | 6.050 ± 0.157 |
| breaking + 1.09 | 2.582 | 3.561 | 4.463 | 6.011 |

Four independent instances, agreeing within errors, with the clearing coefficient
fitted once on the AM biased control and never refitted.

**Half of the written prediction was wrong, recorded.** Before running I wrote this
down as a "**two-regime** claim": an extinction regime at small Ω and a collective
regime at large Ω, with the slope rising toward `(2n−1)/2` past a crossover at a
few times `(2n−1)`. The rise is real and the mechanism is right, but there is **no
crossover** — the slope rises smoothly from the very first pair, and the survivor
count decays smoothly too (4.15 → 6.00 with no knee). It is one continuous
approach, not two regimes. "Regime" was doing rhetorical work that the data does
not support.

**Caveats.** n ≤ 6, and the deviation from `(2n−1)/2` is 10% at n=6 and shrinking
with Ω, so "unbounded in n" is an extrapolation from four points plus the exact
`λ(n) = 1/(2n−1)`. §6's dominance order-parameter convention throughout. The
extinction table is one replicate of 2000 trials per Ω, single seed.

---

### 6.2 Independent verification of §6.1

§6.1 overturns §6, so it was checked with a separate instrument written from the
mechanism rather than by re-running `freezeout_law.py` — ordinary `gillespie_fast`
on `n_winner`, absorption (exactly one committed species left) as the consensus
criterion, symmetric start, 400–500 trials per point.

**The saturation claim is confirmed.** Two arms, identical instrument, differing
only in what is held fixed. Ratio of successive `H*` — approaching 1.0 means
flattening:

| n step | 2→3 | 4→6 | 8→12 | 12→16 |
|---|---|---|---|---|
| fixed Ω = 160 (as §6 did) | 0.886 | 0.913 | **0.969** | **0.955** |
| fixed 40 molecules/species | 0.855 | 0.842 | 0.822 | 0.886 |

At fixed Ω the decline stalls, reproducing §6's floor. At fixed molecules per
species it does not stall anywhere. **§6's saturation was an artifact of holding Ω
fixed while the alphabet grew**, and §6.1 stands.

**The `(2n−1)/2` slope structure is confirmed at small n and unresolved above it.**
Absorption runs one clearing term steeper than the decision, so the testable
structure is `absorption slope − (2n−1)/2 ≈ const ≈ 1`:

| n | 2 | 3 | 4 | 6 |
|---|---|---|---|---|
| absorption slope | 2.523 | 3.600 | 4.146 | 5.868 |
| minus (2n−1)/2 | **1.023** | **1.100** | 0.646 | 0.368 |

At n = 2, 3 that difference matches §5.1's independently measured clearing term of
1.0895 to 1–6%, which is real support for the structure. At n = 4, 6 it falls away
— but those sweeps are 3–4 points over a factor 4–8 in Ω, giving slope errors of
±0.2–0.3, so the deficits are only ~2σ. §6.1 sees the same "approached from below"
pattern with 15 sizes over ×16384 and attributes it to extinction; **this check can
neither confirm nor refute that, and should not be cited as doing either.**

**One error in the checking, recorded because it produced a convincing false
alarm.** The first version compared a sweep in which n and Ω moved *together*
against a law about the Ω-derivative at *fixed* n, ignoring the intercept `B(n)`.
That showed τ* at n=16 as 38.1 against a "predicted" 100.2 — a 2.6× discrepancy
that was entirely an artifact of the comparison, with `B(n)` absorbing all of it.
A law stated as a partial derivative cannot be tested by a sweep that moves both
variables.

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
with the ability to restore.

> **This κ is wrong and §15 corrects it.** It scales the restoring gain with γ
> and leaves the finite-count *diffusion* at its γ=0 value, when the reverse
> reactions add noise along the decision mode: `D₀(γ) = (1+γ)/9`, so
> `κ(γ) = (3/2)(1−2γ)/(1+γ)`. Refitting the 216 cells below with the corrected
> value lifts the pooled collapse from **R² = 0.933 to 0.960** and every per-γ
> slope toward 1. Everything else in this section stands; the numbers quoted
> below are the ones produced with the uncorrected κ and are left as published,
> with the refit tabulated in §15.2. The limits are the two known results:
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

**~~The ceiling's "factor ≈3" is not a prefactor — it is a 7% error in the
exponent.~~ — WITHDRAWN by §76.1: it is precisely a prefactor, the Laplace factor that
`exp(δ*²/2σ²)/4` discards, and the exact ceiling c\*/Φ(−δ*/σ) puts the ratio at
4c\*√(2π)·(δ*/σ) — rising with δ*/σ, not a constant 3. The exponent was never wrong.**
Original reading retained below per rule 7:

**The ceiling's "factor ≈3" is not a prefactor — it is a 7% error in the
exponent.** §12.1 predicted `D_max ≈ exp(δ*²/2σ²)/4` and measured about three times
that, which was filed as a missing Laplace prefactor. Read at a *fixed* Ω the factor
gives two different answers (drifting 3.00 → 2.41 at Ω=64, flat 3.00 → 3.33 at
Ω=128) because the death depth is still converging in Ω, and converging more slowly
the smaller σ is. Extrapolating each σ to Ω→∞ the way this section already did at
σ/δ* = 0.45 (increments shrink geometrically, ratio r):

| σ_ch/δ* | Ω=16 | Ω=32 | Ω=64 | Ω=128 | Ω=192 | D∞ | r | predicted | factor |
|---|---|---|---|---|---|---|---|---|---|
| 0.45 | 6.53 | 7.44 | 8.27 | 8.86 | 9.01 | 9.06 | 0.25 | 2.95 | 3.07 |
| 0.38 | 14.48 | 19.28 | 23.33 | 26.07 | 27.01 | 27.50 | 0.34 | 7.97 | 3.45 |
| 0.32 | 38.32 | 65.19 | 92.18 | 112.48 | 120.56 | 125.91 | 0.40 | 33.00 | 3.82 |
| 0.28 | 90.86 | 205.28 | 354.61 | 488.58 | 547.98 | 595.31 | 0.44 | 147.12 | **4.05** |

The factor drifts monotonically 3.07 → 4.05 across σ, so it is not a σ-independent
prefactor — the same second-axis test that ejected §14 (§14.1). But letting the
coefficient in the exponent float instead fits almost perfectly:

    D_max = 0.663 · exp( 1.0695 · δ*²/2σ² )       R² = 0.999782,  residuals ±2.6%

against the saddle point's `k = 1` and prefactor `1/4`. **The whole residual is a 7%
error in `k`.** Forcing `k = 1` compounds that error across the argument's 2.6×
range (`x = δ*²/2σ²` runs 2.47 → 6.38, a 50× range in the prediction itself) and
re-expresses it as a prefactor that then appears to drift. The apparent factor is
just `4·0.663·exp(0.0695·x)`, which reproduces the measured 3.07 / 3.45 / 3.82 /
4.05 as 3.15 / 3.37 / 3.72 / 4.13 — the "drift" is the exponent error in disguise. The exponential *form* is right, the
ceiling is real and the σ-scaling is right; the coefficient is not exactly 1.

Two honest caveats. The residuals alternate in sign (−, +, +, −), which is a
curvature signature — with four points and two parameters that is suggestive of a
remaining shape error rather than proof of one. And `D∞` is itself an extrapolation
(the geometric ratio `r` grows 0.25 → 0.44 as σ falls, so the σ=0.28 row leans on it
hardest). The drift survives without any extrapolation — the raw Ω=192 factors are
3.05 / 3.39 / 3.65 / 3.72, still 1.22× — so the conclusion does not rest on `D∞`,
though the precise value of `k` does.

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
`c(ε) = (3/2)ε²`, recovered as the n=2 member. **The saturation is genuinely
derived**, and that part is robust: `λ/D₀ = (2n−1)/(2n−3)` involves no δ at all.

But the predicted floor is `δ²/2 = 0.0050` against a measured 0.0022.

### 14.1 The residual is not a prefactor — a correction

An earlier version of this section called the 2.27 offset "a constant, and
therefore a prefactor of the same family as §12's 0.74 and §12.1's ≈3", and
consolidated it into THEORIES Q7 on that basis. **That was wrong**, and the test
that killed it was the one written down with the claim: a genuine prefactor is
δ-independent, so measure the offset at a second δ.

Measured at two more δ, twice, with independent seeds and grids (R² 0.991–0.997):

| n | δ=0.10 | δ=0.16 | δ=0.24 |
|---|---|---|---|
| offset, n=8 | 1.943 | 1.840 | **1.399** |
| offset, n=16 | 2.274 | 2.018 | **1.364** |

The offset **falls ~40% across δ**, so it is not a constant. The reason is direct:
`c = λδ²/(2D₀)` predicts `c ∝ δ²` exactly, and the measured exponent is not 2.
Fitting `c(δ)` over the three δ:

| n | 2 | 8 | 16 | 32 | 64 |
|---|---|---|---|---|---|
| measured exponent p | 2.08 (§2) | 2.27 | 2.48 | 2.53 | **2.40** |

extended to n = 32 and 64 with all three δ measured under a single protocol
(2000 trials, five Ω each, R² 0.986–0.999), so no point is shared with §3.

**Two readings, and only one of them is solid.**

*Solid:* **p > 2 at large n.** The n ≥ 16 values average 2.471 against a
propagated 1σ of 0.21 — **3.9σ above the quadratic**. So `c ∝ δ²` genuinely fails
for large alphabets, the reduction is incomplete in *shape* rather than amplitude,
and §14's residual does not belong in Q7's prefactor cluster.

*Solid:* **p saturates.** n=32 and n=64 differ by 0.45σ, and the n ≥ 16 spread
(0.132) is well inside 1σ. p does **not** keep climbing.

*Not solid, and an earlier version of this section over-read it:* the **monotone
climb** from 2.08 to ≈2.47. With ±0.21 on each exponent that separation is only
**1.9σ**, and n=2→8 alone is 0.9σ. The sequence was presented here as a clean trend
before its uncertainty had been propagated. What the data supports is "p is about
2.5 for n ≥ 8 and about 2.1 at n = 2", not a resolved curve.

**Independent cross-check, and it passes.** Extrapolating each single-protocol fit
back to δ = 0.10 and comparing with §3's separately-measured value: **0.00199 vs
0.00227 (n=32) and 0.00217 vs 0.00224 (n=64)** — 12% and 3%. Two unrelated
protocols agreeing at that level is the main reason to trust the exponents at all.

What survives unchanged: the two closed forms, their exact agreement with the
engine, `D₀(2) = 1/9` recovering `design.md` §9, and the **saturation mechanism**,
which contains no δ and so is untouched by any of this.

**And p saturating is probably the same fact as c(n) saturating.** §3's explanation
— past n ≈ 16 the champion's share has converged to δ and the field is fragmented,
so the contest stops changing — predicts that *nothing* about the escape keeps
evolving, the exponent included. Q9 therefore folds into §14's saturation rather
than standing as a separate mystery.

### 14.2 Is p ≈ 2.5 physics, or band placement? Partly checked

There is a reason to expect the Gaussian picture to fail at large n, and a reason
to suspect the measurement instead. **Both point at the same quantity: molecules
per rival.**

Under a fixed pairwise margin the champion's share tends to δ while each rival's
tends to `(1−δ)/(n−1) → 0`, so at large n the rivals are few-molecule species and
van Kampen's expansion — which is what supplies `p = 2` — is exactly what breaks.
That would make p > 2 real physics.

But `c` is *defined* as the slope of `ln P(error)` in Ω, and each δ had to be
measured in a different Ω band, because a larger δ drives the error down and the
band must keep it fittable. At n=64 that means the δ points sampled different
regimes outright:

| δ | Ω band | molecules per rival |
|---|---|---|
| 0.14 | 300–680 | 4.1 – 9.3 |
| 0.18 | 180–420 | 2.3 – 5.5 |
| 0.24 | 80–190 | **1.0 – 2.3** |

If `c` drifts with Ω, then "c at δ=0.24" and "c at δ=0.14" are different quantities
measured in different regimes and their ratio is not an exponent at all.

**Tested directly** (n=32, δ=0.18, 4000 trials, Ω = 150…540, Poisson-weighted):
`ln P` is **linear**, with the quadratic term giving F = 0.01 on 1,5 dof — not
close to significant. Local slopes bounce (0.0065–0.0133) without drifting, at
about 2× Poisson scatter. So over that sweep `c` is a well-defined constant and
the confound does not bite.

**What that does and does not license.** The sweep covered **4.0 – 14.3 molecules
per rival**. The band that actually worried me — δ=0.24 at n=64, at **1.0 – 2.3** —
is *below* the tested range and remains unchecked. Nor is it easily checked: at
that δ the error rate falls like `exp(−0.0181·Ω)`, so reaching 10 molecules per
rival means Ω ≈ 800, an error rate ~10⁻⁶, and ~10⁷ trials per point. The band
placement was forced by the same constraint that created the confound.

**And the significance is carried entirely by the suspect point.** Refitting with
only the two δ whose bands sit in the verified-safe range:

| fit | exponent | above quadratic |
|---|---|---|
| δ = 0.14, 0.18 only | 2.40 ± 0.32 | **1.2σ — not significant** |
| including δ = 0.24 | 2.47 ± 0.21 | 3.9σ |

The precision comes from the lever arm in `ln δ`, and the lever arm *is* the
δ=0.24 point. Drop it and p is consistent with 2; keep it and p > 2 at 3.9σ. So
the claim does not merely have an untested caveat attached — **its entire
statistical weight sits on the measurement most likely to be biased.**

**Status, stated plainly: p > 2 is unproven.** What is solid is that `c` is a
well-defined constant over 4–14 molecules per rival, that the two protocols
cross-check to 3–12%, and that the *saturation* results (both `c(n)`'s and `p`'s
apparent one) are untouched by any of this, since they involve no δ comparison.

**Why this is parked rather than pushed.** The natural fix — hold molecules per
rival fixed while varying δ — **is not possible**. `c` is *defined* as
`−∂ln P/∂Ω` at fixed (n, δ), and `m = (1−δ)Ω/(n−1)` varies along that derivative
by construction. One cannot differentiate in Ω while holding a quantity
proportional to Ω fixed. Settling this needs either a different observable than
an Ω-slope, or ~10⁷ trials per point to lift the large-δ band into the safe range.
Neither is a good use of effort against a second-order correction to a barrier
coefficient.

**Caveats.** The closed forms are for γ=0 (irreversible), which is what §3 measured;
`breaking_diffusion` computes the general-γ case numerically but it is untested
against a barrier measurement. The comparison inherits §3's fixed-margin convention
(§3.1) and its δ=0.10. The barrier is read from a 1-D reduction along one mode,
which is exact for the *rate* by symmetry but is an approximation for the escape.

---

## 15. The wall coefficient, measured instead of expanded — `quasipotential.py`, `wall_coefficient_exact.py`

§2 derives the restoration wall's coefficient for irreversible AM as the saddle's
restoring gain against the finite-count diffusion, `κ = λ/(2D₀)` with `λ = 1/3`
and `D₀ = 1/9`, giving `κ = 3/2`. §12 carried it to `γ > 0` by scaling the gain:
`κ₁₂(γ) = (9/2)·λ(γ) = (3/2)(1−2γ)`.

**That scales the gain and leaves the noise alone, and the noise depends on γ
too.** The reverses `2X → B+X` and `2Y → B+Y` are extra jumps along the decision
mode `v = (1,−1,0)`, each contributing `γ/9` to `D₀ = Σ_r (v·S_r)²·a_r` at the
symmetric point — the disagreement reaction and its reverse are orthogonal to `v`
and contribute nothing. So

    λ(γ)  = (1−2γ)/3      the gain shrinks
    D₀(γ) = (1+γ)/9       the noise GROWS
    κ(γ)  = λ/(2D₀) = (3/2)·(1−2γ)/(1+γ)

Unchanged at γ = 0, so **§1–2 stand untouched**. The ratio is exactly
`κ/κ₁₂ = 1/(1+γ)`: κ₁₂ is **45% high at γ = 0.45**, or equivalently κ is 31%
below it.

**The correct ingredient was already in the repo.**
`networks/n_winner_reversible.breaking_diffusion(2, γ)` returns exactly `(1+γ)/9`
and has since §13–§14, and `λ_breaking(2,γ)/(2·breaking_diffusion(2,γ))` reproduces
the formula above to machine precision. §12 used a closed form in a different
module with the diffusion baked in at its γ=0 value, and nothing connected them.

### 15.1 Two independent measurements, neither of them an expansion

**Route A — the exact quasipotential.** `W(n) = −(1/Ω)·ln P_ss(n)` is the
*definition*, and `cme.stationary` gives `P_ss` exactly, so the ridge curvature at
the saddle **is** κ. Both limits are taken rather than assumed: fit window `w → 0`
(quartic terms make the curvature drift, exactly as §2's own table drifts
1.586 → 1.809) and `Ω → ∞` (`W` is Ω-independent only to leading WKB order).

| γ | κ measured (Ω→∞) | `(3/2)(1−2γ)/(1+γ)` | ratio | κ/κ₁₂ | `1/(1+γ)` |
|---|---|---|---|---|---|
| 0.35 | 0.333482 | 0.333333 | **1.0004** | 0.7411 | 0.7407 |
| 0.40 | 0.214008 | 0.214286 | **0.9987** | 0.7134 | 0.7143 |
| 0.45 | 0.103570 | 0.103448 | **1.0012** | 0.6905 | 0.6897 |

Agreement to **0.1% at three γ** against a formula with no fitted parameter,
where κ₁₂ overshoots by 35% / 40% / 45%.

**Route B — the original instrument.** Exact first-passage error probability from a
biased start, fit as `−ln P = c·Ω + const`, then `c/ε²` — §1–2's own method, which
touches neither the stationary distribution nor the ridge minimisation, and which
reaches the γ where Route A cannot go. Across the 10 cells with adequate
statistics (γ = 0.25–0.40, ε/δ* = 0.22–0.34):

    mean c/eps^2 divided by (3/2)(1-2g)/(1+g)  =  0.990   (range 0.88 - 1.07)
    mean c/eps^2 divided by (3/2)(1-2g)        =  0.755   (never near 1)

The per-cell ε→0 extrapolations (1.06–1.12) are **not** quoted: three points each,
drifting in the opposite direction from §2's, and the small-ε cells lose
populations to the p < 0.1 cut. The raw cells are the evidence.

**Route A's window runs the wrong way, which is worth stating.** `P_ss` comes from
a double-precision solve, so probabilities more than ~1e-13 below the mode are
round-off; since `W` is a log, the resolvable barrier is capped at `≈30/Ω` — which
**shrinks** as Ω grows. Meanwhile ε is quantised at `1/Ω`, so a narrow window at
small Ω is a parabola through a handful of sites. The two guards pull in opposite
directions and the window between them is **empty for γ ≤ 0.25**: at γ = 0.25 the
lattice needs Ω ≥ 137 and the floor needs Ω ≤ 115. That is a real limit of the
route. Before those guards existed this experiment produced two contaminated
points, one from each side (γ=0.35 at Ω=400, γ=0.45 at Ω=150), both of which
looked like ordinary data. Route B exists partly to cover the gap.

### 15.2 What it does to §12, on §12's own stored data

§12's 216 cells store their raw `p_flip`, so the collapse refits with no rerun:

| | pooled | γ=0.05 | γ=0.15 | γ=0.30 | γ=0.45 |
|---|---|---|---|---|---|
| R², §12's κ | 0.9329 | 0.9894 | 0.9564 | 0.8944 | 0.9688 |
| R², corrected κ | **0.9604** | 0.9916 | 0.9702 | **0.9349** | 0.9741 |
| slope, §12's κ | 0.742 | 0.795 | 0.626 | 0.419 | 0.497 |
| slope, corrected κ | **0.783** | 0.812 | 0.675 | **0.507** | **0.684** |

Every cell improves and every slope moves toward 1. **But the slope does not reach
1 and stays non-monotone in γ** (0.81, 0.68, 0.51, 0.68), so this is part of §12's
residual and not the whole of it — the rest is presumably §12's own second saddle
point, which minimises a sum of two exponents and keeps only the minimum.

**What this explicitly does not do.** It does not explain §12's fitted slopes as a
missing `1/(1+γ)`. Those slopes are non-monotone and no smooth function of γ
reproduces them; that was the first hypothesis, it was checked against §12's
per-γ table, and it failed. The claim here is about the coefficient — verified
against two exact instruments — not about §12's regression.

`information.wall_coefficient` now returns the corrected value;
`wall_coefficient_gain_only` keeps §12's published numbers reproducible.

---

## 16. Tilted landscapes: what asymmetry buys and what it costs — `networks/am_asymmetric.py`

Every network up to here is symmetric under relabelling the symbols, so both
attractors are mirror images and one coefficient describes both. Q4 asked what
happens when they are not. The minimal honest tilt puts a factor on each
autocatalytic branch and keeps every reverse at γ× **its own** forward:

    f2: B + X -> 2X   at k(1+β)      r2: 2X -> B + X   at γk(1+β)
    f3: B + Y -> 2Y   at k(1-β)      r3: 2Y -> B + Y   at γk(1-β)

Every reversible pair keeps ratio `1/γ`, so the **cycle affinity stays `−3 ln γ`
for every β** (checked against the general `cycle_affinity`, not asserted), and
the Wegscheider product is still `γ³`. So β is a clean second axis against γ:
it costs no thermodynamic *force*. It is not free in the dissipation *rate*,
which is force times flux — a distinction this project has got wrong before
(§9.2's threshold, §10.3's control rails) and which is left measured, not assumed.

**The tilt breaks the pitchfork into an imperfect bifurcation.** Raising β
deepens the X basin, shrinks the Y basin, and slides the saddle toward Y until
saddle and Y-attractor annihilate at `β_c(γ)`. Past β_c the network is
monostable: **it answers X no matter what it is shown.** β_c collapses as the
landscape weakens — 0.998 / 0.809 / 0.281 / 0.103 / 0.050 at γ = 0.05 / 0.20 /
0.35 / 0.42 / 0.45 — so near γ_c **a 5% rate mismatch destroys the device
outright**, while at strong drive almost any tilt is survivable. (The γ = 0.05
entry sits against the bisection's upper bound of 0.999 and should be read as
"essentially any tilt", not as a measured number; β_c → 1 as γ → 0 because Y's
autocatalysis only vanishes at β = 1.)

**The bias lives in the saddle, not the attractors — at strong drive.** At
γ = 0.05 the attractors sit at (0.000, 0.952) and (0.952, 0.000) even at
β = 0.5·β_c ≈ 0.50, while the basin boundary has moved to `−0.367·δ*`: X wins from
a starting Y-majority of a third of the landscape width. **Reading a tilt off the
attractor positions would report no tilt at all.** This weakens as γ rises — at
γ = 0.45 the same β/β_c puts the attractors at (0.160, 0.526) and (0.558, 0.134),
visibly unequal — so the "attractors don't move" reading is a strong-drive
statement, not a general one. The boundary shift is the general one: `−0.12·δ*`
to `−0.22·δ*` at β = 0.25·β_c across the whole γ range. This is a *systematic*
error, a different failure mode from the random error the wall protects against,
and the two are easy to confuse because both present as "wrong answer".

### 16.1 A prediction of mine that was wrong

Written before running: *"`c₊ + c₋` is not conserved; it decreases monotonically
with |β|, so symmetric AM maximises total restoration capacity"* — reasoning that
the saddle is driven into the shallow attractor until they annihilate, so `c₋ → 0`
while `c₊` cannot diverge. **The sum goes up, not down.**

Escape barriers from the exact quasipotential, at Ω = 250:

| | β = 0 | β = 0.5·β_c | β = 0.97·β_c |
|---|---|---|---|
| **γ = 0.42** (β_c = 0.1034) | | | |
| ΔW₊ (X basin) | 0.02388 | 0.04361 | 0.06706 |
| ΔW₋ (Y basin) | 0.02388 | 0.00937 | 0.00107 |
| sum | 0.04776 | 0.05298 | **0.06813** |
| **γ = 0.45** (β_c = 0.0500) | | | |
| ΔW₊ | 0.01017 | 0.01790 | 0.02685 |
| ΔW₋ | 0.01017 | 0.00429 | 0.00069 |
| sum | 0.02034 | 0.02219 | **0.02754** |

The same at Ω = 150, so it is not a single-Ω artifact. The favoured basin deepens
**2.8×** while the disfavoured one loses only 0.023: tilting *creates* barrier
height rather than reallocating it. The error in the prediction was assuming `c₊`
is bounded — the saddle rises relative to *both* wells, since it is one bottleneck
and it moves up as it moves sideways.

**But the sum is maximised exactly where the device stops working.** At
β = 0.97·β_c the Y barrier is 0.001 — Y cannot be stored at all. So the total
barrier height is the wrong figure of merit, and that is the useful part of the
refutation: a metric can improve monotonically while the thing it is supposed to
measure collapses.

### 16.2 The figure of merit that does work

Mutual information through the tilted restorer, symmetric source, same |ε| for
both symbols, exact first-passage. Reported as `I/I(β=0)`:

| β/β_c | Ω=120 | Ω=200 | Ω=300 | Ω=400 |
|---|---|---|---|---|
| 0.25 | 0.953 | 0.906 | 0.854 | 0.810 |
| 0.50 | 0.825 | 0.672 | 0.527 | 0.419 |
| 0.75 | 0.650 | 0.408 | 0.230 | 0.132 |
| 0.95 | 0.504 | 0.239 | 0.093 | **0.037** |

γ = 0.42, input ε = 0.10·δ*. **I falls monotonically with |β| at every Ω, so
symmetric AM is optimal for a symmetric source** — β = 0 is a maximum, not merely
a stationary point.

Repeating the sweep at ε = 0.25·δ* reproduces the same shape but **not** the same
numbers: within 1% at β ≤ 0.25·β_c, then diverging to 0.062 against 0.037 at
β = 0.95·β_c, Ω = 400. Both are monotone and both show the Ω-amplification below,
so the conclusion is robust to the input strength while the *rate* of collapse is
not — the deep-tilt cells are two small numbers whose ratio is unstable, and they
should not be quoted as a measured tilt penalty.

**The population makes it worse, not better.** The penalty at fixed β/β_c grows
steadily with Ω: at β = 0.95·β_c a 3.3× population takes the retained information
from 50% down to 3.7%. Molecules buy reliability for the favoured symbol
exponentially and buy nothing for the other, so more molecules *amplify* a design
asymmetry instead of averaging it out. Everywhere else in this file more molecules
help; this is the first place they reliably hurt.

**One harness bug worth recording**, because it produced a clean, plausible,
impossible number. The first run showed `P(ok|X) = 0.638` against
`P(ok|Y) = 0.667` at **β = 0**, where symmetry forces them equal. The cause was
integer division in building the two initial states: `(rest ± d0)//2` gave biases
of +9 and −11 counts when `rest − d0` was odd. Forcing the parity fixes it. The
asymmetry under test was 3% and the artifact was 20% of it.


---

## 17. A design rule: matching the tilt to a biased source — `biased_source.py`

§16 showed `β = 0` maximises information for a 50/50 source, which is the easy
case — symmetry alone makes `β = 0` stationary and the only content was the sign
of the second derivative. **Q4a asks the case with no symmetry to lean on.**

Predictions were written down before running (`biased_source.py`'s docstring
holds them verbatim). With `e₊ = P(err | X sent)`, `e₋ = P(err | Y sent)`, small
errors, and `e± = A·exp(−Ω·c±)`, setting `dL/dβ = 0` on the information deficit
gives `p·h(e₊) = (1−p)·h(e₋)` (**P1**), whose leading form is
`ln(e₋/e₊) = ln(p/(1−p))` (**P2**), hence `β* ∝ ln(p/(1−p))/Ω` (**P3**) and a
deficit falling by `2√(p(1−p))` (**P4**). β* is found by maximising the *exact*
mutual information, and P1–P4 are then checked at the measured β*.

**The headline: β\* > 0 at every prior. Symmetric restoration is not optimal for
a biased source**, and this is the first statement in this project about how to
*build* the chemistry rather than how it behaves.

### 17.1 The form is exact; the coefficient is not 1

At Ω = 200, γ = 0.35, sweeping the prior:

| p | 0.60 | 0.70 | 0.80 | 0.90 | 0.95 |
|---|---|---|---|---|---|
| β\* | 0.00455 | 0.00949 | 0.01546 | 0.02428 | 0.03211 |
| β\*/β_c | 0.016 | 0.034 | 0.055 | 0.086 | 0.114 |
| `ln(e₋/e₊)` measured | 0.319 | 0.664 | 1.083 | 1.702 | 2.253 |
| `ln(p/(1−p))` | 0.406 | 0.847 | 1.386 | 2.197 | 2.944 |
| ratio | 0.786 | 0.784 | 0.782 | 0.775 | 0.765 |

    ln(e-/e+) at beta* = 0.7625 * ln(p/(1-p)) + 0.0178      R² = 0.999867

**P2's shape is confirmed to R² = 0.9999** — the ratio is constant to 2.7% across
a 7.3× range in log-odds, and the intercept is 0.018 against a predicted 0. So the
rule is real and it is exactly the predicted one: **tilt until the log-ratio of
the two error probabilities matches the prior log-odds** — up to a coefficient
that is 0.76 here rather than 1.

Note how *gentle* the optimal tilt is: β\*/β_c runs 0.016 → 0.114, so the tilt
that helps is 1–10% of the tilt that destroys the device (§16's fold). Nothing
about the optimum lives near β_c.

### 17.2 The coefficient rises with Ω, and *this* data cannot say where it stops

At p = 0.80, sweeping Ω:

| Ω | 100 | 150 | 200 | 300 | 400 |
|---|---|---|---|---|---|
| β\* | 0.02259 | 0.01839 | 0.01546 | 0.01175 | 0.00946 |
| `ln(e₋/e₊)` / `ln(p/(1−p))` | 0.644 | 0.711 | 0.781 | 0.844 | 0.879 |

Monotone toward the predicted 1. **But the limit is not determined**, and the
tempting move — fitting `1/Ω` because that is what P3 predicts — is the one to
avoid:

| assumed correction | extrapolated limit | R² |
|---|---|---|
| `1/Ω` | 0.947 | 0.973 |
| `Ω^−0.75` | 1.004 | 0.985 |
| `1/√Ω` | 1.120 | **0.992** |
| `1/lnΩ` | 1.685 | **0.993** |
| free power `c − aΩ^−b` | 1.342 | — |

**The two forms that fit best both overshoot 1.** Quoting the `1/Ω` row alone
would give "extrapolates to 0.947, confirming P2" from the *worst*-fitting ansatz
of the five. Four times in Ω is not enough to resolve this, and the honest
statement is that the coefficient is rising and 1 is inside the plausible range.

**§17.3 resolves it** — with a 10× range and a test that does not require picking
an ansatz. The value of leaving this subsection standing is that the tempting
shortcut and the real answer can be compared: the shortcut would have reached the
right conclusion by the wrong route, which is not the same as being right.

P3 itself shows the same incompleteness directly: `β* ∝ Ω^−x` with x measured at
**0.508 / 0.603 / 0.678 / 0.752** across consecutive Ω pairs — drifting toward the
predicted 1 and nowhere near it yet.

### 17.3 The coefficient is 1 — `tilt_rule_limit.py`

§17.2 could not decide this: over a 4× range the candidate extrapolations gave
0.947 / 1.004 / 1.120 / 1.685 and the best-fitting two overshot. **The obstacle
was cost, not principle, and the measurement was badly posed.** β\* was being
found by *maximising* the exact mutual information — ~40 CME solves per cell. But
`dI/dβ = 0` is a scalar equation in `p`, so the question inverts: instead of "given
p, which β is optimal", ask **"given β, which prior makes it optimal"** — a root
find in `p` needing no solves once the error curve `e±(β)` is known. One sweep of
~13 β values yields the whole line. That is ~20× cheaper, and it reproduces the
direct optimisation to **0.0001 in slope** at Ω=200 (0.7624 against 0.7625).

| Ω | 100 | 150 | 200 | 300 | 400 | 600 | 800 | 1000 |
|---|---|---|---|---|---|---|---|---|
| r | 0.5868 | 0.6760 | 0.7624 | 0.8368 | 0.8752 | 0.9111 | 0.9334 | **0.9459** |
| intercept | 0.055 | 0.033 | 0.019 | 0.007 | 0.003 | 0.0013 | 0.0005 | **0.0002** |
| R² | 0.9974 | 0.9992 | 0.9998 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |

Over a 10× range the extrapolations are still spread — 0.978 (`1/Ω`), 1.029
(`Ω^−0.75`), 1.131 (`1/√Ω`), 1.701 (`1/lnΩ`) — so quoting one is still choosing the
answer. **Two tests that do not require choosing:**

    deficit:      1 - r  =  27.9 * Omega^(-0.902)      R² = 0.9978
    free limit:   r -> 1.037 +- 0.030                  1.0 is 1.23 sigma away

The first is the decisive one. If `r` converged to anything below 1 the deficit
would approach a constant and the log-log plot would flatten; instead it is a
clean power law across the whole 10×, with residuals under 0.06 in ln units.

**So the design rule is exactly parameter-free:**

    ln(e-/e+) at beta*  =  ln( p / (1-p) )

with a finite-population correction of `1 − 28·Ω^{−0.90}`, which is why it reads
0.76 at Ω = 200. The intercept, predicted zero, falls to 0.0002.

### 17.4 The careful refinement was worse than the crude argument

**P1 is refuted.** It was derived *after* P2, as the more careful version keeping
the log factors, and it fails badly — `p·h(e₊) / ((1−p)·h(e₋))` should be 1 and is:

| p | 0.60 | 0.70 | 0.80 | 0.90 | 0.95 |
|---|---|---|---|---|---|
| P1 ratio | 1.21 | 1.48 | 1.90 | 2.80 | **4.05** |

Not only wrong but systematically worse the more extreme the prior, while P2's
ratio stays flat to 2.7% over the same range. The extra structure P1 added is
structure that is not there. Worth recording: *a refinement that fits worse than
the thing it refines is evidence against the refinement, not noise* — the same
logic §5.1 used to reject the three-parameter freeze-out fit.

**P4 is directionally right and quantitatively far off.** The deficit ratio at
Ω = 200 measures 0.995 / 0.979 / 0.947 / 0.888 / 0.837 at p = 0.60 → 0.95 against
the asymptotic `2√(p(1−p))` = 0.980 / 0.917 / 0.800 / 0.600 / 0.436. So the
**realisable** gain from matching the tilt at accessible population is 0.5–16%,
not the asymptotic 2–56%. It is a real gain and it is not a large one; anyone reading
P4 as the payoff would overstate it by 3×.


---

## 18. Which one survives: the relic asymmetry — `relic_asymmetry.py`

AM's disagreement reaction is literally an annihilation — `X + Y → 2B` is
matter + antimatter → two photons, with `2B → X + Y` as pair production — and the
three ingredients Sakharov requires for a matter excess are all knobs already in
the rig: **number violation** (`X−Y` is untouched by annihilation and pair
production; only the recruitment reactions move it), **C/CP violation** (§16's
tilt β), and **departure from equilibrium** (γ < 1, and §5.1's expansion deadline
`1/H`).

**Where it stops being a mapping, stated up front.** Recruitment `B + X → 2X` is
autocatalysis and has no counterpart in the standard picture, where an asymmetry
survives *linearly*: annihilation removes matched pairs and the pre-existing
excess is what is left. Here the excess is *amplified* by an instability. This is
therefore not a model of baryogenesis — it is the question of what changes when an
asymmetry is fed through a **restoring** landscape rather than a passive one. That
difference is the content and it should not be dressed up as a cosmology result.

### 18.1 Dynamical or accidental, and a parameter-free answer

From an exactly symmetric start two things can decide which species survives: the
tilt's deterministic push and shot noise.

    g/λ  =  (2β/3)·(1−γ)/(1−2γ)                     the tilt
    σ    =  sqrt( (D₀/2) / (λΩ) ) = sqrt( (1+γ) / (3(1−2γ)Ω) )    the noise

— the σ uses **§15's corrected `D₀(γ) = (1+γ)/9`**, so this result depends on that
correction. The relic is dynamical when `g/λ > σ`, i.e. above

    β·√Ω  =  (√3/2)·(1−2γ)/(1−γ)          = 0.820 at γ = 0.05

**Prediction, written before running: P(X survives) collapses onto a single curve
in `u = (g/λ)/σ`, and that curve is `Φ(u)`.** Both parts parameter-free. Computed
exactly as a splitting probability from the symmetric start, so no sampling error
can hide a failed collapse.

| u | 0.25 | 0.50 | 1.00 | 1.50 | 2.00 | 3.00 |
|---|---|---|---|---|---|---|
| spread of P(X) across Ω = 60→240 | 0.0013 | 0.0022 | **0.0028** | 0.0019 | 0.0007 | 0.0001 |
| mean P(X) − Φ(u) | −0.0027 | −0.0049 | **−0.0065** | −0.0050 | −0.0025 | −0.0002 |

**Both hold to under 1%.** A 4× population change moves P(X) by at most 0.0028 at
fixed u, and the parameter-free `Φ(u)` is right to 0.0065 at its worst. The
residual is systematically negative and shrinking with Ω (at u=1: 0.0080 / 0.0063 /
0.0052 at Ω = 60 / 120 / 240), consistent with the one approximation in the
derivation — the effective-seed picture treats accumulated noise as a single
Gaussian kick.

So **`β√Ω ≈ 0.82` separates a relic whose sign is set by the chemistry from one
set by a coin flip.** The scaling is the interesting half: the asymmetry needed to
be decisive *falls* as `Ω^{−1/2}`, so a larger system needs a smaller bias, not a
larger one.

### 18.2 The deadline does not just decide whether — it decides which, harder

With expansion, §5.1's exact time change means the expanding SSA is ordinary SSA
stopped at internal time `1/H`, so this is a stopped run, not a second integrator.

**Prediction, written before running:** the deadline decides *whether* a relic
forms, not *which* — `P(X | decided)` stays `Φ(u)` at every H while `P(decided)`
falls with H. **The first half is wrong.**

`P(decided) / P(X | decided)`, γ = 0.05, Ω = 120, 800 trials:

| u | Φ(u) | H=0.02 | H=0.05 | H=0.1 | H=0.2 |
|---|---|---|---|---|---|
| 0.0 | 0.500 | 1.000 / 0.509 | 0.961 / 0.531 | 0.249 / 0.487 | 0.000 / — |
| 0.5 | 0.691 | 1.000 / 0.709 | 0.958 / 0.735 | 0.319 / **0.788** | 0.000 / — |
| 1.0 | 0.841 | 1.000 / 0.848 | 0.975 / 0.863 | 0.448 / **0.958** | 0.001 / — |
| 2.0 | 0.977 | 1.000 / 0.975 | 0.995 / 0.974 | 0.756 / **0.997** | 0.003 / — |

At a generous deadline `P(X | decided) = Φ(u)` as predicted. At `H = 0.1` it runs
far above: 0.958 against 0.841 at u = 1 is **11σ** on 358 decided trials.
Reproduced at Ω = 240 (0.766 against 0.691 at u = 0.5, H = 0.07, 4.4σ), and the
u = 0 row stays at 0.5 throughout, which is the symmetry check that says this is
not an artifact.

**The mechanism is a selection effect.** Beating the deadline requires fast growth
away from the symmetric point; the tilt supplies exactly that; so conditioning on
"a relic formed at all" preferentially keeps tilt-aligned trajectories. The two
effects do not factorise, and the reading is the opposite of the intuition that
rushing the expansion randomises the outcome: **the faster the expansion, the
purer the surviving relic's alignment with the asymmetry** — there is just far
less of it. Past `H ≈ 0.2` nothing decides at all (P(decided) ≤ 0.003).

**One caveat on the table.** The `H = 0.12` cells at Ω = 240 have 10–37 decided
trials out of 900 and are quoted nowhere for that reason; the `H = 0.2` column
above is reported only as "nothing survives", not as a ratio.


---

## 19. Giving the model a temperature — `crnl/cooling.py`, `cooling_relic.py`

**What was missing.** Every AM reaction is 2→2, so under expansion every propensity
scales by the *identical* factor and the ratios never move: `γ`, `δ*(γ)`, `κ(γ)`,
`β` and `γ_c` are all invariant. **The landscape is frozen under expansion; only
the clock slows.** That is why §5.1's reduction to "ordinary SSA stopped at
internal time `1/H`" came out exact — it could not have come out otherwise, and
`expanding.common_order` enforces the uniform-order property that makes it true.

So in this rig, *expanding the volume* and *slowing down time* were literally the
same operation. Real freeze-out is not a clock running out; it is the
**equilibrium moving**. The model had no temperature, and the concrete
consequence is that §18 could measure a relic *sign* but there was no relic
**abundance** to measure — the equilibrium never went anywhere.

**The minimal fix.** Adiabatic expansion gives `T ∝ Ω(t)^{−w}`; the reverses carry
an activation energy, so `γ = exp(−ΔE/T)` and `γ(t) = γ₀^{exp(wHt)}`. Forward
rates are untouched — only the balance moves. On §5.1's internal clock, with
`s = Hτ ∈ [0,1)`:

    γ(s) = γ₀^((1−s)^(−w))

**which does not contain H.** The cooling schedule is universal in `s`; H enters
only as an overall rate. So the same sweep of γ from γ₀ down through γ_c to 0
always happens, and **H decides how many reactions fit inside it** — the
competition (cooling deepens the landscape, dilution starves it) that a fixed-γ
model cannot have. This is genuinely not a time change.

**Anchored on the known result.** At `w = 0` the new integrator reproduces
`gillespie_expanding` **0/300 mismatches**, state-for-state and step-for-step —
the same bit-for-bit check §5.1 used on the time change itself.

### 19.1 A relic abundance, and it is set by the expansion

Two arms differing only in whether the medium cools. Cooling: `γ₀ = 0.55` (just
above `γ_c = ½`, so **there is no landscape at the start** and the pitchfork is
crossed at `s = 0.358` for every H — the symmetry breaking is *driven by the
expansion*). Fixed: `γ = 0.05` held constant, §5.1's model. Ω = 300, 400 trials.
Relic = minority share of the committed population, **conditioned on having
decided**.

| H | 0.005 | 0.010 | 0.020 | 0.050 | 0.100 |
|---|---|---|---|---|---|
| **cooling** relic | 0.00029 | 0.00573 | **0.08373** | 0.26779 | 0.33915 |
| relic / equilibrium(γ_freeze) | 1.9e5 | 3.2e6 | 9.9e7 | 1.5e8 | 4.4e6 |
| undecided | 0.000 | 0.000 | 0.068 | 0.285 | 0.623 |
| **fixed** relic | 0.00013 | 0.00020 | 0.00013 | 0.00950 | 0.17265 |
| relic / equilibrium | 1.0 | 1.5 | 1.0 | 72 | 1312 |
| undecided | 0.000 | 0.000 | 0.000 | 0.013 | 0.158 |

Read only the clean cells (H ≤ 0.02, undecided ≤ 7%), which is a 4× range in H:

* **P1 confirmed.** The cooling relic rises **290×**, 0.00029 → 0.08373. Faster
  expansion freezes the annihilation earlier and leaves more of the minority
  species — the standard cosmological direction.
* **P2 confirmed.** The relic sits `10⁵–10⁸` **above** the equilibrium abundance at
  the γ where it froze. Adiabatic following would give a ratio of 1; this is
  freeze-out.
* **P3 confirmed.** The fixed arm is **flat** — 1.3e-4 / 2.0e-4 / 1.3e-4 across the
  same range — and equal to its equilibrium value (ratio 1.0). Its abundance is
  set by the chemistry and does not know the expansion rate exists.

**Abundance set by expansion versus abundance set by chemistry** is the whole
contrast, and it is the observable the fixed-drive model structurally could not
produce.

**A summary line that argued the opposite, and why it was wrong.** Quoting the
full H range gives "cooling 1163×, fixed 1318×" — no difference between the arms,
i.e. the exact opposite conclusion. At large H *both* arms are dominated by
marginal decisions, and conditioning on "decided" does not rescue a cell that was
62% undecided; the fixed arm's apparent 1318× is entirely its H = 0.10 cell.
`cooling_relic.py` now refuses to summarise cells above 10% undecided. **A
conditional mean is only as good as the thing conditioned on being rare.**

### 19.2 What this does and does not overturn

**It does not overturn `Hc = 0` (§5.1).** The impossibility argument there is that
from an exactly symmetric start the deterministic ODE stays exactly symmetric
forever, so `D(H,Ω) → 0` for every `H > 0`. Cooling does not touch that: the x↔y
symmetry is preserved by *any* `γ(t)` as long as `β = 0`, so the symmetric ODE
still never decides. §5.1 stands as written.

**What it does change** is the scope of §5.1's reduction. "The expanding SSA is
exactly ordinary SSA stopped at internal time 1/H" is a theorem about
**uniform-order kinetics with state-independent rate constants**, not about
restoration under expansion in general. Give the medium a temperature and the
reduction fails immediately — the drive profile `γ(s)` is universal but the
*number of reactions inside it* is not, which is precisely the content of §19.1.

**Still not modelled.** Forward rates have no temperature dependence, so this is
the minimal change that lets the *balance* move rather than a thermochemistry;
the medium is still well mixed; and the drive is still an infinite reservoir that
never depletes — a cascade dissipates but nothing ever runs down, so §12.1's depth
ceiling remains purely noise-limited with no thermodynamic competitor.


---

## 20. A second ceiling: restoration that dies of exhaustion — `networks/am_fueled.py`, `fuel_ceiling.py`

**The last free lunch.** Everywhere above, `γ` is a free parameter held fixed
forever — an infinite reservoir, set once and maintained at no cost. §9 therefore
measures what restoration *dissipates* while nothing ever runs down, and §12.1's
depth ceiling is purely noise-limited with no thermodynamic competitor.

**Fuel as chemistry, not bookkeeping.** A drive `γ < 1` physically *is* a coupling
to fuel hydrolysis, so the fuel becomes a reactant:

    f1: X + Y + F -> 2B + W        r1: 2B + W -> X + Y + F     (× γ∞)
    f2: B + X + F -> 2X + W        r2: 2X + W -> B + X + F     (× γ∞)
    f3: B + Y + F -> 2Y + W        r3: 2Y + W -> B + Y + F     (× γ∞)

The drive the chemistry feels is `γ_eff = γ∞·w/f`, which **rises** as the tank
empties — the mirror of §19's cooling, which drove γ down. Letting the integrator
adjust γ from a running count of firings instead would be the harness doing the
chemistry, the failure mode behind three withdrawn results here.

**`n_F` is a genuinely independent coordinate.** A complete cycle `f1→f2→f3`
returns (X, Y, B) exactly to where it started while consuming three fuel, so the
fixed-γ model is a **projection** that discards a coordinate which must exist.

**Anchored on the known case**, in §19's style: with the tank held fixed, the
(X, Y, B) drift equals `am_reversible(γ∞·w/f)` with time rescaled by `f`, to
**1e-16** over 120 random interior states at four (f, w).

**What it costs structurally.** These reactions are 3→3, so the network leaves the
uniform-order-2 class that made §5.1's reduction exact. It stays *uniform* order 3,
so `expanding.common_order` still accepts it and §19's machinery survives with
`λ = 2H` — but trimolecular steps are a real idealisation, stated rather than hidden.

### 20.1 The two ceilings have different shapes

Bit held at the X attractor; loss is the first state with `n_X ≤ n_Y` (a readout,
never an intervention). Control: the same chemistry at the same drive held forever.

**The control needs its clock fixed, and this is not a detail.** The fueled network
is third order, so its rates carry a factor `f` — at fuel concentration 10 with
`γ_eff(0)=0.30` that is `f₀ = 7.69`. Comparing raw lifetimes across the arms without
scaling the control by `f₀` compares different clocks, which is exactly the
mismatched-control error that produced §10.3's withdrawn result. Unscaled, the
control read 1483 at Ω=30; scaled, 193. Every number below is on a shared clock.

| Ω | 30 | 45 | 60 | 90 | 120 | 180 |
|---|---|---|---|---|---|---|
| fueled lifetime, φ=10 | 7.53 | 7.60 | 7.61 | 7.15 | 7.37 | 7.72 |
| fueled lifetime, φ=30 | 4.23 | 4.32 | 4.59 | 4.68 | 4.91 | 4.85 |
| noise control, φ=10 clock | 193 | — | — | — | — | — |

**P3 confirmed, and it is the headline.** The fuel-limited lifetime is **flat in Ω**
— spread **1.08×** (φ=10) and **1.16×** (φ=30) across a **6× population range**, at
two fuel concentrations, which is the second axis. The noise-limited lifetime over
the same clock is exponential:

    ln(lifetime) = 0.1215·Ω + const        R² = 0.9842

so the two cross at **Ω ≈ 3** (φ=10) and **Ω ≈ 8** (φ=30), and by Ω = 180 the noise
ceiling is **2.3e9×** and **1.2e9×** further away. **Above a population of about
ten, restoration is fuel-limited and more molecules buy nothing** — the exact mirror
of §1's wall, where molecules bought exponential reliability. Both ceilings are
real; only one had ever been measured.

**More fuel gives a *shorter* lifetime** — 7.5 at φ=10 against 4.6 at φ=30 — which
is not a mistake. In waste-fraction units the burn rate is exactly
φ-independent, `dω/dt = (1−ω)(xy+bx+by) − γ∞·ω(b²+x²+y²)`, so **a bigger tank buys
no extra fractional runway**. What it does buy is faster chemistry relative to the
fuel clock, so the state follows the collapsing landscape more adiabatically and
gives the bit up *earlier in the tank's life* (see §20.2). More of the resource does
not help, for the same structural reason more molecules do not.

### 20.2 A prediction of mine that was wrong, in sign

Written before running: *"the bit is lost BEFORE the formal death point, because
γ_eff rises continuously and the barrier degrades all the way up to it — the loss
fraction should approach 1 from below."* **It exceeds 1 everywhere.**

| Ω | 30 | 45 | 60 | 90 | 120 | 180 |
|---|---|---|---|---|---|---|
| w_loss/w_death, φ=10 | 1.240 | 1.287 | 1.307 | 1.328 | 1.350 | **1.367** |
| w_loss/w_death, φ=30 | 1.073 | 1.123 | 1.163 | 1.186 | 1.212 | **1.238** |

**The bit outlives the landscape's death by 7–37%.** What the prediction missed:
once `γ_eff` passes `γ_c` the landscape is gone, but the state still has to
physically *relax* off the old attractor, and that relaxation burns further fuel.
I accounted for the barrier degrading and not for what happens after it vanishes.

The overshoot **grows with Ω and shrinks with fuel concentration**, and both make
sense on the same mechanism: crossing `n_X = n_Y` needs a fluctuation, which is
relatively smaller at large Ω, while a richer tank makes the chemistry fast
compared to the fuel clock so the state tracks the collapse more adiabatically.

**Whether it saturates is undetermined**, and this is flagged rather than resolved:
`0.069·lnΩ + 1.018` fits at R² = 0.969, and a saturating form `c − aΩ^{−b}` fits
comparably with `c = 1.43 ± 0.03`. Six times in Ω cannot separate them — the same
limitation §17.2 hit, and the same refusal to pick the flattering one.

**§12.1's ceiling is untouched but no longer alone.** That ceiling is a depth, set
by channel noise and independent of the drive; this one is a *lifetime*, set by the
budget and independent of the population. A cascade now has two ways to die, and
which one binds is decided by the fuel concentration rather than by either
mechanism on its own.


### 20.3 Which ceiling binds for a cascade — `cascade_fuel_vs_noise.py`

§20.1 gives a fuel-limited *lifetime*; §12.1 gives a noise-limited *depth*. They
are different units and had never been put side by side. Both arms run through one
harness here — same channel, same stage length, same clock — differing only in
whether the drive can run out. The channel kicks the decision coordinate by moving
molecules between X and Y **only**, so `n_X+n_Y`, `n_B` and the tank are untouched;
not resetting `n_B` is deliberate, since a reset would hand the chemistry a fresh
blank pool it did not earn.

**The stage has to be measured in relaxation times, not absolute time.** Fixing
`t_stage = 8` made the first run useless: at fuel concentration 10 the whole tank
lasts ≈7.5 time units (§20.1), so the cascade died at depth 1. The relaxation time
is `1/(λ·f₀)`, which itself shrinks as the tank gets richer, and holding the stage
at a fixed multiple of it is what makes depths comparable across budgets.

Median depth at which the bit is lost, Ω = 40, γ₀ = 0.30, stage = 2 relaxation times:

| σ_ch/δ* | 0.03 | 0.08 | 0.15 | 0.22 | 0.30 | 0.40 |
|---|---|---|---|---|---|---|
| control (drive held forever) | **>400** | **>317** | 60 | 21 | 7 | 7 |
| fuel, φ = 50 | 11 | 10 | 8 | 7 | 5 | 5 |
| fuel, φ = 150 | 23 | 23 | 14 | 12 | 8.5 | 6 |

(The two control cells marked **>** are censored — 62% and 40% of trials reached the
depth cap — so they are lower bounds, not measurements.)

**P2 and P3 confirmed.** The control depth is unbounded as the channel quietens
while the fueled depth is not, so the binding ceiling changes hands: at σ/δ* = 0.40
the two arms agree within a factor 1.2 (noise binds both), and at σ/δ* = 0.03 the
fueled cascade dies **36× earlier** than the control. **Below the crossover a
cascade dies of exhaustion** — a failure mode nothing in §1–§19 could produce.

**P1 refuted, and the reason is the interesting part.** I predicted the fuel-limited
depth would be roughly independent of the channel noise, since the tank is drained
by the restoring chemistry rather than by the channel. It falls **2.2×** (φ=50) and
**3.8×** (φ=150) across the σ range. **The two ceilings are not independent and do
not combine as a `min()`.** Spent fuel raises `γ_eff`, which shrinks both `δ*` and
`κ(γ)` — so a half-empty tank makes the *same* channel noise bite harder. Exhaustion
and noise compound: the drive degrading is itself what lets the channel win.

**P4's premise was wrong, and its conclusion is unsettled.** P4 assumed
`D_fuel ∝ Φ`, which would make the crossover move only as `1/√(ln Φ)` — a
logarithmically weak lever. Measured, `D_fuel` grows **sub-linearly**: 11 → 23 for a
3× budget, i.e. `∝ Φ^0.67`. The crossover does move in the predicted direction (more
fuel ⇒ noise binds sooner: above σ/δ* = 0.40 at φ=50, down to ≈0.26 at φ=150), but
that is a ≥35% shift where the log estimate said ~10%. **Two fuel concentrations and
a six-point σ grid cannot pin the scaling**, and with the ceilings compounding there
is no reason to expect the clean form P4 assumed. Recorded as unsettled rather than
fitted.

> **Settled in §23**, and the two-point reading above was right: five budgets give
> `Φ^0.6498 ± 0.019` at σ/δ* = 0.03. But the exponent **drifts with the channel and
> with Ω**, so there is no single `D_fuel(Φ)` to pin — the numbers above stand for
> their cells and the framing is reworded there.

**What this settles for T10b-ii.** A cascade has two ways to die and the fuel
concentration decides which. But the honest headline is the coupling, not the
competition: *restoration does not run until the fuel is gone and then stop — it
degrades continuously as the fuel goes, and the channel finishes it early.*


---

## 21. What a simulation may throw away — `crnl/approximations.py`, `approximation_hierarchy.py`

CRNL's method is a two-point version of this question — ODE against exact SSA, and
the gap is the subject. This fills in the levels between, so "what can a simulation
discard and still get restoration right?" becomes a measurement.

The observable is the one with an exact answer: `P(error)` from a start biased by
ε, as a CME splitting probability with **no sampling error at all**. Every
approximate level is scored against that reference rather than against each other.

**Why this is not a numerics exercise.** Kurtz's theorem says the density process
converges to the mass-action ODE on finite time intervals, and §5.1 leans on it to
prove `Hc = 0`. It is true, and it does **not** license discarding the molecules for
this observable: restoration lives in tails and long times, where the convergence is
not uniform. The error probability vanishes in the limit while being nonzero at
every finite Ω. **A limit theorem cannot tell you what your simulation may throw
away.**

### 21.1 The error exponent, by level

`c = d(−ln p)/dΩ` fitted over Ω = 40 → 140, γ = 0.30, 6000 trials per cell:

| ε/δ* | CME (exact) | SSA | CLE | τ=0.05 | τ=0.3 |
|---|---|---|---|---|---|
| 0.25 | 0.01503 | 0.01500 | 0.01543 | 0.01612 | 0.01689 |
| 0.40 | 0.03587 | 0.03606 | 0.03673 | 0.03607 | 0.03452 |
| **ODE** | **p = 0** | **p = 0** | **p = 0** | **p = 0** | **p = 0** |

**P1 confirmed. The ODE fails categorically, not quantitatively.** It reports
exactly 0 in all sixteen cells where the truth ranges over 1.5e-3 to 1.6e-1. There
is no refinement parameter that improves it — it is not an inaccurate number, it is
the wrong kind of number.

**P2 confirmed. The SSA is the exact chain.** Exponents agree with the CME to 0.2%
and 0.5%. The per-cell ratios sat 1–3% high across six cells, which is a 1-in-64
sign pattern, so it was rechecked at 60,000 trials: **ratio 1.0085 (z = +0.89) and
1.0100 (z = +0.81)** — sampling noise, not bias. The anchor holds.

### 21.2 The prediction that was wrong, and in sign

P3 argued the CLE should **overestimate** the failure probability: it is the
quadratic truncation of the jump Hamiltonian `Σ a_j(e^{p·S_j} − 1)`, and for a 1-D
birth–death chain the exact barrier `∫ln(a₊/a₋)` exceeds the diffusion
approximation's `∫2(a₊−a₋)/(a₊+a₋)`, since `ln r > 2(r−1)/(r+1)`.

**It underestimates it.** The CLE exponent is **+2.7% and +2.4%** — a *larger*
barrier, hence a smaller error probability. The 1-D intuition does not survive the
move to two dimensions, and the sign of the correction is not something I could
read off the scalar case.

The magnitude prediction fares no better: I expected the discrepancy to **grow with
the barrier**, and it is flat (+2.7% at ε/δ* = 0.25, +2.4% at 0.40). Two ε values
cannot establish flatness either — recorded as unsupported rather than refuted.

**P4 partly confirmed.** τ-leaping degrades monotonically with the window at
ε/δ* = 0.25 (0.01500 → 0.01612 → 0.01689) and does not at 0.40 (0.03607 → 0.03452,
moving the other way). One ε in each direction is not an interpolation law.

### 21.3 The finding all four predictions missed: it is a cliff, not a slope

Every level that keeps **any** noise at all gets the restoration exponent right to
2–12%. The one level that keeps **none** is infinitely wrong. There is no graded
degradation in between — the ODE is not "the coarsest member of a hierarchy", it is
categorically apart, and CLE, τ-leaping and the exact SSA are all in one class.

Put as an answer to the question this section exists to ask: **you can discard the
discreteness (CLE keeps real-valued counts), the exact jump timing (τ-leaping fires
in windows), and the correct jump distribution (Gaussian instead of Poisson) — and
still recover the restoration exponent to a few percent. What you cannot discard is
having noise at all.**

> ⚠ **"Having noise at all" is refined in §24, and this section could not have seen
> it.** Every level on this ladder keeps *every species*, so it retains the noise and
> the coordinates together and cannot separate them. Separating them on one kernel
> (§24) gives the coordinate axis **+0.1000** of exponent against the noise axis's
> **+0.0146** — and once the signal-carrying coordinate is kept, adding noise to a
> bookkeeping one is worth −0.0026. The measurements in §21 stand exactly as printed;
> what changes is which axis they are evidence about. The sharper statement both
> support: **noise matters in the coordinate that carries the signal, and costs
> nothing in a bookkeeping coordinate.** §24 records it as a suspect with a kill test
> aimed at this ladder — **and §24.1 runs it here and it holds**: with the noise
> projected onto the blank pool only, 88% of the variance retained, this section's
> observable reads exactly 0 in all eight cells, the ODE's own failure.

The corollary is about cost. The exponent — the physics of §1–§2, the thing the
whole project is about — is reproduced by a cheap SDE at O(1/dt) per unit time. The
expensive exactness (O(Ω) events for the SSA, O(Ω²) memory for the CME) buys the
*prefactor* and the individual probabilities. So a simulation that only needs to
know **how fast reliability grows with population** can be cheap; one that needs to
know **the actual failure rate** cannot.

**Scope, stated plainly.** One network, one observable, γ = 0.30, Ω ≤ 140, two ε.
Whether the cliff is a general feature of restoring networks or a property of AM is
untested, and the 2.5% CLE offset is small enough that a second network could
plausibly move it either way.


### 21.3a Correction: the CLE arm is not converged, and its small deviations do not stand

Pushing the statistics (a trajectory-batched integrator, `approximations.run_batch`,
which reproduces the reference propensities to 0.0) exposed something the original
run could not see: **the CLE's answer depends on its own step size and on its
negativity policy at the ~10% level** — larger than the deviations §21.1 and §21.4
reported for it.

Two independent symptoms, both at n = 3, ε/δ = 0.25:

| | Ω = 45 | Ω = 60 |
|---|---|---|
| p at dt = 0.02 | 0.0985 | 0.0953 |
| p at dt = 0.005 | 0.1092 | 0.0945 |
| p at dt = 0.001 | 0.1055 | **0.0848** |

and switching the negativity policy alone — halving the step and retrying, versus
rejecting that sweep — moved the CLE/CME ratio from 0.855 to 0.967 at Ω = 45. Euler
–Maruyama drives the state negative often here (thousands of rejections per run),
because the minority species sits near zero, which is exactly where the restoration
observable lives.

**Correction to this correction — the reason above is wrong.** Both symptoms were
themselves measured at 4,000 trials, where the standard error on p is ~5%, so the
"11% drift" was ~2σ. Re-run at **40,000 trials** with the batched integrator, the
CLE is *converged and correct*:

| n=3, Ω=45 | dt=0.02 | 0.01 | 0.005 | 0.002 | 0.001 |
|---|---|---|---|---|---|
| CLE/CME | 0.9916 | 1.0005 | 0.9970 | 1.0271 | 0.9998 |
| z vs exact | −0.58 | +0.04 | −0.21 | +1.85 | −0.01 |

No dt dependence, and agreement with the exact reference at every step size. So the
CLE arm is **not** unconverged — the earlier scatter was sampling noise, and this
subsection's first draft withdrew the right numbers for the wrong reason.

**What this retracts, on the corrected grounds.** §21.2's +2.5% and §21.4's 5.3% are
**withdrawn as sampling noise**, not as discretization artifacts. Re-measured at
40,000 trials the n=2 exponent ratio is **0.980** (ε/δ*=0.25) and **1.001** (0.40),
against the 6,000-trial run's 1.027 and 1.024 — the sign itself flipped. The
refutation of my P3 goes with them: there is no established CLE exponent bias at
either alphabet size, so there is no sign to have been wrong about.

**Followed up over a 16x range in dt, and left unresolved on purpose —
`cle_prefactor.py`.** Three step sizes a factor 4 apart, 60,000 trials per cell.
Euler–Maruyama is weak-order 1, so a discretisation excess decays monotonically
with dt and a prefactor excess is flat; the measurement is neither.

| pooled CLE/CME | dt=0.02 | dt=0.005 | dt=0.00125 |
|---|---|---|---|
| ε/δ* = 0.25 | 1.0036 ± 0.0077 | 1.0298 ± 0.0079 | 1.0175 ± 0.0078 |
| ε/δ* = 0.40 | 1.0092 ± 0.0304 | 1.0305 ± 0.0315 | 1.0096 ± 0.0302 |

The excess is real in the sense that it is positive in 5 of 6 pooled cells and
~+1.7% ± 0.5% at the well-resolved ε, but it is **non-monotone in dt and scatters
by more than its binomial error** — so it is not the clean decay discretisation
would give, nor the flat plateau a prefactor would give. The exponent meanwhile is
untouched at every step size (ratios 1.011 / 0.969 / 0.986 and 1.003 / 0.977 /
1.001), which is P3 confirmed and is the part that matters for §21.3.

**This thread is closed as unresolved rather than pursued further.** It has now run
three rounds — claimed at 2.5%, withdrawn for the wrong reason, withdrawn again for
the right one — and a ~2% effect that will not hold still across dt is at the edge
of what this setup resolves. Separating it would need a positivity-preserving
integrator rather than more samples, which is a different piece of work and buys
nothing for the cliff.

**The original hint, kept for the record.** Across all eight n=2 cells at
40,000 trials the CLE sits *above* the exact p — ratios 1.018 to 1.045, every z
positive, combining to roughly 3σ for a **~+3% uniform excess**. A uniform factor on
p is a *prefactor* effect and would not touch the exponent, which is consistent with
the CLE being the correct diffusion limit with a slightly wrong amplitude. It is one
step size and ~3σ; it is recorded as a hint, not a result.

**What survives untouched.** §21.3's cliff, which is the section's actual result. It
rests on a contrast between *exactly zero* and *within ~10%*, and a 10% wobble in the
CLE cannot bridge a categorical failure. The ODE reports 0 in every cell of both
networks; every level that keeps noise lands within about ten percent of exact. That
statement does not depend on resolving the CLE at the percent level. §21.2's P1 and
P2 are also unaffected — P2 was verified at 60,000 trials against an exact reference,
and the SSA has no step size to converge.

**The lesson, which is the familiar one.** The comparison was set up to measure
differences between approximation levels and was never checked for convergence
*within* a level. An approximation's own numerical parameter is a second axis, and
this file's standing rule — constancy along the axis you happened to sweep is not
constancy — applies to dt exactly as it applies to Ω.

### 21.4 The cliff survives a bigger alphabet — `approximation_hierarchy_nwinner.py`

T11a's kill test: the same ladder on `n_winner_reversible` at n = 3, at the same
fraction of the bifurcation (γ = 0.6·γ_c = 0.121, γ_c(3) = 0.202) so the two
networks are asked the same question.

| Ω | ε realised | CME (exact) | ODE | SSA | CLE | τ=0.05 |
|---|---|---|---|---|---|---|
| 45 | 0.2222 | 1.065e-1 | **0** | 1.111e-1 | 9.106e-2 | 1.111e-1 |
| 60 | 0.2167 | 9.160e-2 | **0** | 8.687e-2 | 8.066e-2 | 8.687e-2 |
| 80 | 0.2250 | 4.120e-2 | **0** | 4.100e-2 | 4.050e-2 | 3.975e-2 |
| 110 | 0.2273 | 2.094e-2 | **0** | 2.262e-2 | 1.937e-2 | 1.988e-2 |

`κ = d(−ln p)/d(ε²Ω)`: **0.4909** (CME) / 0.4714 (SSA) / 0.4650 (CLE) / 0.5098 (τ),
i.e. ratios **0.960 / 0.947 / 1.038**.

**The cliff survives.** The ODE reports exactly 0 in all four cells where the truth
spans 2.1e-2 to 1.1e-1, and every level that keeps noise lands within ~5% of exact.
§21.3's statement is not an AM artifact.

**What this run cannot resolve, and it matters.** The SSA is the exact chain, so
its 4.0% deviation is this measurement's **noise floor**, not a bias — and the
CLE's 5.3% sits barely outside it. At n = 2 the noise floor was 0.2–0.5% and the
CLE's 2.5% was cleanly resolved; here it is not. **So whether the CLE's error grows
with alphabet size is untested by this run**, and quoting 5.3% against 2.5% as a
trend would be reading a difference smaller than the anchor's own scatter.

**Two setup problems, both caught by the reference disagreeing with itself.**
`cme.first_passage` scores its favoured set as `n[0] > n[1]`, which is right at
n = 2 and silently wrong above it — a state where X3 has won can satisfy it —
so `cme.splitting_probability` now takes the predicate (it reproduces
`first_passage` to 0.00e+00 on AM). And the first integer construction let the
realised champion-minus-rival margin overshoot by up to n−1 counts, which made the
exact CME error probability **non-monotone in Ω** (0.115 → 0.064 → 0.092 → 0.026).
A first-passage probability cannot do that, which is how it was caught. The margin
is now pinned exactly, Ω = 30 is excluded (the lattice there is 22% off target),
and the fit is against `ε²Ω` rather than Ω so the residual lattice drift is
absorbed rather than mistaken for physics.


---

## 22. How the barrier dies — T4's kill test — `barrier_near_gamma_c.py`

T4 guessed the restoration barrier vanishes like `(γ_c − γ)`. §12 implied **two**
vanishing factors — `κ·δ*²` with both linear in the gap — hence `(γ_c − γ)²`, and a
population cost diverging like `1/(γ_c−γ)²`. That had only ever been *inferred from
a collapse*, never measured by sweeping γ.

Measured directly with §15's exact quasipotential, whose usable window — large γ,
small barrier — is exactly this region. Each `ΔW` extrapolated to `1/Ω → 0`.

| γ | γ_c−γ | δ* | ΔW (Ω→∞) |
|---|---|---|---|
| 0.400 | 0.1000 | 0.5408 | 3.243e-2 |
| 0.420 | 0.0800 | 0.4925 | 2.106e-2 |
| 0.440 | 0.0600 | 0.4345 | 1.201e-2 |
| 0.455 | 0.0450 | 0.3818 | 6.830e-3 |
| 0.470 | 0.0300 | 0.3165 | 3.069e-3 |
| 0.480 | 0.0200 | 0.2611 | 1.368e-3 |
| 0.487 | 0.0130 | 0.2120 | 5.775e-4 |

    dW = 3.0906 * (gamma_c - gamma)^1.9745        R² = 0.999969

**Exponent 2, confirmed. T4's original guess of 1 is dead.** And the approach is
the right shape: the local slope between consecutive γ climbs **monotonically**
toward 2 as the gap closes — 1.935 / 1.952 / 1.962 / 1.973 / 1.993 / **2.0015** —
which is what an asymptotic normal form must do. Near a pitchfork `W = −a·x²/2 +
b·x⁴/4` with `a ∝ (γ_c−γ)`, so the barrier is `a²/(4b)`; the 1.3% shortfall in the
pooled fit is the non-asymptotic end of the sweep, not a discrepancy.

**This was a 2-for-1 and §15 passes it too.** The corrected `κ = λ/(2D₀)` rests on
`λ = (1−2γ)/3` vanishing *linearly* at γ_c. An exponent of 1 would have restored T4
and falsified that; anything but 2 would have broken the pitchfork picture. THEORIES
said the two stand or fall together — they stand.

### 22.1 What it does to the residuals, which is the point of running it

§12's collapse slope (0.783, not 1) and §12.1's ceiling exponent (`k = 1.0695`, not
1) are both saddle-point residuals, and the obvious shared suspect was that the
saddle gets the barrier's **γ-dependence** wrong. **It does not** — the exact barrier
reproduces `κ·δ*²`'s scaling to 1.3% pooled and to 0.1% asymptotically.

So that hypothesis is eliminated, and the two residuals must live in the *other*
ingredients: §12's **second** saddle point (it minimises a sum of two exponents over
the flip location and keeps only the minimum) and §12.1's σ-dependence. They are not
a shared failure of the barrier.

This is the useful outcome. The four items grouped as "the saddle-point exponent is
inexact" are now **two fewer**, and — unlike the Q7 episode, which consolidated on
constancy and was wrong about two of three members — this narrowing came from
measuring the suspected common cause and finding it innocent, rather than from
noticing that some numbers looked alike.


### 22.2 The second saddle is eliminated too, and what is left

§22.1 cleared the barrier's γ-dependence. The next suspect was §12's **other**
approximation: it minimises `f(δ) = (δ*−δ)²/2σ² + κΩδ²` over the flip location and
keeps only the minimum. The exact version is a one-dimensional convolution of the
channel Gaussian against the wall, and it is analytic:

    p = Phi(-delta*/sigma)  +  [Gaussian prefactor] * exp(-kappa Omega delta*^2 / (1 + 2 kappa Omega sigma^2))

The exponent is *exactly* §12's, so the saddle point loses two things: a prefactor,
and an **Ω-independent** term — the channel crossing the saddle unaided — which no
exponent can represent. Refitting §12's own 216 stored cells with the exact form:

| | pooled slope | pooled R² | per-γ slopes (0.05 / 0.15 / 0.30 / 0.45) |
|---|---|---|---|
| §12's saddle | 0.7830 | 0.9604 | 0.81 / 0.68 / 0.51 / 0.68 |
| exact convolution | 0.7756 | **0.9698** | 0.82 / 0.69 / 0.51 / 0.64 |

**R² improves at every γ, so those two missing pieces are real — and the slope does
not move.** It stays at ~0.78 and stays non-monotone in γ. The second saddle is
eliminated as the explanation of the residual.

**What that leaves, and it is now the only ingredient untested.** §12 uses
`c(δ) = κδ²` — the barrier *quadratic in the displacement*. §22 checked how `κδ*²`
scales with γ, not that `c` is quadratic in δ across the range the channel actually
samples. There is already direct evidence it is not: §2's own table has `c/ε²`
drifting **1.586 → 1.809** over ε ∈ [0.04, 0.20], and §14.1 measured `c ∝ δ^~2.5`
for the n-winner at large n. A barrier that stiffens faster than quadratic away from
the saddle would depress the fitted slope exactly as observed, and would do so
γ-dependently because the sampled range of δ scales with δ*(γ).

**Kill test:** replace `κδ²` with the exact ridge profile `W(0) − W(δ)` from
`quasipotential.ridge_profile` inside the convolution above and refit. The
instrument reaches γ = 0.30 and 0.45 — two of §12's four — which is enough, since
those are the two with the worst slopes (0.51 and 0.64). If the slope goes to 1
there, the residual is the quadratic barrier and §12's formula needs only its
`c(δ)` replaced.

**Three suspects, two down, by elimination rather than consolidation.** This is what
the Q7 episode should have looked like.


### 22.3 The quadratic barrier was the residual — identified, not resolved

The kill test named in §22.2, run: replace `c(δ) = κδ²` with the exact ridge
profile `W(0) − W(δ)` inside the convolution and refit §12's stored cells at the
two γ inside the quasipotential's window — which are also its two worst slopes.

| γ | slope, `c = κδ²` | R² | slope, **exact `c(δ)`** | R² |
|---|---|---|---|---|
| 0.30 | 0.5074 | 0.9349 | **0.9006** | 0.8103 |
| 0.45 | 0.6838 | 0.9741 | **1.2495** | 0.9539 |

**The quadratic barrier was the dominant residual.** The slopes move from 0.51 and
0.68 — roughly 40% low — to 0.90 and 1.25, i.e. from well below 1 to bracketing it.
Nothing else tried moved them at all: §22.1's γ-scaling check cleared the barrier's
γ-dependence, and §22.2's exact convolution shifted the pooled slope by 0.007.

**It is not a clean fix, and the two ways it falls short are different.** R² *drops*
at both γ (0.935 → 0.810, 0.974 → 0.954), so the exact-`c` model fits worse even as
it centres better. And the γ-spread is untouched: the two slopes ran 0.51 → 0.68
(ratio 1.35) and now run 0.90 → 1.25 (ratio 1.39). **The correction moves the centre
to ≈1 and leaves the γ-dependence exactly where it was.**

So §12's residual splits into two things that were being read as one: a **large,
now-explained** deficit from assuming the barrier is quadratic in the displacement,
and a **smaller, still-unexplained** γ-dependence that survives every correction
applied so far. §2 already showed the quadratic is only an ε→0 limit — `c/ε²` drifts
1.586 → 1.809 over ε ∈ [0.04, 0.20] — and §14.1 measured `c ∝ δ^~2.5` at large n, so
the first half was predictable in hindsight and is now measured.

**Why the fit gets worse is the honest open end.** Candidates, none tested: the
ridge profile carries its own Ω-error and was taken at a single Ω; the convolution
still omits the transverse relaxation and the algebraic prefactor; and §12's cascade
channel may not be exactly the Gaussian-on-δ this integral assumes. A correction
that improves the slope while degrading the fit is reporting that it has the right
first-order term and the wrong second-order one.

**Sensitivity to the instrument's own axis, which had to be checked before any of
the above is quotable.** The ridge profile was built at a single Ω, and §15 showed
it converges as 1/Ω — so the profile Ω is a second axis on this measurement:

| γ | profile Ω | slope | R² |
|---|---|---|---|
| 0.30 | 150 / 162 / 174 | 0.894 / 0.902 / 0.909 | 0.8095 / 0.8105 / 0.8114 |
| 0.45 | 300 / 500 / 900 | 1.202 / 1.280 / **1.339** | 0.9548 / 0.9533 / 0.9522 |

Two things, and they point opposite ways.

**R² is flat in the profile Ω at both γ** (fourth decimal). So the fit degradation
reported above is **not** the profile's error — candidate 1 is eliminated, and the
omitted prefactor/transverse relaxation (candidate 2) and the channel model
(candidate 3) survive as the explanation.

**But the γ = 0.45 slope is not converged.** It drifts 11% over a 3× change in
profile Ω and is *still climbing* at Ω = 900. So the `1.2495` quoted above carries
at least that much uncertainty and the true value is ≥ 1.34. The γ = 0.30 slope is
solid (1.6% over its whole available window, which is narrow because the barrier
there is near the floor).

**This weakens the "brackets 1" reading and does not weaken the main result.** The
slope at γ = 0.45 may sit well above 1 rather than just above it, which makes the
overshoot larger and the remaining γ-dependence *worse*, not better — the two
slopes are then 0.90 and ≥1.34, a ratio of ≥1.5 against the quadratic model's 1.35.
The finding that the quadratic barrier is the dominant cause is untouched: it rests
on the move from 0.51/0.68, which is far larger than this uncertainty.

**Scoreboard for §12's slope, three suspects and one session.** Barrier
γ-dependence: cleared. Second saddle: cleared, and worth ~0.01 of slope. Quadratic
barrier: **confirmed as the main cause**, worth ~0.4 of slope. Remaining: the
γ-dependence, which no candidate so far touches.


### 22.4 The convolution framework is not a model of §12, and §22.3's reading is withdrawn

§22.2 and §22.3 fitted a convolution model against §12's stored `p_flip`. Fitting
never tests a model's *absolute* correctness, so the model was finally compared
against the exact same quantity it claims to predict: the single-stage flip
probability, computed directly from `cascade_exact.stage_kernel` and
`channel_matrix` with no fitting at all.

Ratio of model to exact, γ = 0.30 and 0.45, Ω = 30/60/90, σ/δ* = 0.15/0.28/0.45:

| model | ratio to exact |
|---|---|
| convolution with **exact** `c(δ)` | **5 – 3688×** |
| convolution with `κδ²` | 1.0 – 10.4× |

**The framework is wrong by up to three orders of magnitude, and the "improved"
version is far worse than the crude one.** The cause is an approximation §22.2 and
§22.3 never named: the convolution assumes the chemistry **runs to completion**, so
the flip probability from displacement δ is `exp(−Ω·c(δ))`. §12's stage has a
*finite* time, `t_stage = 16`, and a stage that has not finished cannot have
completed an escape. The model therefore overestimates flipping, worst where the
barrier is highest — 3688× at γ = 0.30, Ω = 90, σ/δ* = 0.15.

**What this withdraws.** §22.3's reading — that replacing `κδ²` with the exact
barrier moves the slope toward 1 and therefore identifies the quadratic barrier as
§12's residual — **does not survive**. The slope did move, but a fitted slope from a
model that is 10³ out in absolute terms is not evidence about §12's physics. The
measurement stands; the interpretation does not.

**What it establishes instead, and this is worth more.** The exact barrier is
*shallower* than `κδ²` away from the saddle (that is the direction that makes the
convolution predict *more* flipping). So `κδ²` **overestimates** the barrier, which
suppresses predicted flipping — partially cancelling the framework's overestimate
from assuming completion. **§12's formula fits as well as it does partly by error
cancellation**: two wrong ingredients pulling opposite ways. That is why replacing
one of them made agreement worse, and it is a better explanation of §12's residual
than anything in §22.2 or §22.3.

**§22.1 is untouched.** The `(γ_c−γ)²` barrier scaling was a direct quasipotential
measurement with no convolution anywhere near it, and §15 stands or falls with it as
before.

**The lesson, and it is one of this file's own rules.** §22.2 and §22.3 compared
models to each other and to a *fitted* slope, and never once asked what the model
predicted in absolute terms against a quantity that could be computed exactly — even
though `stage_kernel` has been in the repo since §12 and computes precisely that. A
model that is only ever fitted is a model that is never tested.


### 22.5 The cancellation, measured

§22.4 asserted that §12's formula works partly by error cancellation. Asserting is
not measuring, so the two errors were separated. The exact barrier has *no* barrier
error by construction, so `conv(exact c)/exact` isolates the **framework** error
`F`; the ratio between the two convolutions isolates the **barrier** error `B`.

| | range across γ = 0.30/0.45, Ω = 30/60/90, σ/δ* = 0.15/0.28/0.45 |
|---|---|
| **F** — framework (assumes the chemistry completes) | 5.0 → **3687.6** |
| **B** — barrier (`κδ²` vs exact) | **0.0010** → 0.7840 |
| **net** = F·B — what §12's formula actually costs | 1.03 → 10.4 |

**Two errors spanning about three decades each, cancelling to within one.** That is
the cancellation, and it is larger than §22.4 guessed. It explains both facts that
motivated §22.2–§22.4: why §12's crude formula tracks the exact answer as well as it
does, and why replacing one ingredient with an exact one made agreement *worse*.

The two errors pull opposite ways for reasons that are now clear. `F > 1` because a
stage that has not finished cannot have completed an escape, so assuming completion
overestimates flipping. `B < 1` because `κδ²` is *stiffer* than the true barrier away
from the saddle, so it suppresses flipping. And `B → 1` as the barrier shrinks
(0.78, 0.74, 0.70 at γ = 0.45) exactly as it must, since §15 verified `κ` **is** the
true curvature at the saddle — the quadratic is right there and wrong further out.

**The mechanism's functional form is NOT established, and this is where to stop.**
Finite stage time predicts `F` should grow with barrier height, and it does — but a
pooled log-log fit gives `ln F = 0.969·ln(Ω·c(δ*)) + 2.48` at only **R² = 0.62**, and
the behaviour differs *within* each γ: roughly exponential in the barrier at
γ = 0.30 (F = 113 → 3688 for a 1.5× barrier change) and sublinear at γ = 0.45
(F = 5.0 → 8.7 for a 3× change). Direction confirmed, law not. Given that this
thread has already produced one withdrawn interpretation from over-reading a fit,
the functional form is left unclaimed rather than fitted harder.

**What stands from §22 as a whole.** §22.1: the barrier vanishes as `(γ_c−γ)²`,
exponent measured 1.9745 with the local slope reaching 2.0015 — a direct
quasipotential measurement, and §15 stands with it. §22.4–§22.5: §12's residual is
error cancellation between a framework that assumes completion and a barrier that is
too stiff, both now quantified. §22.2 and §22.3's readings are withdrawn.


## 23. What a fuel budget actually buys — T10b-iii's kill test — `cascade_fuel_vs_noise.py`

§20.3 established two ways for a cascade to die — exhaustion and noise — and that
they **compound** rather than combining as a `min()`. It also left the crossover
argument (§20's P4) resting on an untested premise: that the fuel-limited depth
`D_fuel` is **linear** in the budget `Φ`. Two budgets gave `~Φ^0.67`, which was
recorded as unresolved. T10b-iii named the kill test: *four or more budgets at a
fixed quiet channel; if the exponent drifts with σ, `D_fuel` is not a budget
property and the whole two-ceiling framing needs rewording.*

**Scope.** `am_fueled(γ_inf = 1.0)`, γ₀ = 0.3, stage = 2 relaxation times (so
`t_stage ∝ 1/Φ`, §20.1's fair clock), Ω ∈ {25, 40, 60}, `Φ/Ω` ∈ {25, 50, 100, 200,
400}, σ_ch/δ* ∈ {0.03, 0.15, 0.30}, 100–120 trials per cell, `max_depth` 500. Zero
censoring in every fueled cell. The control is `am_reversible(γ₀, k = f₀)` — the
same drive held forever on the same clock.

### 23.1 The exponent is not 1, and the two-point 0.67 was not an artifact

Depths at Ω = 40 (median over 120 trials):

| σ_ch/δ* | control | Φ/Ω=25 | 50 | 100 | 200 | 400 | exponent (mean depth) |
|---|---|---|---|---|---|---|---|
| 0.03 | 500 (57% censored) | 7 | 12 | 17 | 27 | 44 | **0.6498 ± 0.019** |
| 0.15 | 69 | 6 | 9 | 13 | 17.5 | 28 | **0.4954 ± 0.021** |
| 0.30 | 7 | 5 | 6 | 6 | 7 | 9 | **0.2586 ± 0.016** |

> **Prediction P5, refuted.** I predicted the exponent is **1**, on the argument
> that the cascade dies at a condition on the *burn fraction* `w/(w+f)`, which is
> scale-free in Φ, while fuel burned per stage is set by Ω. Fixed usable fraction
> ÷ fixed per-stage cost = linear. Measured 0.6498 ± 0.019 at the quiet channel,
> with linear **18.4σ** away (the 20.1σ figure belongs to the earlier two-σ run's
> 0.6645 ± 0.0167, not to this one — see the sourcing note in §23.3). I also
> predicted the two-point 0.67 was contaminated
> and that the contamination had the wrong sign to produce it; it was neither
> contaminated nor wrongly signed. The a-priori argument was wrong about which
> factor is scale-free — see §23.3.

### 23.2 `D_fuel` is not a budget property — P6 refuted, and its first pass was under-powered

> **Prediction P6** — the load-bearing one — said the exponent is *the same at
> every σ* within fit error: Φ sets the scaling, σ only the prefactor. On the
> first pass (σ/δ* = 0.03 and 0.08) it **passed** at 1.09σ. That pass was worth
> nothing: both values sit on the same side of §20.3's crossover, both
> exhaustion-bound, with depths differing by 1.1×, so the test never asked the
> question. It is recorded here because a null result from an instrument that
> cannot resolve the effect is exactly the failure this project keeps buying rules
> with (§17.2, T10b-i).

Widened to σ/δ* = 0.30, the exponent **falls monotonically**: 0.6498 → 0.4954 →
0.2586, steps of **5.5σ** and **9.1σ** (15.9σ end to end). And the drift is not an
Ω = 40 accident — measured along the axis it was not swept on (rule 9):

| Ω | exponent @ σ/δ*=0.03 | @ 0.15 | drift |
|---|---|---|---|
| 25 | 0.5473 ± 0.028 | 0.4050 ± 0.020 | −0.1423 ± 0.034 (4.2σ) |
| 40 | 0.6498 ± 0.019 | 0.4954 ± 0.021 | −0.1543 ± 0.028 (5.5σ) |
| 60 | 0.6722 ± 0.028 | 0.5269 ± 0.022 | −0.1453 ± 0.035 (4.1σ) |

**Two things at once.** The drift per σ-step is itself remarkably Ω-independent
(−0.142/−0.154/−0.145 across a 2.4× population range), but the exponent *level*
also climbs with Ω (0.547 → 0.650 → 0.672, apparently saturating). So the budget
exponent is a function of **both** the channel and the population.

**The honest caveat, which cuts the 0.30 column out of the argument.** At
σ/δ* = 0.30 the control depth is 7 and four of the five fueled cells are within
0.8× of it — noise binds, not fuel. In the fully noise-bound limit the depth is
Φ-independent by construction, so an exponent falling toward 0 there is partly
definitional and cannot carry the claim. **The load-bearing comparison is
0.03 → 0.15**, where the fueled arm is the binding ceiling at every budget
(28 vs a control of 69 at the largest tank) and the exponent still moves 4.1–5.5σ
at all three Ω. That is the kill test passing on its own terms.

### 23.3 Where the sublinearity lives, and what it costs §20.3

Decomposing `D = (θ − θ₀)·Φ / c`, with θ the waste fraction at loss, θ₀ = γ₀/(1+γ₀)
the seeded waste, and `c` the fuel burned per stage:

| σ_ch/δ* | slope of usable burn fraction | slope of burn per stage `c` | `c` (molecules, Φ/Ω = 25 → 400) |
|---|---|---|---|
| 0.03 | −0.3630 | −0.0127 | 26.4 → 25.0 |
| 0.15 | −0.5151 | −0.0105 | 27.9 → 28.1 |
| 0.30 | −0.7741 | −0.0326 | 33.7 → 30.5 |

`c` is **flat** — 25–34 molecules per stage across a 16× budget at every σ — which
is the one half of P5's argument that survives, and it was an a-priori prediction:
the burn rate is third order so it goes like `f`, and `stage_time` goes like `1/f`.
The entire effect is that the **usable burn fraction collapses**, 0.198 → 0.069,
and steepens as the channel gets louder.

*(The three-term decomposition sums to the measured exponent to 0.0000 by
construction — `c` was obtained as `(θ−θ₀)Φ/D`, so that is bookkeeping, not a test.
What is measured is that the `c` slope is −0.01 to −0.03 rather than something of
order the others, i.e. 1–3% of the effect.)*

Equivalently, in terms of the drive the cascade dies with: γ_eff at loss is **0.75**
at Φ/Ω = 25 — *past* γ_c = 0.5, so the landscape was already dead and the bit
outlived it — but only **0.43** at Φ/Ω = 400, where the landscape is still alive
when the bit goes.

> **Sourcing correction, found by re-deriving every §23 number from the stored
> results.** The slope column above is from this section's stated scope — the
> 5-budget × 3-σ × 3-Ω sweep in `fuel_depth_scaling.json`. Two figures in the two
> paragraphs above were taken from the *earlier* two-σ run instead, and both are
> left visible with the in-scope value beside them:
>
> | quantity | as printed (2-σ run) | in-scope (3-σ run) |
> |---|---|---|
> | usable burn fraction, Φ/Ω = 25 → 400 | 0.198 → 0.069 | **0.191 → 0.070** |
> | γ_eff at loss, Φ/Ω = 25 | 0.75 | **0.729** |
>
> Neither shifts a conclusion — both runs agree that the fraction collapses by ~2.7×
> and that the smallest tank dies past γ_c — but the section should not quote two
> runs in one paragraph, and §23.4's repetition of γ_eff = 0.75 inherits the same
> slip.

**The mechanism, and the rewording §20.3 owes.** A bigger tank does not buy
proportionally more depth because the extra stages each carry their own per-stage
loss probability. More stages spent in the healthy part of the tank means more
accumulated chances to lose the bit, so **the bit goes while the drive is still
strong** — and a louder channel raises the per-stage loss probability, which is why
the exponent is a function of σ. §20.3's "exhaustion and noise compound" is right
and this is the same fact from the budget side; what has to go is the residual
picture of `D_fuel` as *a* number attached to a tank. There is no fuel-limited depth
in the abstract: there is a fuel-limited depth **for a given channel and
population**, and its budget exponent runs from 0.65 down to 0.26 within the range
tested here.

> ⚠ **The mechanism paragraph immediately above is mostly wrong — corrected in
> §23.4.** It survives only as the explanation of the *σ-drift*. Accumulated hazard
> is not what makes the scaling sublinear, and "more stages spent in the healthy
> part of the tank" is the specific phrase that fails: the healthy-tank per-stage
> loss probability is ≤ 0.00125 at σ/δ* = 0.03. Everything else in §23 stands.

**What §23 does not establish.** No functional form for the exponent's σ- or
Ω-dependence is claimed. Three σ (one of which is disqualified above) and three Ω
are enough to establish the drift and kill P6; they are not enough to fit a law, and
§22.5's lesson about fitting harder than the data supports applies here unchanged.

### 23.4 Testing that mechanism absolutely — and correcting it — `fuel_hazard.py`, `fuel_rail_lag.py`

§23.3's mechanism was read out of a decomposition, which is the failure mode rule 16
exists for. Tested as an **absolute prediction with no free parameter**: measure the
per-stage loss probability `q(θ)` and waste made per stage `c(θ)` in single stages
from the attractor of each burn fraction θ, then integrate the survival product
`θ_{d+1} = θ_d + c(θ_d)/Φ`, `S_{d+1} = S_d(1 − q(θ_d))` and read off where `S` crosses
½. Ω = 40, γ₀ = 0.3, 800 trials per (θ, Φ) cell, 14 θ from 0.2308 to 0.48.

**The absolute depths come out right — the exponent does not.** Predicted ÷ measured
median depth, over the five budgets:

| σ_ch/δ* | Φ/Ω=25 | 50 | 100 | 200 | 400 |
|---|---|---|---|---|---|
| 0.03 | 0.571 | 0.667 | 0.824 | **1.000** | **1.091** |
| 0.15 | 0.667 | 0.667 | 0.923 | **1.200** | **1.107** |

Parameter-free agreement to 1.09–1.20× at the large budgets is a real success for
the picture. But the ratio *rises monotonically with Φ*, so the predicted exponent
must be steeper than the measured one, and it is:

| | σ/δ* = 0.03 | σ/δ* = 0.15 | σ-drift |
|---|---|---|---|
| hazard integral | 0.8925 ± 0.017 | 0.7716 ± 0.042 | **−0.121** |
| same integral, hazard forced to 0 | 0.9156 ± 0.015 | 0.8857 ± 0.041 | −0.030 |
| **measured (§23.1)** | **0.6474 ± 0.022** | **0.5404 ± 0.020** | **−0.107** |

> **§23.3's mechanism, corrected.** Switching the hazard off changes the predicted
> depth by **0.0% at four of the five budgets** (7.7% at the largest, quiet
> channel). Accumulated hazard therefore cannot be what makes the scaling
> sublinear — with it switched off the exponent is 0.916, near the linear value P5
> argued for. What the hazard *does* explain is the **σ-drift**: −0.121 predicted
> against −0.107 measured, agreeing within the fit errors, while the no-hazard
> version drifts only −0.030. So the σ-dependence of the exponent is accumulated
> hazard, exactly as §23.3 said; the *sublinearity itself* is not.

**Which leaves a 0.25 gap, and it is not the big tanks under-performing.** The
integral matches the large budgets and under-predicts the small ones by up to
1.75×. The sublinearity is mostly **small tanks over-performing** — surviving
longer than a quasi-static account permits. §23.3 recorded the symptom without
reading it: the smallest tank loses the bit at γ_eff = 0.75, past γ_c = 0.5, where
δ* does not exist and **the landscape is monostable**.

Measured directly (200 trials per cell, tracking γ_eff and δ every stage):

| σ/δ* = 0.03 | Φ/Ω=25 | 50 | 100 | 200 | 400 |
|---|---|---|---|---|---|
| median stages held at γ_eff ≥ γ_c | **3 of 7** | 3 of 11 | 1.5 of 16.5 | 0 of 27 | **0 of 46** |
| trials reaching γ_eff ≥ γ_c | **97.0%** | 91.0% | 66.5% | 41.0% | **9.5%** |
| median δ while past γ_c | 0.188 | 0.200 | 0.225 | 0.213 | — |

The smallest tank spends a median **3 of its 7 stages holding a bit in a landscape
with no stable rail**, and the bit is a real one — median separation 0.19, about 8
molecules at Ω = 40 — not a marginal one. The largest tank essentially never enters
that regime. Restoration has stopped being available and the bit persists on
kinetics: the stage is two relaxation times at γ₀, but `λ(γ_eff) → 0` at γ_c, so a
tank that drains in a handful of stages leaves the state no time to follow the
collapsing landscape. This is why `fuel_hazard.py` matches the large budgets (P3:
quasi-static, loss at γ_eff = 0.42 with the rail intact) and under-predicts the
small ones.

> **Prediction P1 of `fuel_rail_lag.py` is NOT supported, and its statistic was
> badly chosen.** I predicted the carried separation would exceed the rail by more
> at small Φ. Measured one stage before loss the ratio is confounded — δ is
> necessarily small there, so the ratio tracks proximity to loss, not lag. Measured
> instead at *matched* γ_eff, it scatters 0.70–1.65 with **no ordering in Φ** on the
> 12 stored traces per cell. A general "the state lags the rail" claim is therefore
> unestablished; what is established is the narrower and larger-sample fact in the
> table — bits held past γ_c, 97% of trials to 9.5% across the budget range. P4,
> which would have forced outright withdrawal, did not fire.

**What §23.4 leaves open.** The 0.25 exponent gap is *attributed* to holding past
γ_c but not *accounted for* — no calculation here predicts the measured 0.647 from
the quasi-static 0.893 plus a kinetic term. That is T10b-iii-b.

> ⚠ **The kinetic attribution above is withdrawn in §23.5.** The *phenomenon* — bits
> held past γ_c, 97% of trials to 9.5% across the budget range — stands as measured.
> The *explanation* offered for it, that the state cannot follow the collapsing
> landscape because the stage is fixed at γ₀'s relaxation time, is not supported.

### 23.5 The kinetic attribution fails both ways — `fuel_quasistatic.py`

Two arms, Ω = 40, σ_ch/δ* = 0.03, 80 trials, five budgets. Both aimed at §23.4's
claim that the small-tank excess is the state failing to follow the rail.

**Arm A — impose quasi-staticity by hand.** Re-seed (X, Y, B) onto the attractor of
the current γ_eff at every stage, leaving the tank exactly as the chemistry left it.

> **This arm contains a free restoring element and is not a physical result.** It is
> the "harness doing work the chemistry cannot" that cost this project three
> withdrawn results (rule 10), used deliberately as an instrument so the simulation
> obeys the same assumption the hazard integral does. It re-seeds (X, Y, B) only —
> `n_X+n_Y+n_B = Ω` holds and the tank is untouched, so it is not also a reset blank
> pool. Nothing about a real network may be read off it.

| | exponent in Φ |
|---|---|
| hazard integral (§23.4) | 0.8925 ± 0.017 |
| **Arm A, rail-reseeded every stage** | **0.7077 ± 0.030** |
| plain simulation (§23.1) | 0.6474 ± 0.022 |

> **Prediction P1 refuted.** I predicted Arm A would reproduce the integral at ~0.89.
> It lands at 0.7077 — **5.4σ from the integral** and only **1.6σ from the plain
> simulation**, i.e. statistically indistinguishable from doing nothing. Re-seeding
> onto the rail on a median 52 of 52 stages closes essentially none of the 0.25 gap.
> **The quasi-static state assumption is therefore not what separates the integral
> from the simulation**, which is the opposite of what §23.4 assumed, and it also
> means the integral has an error somewhere other than the state.

**Arm B — give the state time to follow.** Stage time set from the current state,
`t_stage = 2/(λ(γ_eff)·f)`, capped at `cap ×` the baseline and the cap swept
(rule 13):

| cap | 1 | 3 | 10 | 30 |
|---|---|---|---|---|
| exponent | 0.6206 ± 0.012 | 0.6437 ± 0.030 | 0.6408 ± 0.021 | 0.6451 ± 0.053 |

Flat across a 30× sweep. **The `cap = 1` cell is a built-in control on the whole
harness** — it is §23's own configuration, and it returns 0.6206 ± 0.012 against
§23.1's 0.6474 ± 0.022 (1.1σ) with depths 7/11/17/28/42 against 7/12/17/27/44. The
harness reproduces §23.

> **Arm B is inconclusive, not evidence, and the reason is a design fault I did not
> anticipate.** The adaptive prescription is within 2.4× of the fixed stage until
> θ > 0.29 and only diverges in the final stage or two (ratio 1.00 → 1.40 → 2.37 →
> 4.40 → 7.69 → 30.8 as θ goes 0.231 → 0.33). So it never tests the bulk of the
> cascade. Worse, longer stages simultaneously let the state follow the rail *and*
> drain the tank faster — the depth drop from cap 1 to cap 3 (7/11/17/28/42 →
> 5/8/12/20/33) is that extra burn — so a flat exponent could be two effects
> cancelling. Arm B neither supports nor refutes the lag; it is recorded because a
> null from an instrument that cannot resolve the effect is exactly what §23.2
> already had to disown once.

**Where this leaves §23 as a whole.** §23.1–§23.2 stand: the exponent is 0.6498 ±
0.019 at the quiet channel, it drifts with σ and with Ω, and `D_fuel` is not a budget
property. §23.3's mechanism survives only as the account of the σ-drift, which the
hazard integral reproduces (−0.121 vs −0.107). §23.4's *measurement* stands — small
tanks hold real bits past γ_c and large tanks do not — but §23.4's *explanation* of
it is withdrawn. **The 0.25 exponent gap between the quasi-static integral and the
measurement is now unexplained**, with the state assumption eliminated as its source.

The surviving suspect is the integral's own hard stop: it kills the bit the moment
θ crosses 1/3, because `initial_counts` gives zero separation past γ_c and so the
measured `q(θ)` is 1 by construction there — while the simulation demonstrably holds
a δ ≈ 0.19 bit through a median 3 of 7 stages in that region. **How to kill:** extend
the `q(θ)` measurement past γ_c starting from the separation the cascade actually
carries there rather than from the (nonexistent) rail, and re-integrate. Note this
costs the integral its parameter-free status — the starting δ becomes an empirical
input — so the result is a one-parameter model and must be reported as one, not as
the absolute test §23.4 ran.

### 23.6 The hard stop is 40% of the error — `fuel_hazard_pastgc.py`

Past γ_c the state is seeded at an **imposed** separation `δ_past` instead of the
nonexistent rail, then the stage runs unchanged. `δ_past` is an empirical input read
off §23.4, so this is a **one-parameter model, not the absolute test §23.4 ran**, and
it is swept rather than fitted. Seeding a state the simulation was *observed* to
occupy is not a free restoring element — no separation is created that the chemistry
was not already carrying, and the conditional hazard measured from it is a real one.

**First, an instrument fault that had to be fixed before any of it could be read.**
`predict_depth` counts whole stages. At Φ/Ω = 25 the predicted depth is ~5 stages, so
one stage is 20% of it — and the four swept models returned **byte-identical
exponents to four decimals**, which is what exposed it. The sweep was working: the
measured past-γ_c hazard falls from `q = 0.429` to `0.302` as δ_past goes 0.10 → 0.25.
The integer depth simply could not move. Interpolating the survival crossing in `ln S`
removes the integrator's own discreteness (rule 13) and makes the dependence visible.
**No published number changes:** the continuous integrator puts §23.4's hard-stop
exponent at 0.8867 ± 0.016 against the 0.8925 ± 0.017 printed there, 0.25σ apart.

| model | exponent in Φ |
|---|---|
| hard stop at θ = 1/3 (§23.4), continuous | 0.8867 ± 0.016 |
| past-γ_c seeded, δ_past = 0.10 | 0.7987 ± 0.022 |
| δ_past = 0.15 | 0.7960 ± 0.022 |
| δ_past = 0.20 | 0.7919 ± 0.022 |
| δ_past = 0.25 | 0.7884 ± 0.022 |
| **measured (§23.1)** | **0.6474 ± 0.022** |

**The hard stop was 39.6% of the error and no more.** Removing it moves the exponent
0.8867 → 0.7919 at the δ_past matching §23.4's observed 0.19–0.23, leaving **+0.1445,
4.6σ** still unexplained.

> **Prediction P1, approximately right.** I predicted ~0.76 from §23.4's stage counts
> (3 extra stages of 7 at the small end, 0 of 46 at the large end, so ~0.13 off a
> log-log fit spanning ln 16). Measured 0.7919 — the estimate was 0.03 low because
> the small-tank depth rose ×1.35 rather than the ×1.43 I assumed. **P2** predicted
> "about half the gap, ~0.11 left"; actual is 40% with 0.145 left. **P3** predicted
> monotone decrease in δ_past: confirmed, 0.7987 → 0.7884 — but the whole effect is
> **0.0103 across a 2.5× sweep**, so the one imported parameter barely matters and
> the model is nearly parameter-free in practice. **P4** did not fire.

**What remains is a shape error, not a scale error.** Predicted ÷ measured depth at
δ_past = 0.20 runs 0.78, 0.72, 0.88, 1.00, 1.10 across the five budgets — still
under-predicting small tanks and now over-predicting large ones. No overall
normalisation fixes both ends.

> **Suspect for the remaining 0.145, named as a suspect** (rule 17 — this thread has
> already withdrawn two mechanisms stated in the confident voice of their
> measurements). The survival product is mean-field: it applies the *unconditional*
> hazard `q(θ)` at every stage. But δ_past is the median separation observed *while
> past γ_c*, and among trials that survive several such stages the separation is
> selected upward — survivors are the ones that happened to keep a big bit, and they
> face a lower hazard than `q(θ)` says. That bias grows with the number of past-γ_c
> stages, which is 3 of 7 at the smallest tank and 0 of 46 at the largest, so it acts
> exactly where the integral under-predicts. **How to kill:** measure `q` conditioned
> on having already survived `k` past-γ_c stages and check whether it falls with `k`.
> If it is flat in `k`, survivor selection is not the residue and the survival
> product's independence assumption is the next thing to test. Either way the fix
> makes the model non-Markovian in θ, which is a different object from the integral
> §23.4 set out to build.

### 23.7 Survivor selection is not the residue — the bit is ground down, not selected — `fuel_survivor_bias.py`

Direct measurement, 3000 trials at each of Φ/Ω = 25, 50, 100, σ_ch/δ* = 0.03. Stages
with γ_eff ≥ γ_c are indexed k = 1, 2, …; the conditional hazard at k is P(lost during
the k-th such stage | survived k−1), and the control for the fact that γ_eff is *also*
rising with k is to divide by the unconditional `q(θ)` that §23.6's integral actually
used at the same θ.

| Φ/Ω = 25 | k=1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| trials entering | 2465 | 1741 | 1048 | 560 | 281 | 148 | 69 |
| conditional hazard | 0.293 | 0.398 | 0.466 | 0.498 | 0.473 | 0.534 | 0.507 |
| ÷ unconditional `q(θ)` | 0.798 | 0.899 | 0.983 | 0.998 | 0.918 | 1.061 | 1.014 |
| mean carried δ | 0.317 | 0.231 | 0.172 | 0.144 | 0.121 | 0.119 | 0.135 |

> **Predictions P1, P2 and P3 all refuted; P4 fires.** I predicted the hazard ratio
> would sit below 1 and *fall* with k, and that the carried separation would *rise* —
> the signature of survivors being the ones that kept a big bit. Both go the other
> way at all three budgets: the ratio climbs from 0.80 to ≈1 and then hovers there,
> and δ **decays monotonically** (Φ/Ω = 25: 0.317 → 0.117; 50: 0.309 → 0.117; 100:
> 0.303 → 0.176). **Survivors are not selected upward, they are ground down** — which
> is what a monostable landscape should do to a bit, and I had the sign backwards.
> Survivor selection is eliminated as the residue.

**What the same data does establish, and it is better news than the suspect it
killed.** The decay is a function of **θ, not of survival history**: pooling all 22
cells with ≥ 50 trials across three budgets and eight survival indices, `δ ~ θ` alone
gives **weighted R² = 0.9177**, and adding `ln Φ` lifts it only to 0.9675 with
coefficient −0.028 (≈23% of the δ range, and the medians are lattice-quantised at
1/Ω = 0.025 which limits how much better either fit could look). At matched θ ≈ 0.35
every budget carries δ ≈ 0.30 regardless of how many stages it took to get there.

So §23.6's closing note — that any repair would make the model non-Markovian in θ —
is **wrong in the useful direction.** The constant `δ_past` was simply the wrong
constant, and the fix is a measured `δ_past(θ)`, which keeps the integral Markovian.

> ⚠ **That last clause is a claim about a fix, not a measurement, and §23.8 refutes
> it.** Imposing the measured `δ_past(θ)` moves the exponent 0.7919 → 0.7880 and
> closes 2.7% of the residue. The δ(θ) *collapse itself stands* — it is a fact about
> the data — but the inference that it repairs the integral does not. This is rule 17
> being violated in the very paragraph that reported rule 17 working: I let a fix
> inherit the credibility of the measurement it sat next to.

**And the residue splits into two errors in two regions**, which reframes what is
left to explain:

| | predicted (δ_past = 0.20) | measured | error |
|---|---|---|---|
| Φ/Ω = 25 | 5.47 | 7 | **−22%**, past-γ_c region |
| Φ/Ω = 400 | 48.23 | 44 | **+10%**, pre-γ_c region (9.5% of trials ever reach γ_c) |

No single correction touches both. `δ_past(θ)` can only act on the small-tank end.

### 23.8 The integral is insensitive to the imposed separation — closing the thread — `fuel_hazard_deltacurve.py`

Imposing §23.7's measured `δ_past(θ)` — 0.30 at θ = 0.34 falling to 0.10 at θ = 0.49,
22 pooled points — in place of the constant 0.20:

| | Φ/Ω=25 | 50 | 100 | 200 | 400 | exponent |
|---|---|---|---|---|---|---|
| δ_past(θ) measured | 5.50 | 8.82 | 14.97 | 26.97 | 48.23 | **0.7880 ± 0.020** |
| δ_past = 0.20 constant | 5.47 | 8.66 | 14.94 | 26.97 | 48.23 | 0.7919 ± 0.022 |
| measured | 7 | 12 | 17 | 27 | 44 | 0.6474 ± 0.022 |

> **P1 and P3 refuted, P2 confirmed, P5 fires.** I predicted the small-tank depth
> would rise to 6–6.5 and the exponent to ~0.75. The depth moved 5.47 → 5.50 and the
> exponent 0.7919 → 0.7880, closing **2.7%** of the residue. P2 was right that large
> tanks cannot move (48.23 → 48.23, unchanged to two decimals).

**Why it cancels, and what that says about the whole repair programme.** The measured
curve raises the imposed separation at low θ (0.20 → 0.30, lowering the hazard) and
lowers it at high θ (0.20 → 0.10, raising it). The two corrections very nearly cancel
across the θ range the integral traverses. This was already visible in §23.6 and I
under-weighted it: a **2.5× sweep of the constant moved the exponent by 0.0103**. The
integral is simply not sensitive to the separation it imposes past γ_c. So §23.6's
39.6% came from *permitting past-γ_c stages at all*, not from what separation they
were given — a structural correction, not a parametric one.

**Three repairs tested, and the thread closes here.**

| repair | region | worth |
|---|---|---|
| quasi-static state imposed by hand (§23.5) | pre- and past-γ_c | 1.6σ from doing nothing |
| removing the hard stop at θ = 1/3 (§23.6) | past-γ_c | **39.6%** of the error |
| correcting the imposed separation (§23.8) | past-γ_c | 2.7% |

What survives is **+0.14 in the exponent, split as −22% at the smallest tank and +10%
at the largest** — two errors in two regions, neither of which is about the imposed
separation, and the larger-tank one sitting in a region where §23.5 already showed the
state assumption is not the problem. Following §22.5's precedent, the residue is left
**named and unexplained** rather than fitted further, and no fourth repair is attempted
in this session.

**What §23 leaves standing.** The measured results need none of this: the budget
exponent is 0.6498 ± 0.019 at the quiet channel (§23.1); it drifts with σ and Ω, so
`D_fuel` is not a budget property (§23.2); burn per stage is flat in Φ while the usable
burn fraction collapses (§23.3); small tanks hold real bits past γ_c and large ones do
not (§23.4); the carried separation decays as a function of θ, R² = 0.9177, not of
survival history (§23.7). The integral's own contribution is that it reproduced the
**σ-drift** of the exponent, −0.121 against −0.107 measured — the single part of
§23.3's proposed mechanism that survived contact with a test — and that it bounded
three candidate errors well enough to eliminate two of them.

### 23.9 Testing the integral's *form* instead of its numbers — `fuel_ensemble_integral.py`

§23.8 closed the thread against further **parametric** repair. This reopens it on a
different axis, because the two surviving errors have the signs and Φ-scalings of the
integral's two **structural** approximations:

```
θ_{d+1} = θ_d + c(θ_d)/Φ        ← deterministic burn
S_{d+1} = S_d · (1 − q(θ_d))    ← independent stages
```

**(a)** Real burn is stochastic, so at stage *d* there is a *distribution* of θ. Fast
burners reach high θ and die, so survivors are biased below the mean; `q(θ)` rises
steeply, so evaluating at the mean overstates the hazard among survivors and the
integral under-predicts. Relative spread of accumulated burn ~ 1/√Φ, so this is
strongest at small tanks — the sign and scaling of the −22%. **(b)** δ carries memory
across stages, so correlated trials die sooner than an independent product — the +10%.

This tests **(a) only**, and changes **no hazard value**: same measured `q(θ)`, same
measured burn. The single change is that θ is propagated as an ensemble of stochastic
trajectories, with per-stage burn drawn from its *measured* mean and spread.

| Φ/Ω | deterministic | ensemble | measured | det/meas | ens/meas |
|---|---|---|---|---|---|
| 25 | 5.47 | 5.72 | 7 | 0.782 | 0.817 |
| 50 | 8.66 | 8.87 | 12 | 0.722 | 0.739 |
| 100 | 14.94 | 14.92 | 17 | 0.879 | 0.878 |
| 200 | 26.97 | 26.49 | 27 | 0.999 | 0.981 |
| 400 | 48.23 | 46.75 | 44 | 1.096 | 1.063 |
| **exponent** | **0.7919 ± 0.022** | **0.7641 ± 0.023** | **0.6474 ± 0.022** | | |

> **The same instrument fault as §23.6, reintroduced by me two sections later.** The
> first pass took `np.median` of *integer* depths for the ensemble while comparing it
> against the deterministic *continuous* crossing. That inflated the headline from
> **19.2% to 27.3%**. Both numbers are recorded; 19.2% is the real one. Catching this
> twice in one thread is the argument for making the continuous crossing the default
> rather than a thing each new script re-derives.

**Verdict: the effect is real, systematic, and small.**

> **P1 refuted, and in the helpful direction.** I predicted ensemble ≥ deterministic
> at every budget. It is higher only at the two smallest tanks and *lower* at the
> three largest (48.23 → 46.75), so both ends move toward the measurement rather than
> one. Caveat against over-reading that: part of the large-tank decrease is
> trajectories fluctuating across the θ_max cutoff, which is a model boundary, not
> physics.
>
> **P2 confirmed, but its escape clause was wrong.** I predicted < 50% of the
> small-tank gap would close *unless* burn is much burstier than Poisson. Burn **is**
> super-Poisson — measured sd/Poisson sd = 2.60, 2.66, 2.85, 3.00, 2.95 — and only
> **16.2%** of the small-tank gap closed anyway. So the dispersion is ~3× larger than
> my estimate assumed and still does not carry the error.
>
> **P3 approximately confirmed:** predicted ~0.72–0.75, measured 0.7641.
>
> **P4 fires on its own pre-registered criterion, and I am not going to argue it
> away.** (§23.10 later shows the same cross-run baseline problem it corrects in
> §23.5 applies to this comparison too — paired against the same run, θ-dispersion is
> 17% of the integral's error.) I wrote that if the exponent moved by less than its fit error, θ-dispersion
> is not the small-tank error. The move is 0.0278 against a combined SE of 0.032 —
> **0.87σ**. Against that: the two integrators share identical `q(θ)` and `c(θ)`, so
> the *paired* difference is far better determined than independent fits imply, and
> the per-budget ratios are monotone across all five (1.046, 1.024, 0.999, 0.982,
> 0.969). Both readings are stated because picking the favourable one is exactly what
> §23 has had to withdraw three times. **The honest position: (a) is a real
> contributor of the right sign and scaling, worth ~19% of the residue, and not
> significant by the test I set myself in advance.**

**Where the 0.2393 gap now stands, measured from the hard-stop integral:**

| | closes |
|---|---|
| removing the hard stop at θ = 1/3 (§23.6, structural) | 39.6% |
| θ-dispersion (§23.9, structural) | 11.6% |
| imposed separation (§23.8, parametric) | 2.7% |
| quasi-static state (§23.5, parametric) | ~0% |
| **unexplained** | **~48.8%** |

The two structural repairs are worth 51%; the two parametric ones essentially nothing.
That is the pattern the whole thread has been converging on, and it points the
remainder at **(b)**, the independence assumption — untested here, and the only
structural approximation left.

### 23.10 A defect in the instrument that does not matter, and one in my reading that does — `fuel_burn_conditioning.py`

An independent analysis of §23's stored data (a separate model, given the numbers and
the eliminations but not my conclusions) raised two objections. Both were checked
against the code and the stored results rather than adopted. **One is a real defect
with no consequence; the other is a real error in how I read my own experiment.**

**Objection 1 — `c(θ)` is a conditional-mean error.** Verified in the code:
`hazard_at` runs `waste.append(...)` for *every* trial including the ones that lost
the bit, and the integral spends that pooled mean on survivors. The stoichiometric
argument is verified too — `f1: X + Y + F → 2B + W` has propensity ∝ n_X·n_Y·n_F,
maximal at δ = 0, so a losing trial rides δ≈0 and burns hard. **This is rule 12
inside the instrument**, and §23.6 and §23.8 both swept δ's effect on `q` while nobody
swept its effect on `c`.

Measured, 2000 trials per cell, burn split by outcome (Φ/Ω = 25):

| θ | 0.231 | 0.269 | 0.308 | 0.327 | 0.384 | 0.423 | 0.461 |
|---|---|---|---|---|---|---|---|
| `q` | 0.001 | 0.012 | 0.129 | 0.301 | 0.419 | 0.467 | 0.494 |
| `c_lost/c_surv` | **4.28** | 2.25 | 1.33 | 1.16 | 1.04 | 1.04 | 1.01 |
| contamination `q·Δc`, as % of `c` | 0.3% | 1.4% | 4.1% | **4.5%** | 1.8% | 1.9% | 0.5% |

> **The premise is confirmed and the consequence is refuted, by an anti-correlation
> neither I nor the analysis anticipated.** Losers really do burn up to **4.3×** more.
> But the pooled mean is contaminated by `q·(c_lost − c_surv)`, and `q` is tiny
> exactly where the ratio is large — at θ = 0.231 the two losers in 2000 burn 4.3×
> more and shift the mean by 0.3%. Where `q` is large the two populations burn
> almost identically. Peak contamination is **4.5%**, and re-integrating with
> `c_surv(θ)` moves the exponent 0.7914 → **0.7918**. My P3 (0.68–0.74) is refuted;
> the independent analysis's sharper prediction (depths 7/10/16/28/52, exponent ~0.71)
> is refuted harder — measured 5.49/8.89/15.02/26.89/49.11. **P5 fires: the burn side
> is worth −0.3%, against θ-dispersion's 19.2%.** A genuine defect worth fixing on
> principle, with no measurable consequence — and another instance of §22.4's
> cancellation motif.

**Objection 2 — §23.5's "1.6σ from doing nothing" was a bad read of a real effect.
This one is upheld and the published conclusion is corrected.**

§23.5 compared Arm A (0.7077 ± 0.030) against **§23.1's separate run** (0.6474 ±
0.022) rather than against the plain arm measured in *its own* run. Re-run at 400
trials with both arms paired:

| | Φ/Ω=25 | 50 | 100 | 200 | 400 | exponent |
|---|---|---|---|---|---|---|
| Arm A (memory removed) | 7 | 10 | 17 | 29 | **52** | 0.7217 ± 0.024 |
| plain, same run | 7 | 11 | 17 | 27 | **45.5** | 0.6325 ± 0.016 |

The paired difference is **+0.0893 ± 0.0286 = 3.1σ**, and it is localised exactly
where the integral over-predicts: nothing at the two smallest budgets, +7% and +14%
at the two largest. §23.5's conclusion that re-seeding "closes essentially none of
the gap" is **withdrawn** — it closes none of the *small-tank* gap and most of the
*large-tank* one, and an exponent-level summary averaged that structure away. This is
rule 9 applied one level up: the axis was the budget, and the statistic was a fit
*across* it.

**The decomposition this yields**, all paired within the 400-trial run:

| | exponent | share of the integral's error |
|---|---|---|
| deterministic integral | 0.7919 | — |
| Arm A — simulation with memory removed | 0.7217 | integral vs memoryless: **+0.0702 (44%)** |
| plain simulation | 0.6325 | inter-stage memory: **+0.0893 (56%)** |

So **inter-stage memory is the single largest term at 56%**, θ-dispersion is 17%, and
the burn conditioning is 0%. Hypothesis (b) of §23.9 is now supported by direct
measurement rather than by elimination.

> **Left unexplained, and flagged rather than smoothed:** Φ/Ω = 50 moves the *wrong*
> way under Arm A (10 against plain 11) and does so at both 80 and 400 trials, so it
> is not sampling noise. Every other budget moves in the predicted direction or not
> at all. No account of it is offered here.


## 24. Two axes of coarse-graining — which coordinate, not how much noise — `coarse_grain_axes.py`

§21 is this project's sharpest result: against an exact CME reference, every
approximation that keeps **any** noise recovers the restoration error exponent to
2–12%, while the deterministic ODE is categorically wrong. The lesson written down
was *noise is what matters*.

§23 then built a stage-level approximation and spent seven subsections failing to
repair it, and the decomposition (§23.10) put its error at 56% inter-stage memory and
17% θ-dispersion. But that model is not merely deterministic — it also **collapses the
state**, keeping the burn fraction θ and discarding δ, the coordinate the bit lives
in. §21's levels (CLE, tau-leaping) keep every species, so they retain the noise *and*
the coordinates, and could never separate the two. That makes the lesson ambiguous in
a way nobody had noticed, and testable.

**Design.** One measured kernel `K(θ, δ)` — per-stage loss probability, outgoing-δ
distribution, and burn — propagated four ways, so nothing differs between cells except
what the propagation retains. Ω = 40, σ_ch/δ* = 0.03, 12 θ × 8 δ × 600 trials,
10⁵ trajectories per cell.

| | δ kept | δ collapsed |
|---|---|---|
| **θ stochastic** | **0.7026 ± 0.019** | 0.7855 ± 0.025 |
| **θ deterministic** | **0.7001 ± 0.022** | 0.8001 ± 0.025 |

*(plain simulation, paired: 0.6325 ± 0.016)*

**From the fully collapsed corner, the two axes are worth wildly different amounts:**

| retained | buys |
|---|---|
| the coordinate δ only | **+0.1000** |
| the noise in θ only | +0.0146 |
| both | +0.0975 |

> **P1 confirmed.** Keeping the coordinate while propagating θ *deterministically*
> (0.7001) beats keeping the noise while discarding the coordinate (0.7855) — 2.6σ
> unpaired, and the paired per-budget ratios are monotone across all five budgets
> (0.991, 0.960, 0.879, 0.840, 0.788), which by rule 18 is the reading that counts
> since all four cells share one kernel. **δ is worth about 7× what θ-noise is
> worth.**
>
> **And once δ is kept, θ-noise is worth nothing at all** — 0.7001 → 0.7026, the
> wrong sign, with per-budget ratios 0.969–0.985 that are systematically *below* 1.
> Adding noise to the bookkeeping coordinate slightly *hurts*.
>
> **P3 confirmed:** not additive. The two singles sum to +0.1146 and keeping both buys
> +0.0975 — another instance of the cancellation motif running through §22.4, §23.8
> and §23.10.

**The refinement this forces, and it reconciles rather than overturns §21.** "Noise is
what matters" cannot be right as stated, because here the noise axis buys +0.0146
against the coordinate axis's +0.1000. But §21's ODE is deterministic in **every**
coordinate, *including δ*. The `δ kept / θ deterministic` cell is deterministic only in
θ, while δ still fluctuates through its measured transition histogram. In both cells
that work, the bit-carrying coordinate keeps its noise; in both that fail, it has been
removed or discarded.

> **So the statement both results support is: noise matters in the coordinate that
> carries the signal, and costs nothing in a bookkeeping coordinate.** θ is
> bookkeeping — it sets how hostile the landscape is, but no decision happens in it.
> δ is where the bit is. §21's cliff removed noise from δ; §23's integral removed δ
> outright. Same failure, seen twice, and it explains the otherwise awkward 17% that
> §23.9's θ-dispersion bought.
>
> **Stated as a suspect, per rule 17.** This is one system, one channel, one Ω, and the
> signal/bookkeeping distinction is clean here in a way it may not be elsewhere. **How
> to kill:** run §21's ladder with a level that keeps full noise but *collapses* a
> signal-carrying species, and one that is deterministic in a bookkeeping species
> only. If the cliff tracks the coordinate rather than the noise there too, it
> generalises; if §21's levels all fail the moment any species is collapsed
> regardless of role, the signal/bookkeeping distinction is an artifact of the fuel
> system's peculiar structure.

**What this does not establish.** The full cell reaches 0.7026 against a measured
0.6325, so **P2 is refuted** — (θ, δ) is not the whole state either, and 0.0701 of
exponent remains outside both coordinates. And the collapsed cells here are **not**
§23.4's integral: this collapse forces δ onto the mean the full model occupies, giving
60.7 stages at Φ/Ω = 400 where §23.4's rail-seeded integral gave 48.2. Absolute depths
do not compare across sections; only the ordering within this one kernel does. A first
draft of the experiment advertised those cells as a built-in control on §23.4, which
a smoke run refuted before the real one was launched.

### 24.1 The kill test, on §21's own ladder — `noise_placement.py`

§24's claim was recorded as a suspect on one system. Its kill test was to run the
coordinate-vs-noise question in §21's ladder, with §21's observable (P(error) from a
start biased by ε, decided at |n_X − n_Y| ≥ 0.8·δ*·Ω), against §21's exact CME
reference. `am_reversible` conserves n_X + n_Y + n_B, so the CLE noise increment sums
to zero and splits by role: **δ = n_X − n_Y is the signal**, **s = n_X + n_Y is
bookkeeping** (the blank pool, `Ω − n_B`). Project the noise, keep the drift full:

| ε | Ω | CME (exact) | CLE full | δ-only | s-only | uniform 11% |
|---|---|---|---|---|---|---|
| 0.20 | 40 | 0.23385 | 0.23728 | 0.24440 | **0** | 0.01405 |
| 0.20 | 60 | 0.17307 | 0.17340 | 0.17608 | **0** | 0.00235 |
| 0.20 | 80 | 0.15607 | 0.15770 | 0.16105 | **0** | 0.00112 |
| 0.20 | 100 | 0.11997 | 0.11975 | 0.12172 | **0** | 0.00015 |
| 0.35 | 40 | 0.09833 | 0.09940 | 0.11115 | **0** | 0.00007 |
| 0.35 | 60 | 0.05195 | 0.05003 | 0.05720 | **0** | 0 |
| 0.35 | 80 | 0.02347 | 0.02228 | 0.02633 | **0** | 0 |
| 0.35 | 100 | 0.01312 | 0.01315 | 0.01553 | **0** | 0 |

δ-only keeps **11%** of the total noise variance; s-only keeps **88%**.
40,000 trials per cell.

> **P1 confirmed, and not as a matter of degree.** `s-only` — 88% of the noise
> variance retained, all of it in the blank pool — reports **exactly 0 in all eight
> cells**, which is the ODE's signature failure. **P4 did not fire:** there is no
> graded middle to restate as a continuum. A model can keep seven-eighths of the
> noise and be as categorically wrong as one that keeps none.
>
> **P2 confirmed.** `δ-only` recovers the exponent — −0.7122 against the exact
> −0.6867 at ε = 0.20, and −2.1719 against −2.2142 at ε = 0.35 — on 11% of the noise.
> **P3 confirmed**, the harness control: the full CLE tracks the CME to 0.2–5%.

**And the arm I added expecting it to pass is the one that makes the argument.** I
built `uniform 11%` — the same *total* variance as δ-only but spread over both
coordinates — to show that amplitude was not the driver, writing in the code that "if
11% everywhere works while 88% in the wrong place fails, the total is not what the
observable is sensitive to." **It does not work.** It is wrong by factors of 17 to 770
at ε = 0.20 and categorically zero at ε = 0.35.

So the sharp statement is neither "noise matters" nor "placement matters", but both
together, and the three arms pin it exactly:

| | δ-noise | s-noise | result |
|---|---|---|---|
| full CLE | 100% | 100% | correct to 0.2–5% |
| **δ-only** | **100%** | 0% | **correct to 2–18%** |
| uniform 11% | 11% | 11% | wrong by 17–770×, or categorically 0 |
| **s-only** | **0%** | **88%** | **categorically 0, like the ODE** |

**The observable is sensitive to the noise in the signal coordinate at its correct
amplitude, and essentially indifferent to everything else.** Removing 88% of the
noise costs 2–18%; removing the remaining 11% costs everything; keeping 11% of the
*signal's own* noise is as fatal as keeping none, because barrier crossing depends
exponentially on the noise amplitude in the crossing direction.

> **What it costs to drop the bookkeeping noise, stated honestly rather than
> rounded.** δ-only is not free: it runs **+1.5% to +4.5% high at ε = 0.20** and
> **+10% to +18% at ε = 0.35**, systematically over-estimating and worsening with the
> bias. That is inside §21's own 2–12% band at the smaller bias and outside it at the
> larger. Blank-pool noise is *mostly* discardable for this observable, not exactly
> discardable, and the error grows with the barrier.

**§24 generalises.** The signal/bookkeeping distinction was not an artifact of the
fuel network: it reproduces in the system §21 was measured on, with §21's observable
and an exact reference, and more sharply there than in §23 — categorical rather than
7×. §21's measurements stand as printed; **"having noise at all" is superseded by
"having the signal coordinate's noise, at its own amplitude."**

> ⚠ **What the projection arms are, and are not — a scoping note added after reading
> the coarse-graining literature.** These arms zero a noise component and leave the
> drift untouched. That is **not** a valid reduced model in the Mori–Zwanzig sense.
> The rigorous construction (Zwanzig projection; see
> [arXiv:2512.03706](https://arxiv.org/abs/2512.03706), which derives it for
> underdamped Langevin dynamics) returns a reduced system with a **generalized force**
> and a **state-dependent diffusion** — the eliminated degrees of freedom come back as
> modified drift and memory, not as deletion.
>
> So §24 and §25 measure **what naive discarding costs**, which is what a practitioner
> actually does when they drop a fluctuation term. They do **not** measure what a
> correctly reduced model achieves, and the two claims must not be conflated. §24.1b
> tests the difference. This also corroborates, from a published direction, the
> thermodynamic-consistency objection that led §25 to drop entropy production as
> ill-posed for a projected CLE.

### 24.1b The residual IS a generalized force — measured, but not in closed form — `zwanzig_correction.py`

If §24's arms are naive discarding, the Zwanzig picture predicts the cost should be a
**generalized force**: the eliminated pool comes back as modified drift. §24.1a's
residual is the target — `δ-only` overshoots P(error) by +13.2% at γ = 0.05 falling
monotonically to +0.4% at γ = 0.45.

**It cannot be a diffusion correction, and that is settled by construction.** The
`δ-only` projection preserves `a − b` exactly, so **δ's own diffusion is unchanged**.
Whatever is missing must act on the drift.

**Measured directly, rather than derived.** Running full and `δ-only` from an
identical state over a window of a few pool-correlation times and differencing the
mean δ-increment gives the missing force with no formula and no sign ambiguity:

| γ | missing force (window 2.0) | §24.1a residual |
|---|---|---|
| 0.05 | **+0.0483** | +13.2% |
| 0.30 | −0.0015 | +3.1% |
| 0.45 | −0.0041 | +0.4% |

It is **positive and large exactly where the residual is large**, and ~0 where the
residual vanishes. It also *builds over ~1/λ_s* (+0.018 at window 0.5 → +0.048 at
2.0), which is what a memory term does and a Markovian one does not.

> **So the Zwanzig reading is confirmed in kind: what naive discarding costs is a
> drift term, not a fluctuation term.** That is the useful half.

**Both closed forms for it fail, and rule 16 is what caught them.**

> **P1 refuted, structurally.** I predicted the curvature term
> `½·∂²b_δ/∂s²·Var(s)`. It is **identically zero**: `b_δ` is *exactly linear in s* —
> measured second derivative 4.4×10⁻¹⁶, and a sweep across s±3 is perfectly straight.
> The AM drift is bilinear in (δ, s), so this term vanishes by the network's own
> structure. The `+force` arm was bit-identical to `δ-only`, which is how it was
> caught.
>
> **The cross-correlation term is worse.** `∂b_δ/∂s · D_δs/λ_s` is the natural next
> candidate, and `D_δs = +2.25` is genuinely nonzero. Against the measured force it
> is **wrong in sign at γ = 0.05 (−0.085 predicted, +0.048 measured)** and **11× to
> 140× too large at γ = 0.30 and 0.45**. Had I fitted a coefficient instead of
> computing the absolute prediction, it would have looked like a success at two of
> three γ.

**Where this leaves the naive-vs-correct question.** The direction is established —
replace the deleted noise with a generalized force, not with nothing — but **I cannot
supply that force in closed form**, and per rule 16 I am not going to fit one. That
the measured force *grows with the correlation window* points at my own P4 fallback:
a correct reduction here may need the full memory kernel rather than any Markovian
drift correction, which this experiment cannot construct.

> **Scope, stated because a positive result here would have been over-read.** None of
> this touches the *categorical* failure. `s-only` deletes the noise in the coordinate
> the observable is defined on, and no generalized force on a deterministic coordinate
> produces barrier crossings. A proper projection may repair the few-percent residual;
> nothing suggests it repairs an exact zero.

### 24.1c Learning the closure — and what a training R² of 0.99 is worth in the tail — `learned_closure.py`

§24.1b measured the missing force but could not write it down. So it was **learned**,
following the coarse-graining literature's method (arXiv:2512.03706 fits reduced drift
and diffusion with random features and shallow networks) — with the one thing that
literature usually lacks: **an exact reference to score against.**

The design exists to stop the network being tuned into looking right:

- **Trained** on the local missing force, measured by paired full-vs-projected runs
  under *common random numbers*, on 500 states sampled from real trajectories.
- **Scored** on exact P(error) — a tail probability, **never in the loss**.
- **Controls:** a linear model in (δ, s) alongside the MLP, plus the uncorrected
  `δ-only` and `full` arms. γ = 0.05, Ω ∈ {40, 60, 80}, 40,000 trials.

**At the longest window tested (τ = 4 ≈ 4.4 pool-correlation times):**

| Ω | exact | full | δ-only | linear | **MLP** | MLP train R² |
|---|---|---|---|---|---|---|
| 40 | 0.027095 | +2.4% | +4.3% | +2.1% | **+1.6%** | — |
| 60 | 0.016144 | +7.2% | +16.9% | +16.1% | **+9.5%** | 0.992 |
| 80 | 0.009425 | +6.1% | +9.6% | +8.0% | **+4.0%** | 0.990 |

> **A learned closure does beat naive deletion.** At Ω = 80 the MLP takes δ-only's
> +9.6% down to +4.0%; averaged over all nine cells, δ-only is +14.1% and the MLP
> +10.2%. So the §24.1b direction holds: replacing the deleted noise with a force
> helps, and it helps more than either closed form did.
>
> **P4 refuted — the force is genuinely nonlinear.** I predicted linear would do as
> well. It does not: MLP train R² is 0.98–0.99 while linear's swings from **0.064**
> to 0.91, and averaged over cells linear (+22.2%) is *worse than doing nothing*
> (+14.1%). Reporting only the MLP would have been dressing up a straight line; here
> the straight line genuinely fails.
>
> **P3 confirmed, and it is the point of the section.** The MLP fits its training
> target at **R² = 0.99** and still misses the tail by 4–10%. A closure that
> reproduces the local force to one part in a hundred buys roughly *half* the
> residual, not 99% of it. **Local accuracy is not tail accuracy** — which is rule 16
> restated for universal approximators, and the reason the scoring target was fixed
> in advance and kept out of the loss.

**A′ — the Markovian consistency test, and it does not pass.** A well-defined
Markovian closure must give the same answer whichever window it was estimated at:

| | τ=1 | τ=2 | τ=4 |
|---|---|---|---|
| Ω=40, MLP | +18.1% | +13.2% | +1.6% |
| Ω=60, MLP | +20.5% | +9.2% | +9.5% |
| Ω=80, MLP | +7.2% | +9.0% | +4.0% |

The fitted closure depends strongly on its window, improving monotonically toward
τ = 4 in five of six model×Ω pairs. **T13-a is therefore not closed.** Two readings
survive and this run cannot separate them: either no Markovian closure exists (the
τ-dependence is real), or τ = 4 is simply the first adequately-converged window and
shorter ones undersample an effect that needs several correlation times to develop.
**τ = 8 was not run for the closure**, and that is the missing cell — the jump from
+35.4% to +2.1% between τ = 2 and τ = 4 at Ω = 40 is large for a quantity that is
supposed to be converging.

> **Sampling floor, stated because several comparisons here sit near it.** At 40,000
> trials and p ≈ 0.0094, the relative standard error is ≈ 5%. So single-cell
> differences of a few percent — including "MLP beats the full CLE" at Ω = 40 and 80
> — are **not** individually significant. The arms share seeds, so their differences
> are better determined than independent errors imply, but the claims above rest on
> the *pattern across nine cells*, not on any one of them.

### 24.1a The competing explanation I failed to name, and what the γ sweep does to it

An independent review of §24.1 raised an account of these results that is more
standard than mine and that I had not written down. `am_reversible`'s saddle has
`λ_antisym = (1−2γ)/3 = +0.133` at γ = 0.30 (δ, unstable) against
`λ_sym = −(1+2γ) = −1.60` (s, strongly stable) — a timescale ratio of **12** — and
`_setup` starts the system *on* the s-nullcline `b* = γ/(1+γ)`. So s may simply be a
fast stable variable **slaved** to δ, and §24.1 may be the textbook large-deviation
result that an escape rate is set by diffusion along the unstable direction, with
transverse diffusion entering at subleading order. That is a property of the saddle's
geometry — of the *system* — and it is exactly the collapse §24 named as its own
failure mode. **No observable-swap test touches it.** The discriminating axis is γ,
because `|λ_sym/λ_antisym| = 3(1+2γ)/(1−2γ)` spans a factor of 15 over the usable
range.

| γ | timescale ratio | `s-only` | `δ-only` cost |
|---|---|---|---|
| 0.05 | 3.7 | **exactly 0 in 3/3** | +13.2% |
| 0.15 | 5.6 | **exactly 0 in 3/3** | +8.1% |
| 0.30 | 12.0 | **exactly 0 in 3/3** | +3.1% |
| 0.45 | 57.0 | **exactly 0 in 3/3** | +0.4% |

**The review is right about one half and wrong about the other, and the split is
clean.**

> **The cliff is not slaving.** `s-only` is categorically 0 at *every* γ, including
> γ = 0.05 where the timescale separation is only 3.7 and the pool holds
> `n_B ≈ 1.9` molecules — a near-empty pool that gates the recruitment reactions
> multiplicatively and cannot be Gaussian-slaved at all. Adiabatic elimination
> predicts `s-only` should recover as γ → 0. It does not, over a 15× range.
>
> **But the residual is slaving, and that is a better account than I had.** `δ-only`'s
> cost — the price of throwing the pool's noise away — falls **monotonically with the
> timescale separation**, 13.2% → 8.1% → 3.1% → 0.4%. That is precisely
> subleading-order transverse diffusion, shrinking as the separation grows. §24.1
> reported this residual as "2–18%, worsening with the barrier" and left it
> unexplained; it now has a mechanism, and it is the reviewer's, not mine.

So the two effects are different things: **the categorical failure is about which
subspace carries the observable; the few-percent residual is about timescale
separation.** §24.1's claim survives the sharpest attack available to it, and comes
back with its error term explained.

> **Also checked, from the same review: a real bug with no consequence here.**
> `run_projected` advances the clock only for accepted steps (`t[upd] += dt`), so a
> step rejected for negativity freezes that trajectory's time — an arm-dependent clock
> distortion, and a rule-10 free-restoring-element instance (rejection acts as a
> reflecting wall the chemistry lacks). Measured, the rejection rate is ≤ 7.7×10⁻⁴ and
> exactly 0 for two of the four arms, so it cannot account for effects of a few
> percent, let alone the cliff. **The mechanism is real and would matter in a regime
> with more boundary contact**; it is recorded rather than fixed because fixing it
> would change §21's published numbers, and it is bounded here.

> **§24.2 sharpens "coordinate" to "subspace".** At n = 2 the difference subspace is
> one-dimensional, so "the signal coordinate" and "the signal subspace" are the same
> object and this section cannot tell them apart. At n = 3 they differ, and the
> subspace-level split is what survives.

### 24.2 At n = 3: the split is by SUBSPACE, and one arm was invalid — `noise_placement_nwinner.py`

n = 2 has exactly one signal direction and one bookkeeping direction, so §24.1's
split had nowhere to hide. At n = 3 the noise decomposes into three orthogonal free
directions: **d1 = (1, −½, −½)** champion-vs-rivals, **d2 = (0, 1, −1)**
rival-vs-rival, and the committed-vs-blank pool. `n_winner_reversible`, γ = 0.60·γ_c(3)
= 0.1214 matching §21.4's convention, exact CME reference via
`cme.splitting_probability`, 40,000 trials per cell.

| ε | Ω | CME | full | signal-only | bookkeeping-only | rivals-only |
|---|---|---|---|---|---|---|
| 0.25 | 30 | 0.11530 | 0.11918 | 0.14780 | **0** | **0** |
| 0.25 | 45 | 0.10654 | 0.10723 | 0.12392 | **0** | **0** |
| 0.25 | 60 | 0.09160 | 0.09070 | 0.10378 | **0** | **0** |
| 0.25 | 80 | 0.04120 | 0.04200 | 0.05025 | **0** | **0** |
| 0.40 | 30 | 0.02858 | 0.03195 | 0.05378 | **0** | **0** |
| 0.40 | 80 | 0.00075 | 0.00068 | 0.00193 | **0** | **0** |

Variance retained: signal-only 19–20%, bookkeeping-only 80%, rivals-only 6%.
Nothing unfinished in any cell.

> **P1 confirmed.** `bookkeeping-only` — 80% of the noise, all in the blank pool —
> is **categorically 0 in all eight cells**, reproducing §24.1's `s-only` at n = 3.
> **P2 confirmed on the exponent, with a caveat on the probabilities:** `signal-only`
> gets the exponent to 6.6% (−1.0239 against the exact −0.9608) and 8.5% (−3.3839
> against −3.6995), but the probabilities themselves run **13–28% high at ε = 0.25
> and 88–157% high at ε = 0.40** — substantially worse than n = 2's 2–18%, and
> worsening with the barrier as it did there. Dropping the pool's noise costs more
> when there are more rivals.

> ⚠ **`decision-only` is INVALID as constructed and P3/P4 are unresolved.** Its
> values were erratic and non-monotone in Ω — 0, 0.0537, 0, 0.0194 — and the cause is
> not sampling. With d2's noise zeroed, X2 − X3 evolves *deterministically* from its
> initial value, so the arm's entire behaviour is set by an integer-rounding parity
> in the start state:
>
> | Ω | 30 | 45 | 60 | 80 |
> |---|---|---|---|---|
> | rivals at start | tied | differ by 1 | tied | tied |
> | `decision-only` | 0.00000 | 0.05370 | 0.00000 | 0.00003 |
>
> When the rivals start tied they stay tied forever, no rival can break away, and X1
> wins by construction. **This is the same class as the floor-division artifacts
> already catalogued here** (§12's non-monotone CME error, §17's MI integer bias) and
> no number of trials fixes it. **How to kill properly:** seed the rivals with a
> deliberate fixed asymmetry independent of Ω parity, so the arm is not decided by
> rounding. Until then, whether the signal subspace decomposes *within itself* is
> untested. `rivals-only` is valid but near-definitional — with d1's noise removed,
> X1 versus the rival mean is deterministic and X1 never loses.

> **A second error caught before it produced anything, and worth recording.** The
> first n = 3 run copied γ = 0.30 from §21's AM defaults. But γ_c falls with n —
> γ_c(3) = 0.2023 — so that run measured a network with **no landscape at all**:
> `landscape_width` = 0, the threshold collapsed to its floor of 2, the
> bookkeeping-only state sat at a perfectly symmetric [9.2, 9.2, 9.2], and the "CME
> reference" was the probability of a 2-count lead in a system with no attractor to
> restore toward. It would have read as a clean set of results. The experiment now
> takes `--gamma-frac` of γ_c like §21.4 and refuses to run if the width is zero.

**What §24.2 settles.** The signal/bookkeeping split holds at n = 3: 80% of the noise
is in the pool and removing it is categorical, while the 19% in the difference
subspace carries the exponent. So §24's claim generalises to a **multi-dimensional
signal subspace** — and its wording should be "the signal subspace", since n = 2
could not distinguish that from a single coordinate. Whether the subspace decomposes
further is open, with the arm that would test it named and currently broken.

### 24.3 The signal subspace does not decompose — `--rival-skew` — T11-REFINED-b

§24.2's `decision-only` arm was decided by an Ω parity. The repair forces
`max(rivals) − min(rivals) = skew` exactly in the start state, independent of Ω, and
applies it to every arm and to the CME reference alike. Run at skew = 2 and skew = 4,
40,000 trials per cell.

**Ratios to the exact CME, paired on an identical start state** (rule 18 — these are
matched-start comparisons, which is what the following claim rests on):

| ε | Ω | full | signal-only | decision-only (skew 2) | decision-only (skew 4) |
|---|---|---|---|---|---|
| 0.25 | 30 | 1.026 | 1.270 | 0.612 | 0.705 |
| 0.25 | 45 | 1.003 | 1.256 | 0.589 | 0.657 |
| 0.25 | 60 | 1.002 | 1.135 | 0.550 | 0.596 |
| 0.25 | 80 | 1.043 | 1.296 | 0.544 | 0.590 |
| 0.40 | 45 | 0.985 | 2.482 | 0.799 | 0.838 |
| 0.40 | 60 | 1.091 | 1.975 | 0.700 | 0.833 |
| 0.40 | 80 | 1.047 | 2.454 | 0.818 | 0.793 |

> **P6 confirmed.** The parity is gone: `decision-only` is smooth and finite at every
> Ω, and `full` and `signal-only` barely move, so the skew is a neutral repair rather
> than a different experiment.
>
> **P7 confirmed in direction.** `decision-only` **never recovers** — it
> under-estimates in all fifteen surviving cells across both skews, averaging
> **0.666** of the exact answer at skew 2 and **0.716** at skew 4. Keeping only the
> champion-vs-rivals direction loses roughly a third of the failure probability. The
> order-statistic mechanism I predicted is *consistent* with this — the champion
> loses to the best of two noisy rivals, which sits higher than the best of two whose
> separation is set by drift alone — but nothing here isolates that mechanism from
> other consequences of deleting d2, so it stays a suspect (rule 17).
>
> **P8 partially fails, and the magnitude is therefore not pinned.** The ratio drifts
> systematically upward with the skew — +0.063 averaged over the four cells present
> in both runs. The direction of P7 is robust to that; its size is not.

**Two limits on these runs, stated rather than worked around.**

> **Exponents are not fittable here.** `setup_skewed` absorbs its divisibility
> remainder into the *margin*, so realised ε wobbles ±11% around target (0.912–1.115)
> and the CME reference is non-monotone in Ω — precisely the corruption
> `approximation_hierarchy_nwinner.setup`'s docstring was written to prevent, which I
> reintroduced while fixing a different artifact. Only the paired per-cell ratios
> above are used. **The fix** is to absorb the remainder into the *skew* instead,
> letting it vary 2–3 while the margin stays exact: P8 shows the skew dependence is
> weak and §24.2 shows the margin dependence is not.
>
> **One cell is excluded on a stated criterion, not by inspection.** At ε = 0.40,
> Ω = 30, skew = 4 the realised ε is 14% high and **the full-CLE control is itself
> 1.53× the CME** — the level is outside its own validity there, independently of any
> projection. The gate applied is |full/CME − 1| < 0.25, which drops that one cell and
> keeps 15 of 16.

**Where this leaves §24.** The split is by **subspace**, and the subspace is not
further decomposable: removing the pool's noise (80% of the variance) is categorical,
removing one of the two signal directions is not categorical but costs a third of the
answer, and only the full difference subspace reproduces it. For a simulation that
means the saving is bounded — you may discard the bookkeeping noise entirely and get
the exponent, but there is no further cheap truncation *inside* the signal directions.


## 25. The required subspace belongs to the OBSERVABLE, not the system — `observable_dependence.py`

§24 asserted that the required noise subspace is a property of the *(system,
observable)* pair. That was an **inference**, and a load-bearing one: in
twenty-four sections every observable this project had ever tested was the
restoration error probability. If the required subspace were the same for every
question, §24 would collapse to "AM has one stiff direction" — a fact about this
network, not about simulation.

**Design.** Identical machinery to §24.1 — same projected-noise CLE arms, same drift,
same network, same start states — with **two observables read off the same
trajectories**:

- **P(error)** — a *tail* event: δ must cross against its drift.
- **MFPT** — mean time to first passage at |n_X − n_Y| ≥ thr, a *bulk* quantity
  carried by the drift.

Both come out of **one `first_passage` solve**, so the references are exact and paired
by construction. Projecting the noise does not rescale time, so the arms share a clock
(rule 11). And this is not a definitional trap: under `s-only` the signal still
reaches the threshold, carried by the drift every arm retains, so MFPT is an ordinary
finite number whose correctness the run decides.

| arm | variance | **P(error)** | **MFPT** |
|---|---|---|---|
| full CLE | 100% | −0.6% mean | +1.4% mean |
| δ-only | 11% | +8.1% mean | +4.6% mean |
| **s-only** | **88%** | **exactly 0 in 8/8 cells** | **within 6.1% everywhere** |
| **uniform 11%** | **11%** | **−99.0%, exactly 0 in 3/8** | **+0.1% mean, within 15.1%** |

Ω ∈ {40, 60, 80, 100}, ε/δ* ∈ {0.20, 0.35}, 40,000 trials per cell.

> **P1, P2 and P3 all confirmed, and P2 is the point.** The **same arm**, on the
> **same trajectories**, is categorically wrong about one observable and correct to a
> few percent about another. `uniform 11%` is sharper still: wrong by two orders of
> magnitude on P(error) while reproducing MFPT to **+0.1% on average**. The asymmetry
> is one-sided as predicted — every arm is adequate for the bulk quantity, only some
> for the tail — which is a difference in *requirement*, not a swap.

> **P5's confound, controlled rather than mentioned.** The unconditional MFPT averages
> over both outcomes and wrong-way trajectories take longer, so `s-only` — which
> produces no wrong-way trajectories at all — could inherit its MFPT agreement from
> its P(error) failure. Conditioning on the correct outcome puts the arms on matched
> populations: `s-only / full` is **1.041 mean (0.987–1.167)** and `uniform / full` is
> **1.064 (0.983–1.248)**. There is no exact reference for the conditional quantity —
> it needs a Doob h-transform that is not implemented here — so it is an arm-to-arm
> comparison only. **The conclusion survives both readings.**

**What this settles.** The required noise subspace is **not a property of the
system**. The same chemistry, at the same parameters, demands the signal subspace for
one question and essentially nothing for another. §24's framing is now measured rather
than inferred, and the statement that survives is:

> **Which degrees of freedom a simulation may discard is a property of the question
> being asked of it, not of the system being simulated.** A model can be certified
> against one observable, retain 88% of the system's fluctuations, and still be
> categorically wrong about a different observable of the same run — and no limit
> theorem licensing the approximation will say so, because Kurtz's theorem is true and
> covers neither case (§21).

**Scope, and it is the WEAK form of the claim.** An independent review of this design
made the sharp version of the objection, and it stands: the only outcome that would
establish observable-dependence *strongly* is an observable requiring **span(s)** or
**span(δ)⊕span(s)**. An observable requiring the **empty** subspace — any
drift-dominated mean — is not evidence that requirements differ by question so much as
evidence that means are means. MFPT is drift-dominated, so §25 demonstrates:

- ✅ the required subspace is **not** a fixed property of the system (P(error) needs
  span(δ); MFPT needs nothing) — this is measured and it is real;
- ❌ **not** that some observable of this system needs the *pool's* noise. No
  observable tested here does.

> ⚠ **Superseded by §25.2:** the requirement varies in **both** directions. The
> two-target race is an observable of this same network whose requirement points at
> `span(s)`, with `span(δ)` alone giving exactly 0 in 4/4 cells — the mirror of the
> table above.

So the honest statement is **"the requirement varies by observable, downward"**. The
strong form needs an observable whose requirement points at s, and the review names
two candidates with exact references and no definitional shortcut: **Var(T)** at first
passage (the ODE gives 0, so it is a pure-noise quantity like P(error), and the two
competing mechanisms — jitter from crossing-direction diffusion versus jitter from the
*rate of advance*, which n_B sets through the recruitment propensities — are both
credible and neither is guessable), and a **two-target race** absorbing on
`|n_X − n_Y| ≥ thr` **or** `n_B ≤ m`, tuned to ~50/50, which `cme.splitting_probability`
already supports with no new code. That is **T12**.

Two further limits recorded rather than smoothed:

> **`Var(n_B)` is deliberately not among the observables.** Under `s-only` the noise
> *is* n_B's noise, so the arm is right by construction and the answer can be written
> down without running anything — the exact mirror of §24.2's `decision-only` trap. It
> has one legitimate use, as a **unit test of the projection code**, and none as
> evidence.
>
> **Entropy production was considered and dropped as ill-posed, not merely trappy.**
> For `dn = b dt + σ dW` with `D = σσᵀ` singular and drift components outside
> `range(D)`, the process is not absolutely continuous with respect to its time
> reversal and the EP rate is formally infinite. δ and s are both molecule counts and
> both even under time reversal, so there is no parity trick to rescue it. Anything
> computed by accumulating `ln(a_f/a_r)` along a projected path would be the EP the
> *jump* process would have paid — a state functional that every arm reproduces
> because every arm keeps the drift, i.e. right for reasons unrelated to noise
> placement.

### 25.1 T12's best candidate, and the strong form is not found — `timing_jitter.py`

`Var(T)` at first passage was the sharpest candidate for an observable requiring the
pool's noise: it is a **pure-noise** quantity (the deterministic limit gives exactly
zero, so there is no drift-dominated leading order to dilute against), and it has two
credible mechanisms neither of which is guessable — jitter from diffusion along the
**crossing direction** (→ `span(δ)`), or jitter from fluctuations in the **rate of
advance**, since every recruitment propensity carries `n_B` and the pool holds only
~9 molecules at Ω = 40 (→ `span(s)`). Exact reference added as
`cme.first_passage_moments` (`Qtt m2 = −2T`), pinned against sampled jump
trajectories by `test_first_passage_variance_matches_the_ssa` — necessary because a
wrong factor in that recursion still yields a positive, plausible variance.

**Verdict bands were fixed before looking at output** (recovers ≤ 25% off, fails
> 3× off), because §24.1's exact zeros will not recur and "categorical" must not be
decided by eye.

| arm | Var(T) ratio | verdicts | P(error), same trajectories |
|---|---|---|---|
| full | 1.013 [0.977–1.039] | 8/8 recovers | fine |
| **δ-only** | **1.183 [1.116–1.241]** | **8/8 recovers** | +8.1% |
| **s-only** | **0.160 [0.109–0.199]** | **0/8, 8 fails** | **exactly 0 in 8/8** |
| uniform 11% | 0.148 [0.075–0.380] | 7 fails, 1 partial | 0 in 3/8 |

> **P1 and P2 refuted, P3 confirmed. The strong form is not found.** I predicted
> neither arm would recover alone, making Var(T) an observable requiring
> `span(δ)⊕span(s)`. Instead **`δ-only` recovers it alone** and `s-only` captures only
> 16% of the variance. The P4 control rules out bimodality as the explanation:
> `Var(T | correct)` for `s-only` is 0.198 of full, so the shortfall is genuine
> missing jitter, not the absence of slow error paths.
>
> **So `Var(T)` is a third δ-observable**, and the honest reading is the one the
> review pre-supplied: three observables tested, all of them saddle- or
> drift-dominated. That does **not** establish "one stiff direction" — it establishes
> that I have now asked the same question three ways.

**One thing worth keeping from a negative result.** The pool's noise is not
*irrelevant* to Var(T), just insufficient: `s-only` captures 16–20% of it, and
`δ-only` overshoots by 18% where it overshoots P(error) by only 3.1% at the same γ
(§24.1a). **So the pool's contribution is roughly 6× more important to the timing
jitter than to the error probability** — observable-dependence in the *magnitude*,
even where the categorical requirement is the same. That is a weaker claim than T12
set out to make and it is the one the data supports.

**T12's remaining candidate, with a caveat the review did not raise.** The two-target
race — absorbing on `|n_X − n_Y| ≥ thr` **or** `n_B ≤ m`, tuned to ~50/50 — is still
untested and `cme.splitting_probability` supports it with no new code. But it may sit
closer to the definitional line than it appears: under `δ-only` the pool carries no
noise of its own, so whether it reaches `m` is driven only by δ's fluctuations feeding
the drift. That is not *forced* the way `Var(n_B)` is, but it is one step removed, and
the pre-committed criterion should be written before it runs.

### 25.2 The race, and the requirement reverses — `two_target_race.py`

Two absorbing targets, one in each coordinate: **decision** at `|n_X − n_Y| ≥ thr`,
**pool** at `n_B ≤ m`, with `m` chosen per cell to minimise |P(pool first) − ½| against
the exact CME. At Ω = 40 that is m = 5 against a starting pool of 9 — **the pool must
fluctuate down by 4 to win**, so it is a genuine fluctuation event, not a formality.
Criterion fixed before any output: recovers |Δp| ≤ 0.05, fails |Δp| > 0.20.

| Ω | m | n_B(0) | P(pool) exact | full | δ-only | s-only | uniform 11% |
|---|---|---|---|---|---|---|---|
| 40 | 5 | 9 | 0.4396 | 0.4970 | **0.0000** | 0.5727 | **0.0000** |
| 50 | 7 | 12 | 0.4130 | 0.4522 | **0.0000** | 0.5216 | **0.0000** |
| 60 | 10 | 14 | 0.5993 | 0.6099 | **0.0000** | 0.6648 | **0.0000** |
| 70 | 12 | 16 | 0.5595 | 0.5651 | **0.0000** | 0.6229 | **0.0000** |

Verdicts: full 3 recovers / 1 partial (P3 control passes); **δ-only 4/4 fails,
categorically zero**; s-only 4/4 partial at +0.063 to +0.133; uniform 11% 4/4
categorically zero. Simultaneous hits on both targets: 612 in ~640,000 trajectories,
≈0.1%, assigned to the decision target.

> **P1 confirmed and it carries no weight, exactly as declared in advance.** `δ-only`
> reports the pool never winning. But the pool target *is* a pool-fluctuation event
> and `δ-only` has no pool fluctuations, so that arm answers by construction. It is
> the mirror of §24.1's `δ-only` anchor and I am not counting it as evidence.
>
> **P2 not met on its own terms:** `s-only` lands *partial*, not *recovers*, in all
> four cells. The band was fixed in advance and it does not get moved now.
>
> **What is informative is the bias, and it is 36 standard errors.** `s-only`
> over-predicts the pool winning by **+0.06 to +0.13** (SE ≈ 0.0025 at 40,000
> trials). Without δ's fluctuations the decision target is reached less readily, so
> the pool wins too often. That is a real, quantitative contribution from the signal
> noise to an observable whose *necessary* ingredient is the pool noise — and it is
> not definitional in either direction.

**The mirror, which is the actual result.** Set the two observables of the same
network, same arms, side by side:

| | needs | other subspace alone | cost of dropping the other |
|---|---|---|---|
| **P(error)** | `span(δ)` | `span(s)` → **exactly 0**, 8/8 cells | +3.1% (γ = 0.30, §24.1a) |
| **race** | `span(s)` | `span(δ)` → **exactly 0**, 4/4 cells | +9% to +13% |

**T12's strong form is established, with its caveat attached.** There is an
observable of this system whose requirement points at the pool, and it is the
mirror-image of the one that points at the signal. §25's conclusion — that the
requirement varies by observable "downward" only — is superseded: **it varies in both
directions.** The honest qualification is that each arm's categorical zero is
semi-definitional in its own observable, so the weight rests on the *quantitative*
halves: dropping the pool's noise costs P(error) 3.1%, and dropping the signal's costs
the race 9–13%. Those are the two numbers no construction forces.

> **And the earlier negatives now read differently.** §25.1 concluded that three
> observables all living in `span(δ)` meant I had asked one question three ways. That
> was right, and the reason is visible here: P(error), MFPT and Var(T) are all
> first-passage functionals *of δ*. The race is the first observable in this project
> whose absorbing set is defined on a different coordinate — which is what it took,
> and which supports the §25.1 worry that the limitation was in the **instrument**
> rather than the chemistry.

### 25.3 Stress-testing the only usable rule the arc produced

Across §25–§25.2 one predictive summary fits every result:

> **The required subspace is the one the observable is a functional of.**

P(error), MFPT and Var(T) are first-passage functionals of δ and all three need
`span(δ)`; the race's pool target is a functional of `n_B` and needs `span(s)`. If
that holds it is the arc's practical payoff — a simulation designer could apply it by
inspection, without running a projection at all. So it is worth attacking.

**Where it should break.** Every recruitment propensity carries `n_B`
multiplicatively, and the pool is `b* = γ/(1+γ)`: ~9 molecules at γ = 0.30 but **~2 at
γ = 0.05**. At two molecules the pool's relative fluctuation is order one and it gates
the rate of advance directly, so a δ-functional like Var(T) might be forced to require
`span(s)` — making the required subspace a function of the *parameters*, not just the
observable's definition. Predicted in advance; measured across γ:

| γ | pool at Ω=40 | `δ-only` on Var(T) | `s-only` share |
|---|---|---|---|
| 0.05 | ~2 | **1.160, 3/3 recovers** | 0.140 |
| 0.15 | ~5 | **1.168, 3/3 recovers** | 0.129 |
| 0.30 | ~9 | **1.183, 8/8 recovers** | 0.160 |
| 0.45 | ~12 | **1.063, 3/3 recovers** | 2.362 |

> **P5 refuted, and the rule survives its sharpest available stress.** Shrinking the
> pool to two molecules does *not* make its noise necessary for a δ-functional:
> `s-only`'s share of Var(T) is 0.140 at γ = 0.05 against 0.160 at γ = 0.30 —
> unchanged, and if anything slightly lower. **`δ-only` recovers Var(T) at every γ
> from 0.05 to 0.45.** The rule is not an artifact of the one γ it was found at.

**An unpredicted finding worth more than the confirmation: `s-only`'s error changes
SIGN with γ.** It *under*-predicts Var(T) by ~7× at γ ≤ 0.30 and *over*-predicts by
2.4× on average at γ = 0.45 (up to 4.1× at Ω = 40), crossing somewhere between. The
full-CLE control recovers in all cells at every γ, so this is not the instrument
failing.

> *Suspect, stated as one (rule 17).* Near γ_c = 0.5 the landscape flattens, so with
> δ deterministic the drift to the threshold becomes weak and slow, and pool
> fluctuations modulating a weak drift produce enormous timing variance — whereas the
> true process has δ's own noise, which lets it cross sooner and more consistently. On
> that reading, removing the signal's noise *inflates* timing jitter near the critical
> point rather than deflating it. **How to kill:** measure the mean advance rate and
> its variance as a function of γ under `s-only`; if the inflation tracks the
> flattening of the drift it is this, and if it does not, something else changes sign
> between γ = 0.30 and 0.45.

**Scope.** One network, one family of observables, four γ. The rule is stated as a
finding rather than a summary because it was tested where it was predicted to fail and
did not — but it has never been tried on a system where the "coordinate the observable
is a functional of" is ambiguous, which is the obvious next attack.


## Open questions

1. ~~**Universality class of the freeze-out transition** (§5). Is a = 0.38 really
   1/3 or 2/5, and can the quasipotential of §2 predict it?~~ **VOID** → §5.1. The
   quasipotential does predict it, and what it predicts is that there is no
   transition: `Hc = 0`, `1/H* = (3/2)lnΩ`, and `a` is a parameter of the wrong
   functional form. Still open (small): the intercept `B`, and the 9% excess in the
   loser-clearing coefficient (measured 1.0895 ± 0.0176 against a predicted 1).
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


## 26. An MLRift SSA that agrees with the exact CME — `mlrift/am_ssa.mlr`

Every tail claim in this project is capped by sampling, not by physics. §24.1c's
conclusions rest on the pattern across nine cells precisely because the per-cell
sampling error (~5% relative at p ≈ 0.0094) is comparable to the effects being
compared. And the project's founding claim concerns a switch that errs at 10⁻¹⁵,
while every number measured here sits between 10⁻¹ and 10⁻².

The two instruments are blocked in opposite directions: the **exact CME** reaches
arbitrarily small probabilities but its state space grows as ~Ω²/2, and **sampling**
handles any Ω but has a floor at 1/N. Large Ω *and* small probability is reachable by
neither.

**This is a probe of whether MLRift can move the sampling floor** — an exact Gillespie
SSA for reversible AM, written in MLRift (Pantelis Christou's self-hosted compiler,
which AOT-compiles to native machine code with no Python, LLVM or ROCm in the path).

> **The gate, fixed before writing a line: it must reproduce the Python/scipy exact
> CME.** Until it does, nothing it computes where no exact reference exists is
> admissible. The propensity conventions were copied from `crnl/reactions.py` rather
> than rederived — heterobimolecular `c = k/Ω`, homodimer `c = 2k/Ω` with `n(n−1)/2` —
> because getting those wrong is the likeliest way for a probe to fail for a reason
> that has nothing to do with the language.

At γ = 0.30, Ω = 40, ε/δ* = 0.20, θ = 0.80 (start X=18, Y=13, B=9; thr = 23),
200,000 trajectories:

| | MLRift SSA | exact CME | deviation | sampling error |
|---|---|---|---|---|
| P(error) | 0.232450 | **0.233847** | −0.60% | ±0.41% (1.5σ) |
| MFPT | 13.42224 | **13.41524** | +0.052% | ±0.14% (0.37σ) |

Both agree within sampling error, on the first run, and the file compiled with zero
errors on the first attempt. The RNG is **counter-based splitmix64 rather than a
stateful stream**, so the same draws are produced whether trajectories run in sequence
on a CPU or in parallel lanes on a GPU — which is what will make a later GPU port
checkable against this file rather than merely faster than it.

**What it buys, measured on the identical cell rather than estimated:**

| | trajectories/sec | time for 10⁸ (p ~ 10⁻⁶) |
|---|---|---|
| Python exact SSA (`gillespie_instrumented`) | 351 | 79 hours |
| MLRift, one core | **6,211** | **4.5 hours** |
| MLRift, 24 cores (estimate, not yet built) | ~149,000 | ~11 min |

**18× on a single core**, which puts the deep-tail regime within reach on CPU alone —
the GPU is not required for the first result, and that ordering is deliberate: a CPU
implementation validated against an exact reference is the thing a GPU port must later
be diffed against.

> **Scope.** One cell, one γ, one Ω. The 24-core figure is an *estimate* — no
> multithreading is implemented yet. And this validates the sampler, not any physics:
> its only purpose is to license deep-tail measurements that no exact reference can
> reach, and each of those will still need its own scope stated.

### 26.1 Threading: measured, and bit-exact — `mlrift/am_ssa_mt.mlr`

`std/thread.mlr` initially segfaulted on this machine (its raw `clone`/`futex` path
predates Linux 7); after Pantelis fixed it, the pool works and the estimate above can
be replaced by a measurement. The estimate stands as printed, per rule 7.

> **The acceptance test is bit-exactness, not statistical agreement.** Because the RNG
> is counter-based on the *global* trajectory index, trajectory k draws the same
> numbers whichever worker runs it. So a threaded run must reproduce the
> single-threaded result exactly. At NW = 1, 4 and 12 on the same 200,000
> trajectories: `wrong = 46490` and `mfpt = 13.42224` **identically at every thread
> count**, matching the single-threaded reference. That is a stronger check than the
> CME gate, which allows sampling slack.

Fixed work of 1.2M trajectories at the gate cell (Ryzen 9 7900X, 12 physical cores /
24 threads):

| NW | wall | traj/sec | speedup | efficiency |
|---|---|---|---|---|
| 1 | 193.6 s | 6,198 | 1.00× | — |
| 2 | 97.8 s | 12,274 | 1.98× | 99% |
| 4 | 48.7 s | 24,620 | 3.97× | 99% |
| 8 | 24.4 s | 49,241 | 7.94× | 99% |
| 12 | 16.5 s | 72,816 | 11.75× | **98%** |
| 16 | 13.5 s | 88,889 | 14.34× | 90% |
| 24 | 9.7 s | 123,967 | 20.00× | 83% (SMT) |

**98% efficiency to 12 physical cores** — the per-worker accumulator design has no
contention, and NW = 1 reproduces the standalone single-threaded rate (6,198 vs
6,211) so the pool adds no overhead.

> **Checked for thermal throttling, because the scaling points above ran only
> 10–194 s and NW = 24 ran just 9.7 s — inside the boost window.** Re-run sustained
> at 15M trajectories (~2–3 min each): NW = 16 gives 95,063 traj/s against its 88,889
> burst (107% retention) and NW = 24 gives 127,866 against 123,967 (103%). **Both are
> faster sustained than in burst** — the short runs were slightly pessimistic from
> startup overhead amortised over less work — so there is no throttling and NW = 24
> stays 34% ahead. Worth noting this workload is atypical: three integers and a
> counter per trajectory, entirely in registers with no memory traffic, which is the
> best case for SMT. A bandwidth-bound workload on the same machine would peak near
> the physical core count instead. Against Python's 351 traj/sec this is **353×**,
which puts 10⁸ trajectories at ~13 minutes.

> **What this does and does not change.** §27 found the exact CME reaches 2.9×10⁻⁸
> unaided, so this throughput is *not* what unlocks the two-species deep tail — that
> claim in §26 was wrong and is corrected there. It matters where the exact state
> space explodes: `~Ω^(n−1)/(n−1)!` in species count, i.e. multi-species networks and
> large Ω.


## 27. The restoration collapse holds over 6.5 decades — exactly — T14-a

§12 and §15 measured `P(error) ~ exp(−Ω·ΔW)` and §22 the barrier dying as
`(γ_c−γ)^1.9745`, both across roughly **one decade** of probability. The transistor
analogy that motivates this project lives at 10⁻¹⁵. T14-a asked whether the collapse
survives further down, or whether one decade of agreement is being extrapolated past
its evidence.

> **§26's stated motivation was wrong, and the correction matters more than the
> result.** §26 argued the deep tail needed a faster sampler because the exact CME's
> state space "grows as ~Ω²/2". It does — but Ω²/2 is *small*. At Ω = 620 that is
> 193,131 states and scipy solves it in **13.6 seconds**. The deep tail was reachable
> exactly the whole time; I assumed a limit instead of measuring one, and built a
> sampler partly on that assumption. What MLRift is actually for is stated at the end
> of this section.

**Measured exactly, γ = 0.30, ε/δ* = 0.35, θ = 0.80:**

| Ω | states | P(error), exact | solve |
|---|---|---|---|
| 40 | 861 | 9.833×10⁻² | 0.0 s |
| 120 | 7,381 | 7.430×10⁻³ | 0.4 s |
| 200 | 20,301 | 8.800×10⁻⁴ | 1.1 s |
| 340 | 58,311 | 2.875×10⁻⁵ | 3.6 s |
| 500 | 125,751 | 6.004×10⁻⁷ | 8.1 s |
| **620** | **193,131** | **2.927×10⁻⁸** | **13.6 s** |

**6.53 decades**, no sampling anywhere, ~35 s of solving in total.

**The raw local slopes bounce, and §24.3 already taught me not to trust that.** They
run −0.0358, −0.0288, −0.0218, −0.0315, −0.0229 — non-monotone by ±40%. The cause is
the same integer-lattice artifact §24.3 found: the start-state margin `d0` is an
integer, so **realised ε swings 10.2%** across the sweep (0.909 to 1.010 of target,
worst at Ω = 40). The threshold's own rounding contributes nothing — adding it to the
fit moves R² from 0.999533 to 0.999534.

| fit | R² | max residual |
|---|---|---|
| `ln P ~ Ω` | 0.998639 | 0.4063 |
| **`ln P ~ Ω + Ω·(ε deviation)`** | **0.999533** | **0.1622** |

**With the lattice artifact removed:**

```
ln P = −0.024904·Ω − 0.2145·Ω·(ε−ε̄) − 1.7280        R² = 0.999533
```

ε-corrected local slopes have **mean −0.025283 and 7.1% scatter**, against raw slopes
that swung ±40%. First half −0.026064 versus second half −0.024212 — a **+7.3%
difference against 7.1% scatter, i.e. within noise.**

> **T14-a answered: the exponential collapse holds across 6.5 decades with a local
> slope stable to ~7%, and no significant drift.** This is the widest-range
> quantitative claim in the project and it required no sampling. The extrapolation
> §12 and §15 made from one decade is supported, at this γ and ε.

**MLRift validated 265× deeper than its gate.** §26 checked the SSA at Ω = 40 where
P = 0.234. Re-checked at Ω = 200 where the exact answer is 8.800×10⁻⁴: 2.4M
trajectories over 12 processes gave **2097 errors = 8.7375×10⁻⁴, −0.71% against a
2.18% sampling error (0.33σ)**, in 118 s wall. Process sharding is bit-exact — four
shards summed to 46490 errors against the single process's 46490 — because the RNG is
counter-based on the global trajectory index. (`std/thread.mlr` segfaults on Linux 7;
in-process threading turned out to be unnecessary.)

**So what is the sampler actually for?** Not deep tails in two-species AM — the CME
owns those. The exact state space is `~Ω^(n−1)/(n−1)!` in the number of species: fine
at n = 3, but for `n_winner_reversible` at n = 4 it is ~Ω⁴/24, which is 6.7×10⁷ states
at Ω = 200 and out of reach. **Multi-species networks and large Ω are where sampling
is the only instrument**, and that is where §26's 18×-per-core belongs — not where I
first pointed it.


## 28. Predicting the collapse slope from closed forms — T14-c — `collapse_slope_absolute.py`

§27's `ln P = −0.024904·Ω + …` is a **fitted** slope, and rule 16 says a law only ever
fitted is never tested. §15 gives `κ(γ) = (3/2)(1−2γ)/(1+γ)` and `δ*(γ)` in closed
form and `breaking_diffusion` gives `D₀`, so the slope is computable with **no free
parameter**.

For a 1-D diffusion `dδ = μ dt + √(D/Ω) dW` between absorbing barriers, the scale
function is `exp(−2Ω V(x))` with `V(x) = ∫₀ˣ μ/D`. By Laplace the numerator is
dominated by the start and the denominator by the saddle, giving

```
d(ln P)/dΩ = −2·V(x₀)          exact         V from the slaved 1-D reduction
d(ln P)/dΩ = −κ·x₀²            quadratic     μ ≈ λδ, D ≈ D(0) near the saddle
```

> **The normalisation is the trap, and it is worth stating because it would have read
> as physics.** `breaking_mode` is a *unit* vector, so `D₀` is computed with
> `(1,−1,0)/√2` while the coordinate `δ = x−y` corresponds to the *unnormalised*
> `(1,−1,0)`. Hence `D_δ = 2D₀`, and the near-saddle exponent is `κx₀²/2`, giving
> `−κx₀²` and **not** `−2κx₀²`. The wrong factor turns a 12% agreement into a 2.3×
> failure. A numerical guard is built in — the exact integral's `x→0` limit divided by
> the quadratic — and it returns **1.0000 at every γ**.

| γ | measured | quadratic | ratio | **exact integral** | **ratio** | decades | R² |
|---|---|---|---|---|---|---|---|
| 0.15 | −0.062353 | −0.083230 | 1.335 | −0.077306 | **1.240** | 12.39 | 0.9907 |
| **0.30** | **−0.025108** | −0.028293 | 1.127 | **−0.026132** | **1.041** | 5.21 | 0.99959 |
| ~~0.45~~ | ~~−0.002196~~ | ~~−0.002033~~ | ~~0.926~~ | ~~−0.001895~~ | ~~0.863~~ | **0.40** | 0.9960 |

> **γ = 0.45 is excluded on a stated criterion, not by eye.** Its collapse spans
> **0.40 decades** over Ω = 40–500; reaching three decades needs
> `ΔΩ ≈ 3·ln10/0.0022 ≈ 3140`, i.e. Ω ~ 3000 and ~4.5M states. The slope is
> under-determined and no conclusion rests on it.

**P1 confirmed:** the exact integral beats the quadratic at every γ, as it must —
the quadratic linearises a drift evaluated 35% of the way to the attractor.
**P2 confirmed on the valid cells:** the quadratic over-estimates the slope magnitude,
consistent with §22.4's finding that `κδ²` is stiffer than the exact barrier away from
the saddle.

> **The genuine success: at γ = 0.30, over 5.21 decades, the parameter-free prediction
> is 4.1% from the measurement.** That is the project's central law tested in absolute
> terms rather than fitted — the thing rule 16 exists to demand, on the claim §12 and
> §15 rest on.

**P4 fires: the ratio drifts with γ, and it is not finite-Ω.** 1.041 at γ = 0.30
against 1.240 at γ = 0.15. The prediction is a large-Ω Laplace asymptotic, so the
obvious defence is that small-Ω cells pollute the fit — but refitting on the upper
half moves the ratios **further from 1**, not closer:

| γ | Ω 40–500 | Ω 200–500 | Ω 340–620 |
|---|---|---|---|
| 0.15 | 1.240 | **1.362** | — |
| 0.30 | 1.041 | 1.053 | **1.063** |

So the closed form is good to ~4–6% at γ = 0.30 and off by 24–36% at γ = 0.15, and
the degradation is real.

> **Suspect, named as one (rule 17).** The prediction uses a **1-D reduction with the
> pool slaved to its nullcline**, and §24.1a measured that reduction's error shrinking
> monotonically with the timescale separation `3(1+2γ)/(1−2γ)` — which is **5.6 at
> γ = 0.15 against 12.0 at γ = 0.30**. Worse reduction where the prediction is worse,
> in the right direction. **How to kill:** measure at γ ∈ {0.20, 0.25, 0.35} and test
> whether the discrepancy collapses against the separation. Two points fit anything;
> five would make it a curve or kill it. If it does not collapse, the γ-dependence is
> in `κδ*²` itself and §15's closed form does not survive the absolute test.

### 28.1 The collapse test — and §28's drift was partly the instrument — T14-c-i

§28 predicted the collapse slope from closed forms and found 4.1% agreement at
γ = 0.30 but 24–36% at γ = 0.15, with the suspect being the 1-D slaved reduction
whose error §24.1a showed shrinking with the timescale separation. Five γ, with
criteria fixed before any ratio was computed: **≥ 2 decades**, **P ≥ 10⁻¹²**, and
eps-controlled fits.

> **The precision floor is the new criterion and it matters.** `p_cme` returns
> `1 − split`, so a tail at 10⁻¹⁷ is pure cancellation noise in double precision.
> §28's γ = 0.15 sweep ran to Ω = 500, which at that slope reaches ~10⁻¹⁷ — **below
> where the reference means anything.** γ = 0.15 stays monotone to 2.0×10⁻¹³ when
> checked, so the floor is set at 10⁻¹² — above where it was verified, not at it.

| γ | sep | decades | P range | measured | predicted | **ratio** | R² |
|---|---|---|---|---|---|---|---|
| 0.15 | 5.57 | 9.91 | 1.5e-2 → 1.9e-12 | −0.066953 | −0.077306 | **1.155** | 0.99983 |
| 0.20 | 7.00 | 8.94 | 2.4e-2 → 2.7e-11 | −0.050710 | −0.056682 | **1.118** | 0.99982 |
| 0.25 | 9.00 | 7.77 | 5.3e-2 → 8.8e-10 | −0.041823 | −0.039834 | **0.952** | 0.99967 |
| 0.30 | 12.00 | 7.10 | 9.8e-2 → 7.9e-9 | −0.024630 | −0.026132 | **1.061** | 0.99973 |
| 0.35 | 17.00 | 5.47 | 1.6e-1 → 5.4e-7 | −0.014605 | −0.015255 | **1.044** | 0.99938 |

**§15's closed form predicts the collapse slope to within ±15% at every γ, across
5.5–9.9 decades each, with no free parameter.** Mean ratio 1.066, spread 0.952–1.155.

> **§28's γ-drift was partly the instrument.** It reported γ = 0.15 at 1.240; with the
> precision floor applied that becomes **1.155**, because §28's deepest cells were
> below double-precision cancellation. §28's numbers stand as printed (rule 7) but its
> "24–36% off" reading overstates the failure.

> ⚠ **My own analysis script produced a false result and I am recording it rather
> than the result.** It printed `excess = 1.069·sep^(−1.134), R² = 0.9966` and the
> line "a clean power law supports P1". That fit uses **4 of 5 points**: γ = 0.25 has
> ratio 0.952, so its excess is negative, so `log(excess)` is undefined, so the code
> dropped it and fitted the rest. The dropped point is precisely the one contradicting
> the pattern. **No power law is claimed.**

**With all five points, P1 fails and P3 fails.** The excesses run +0.155, +0.118,
**−0.048**, +0.061, +0.044 — **not monotone in γ**, and the sign flips at γ = 0.25,
which no account here predicts (P3 said the prediction should be too steep
everywhere). Spearman rank correlation with the separation is −0.70 on five points,
which is not evidence of anything.

> **So the honest position is better for §15 and worse for the mechanism than §28
> suggested.** The closed form is good to ±15% across five γ and nine decades — that
> is the parameter-free test of the project's central law, and it passes. But the
> residual is **scatter, not drift**, so the 1-D-reduction story §28 proposed is
> unsupported: a reduction error that tracks timescale separation cannot change sign
> at γ = 0.25.
>
> **A confound I introduced and should have avoided:** each γ used a *different* Ω
> list, so the cells are not matched across γ and each carries its own realised-eps
> pattern. Since that wobble alone moved §27's raw local slopes by 40%, ±10% scatter
> from unmatched Ω grids is entirely plausible. **How to settle it:** re-run with Ω
> chosen by one rule for every γ — same decade span and same P range — so the only
> difference between cells is γ.

### 28.2 Matched grids revive the mechanism §28.1 killed — T14-c-ii

§28.1 concluded the residual was **scatter, not drift**, and that §28's
1-D-reduction story was unsupported — on the strength of γ = 0.25 coming in at
0.952, below 1, which no reduction-error account can produce. But §28.1 also flagged
its own confound: each γ used a **different Ω list**, so cells were unmatched and each
carried its own realised-ε rounding pattern.

**One rule for every γ:** span the same probability window, 10⁻² → 10⁻⁶, with ln P
equally spaced (hence Ω equally spaced), 12 cells each. Endpoints taken from §28.1's
own measured slopes, so the grid is set by data rather than by hand. Same dynamic
range, same cell count, same P window — the only difference between two γ is γ.

| γ | sep | Ω range | decades | measured | ±se | predicted | **ratio** | R² |
|---|---|---|---|---|---|---|---|---|
| 0.15 | 5.57 | 46–184 | 4.28 | −0.068360 | 0.000281 | −0.077306 | **1.131** | 0.99988 |
| 0.20 | 7.00 | 57–239 | 4.16 | −0.049833 | 0.000278 | −0.056682 | **1.137** | 0.99985 |
| 0.25 | 9.00 | 80–300 | 3.75 | −0.036879 | 0.000150 | −0.039834 | **1.080** | 0.99990 |
| 0.30 | 12.00 | 133–507 | 4.24 | −0.024775 | 0.000119 | −0.026132 | **1.055** | 0.99990 |
| 0.35 | 17.00 | 229–860 | 4.21 | −0.015051 | 0.000054 | −0.015255 | **1.013** | 0.99989 |

> **§28.1's central claim is withdrawn: the sign flip was a grid artifact.** With
> matched grids **every ratio is above 1**, γ = 0.25 included, and the scatter falls
> (sd 0.077 → 0.052). What is left is not scatter but a **monotone drift** — 1.131,
> 1.137, 1.080, 1.055, 1.013, declining toward 1 as γ (and the timescale separation)
> rises. The only break in monotonicity is a 0.006 tie between the first two points.
>
> **The drift is 11.7% against ratio uncertainties of ~0.48% — 24× the measurement
> error.** It is real, and §28.1 could not see it because an unmatched grid had pushed
> one point across 1.

**So §28's 1-D-reduction story is back, and now on all five points:**

```
excess = 5.987 · sep^(−2.031)        R² = 0.887, nothing dropped
```

The reduction slaves the pool to its nullcline, and §24.1a independently measured that
reduction's error shrinking with the timescale separation `3(1+2γ)/(1−2γ)`. Here the
excess falls as roughly **sep⁻²** — and unlike §28.1's discarded fit, this one uses
**all five points**, because with matched grids none has negative excess to drop.

> **What this leaves standing.** §15's closed form predicts the collapse slope with no
> free parameter to **1.3% at γ = 0.35** and within 14% everywhere, across ~4 decades
> per γ. The residual is not a defect in `κδ*²` but in the **1-D reduction used to
> evaluate it**, and it vanishes in the direction where that reduction is exact.
> **§15 survives the absolute test.**

> **Scope, and one honest caveat about R².** The power law's R² = 0.887 is carried by
> five points spanning a factor of 3 in `sep`, and the exponent −2.03 should be read as
> "roughly inverse-square", not as a measured constant. The first two γ are within
> 0.006 of each other, which a genuine `sep⁻²` law would not produce (7.00/5.57 = 1.26,
> so ~37% of excess is expected between them). Something beyond a clean power law is
> present at small γ, and this run cannot say what.


## 29. §24.1's categorical zero is a theorem — and that bounds what it means — `mlrift/am_cle_proj.mlr`

§24.1 reported `s-only` giving **exactly 0 in 8/8 cells** and called the failure
categorical. But that rested on 40,000 trials, so the defensible statement was
`P < 2.5×10⁻⁵` against a truth of 2.3×10⁻¹ — a bound, not a zero. §24's language
asserted the stronger thing on evidence for the weaker one.

> **Going deeper in Ω cannot settle it, and that was my first instinct.** Deeper cells
> make the *truth* fall faster than the sampling floor, so the gap between bound and
> truth shrinks and the claim gets **weaker**. The test that settles it is more
> trajectories at a fixed shallow cell.

The projected CLE was ported to MLRift and gated against the Python arms at §24.1's
cell (γ=0.30, Ω=40, ε/δ*=0.20, 40,000 trials):

| arm | MLRift | Python §24.1 |
|---|---|---|
| full | 0.23820 | 0.237275 |
| δ-only | 0.24813 | 0.244400 |
| s-only | **0** | 0.000000 |
| uniform 11% | 0.01498 | 0.014050 |

At **2,000,000 trajectories** `s-only` is still exactly 0 — and that prompted the
right question, which turns out to have an algebraic answer rather than a statistical
one.

**The drift of the signal coordinate is exactly proportional to the signal.** For
`am_reversible`, collecting the six reactions' contributions to `b_δ = ẋ − ẏ`:

- `f1: X+Y→2B` and `r1: 2B→X+Y` change X and Y **equally**, contributing **0**;
- `f2/f3` give `c_het·B·(X−Y) = c_het·B·δ`;
- `r2/r3` give `−c_hom·[X(X−1) − Y(Y−1)] = −c_hom·δ·(s−1)`.

```
b_δ = δ · [ c_het·B − c_hom·(s−1) ]        no additive term
```

Verified numerically: `b_δ/δ` is constant to **4.6×10⁻¹⁶** across a 25× range in δ at
three pool states. So under the `s-only` projection, where δ carries no noise of its
own, δ obeys `dδ/dt = δ·g(t)` with g depending only on the (stochastically evolving)
pool — giving `δ(t) = δ₀·exp(∫g)`, which **cannot change sign in finite time**.

> **`P(error) = 0` for `s-only` is exact and structural, not a sampling result.**
> §24.1's "categorical" was right, and is now right for a stated reason. No number of
> trajectories was ever going to find a crossing.

**But the reason is not the one §24 gave, and that bounds the claim.** §24 read the
zero as *noise in the wrong subspace is worthless for a barrier-crossing observable*.
The actual mechanism is that AM's drift is exactly proportional to δ, making
**sign(δ) a conserved quantity** once δ's own noise is removed. That is a property of
this network's bilinear structure, not a general fact about coordinate roles.

> **What survives and what narrows.** The *quantitative* half of §24.1 is untouched:
> `δ-only` keeps 11% of the variance and recovers the answer to 2–18%, and the
> `uniform 11%` arm is wrong by 17–770× while retaining δ-noise at reduced amplitude —
> both real measurements about placement and amplitude. What narrows is the
> *categorical* half: **a network whose `b_δ` carries an additive term would give a
> small but nonzero `s-only`, not a zero**, and §24's "a model can keep seven-eighths
> of the noise and be as categorically wrong as one that keeps none" is then a
> statement about AM rather than about coarse-graining.

**Open, and it is the obvious next check:** does the same identity hold for
`n_winner_reversible`? §24.2's `rivals-only` arm also returned categorical zeros at
n = 3, which the same structure would explain — and if it does, §24.2's subspace
result needs the same qualification as §24.1's.

> **→ Answered in §30, and the guess about which arm was wrong.** The identity does
> generalise — to every n, every γ, and every *pair* of committed species — and it
> makes `bookkeeping-only`'s zero a theorem. But it does **not** cover `rivals-only`,
> which starves no pairwise difference direction; that arm's zero is barrier height
> (§30.2), and it is not categorical at all.

### 28.3 The ε axis confirms the attribution; the power law does not survive — T14-c-iii

§28.2 attributed the residual to the 1-D slaved reduction and described it as
`excess = 5.987·sep⁻²·⁰³`. Two things were untested: whether the excess keeps growing
at small γ as that law demands, and whether it depends on **ε**, which the closed form
had never been tested against. Self-calibrating grids — every cell bisects Ω to span
P = 10⁻² → 10⁻⁶ with 12 points, one rule everywhere — over γ ∈ [0.05, 0.35] and
ε ∈ {0.35, 0.50}. **Four cells exceeded the Ω ≤ 900 cap and are reported, not
dropped:** (0.05, 0.35), (0.10, 0.50), (0.15, 0.50), (0.30, 0.50).

**P2 is the discriminating test, and it passes decisively.** If the residual is a
property of the slaved *manifold*, it should not care where on the manifold the
trajectory starts:

| γ | ε = 0.35 | ε = 0.50 | spread |
|---|---|---|---|
| 0.20 | 1.110 | 1.079 | **2.8%** |
| 0.25 | 1.063 | 1.036 | **2.6%** |
| 0.35 | 1.004 | 0.975 | **3.0%** |

**~3% across a 43% change in ε, against a 16-point drift across γ.** The residual is a
manifold property, not a start-point one — which is what the 1-D-reduction attribution
predicts and what a defect in `κδ*²` would not.

**P1 fails, and the power law goes with it.** The excess does keep growing toward
small γ — 0.0044, 0.0339, 0.0634, 0.1098, 0.1255, **0.1598** at γ = 0.35 → 0.10, so
the "saturation" alternative is also dead — but not as `sep⁻²`:

| γ | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.35 |
|---|---|---|---|---|---|---|
| excess | 0.1598 | 0.1255 | 0.1098 | 0.0634 | 0.0339 | **0.0044** |
| §28.2 law predicts | 0.2822 | 0.1829 | 0.1150 | 0.0690 | 0.0385 | 0.0190 |
| ratio | 0.57 | 0.69 | 0.95 | 0.92 | 0.88 | **0.23** |

> **§28.2's exponent is withdrawn as a description.** Refitted over this wider range
> it is **−2.53** (R² = 0.87), not −2.03; an exponent that moves with the fitting
> window is not an exponent. And the decisive objection is structural: **the excess
> reaches zero** (0.0044 at γ = 0.35), which no power law in `sep` can do. A straight
> line in γ fits far better — `excess = 0.2240 − 0.6276·γ`, **R² = 0.9905**, crossing
> zero at γ ≈ 0.357 — but that is six points and a two-parameter fit, and is offered
> as a description rather than a law.

> ⚠ **§36 supersedes the attribution below.** The residual is neither `κδ*²` nor the
> *inaccuracy* of the slaved reduction — it is that the exact runs were started off the
> slaved manifold entirely. The ε-independence measured here is consistent with that
> (the effect is a manifold property), but the mechanism named below is wrong.

> **What stands: the attribution, not the formula.** The residual is ε-independent to
> 3% and γ-dependent by 16 points, which is exactly the signature of the slaved
> reduction and not of `κδ*²`. §15's closed form survives the absolute test, and at
> γ = 0.35, ε = 0.35 the parameter-free prediction is within **0.4%** — the closest
> agreement anywhere in this project.

> ⚠ **§35 QUALIFIES THIS, and the 0.4% is withdrawn as an asymptotic claim.** Every
> slope above is an *effective* slope on the P = 10⁻² → 10⁻⁶ window. §35 solves the
> collapse to 10⁻³³ and finds the local slope drifts, so the asymptotic rate differs:
> at γ = 0.35 the excess is **7.5%, not 0.4%** — an 18.7× correction — and across
> γ = 0.20–0.35 the asymptotic excess falls only 0.155 → 0.075 rather than 0.110 →
> 0.004. **The numbers here stand as measured and as effective slopes; what is withdrawn
> is reading them as the asymptotic rate.** P2's ε-independence is untouched, being a
> comparison at fixed window.

> **One ratio sits below 1** (γ = 0.35, ε = 0.50: 0.975). With the ε-spread at 3%, the
> excess at γ = 0.35 is consistent with zero and **its sign is not resolved**. So the
> honest statement is that the excess declines monotonically from 0.16 at γ = 0.10 to
> ~0 near γ = 0.35, and beyond that this instrument cannot say.

### 30 The identity is general, and it sorts §24.2's arms — `pairwise_identity.py` — T15-a

§29 found that §24.1's `s-only` zero at n = 2 is a theorem: AM's signal drift carries
no additive term, so removing that coordinate's noise conserves sign(δ) exactly.
**T15-a asked whether the same structure explains §24.2's zeros at n = 3.** It does,
in a stronger form, and the general statement was derived by hand from the
count-level propensities of `n_winner_reversible` before anything was run:

> **For every n, every γ, and every pair (i, j) of committed species,**
>
> **d(nᵢ − nⱼ)/dt = (nᵢ − nⱼ) · (k/Ω) · [ n_B − Σ_{l≠i,j} n_l − γ·(nᵢ + nⱼ − 1) ]**

Three cancellations produce it: the pair's own disagreement reaction consumes both
equally; its reverse `2B → Xᵢ + Xⱼ` reaches Xᵢ and Xⱼ through the same (n−1) pairs, so
the entire γ·n_B² term drops; and everything left is bilinear in (nᵢ − nⱼ). At n = 2
the middle sum is empty and this collapses to §29's `b_δ = δ·[c_het·B − c_hom·(s−1)]`
exactly.

**Verified against the network's own stoichiometry and propensities** — not a
hand-rolled copy of them, since the question is whether the shipped conventions
satisfy it — at n = 2, 3, 4, 5, 6 × γ ∈ {0.10, 0.40, 0.60, 0.90}·γ_c(n), 40 random
states per cell, 4,600 pairs in total:

| measure | result |
|---|---|
| residual ÷ total absolute traffic that must cancel | worst **4.4×10⁻¹⁶** |
| strict residual ÷ max(\|lhs\|, \|rhs\|), ties excluded | median **2.1×10⁻¹⁶**, worst 2.9×10⁻¹³ |
| \|bᵢ − bⱼ\| at nᵢ = nⱼ, where the identity forces exactly 0 | worst 1.3×10⁻¹⁵ |
| spread of bᵢⱼ/(nᵢ − nⱼ) as the split varies at fixed nᵢ + nⱼ | worst 5.7×10⁻¹⁴ |

**Both normalisations are reported, not the flattering one.** The strict worst case of
2.9×10⁻¹³ occurs where |nᵢ − nⱼ| = 1 *and* the bracket is near zero, so the left side
is itself near zero and the denominator collapses; the median of 2.1×10⁻¹⁶ is one ulp.
The `traffic` denominator is the honest one for a cancellation claim, and
`max(|lhs|,|rhs|)` is degenerate at nᵢ = nⱼ where both sides vanish — the first pass of
this test reported "P1 FAILS, relative residual 1.0" purely from that 0/0.

**What the identity conserves is not "the signal coordinate" but sign(nᵢ − nⱼ) for
every pair independently.** A projection is covered by the theorem exactly when it
leaves some pairwise difference direction noise-free — and that sorts §24.2's arms,
which had been read three different ways:

| arm | pairwise directions starved | consequence |
|---|---|---|
| `s-only` (n = 2) | the only one | §29's theorem |
| `bookkeeping-only` | **all C(n,2)** | the entire ordering is frozen; P = 0 is a theorem |
| `decision-only` | d2 = (0,1,−1) only | sign(n₂ − n₃) frozen — **this is §24.2's Ω-parity trap** |
| `rivals-only` | **none** | not covered; see §30.1 |

Dynamical confirmation at n = 3, γ = 0.60·γ_c(3) = 0.1213, Ω = 30, ε = 0.25, skew = 2,
40,000 trajectories per arm, all arms and the CME reference on an **identical start
state** [12, 5, 3, 10] (rule 18). Champion-vs-rival flips and rival-vs-rival flips are
counted **separately**, because pooling them makes the counter meaningless — only the
first kind can make the observable wrong:

| arm | var kept | P(error) | champ flips | rival flips | closest champ approach |
|---|---|---|---|---|---|
| exact CME | — | 0.12017 | — | — | — |
| `full` | 1.000 | 0.12350 | 8419 | 36367 | −3.57 |
| `bookkeeping-only` | 0.806 | **0** | **0** | **0** | **+0.4696** |
| `rivals-only` | 0.066 | **0** | **0** | 38913 | **+0.5629** |
| `decision-only` | 0.114 | 0.07310 | 6331 | **0** | −3.23 |

> **P2 confirmed, in the form the theorem actually predicts.** §24.2 measured "no
> trajectory finishes wrong". The stronger statement holds: under `bookkeeping-only`
> **no pairwise sign flips at all**, and the closest any pair came to crossing was
> **9.1×10⁻⁶ of its initial gap**. It approaches zero exponentially and cannot reach
> it — which is a conservation law, not distance from the boundary.

> **P2b: §24.2's parity trap is this theorem.** `decision-only` froze the rival pair
> and only the rival pair — 0 rival flips against 6331 champion flips. §24.2
> diagnosed the trap empirically and filed it with "the floor-division artifacts",
> blaming integer rounding. **The rounding only chose the initial condition; the
> conservation law is what made it fatal.** Rivals that start tied stay tied because
> nothing can move sign(n₂ − n₃), not because 30 is even.

> ⚠ **P3 FAILED, and it is the informative failure.** `rivals-only` keeps d2, which
> has a nonzero component along *every* pairwise difference (δ₁₂ gets −h, δ₁₃ gets
> +h, δ₂₃ gets 2h), so no sign is conserved and I predicted it would turn nonzero once
> the trial count rose far above §24.2's 40,000. At **440,000 trajectories the
> champion's margin never flipped once**, and never came closer than **0.51 of its
> initial gap**. That is a barrier, not rarity. The prediction was wrong and the
> mechanism is elsewhere — §30.1.

**§24.2's stated reason for the `rivals-only` zero is also wrong**, and in a way worth
recording: "with d1's noise removed, X1 versus the rival mean is deterministic and X1
never loses" is a non-sequitur. X1 does not lose to the rival mean, it loses to the
best rival, and δ₂₃ is exactly the direction being driven. The zero is real; that
reason does not establish it.

### 30.1 The champion's margin has a sink, quadratic in rival spread — `rival_erosion.py` — T15-a-i

§30's P3 failed: `rivals-only` never flipped a champion-vs-rival sign in 440,000
trajectories despite no sign being conserved. The mechanism follows from the same
identity rather than from a guess. Under `rivals-only` the noise is (0, +h, −h, 0), so
n₁, n₂+n₃ and n_B are all noise-free and so is the champion's **mean** margin
u = n₁ − (n₂+n₃)/2. Differencing two brackets gives
G₁₃ − G₁₂ = −(1−γ)·δ₂₃/Ω, and therefore

> **du/dt = u·Ḡ − (1 − γ)·δ₂₃² / (4Ω)**, with Ḡ = (G₁₂ + G₁₃)/2

**u carries an additive term, so it is not sign-conserved and the champion is not
protected by a conservation law.** The term is negative for every γ < 1: rival spread
erodes the champion's mean margin at a rate **quadratic in that spread**. This is
§24.2's P7 order-statistic intuition — "the champion loses to the best rival, not to
the rival mean" — as an exact identity. Verified to a worst relative residual of
**1.3×10⁻¹⁵** at γ ∈ {0.10, 0.40, 0.60, 0.90}·γ_c(3), and pinned in the suite.

By §30's identity δ₂₃ grows exactly when G₂₃ = (1/Ω)[n_B − n₁ − γ(n₂+n₃−1)] > 0, which
is negative whenever the champion is well ahead of the blank pool — the case in every
cell §24.2 and §30 happened to run. Sweeping ε at n = 3, γ = 0.1213, skew = 2, 40,000
trajectories per arm, all arms paired on an identical start state per cell:

| Ω | ε | start | G₂₃·Ω | CME | `full` | `bookkeeping` | `rivals` | `rivals`/`full` |
|---|---|---|---|---|---|---|---|---|
| 90 | 0.03 | [21, 20, 18, 31] | +5.51 | 0.59744 | 0.59900 | **0** | 0.93470 | 1.5604 |
| 90 | 0.05 | [23, 19, 17, 31] | +3.75 | 0.47804 | 0.47613 | **0** | 0.42568 | 0.8940 |
| 90 | 0.08 | [25, 18, 16, 31] | +2.00 | 0.35510 | 0.35215 | **0** | 0.10912 | 0.3099 |
| 90 | 0.12 | [27, 17, 15, 31] | +0.24 | 0.24140 | 0.23797 | **0** | 0.01183 | 0.0497 |
| 90 | 0.18 | [29, 16, 14, 31] | −1.52 | 0.14806 | 0.14547 | **0** | 0.00020 | 0.0014 |
| 90 | 0.25 | [33, 14, 12, 31] | −5.03 | 0.03878 | 0.03910 | **0** | 0 | 0 |
| 60 | 0.03 | [16, 13, 11, 20] | +1.21 | 0.48003 | 0.47658 | **0** | 0.42455 | 0.8908 |
| 60 | 0.25 | [22, 10, 8, 20] | −4.06 | 0.09345 | 0.09190 | **0** | 0 | 0 |

> **P2 confirmed in the strongest available form, and this is the durable result of
> the section.** At Ω = 90, ε = 0.03 the champion leads by a **single count** (21
> against 20), the exact CME error is **0.597**, and full noise fails 59.9% of the
> time — and `bookkeeping-only` *still* returns exactly 0 with zero pairwise flips in
> 40,000 trajectories. **A conservation law does not care how low the barrier is.**
> §24.2's zero could have been barrier height; this one cannot be.

> **P4 note, reported rather than clipped:** at that same cell `rivals-only` fails
> *more* than full noise (0.935 against 0.599, ratio 1.56). Removing d1's noise removes
> restoring fluctuations along with destroying ones, so an arm with 6.6% of the
> variance can be **worse** than the full model, not merely degraded. `uniform-11pct`
> in §24.1 was wrong by 17–770× in the same spirit, but always in the safe direction.

> ⚠ **P3 was scored "FAILS as stated" and the honest reading is worse than that: the
> sweep cannot decide the question at all.** The scoring is that 5/5 cells with
> G₂₃ > 0 produced champion flips while 1/3 with G₂₃ < 0 did too (Ω = 90, ε = 0.18:
> 14 flips at G₂₃ = −1.52, just barely negative, where the bracket is state-dependent
> and can turn positive along the way). But **ε sets the champion's margin, which sets
> both the barrier and G₂₃**, so the beautiful four-decade ordering by G₂₃ is equally
> consistent with "`rivals-only` fails when failure is easy". That is rule 9, and I
> read the ordering as support before noticing it.

> **The one accidental matched-barrier pair already argues against the mechanism.**
> Ω = 60/ε = 0.03 and Ω = 90/ε = 0.05 have barriers matched to **0.42%** (CME 0.48003
> against 0.47804) while G₂₃ differs by **3.1×** (+1.21 against +3.75) — and the paired
> ratio is **0.8908 against 0.8940**, a difference of 0.36%. If G₂₃ set the rate, three
> times the bracket should not leave the ratio unchanged to a third of a percent. One
> pair is not a test and it also varies Ω, so §30.2 separates the axes deliberately.

### 30.2 The rival bracket does not gate `rivals-only` — `rival_bracket_scan.py` — T15-a-ii

§30.1's mechanism predicted that `rivals-only` can only fail when G₂₃ > 0, and the ε
sweep ordered it over four decades in exactly that sequence. But ε sets the champion's
margin, which sets **both** the barrier and G₂₃, so that sweep cannot separate them.
G₂₃ depends on n_B, which the champion's margin does not: holding the margin fixed at
7 and the skew at 2, the start state is fixed by n_B alone through
3R = Ω − n_B − m + skew, so raising n_B lowers R and drives G₂₃ across a range **four
times wider than the whole ε sweep produced, with the margin never changing**. Ω = 90,
γ = 0.1213, threshold 63, 40,000 trajectories per arm, paired within each start state.

| n_B | start | G₂₃·Ω | CME | `full` | `bookkeeping` | `rivals` | `rivals`/`full` |
|---|---|---|---|---|---|---|---|
| 13 | [31, 24, 22, 13] | **−23.46** | 0.39641 | 0.38843 | **0** | **0.18222** | 0.4691 |
| 22 | [28, 21, 19, 22] | −10.73 | 0.37678 | 0.37345 | **0** | 0.14427 | 0.3863 |
| 31 | [25, 18, 16, 31] | +2.00 | 0.35510 | 0.35830 | **0** | 0.10500 | 0.2931 |
| 40 | [22, 15, 13, 40] | +14.72 | 0.33118 | 0.32620 | **0** | 0.07185 | 0.2203 |
| 49 | [19, 12, 10, 49] | +27.45 | 0.30484 | 0.29983 | **0** | 0.04530 | 0.1511 |
| 58 | [16, 9, 7, 58] | **+40.18** | 0.27601 | 0.27535 | **0** | 0.02755 | 0.1001 |

> ⚠ **§30.1's mechanism is WITHDRAWN. The rival bracket does not gate this arm, and
> the ε sweep's four-decade ordering by G₂₃ was the barrier all along.** At
> G₂₃·Ω = −23.46 — deeply negative, where §30.1 said δ₂₃ cannot grow and `rivals-only`
> cannot fail — the arm fails **18.2%** of the time. Across the scan
> **corr(ratio, G₂₃) = −0.9957**: not merely absent, but near-perfectly *opposite* to
> the prediction, while corr(ratio, CME) = +0.9870.

> **What breaks the confound is that the two sweeps confound it in opposite
> directions.** In the ε sweep G₂₃ rises *with* the barrier; in the n_B scan it falls
> *against* it. The paired ratio tracks the barrier in both, so the barrier is the
> consistent explanator and G₂₃ is not. Neither sweep alone could have shown this —
> which is the whole content of rule 9.

> **The decisive single number.** Interpolating the ε sweep's ratio-vs-barrier curve
> to n_B = 13's barrier predicts 0.4424; the measured value is **0.4691**, a ratio of
> **1.060** — while G₂₃ differs from the interpolating cells by **25.5 units**. Six
> percent of residual for a quantity §30.1 claimed was gated by that bracket's *sign*.

> **P3 confirmed again, and the control is now the load-bearing result.**
> `bookkeeping-only` returned exactly 0 with zero pairwise flips in all six cells, over
> a range of n_B that moves G₂₃ by 64 units and the barrier by 44%. The theorem does
> not care about either.

> **Two checks that make the table admissible.** *(i) Rule 12:* **zero unfinished
> trajectories in every cell of every arm** — 40,000 of 40,000 absorbed — so no ratio
> here is a conditional mean over a censored population. *(ii) Reproducibility:*
> n_B = 31 and §30.1's ε = 0.08 are the **same start state** [25, 18, 16, 31] reached
> by two different constructions with independent seeds, and they agree at **1.8σ**
> (`full`) and **1.9σ** (`rivals-only`). Rejections held at 0.0014 of steps across the
> scan and do not track the effect, so the harness is not producing it (rule 10).

**What §30–§30.2 leave standing.** The pairwise identity is exact and general, and it
makes `bookkeeping-only`'s categorical zero a theorem at every n — that is the result.
The `du/dt` erosion identity is also exact and is pinned in the suite, but the physical
reading built on it lasted one experiment. **Three mechanisms proposed in this arc,
two withdrawn within the session that proposed them** (rule 17): the sign-conservation
account of `rivals-only`, then the G₂₃ gate. The measurements survive; the stories
about them keep not surviving, and the one that did survive is the one that was derived
algebraically and checked to 4×10⁻¹⁶ rather than inferred from a monotone table.

### 31 An additive drift term breaks the categorical zero — `additive_term.py` — T15-b

§29 predicted that "a network whose `b_δ` carries an additive term would give a small
but nonzero `s-only`, not a zero." That prediction had never been tested, which left
"the identity is *why* the zero happens" as an explanation with no independent
confirmation — rule 17's exact situation. `am_asymmetric` supplies the term. Derived by
hand from its propensities, with s = n_X + n_Y and δ = n_X − n_Y:

> **dδ/dt = δ·(k/Ω)·[ n_B − γ(s−1) ] + (kβ/Ω)·[ n_B·s − γ((s²+δ²)/2 − s) ]**

The first term is §29's identity verbatim; the second exists for every β ≠ 0 and
vanishes identically at β = 0. **P1: verified against the network's own propensities at
six β, worst relative residual 1.8×10⁻¹⁴.** So β tunes precisely the thing §29 says the
zero depends on, in a network already in the repo.

The champion is **Y**, the *disfavoured* symbol, so the additive term erodes its lead
rather than protecting it. γ = 0.25, Ω = 80, threshold 49, ε = 0.06, 100,000
trajectories per arm, `full`/`delta-only`/`s-only` paired on one start state per cell.

**Two start rules, because β lowers the barrier as well as adding the term** — the
§30.2 lesson applied before the run instead of after. `fixed` holds the start state
across β; `matched` holds a fixed distance from the saddle, so the barrier is held.

**The `matched` rule, with the barrier held to 7.1% across the whole sweep:**

| β | saddle | start | CME | `full` | `delta-only` | **`s-only`** | δ flips | closest |
|---|---|---|---|---|---|---|---|---|
| 0.00 | 0.0000 | [30, 34, 16] | 0.32291 | 0.32412 | 0.32568 | **0** (theorem) | 0 | +0.875 |
| 0.02 | −0.0200 | [29, 35, 16] | 0.30801 | 0.31020 | 0.31397 | **< 3.0×10⁻⁵** | 0 | +0.826 |
| 0.05 | −0.0500 | [28, 36, 16] | 0.32633 | 0.32830 | 0.33648 | **< 3.0×10⁻⁵** | 0 | +0.764 |
| 0.10 | −0.0999 | [26, 38, 16] | 0.32918 | 0.32989 | 0.34424 | **0.00001** | 1 | −4.084 |
| 0.15 | −0.1496 | [24, 40, 16] | 0.33107 | 0.33329 | 0.35716 | **0.00069** | 69 | −3.075 |
| 0.30 | −0.2971 | [18, 46, 16] | 0.32672 | 0.32770 | 0.37448 | **0.03210** | 3210 | −1.761 |

> **P3 confirmed under the barrier-controlled rule: §29's prediction holds.** At an
> essentially constant barrier — CME 0.323 → 0.327, spread 7.1%, `full` flat at
> 0.324–0.333 — `s-only` goes from an **exact zero** to **3.2%**. The claim under test
> is categorical, and no change in barrier height can turn an exact zero into a
> nonzero: a theorem does not weaken, it applies or it does not.

> **The control that makes it airtight.** `s-only` retains **0.891 of the variance at
> every β, identical to three decimals** — the noise amplitude is constant across the
> entire sweep and only the *drift structure* changes. So this cannot be an amplitude
> effect, which is the one alternative §24.1's `uniform-11pct` arm showed matters.
> Rejections held at 5–10×10⁻⁵ of steps, flat in β and not tracking the effect (rule
> 10); zero unfinished trajectories anywhere (rule 12).

> ⚠ **β = 0.02 and 0.05 are BOUNDS, not zeros — only β = 0 is a theorem zero.** Zero
> of 100,000 gives a 95% upper bound of 3.0×10⁻⁵, and at β = 0.10 the measured value is
> 1×10⁻⁵, i.e. just at the resolution. So the transition is *smooth and turns on below
> the instrument's floor*, not sharp at some β. Writing "the zero survives to β = 0.05"
> would repeat exactly the bound-versus-zero error §24.1 made and §29 corrected.

> ⚠ **P3 as written FAILED, and the scoring is reported as it came out.** P3 demanded
> nonzero under *both* start rules. Under `fixed` it is not, because **4 of that rule's
> 6 cells are inadmissible on P4** — at β ≥ 0.05 the saddle slides past the fixed start
> and the noiseless ODE no longer reaches Y's attractor. Its only admissible β > 0 cell
> is β = 0.02, which sits below the resolution floor. So `fixed` does not contradict
> P3; **it cannot test it.** The honest statement is: confirmed under `matched`,
> untestable under `fixed`, and the two-rule design was worth having anyway for the
> reason below.

> **P4 was load-bearing, not ceremonial.** In the four inadmissible `fixed` cells
> `s-only` reads **0.597, 1.000, 1.000, 1.000**. Reported without the ODE control that
> would have looked like a spectacular confirmation — "the additive term drives
> `s-only` from 0 to certainty". It is nothing of the kind: the champion is losing
> *deterministically*, which is `am_asymmetric`'s systematic tilt and not a restoration
> failure at all. That module's own docstring warns the bias and the random error "both
> show up as a wrong answer", and that warning was worth what it cost.

> **P6, unpredicted and recorded:** `delta-only`/`full` runs 1.005 → 1.143 as β rises,
> so the signal arm stays inside §24.1's 2–18% band even with the tilt, and degrades
> smoothly rather than breaking.

**What §31 settles.** The mechanism behind §24.1, §29 and §30 now has the independent
test it was missing: the categorical zero tracks the *presence of an additive drift
term*, not the amount of noise removed. **This is the first mechanism in this arc to
survive a test aimed at it** — §30's two did not — and it survived because it was
derived algebraically and checked in absolute terms first, then given a knob that turns
it off and on at fixed barrier. Rule 16, one level up.

### 32 Concatenation works and is a waste — `concatenation.py` — T16

FINDINGS 1's wall says a single AM stage has a fidelity ceiling; the threshold theorem
says concatenation beats the ceiling for a restoring code. Nothing here had ever
concatenated. **The construction is a physical POOL MERGE** — k independent AM tanks of
size Ω run to commitment, their contents are combined into one tank of size kΩ, and
that tank runs AM itself. The merged margin is the sum of the k committed margins,
positive iff a majority answered correctly, and the combining stage carries its own
noise because it *is* an AM tank. **No `sign()`, no free comparison** — the numpy-majority
version would insert a noiseless infinitely-reliable gate, which is the class of error
that has cost this project three withdrawn results (rule 10).

Everything is exact — no Monte Carlo anywhere:
`p₁ = Σⱼ C(k,j)·p₀ʲ·(1−p₀)^(k−j)·p_merge(j)`, with every term an exact CME solve.
γ = 0.25, k = 3, θ = 0.80.

**P4 first, because it is the substantive check and it passes.** The merged stage is
reliable, so the vote is real and the merged tank has no error floor of its own:

| j wrong of 3 | merged start | p_merge | weight |
|---|---|---|---|
| 0 | [93, 3, 24] | **0** | 0.739 |
| 1 | [63, 33, 24] | 0.0030 | 0.235 |
| 2 | [33, 63, 24] | 0.9970 | 0.0249 |
| 3 | [3, 93, 24] | 1.0000 | 0.00088 |

**P1 confirmed and flagged as near-definitional**, exactly as §24.2 flagged its own:
`p₁ ~ p₀^1.85`, R² = 0.9997 (attractor readout) and `p₀^1.78`, R² = 0.9999 (threshold).
Given p_merge ≈ 0 for j ≤ 1 and ≈ 1 for j ≥ 2, `p₁ = 3p₀² − 2p₀³` follows identically,
so the exponent near 2 is a check on the machinery rather than a finding.

**P2 is the result: voting loses to pooling, and by a margin that grows exponentially.**
The control is the same kΩ molecules in ONE tank at the same relative margin:

| Ω | p₀ | p₁ (vote) | p_pool | p₁/p_pool |
|---|---|---|---|---|
| 20 | 1.674×10⁻¹ | 8.39×10⁻² | 9.29×10⁻² | **0.90** |
| 32 | 1.360×10⁻¹ | 5.56×10⁻² | 3.91×10⁻² | 1.42 |
| 44 | 1.370×10⁻¹ | 5.17×10⁻² | 1.74×10⁻² | 2.97 |
| 62 | 6.045×10⁻² | 1.09×10⁻² | 4.84×10⁻³ | 2.25 |
| 74 | 6.130×10⁻² | 1.10×10⁻² | 2.28×10⁻³ | **4.83** |

Voting beats pooling in **1 of 10 cells** — Ω = 20, at a ratio of 0.90, essentially a
tie at the smallest tank where the asymptotic argument is weakest. Everywhere else the
same molecules do better undivided.

> **P3, and this is why the section is more than a trend.** The reason is an exponent
> count that can be predicted in absolute terms: `p₀ ~ exp(−Ωc)`, so voting gives
> `3p₀² ~ exp(−2Ωc)` while pooling gives `exp(−3Ωc)`. **Voting squares the error;
> pooling cubes the exponent.** Hence `ln(p₁/p_pool)` must be linear in Ω with slope
> **c — the collapse rate this project already measured independently.**
>
> | readout | ε-controlled slope | R² | `−2·V_exact` | ratio |
> |---|---|---|---|---|
> | attractor | 0.023059 | 0.975 | −0.021134 | **1.091** |
> | threshold | 0.020844 | 0.970 | −0.021134 | **0.986** |
>
> The prediction comes from §15's closed forms with no free parameter, computed for an
> entirely different purpose (§28's collapse slope), and it lands within **1.4%** on one
> readout convention and **9.1%** on the other. That is rule 16's absolute test, passed.

> **The raw fits were R² = 0.65–0.67 and are reported beside the controlled ones.** p₀
> is non-monotone in Ω (0.167, 0.152, 0.136, 0.155, …) because realised ε wobbles
> 0.206–0.261 on the integer lattice — §27's effect exactly, which it measured as
> bouncing raw local slopes by 40%. ε-controlled fitting lifts R² to 0.97. The
> correction is the established one here, and the uncorrected number stays visible.

> **P5 (rule 13): the two readout conventions give 0.0231 and 0.0208 — a 10% spread
> that brackets the prediction of 0.0211.** So the answer does not depend on the
> convention, and the residual uncertainty in the slope is comparable to the convention
> dependence itself. Neither number is quoted alone.

> **A lattice artifact worth recording:** at Ω = 40 the nominal ε = 0.25 and ε = 0.30
> produce the *identical* start state (both d₀ = 8 after the parity fix), so the ε sweep
> has 5 distinct points, not 6, and the P1 fit contains a duplicated point.

**What §32 settles, and what it does not.** The ceiling survives concatenation. Voting
does suppress error — the merged stage is reliable and the exponent is right — but it
is strictly the worse use of the same molecules, by a factor growing as exp(Ωc). **The
structural difference from quantum error correction is now nameable: in QEC the physical
error rate is fixed and cannot be lowered by using more of the same qubit, so
concatenation is the only lever there. Here error falls exponentially in Ω, so a bigger
tank is a lever, and a better one. Chemistry does not need the code because it has a
cheaper knob.**

> ⚠ **That last sentence was too strong, and §33 corrects it. It is left standing as
> written.** Time-extended re-merging *does* beat the single-tank hold — by up to 5×,
> while burning 29% less — below a crossover in Ω and when cycling is fast. So "does not
> need the code" holds only *above* that crossover. What survives untouched is the
> exponent count and the asymptotic conclusion. The scope note immediately below named
> this as untested, and naming it is what made the correction cheap.

⚠ **Scope, stated rather than left implicit.** This is *one-shot* voting: each tank runs
once to commitment. It is not the time-extended error correction QEC actually performs,
where fresh ancillas repeatedly remove errors that accumulate during storage. §12.1's
depth ceiling is about a bit *held over time*, and the analogue there would refresh
periodically. **T16-a, open: does periodic re-merging beat the single-tank hold?** That
is the comparison where QEC's advantage genuinely lives, and §32 does not touch it.

### 33 Re-merging wins, but only in a bounded window — `remerge_hold.py` — T16-a

§32 compared one-shot voting against pooling and concluded the ceiling survives. But
QEC is *time-extended* — the threshold theorem is about repetition, not a single vote —
so T16-a asks the version §32 could not. **HOLD**: one tank of N = kΩ sits committed
until it spontaneously crosses. **RE-MERGE**: k tanks of Ω, combined into one kΩ tank
every τ and immediately re-split into k portions (mixing and aliquoting, both physical);
a flipped minority is outvoted and the split hands every portion back the corrected
state. Exact MFPTs throughout, γ = 0.30, θ = 0.80.

**P2, the renewal model's own kill test, passes.** The two-state reduction is only legal
if first passage is near-exponential, and `first_passage_moments` settles it exactly —
for an exponential law std = mean:

| N | 12 | 20 | 32 | 48 | 60 | 72 |
|---|---|---|---|---|---|---|
| std/mean | 0.9890 | 0.9891 | 0.9960 | 0.9994 | 0.9998 | 1.0000 |

It is near-exponential everywhere and becomes *more* so as the barrier deepens, which
is Kramers' prediction. **P3 also passes**: L_remerge ∝ 1/τ to 0.6% at small τ (×1.987
against ×2, ×4.867 against ×5), deviating only at large τ where q is no longer small.

> ⚠ **P6 FIRES, and it partly reverses §32.** Re-merge beats hold on lifetime in **10 of
> 13 cells**. At (k = 3, Ω = 8) it lives **5.0× longer while burning 29% slower**; at
> (k = 5, Ω = 14) it lives **3.7× longer**. §32's one-shot result did not generalise to
> the time-extended protocol, which is exactly what §32's own scope note flagged as
> untested.

| k | Ω | N | L_hold | L_remerge | hold/remerge | burn-rate ratio |
|---|---|---|---|---|---|---|
| 3 | 8 | 24 | 2.20×10³ | 1.10×10⁴ | **0.200** | 0.706 |
| 3 | 14 | 42 | 2.03×10⁴ | 4.21×10⁴ | **0.482** | 1.024 |
| 3 | 18 | 54 | 9.02×10⁴ | 1.11×10⁵ | **0.809** | 1.087 |
| 3 | 20 | 60 | 1.91×10⁵ | 1.81×10⁵ | 1.051 | 1.098 |
| 3 | 24 | 72 | 8.54×10⁵ | 4.87×10⁵ | 1.752 | 1.101 |
| 5 | 8 | 40 | 1.58×10⁴ | 3.30×10⁵ | **0.048** | 0.745 |
| 5 | 14 | 70 | 6.65×10⁵ | 2.47×10⁶ | **0.269** | 1.062 |

> **Dissipation is not what drives this, and that is measurable rather than assumed.**
> **P1 FAILED** — `ep_rate/N` is *not* size-independent, running 0.02819 → 0.02534 over
> a 4× range in N (10.8%). So dissipation is accounted explicitly instead of waved
> through. The burn-rate ratio spans **0.71–1.10** across every cell, i.e. within ±30%,
> against lifetime ratios of 5× to 20×. **Dissipation cannot explain a difference an
> order of magnitude larger than itself**, and in the cells where re-merge wins biggest
> it is also the *cheaper* protocol.

> ⚠ **THE RESULT IS CONDITIONAL ON CYCLE SPEED, and this is the dominant sensitivity,
> not a footnote.** Because L_remerge ∝ 1/τ exactly (P3), the winning region depends
> entirely on how fast the merge can be repeated. Largest Ω at which re-merge still
> wins:
>
> | τ / t_relax | 0.5 | 1 | 2 | 5 | 10 | 30 |
> |---|---|---|---|---|---|---|
> | k = 3 | 24 | 18 | 14 | 8 | **never** | **never** |
> | k = 5 | 14 | 14 | 12 | never | **never** | **never** |
>
> Re-merge's advantage exists only if the tanks can be cycled within a few relaxation
> times. At τ ≥ 10·t_relax it never wins at any Ω tested. Quoting the win without the
> cycle time would be quoting half the result.

> **P5, the integer test, as it came out.** Predicted slope ratio k − ⌈(k+1)/2⌉ = 1 at
> k = 3 and 2 at k = 5. Measured **1.0755** (MATCH under the ±0.15 tolerance fixed in
> advance) and **2.3139** (**MISMATCH** — 0.314 over, twice the tolerance). Both run
> high, and *in proportion* — consistent with the Kramers power-law prefactor that the
> pure-exponential ansatz omits and that a straight-line fit to ln T(N) partly absorbs.
> **Comparing the two measured slopes to each other cancels that shared bias:
> 0.285896/0.132891 = 2.1514 against a predicted 2.0000, agreeing to 7.6%.** So the
> integer structure across k is real; the absolute values are not clean enough to claim
> it exactly, and P5 is recorded as failed-as-stated with that diagnosis rather than
> restated to fit.

**What §33 settles, and the correction it forces on §32.** §32 closed with "chemistry
does not need the code because it has a cheaper knob." **That was too strong and is
corrected here:** below a crossover in Ω, and only when cycling is fast, the code *does*
help — substantially. What survives is the asymptotic statement, and it survives with
its exponent count intact: hold/remerge grows as exp((k−m)·c·Ω), so hold wins for large
enough Ω at any fixed τ.

**The sharp contrast with QEC is therefore not "concatenation fails here" but something
better.** In QEC, below threshold, concatenation's advantage **grows without bound with
level**, because the physical error rate is fixed. Here re-merging's advantage occupies
a **bounded window** and then *reverses*, because the physical error rate itself falls
exponentially with Ω — growing the tank eventually outruns the code. **Chemistry has a
knob QEC lacks, and the code wins only until that knob is turned far enough.**

### 34 A closed form for §33's crossover — `crossover_law.py` — T16-b

§33 located the re-merge/hold crossover empirically. The exponent count it validated
determines it. With `m = ⌈(k+1)/2⌉` and `ln T(N) = c·N + a`, setting
`T(kΩ) = T(Ω)^m / (C(k,m)·τ^(m−1))` gives

> **Ω× = [ (m−1)(a − ln τ) − ln C(k,m) ] / [ c·(k−m) ]**

**For every odd k, m−1 = k−m = (k−1)/2 exactly**, so the ratio is 1 and

> **Ω× = (a − ln τ)/c − ln C(k,m)/(c·(k−1)/2)** — the leading term contains no k.

`c` and `a` come from a straight-line fit to `ln T(N)` on the **hold protocol alone**.
**No crossover measurement enters the prediction anywhere.** The crossover is measured
as the *continuous* zero of `ln(L_hold/L_remerge)`, never as the largest integer Ω that
still wins — that quantised form is failure pattern 2 in §4, and §33's own table
reports exactly that quantised version.

**P1 CONFIRMED, and it is the structurally surprising one.** The predicted crossover is
nearly k-independent — spread across k = 3, 5, 7 is **3.00% / 3.61% / 4.39%** at
γ = 0.25 / 0.30 / 0.35 — *even though the win margin at fixed Ω differs by more than 2×
between k = 3 and k = 5* (§33: 0.482 against 0.269 at Ω = 14). The margins differ; the
crossings nearly coincide, because k cancels out of the leading term.

**P2, the absolute test, at the protocol's natural cycle time τ = t_relax:**

| γ | c (hold-only) | a | predicted Ω× | measured Ω× | pred/meas |
|---|---|---|---|---|---|
| 0.25 | 0.189509 | 4.97194 | 15.47 | 14.98 | **1.0326** |
| 0.30 | 0.123526 | 4.74104 | 19.76 | 19.62 | **1.0070** |
| 0.35 | 0.071578 | 4.63967 | 28.18 | 29.15 | **0.9668** |

Across three γ, a 1.9× range in Ω× and a 2.6× range in c, the closed form lands within
**3.3%** with nothing about crossovers ever fitted.

> ⚠ **P2 as stated FAILED over the full grid**, and the failure is patterned. Over all
> 31 reachable cells the linearized form gives pred/meas **0.9191 ± 0.0976**, range
> 0.685–1.033 — well outside the ±10% predicted. The deviation **tracks the
> derivation's own assumption**: correlation with τ/T(Ω×) is **−0.696**, and splitting
> at τ/T = 0.02 gives 0.961 below against 0.880 above. The form was derived for τ ≪ T
> and I did not state that domain in the predictions.

> **Solving the crossover condition EXACTLY — same c and a, still no crossover data —
> removes the systematic:** mean **0.9191 → 0.9864** (8.1% low → 1.4% low), cells within
> 2% go **7 → 11 of 31**. So the linearization was the bias and the exponent count
> underneath it is sound. **Scatter is only reduced 1.14×** and the range widens to
> 0.773–1.196: residual disagreement persists at large τ and small Ω×, where the
> crossover lands near Ω ≈ 7 and both renewal assumptions are marginal — a tank of
> seven molecules, and τ no longer short against T. **I have not isolated that residual
> and am not naming a cause for it** (rule 17).

> **P2's predicted bias direction was also wrong.** I predicted a consistent-sign bias
> from the Kramers prefactor, as §33 saw. 26 of 31 cells fall below 1, but not all, so
> "consistent in sign" fails literally — and the dominant systematic turned out to be
> the small-τ linearization, not the prefactor.

**P3: d(Ω×)/d(ln τ) = −1/c exactly**, with no k and no combinatorics. Measured against
`−1/c` from the same hold-only fit: **0.78 (γ = 0.25), 0.87 (γ = 0.30), 0.996
(γ = 0.35)**. The agreement improves with γ because the low-γ τ-sweeps are dominated by
the large-τ cells where the linearization fails; at γ = 0.35, where the reachable range
is widest, it is exact to **0.4%**.

**P4 holds everywhere used:** std/mean over the fitted sizes stays in 0.950–1.033, so
the renewal reduction is legal at every cell quoted.

> ⚠ **P5's direction was BACKWARDS and it is worth recording why.** I predicted the
> high-γ end would fall out of reach first, reasoning that shallower barriers give
> smaller c and hence larger Ω×. That is true — Ω× runs 15.5 → 28.2 from γ = 0.25 to
> 0.35 — but it is the wrong effect. **The MFPT validity ceiling moves faster than Ω×
> does**: a shallower barrier makes T(N) grow slowly, so exact solves stay trustworthy
> to N = 126 at γ = 0.35 against N = 72 at γ = 0.30 and far less at γ = 0.15. Result:
> γ = 0.35 had **10** reachable cells and γ = 0.15 only **2**. The instrument's reach
> and the physics move the same way and the instrument wins.

**What §34 settles.** §33's crossover is not an empirical boundary but a consequence of
the exponent count, predictable from two numbers measured on the hold protocol alone —
to 3.3% at the natural cycle time, and to 1.4% on average across the whole grid once
the linearization is removed. The k-independence is the load-bearing check, because it
is a structural prediction that could easily have failed and that no fit to the
crossover data would have suggested.

### 34.1 The residual is the Kramers prefactor — T16-c

§34 left scatter of ±8.5% after the linearization was removed, concentrated where
Ω× ≈ 7 and τ/T > 0.05, with three candidate causes and none preferred: (i) the
committed state is meaningless at ~7 molecules, (ii) τ no longer long against t_relax
so portions do not re-amplify, (iii) curvature in `ln T(N)` — the Kramers prefactor —
which a straight-line fit absorbs into `a`. THEORIES named (iii) as the one to kill
first, because it needs no new machinery and **the cheap instrumental explanation has
to be eliminated before a physical mechanism is proposed** (the §28.1/§28.2 lesson).

Refitting the hold data as `ln T = c·N + b·ln N + a` — still hold-only, still no
crossover measurement anywhere — and re-running the identical absolute test:

| ansatz | mean pred/meas | sd | range | within 2% |
|---|---|---|---|---|
| `ln T = c·N + a` | 0.9864 | 0.0854 | 0.773–1.196 | 11/31 |
| `ln T = c·N + b·ln N + a` | **0.9988** | **0.0221** | 0.951–1.052 | **22/31** |

**Scatter falls 3.86× and the mean lands within 0.12% of unity.** (iii) is the cause;
(i) and (ii) are not needed and are withdrawn as candidates rather than left standing.

| γ | 0.15 | 0.20 | 0.25 | 0.30 | 0.35 |
|---|---|---|---|---|---|
| b | −0.653 | −0.459 | −0.285 | −0.069 | +0.074 |
| R² | 0.999995 | 0.999989 | 0.999983 | 0.999991 | 0.999968 |

**The prefactor exponent is γ-dependent and crosses zero near γ ≈ 0.32.** That is
recorded as measured and **no mechanism is proposed for it** (rule 17) — it is the
obvious next thing to explain and explaining it is not this section's job.

> **What this is and is not.** It is 31 cells across 5 γ, 3 k and a range of τ,
> predicted to **0.12% mean with 2.2% scatter, with nothing about crossovers ever
> fitted** — the broadest absolute test in the project. It is **not** parameter-free in
> §28.3's sense: that section's 0.4% came from closed forms with no fitted quantity at
> all, whereas this uses three parameters fitted to the *hold* protocol and then
> extrapolated to a different protocol. Weaker in kind, broader in reach, and the two
> should not be quoted as if they were the same claim.

### 35 The probability floor was an artifact, and the collapse slope is a window property — `deep_tail.py` — T14

THEORIES T14 has said since §21 that the binding constraint is instrumental: *"the
founding claim concerns a switch that errs at 1e-15 while every measured number sits
between 1e-1 and 1e-2 … Large Ω AND small probability is reachable by neither."*
**The probability half of that was an implementation artifact.**

`p_cme` computes the error as `1 − split` — a difference of two numbers near 1 — so it
dies to catastrophic cancellation near 1e-12, and §28 lost its γ = 0.15 cells to exactly
that. But `splitting_probability` takes the favoured-set predicate, so **naming the
wrong outcome as favoured solves for the small number directly, with no subtraction
anywhere.**

| Ω | states | P(error), direct | sec |
|---|---|---|---|
| 400 | 80,601 | 6.663482×10⁻⁸ | 3.0 |
| 1000 | 501,501 | 1.960333×10⁻¹⁷ | 22.3 |
| 1500 | 1,127,251 | 5.600480×10⁻³¹ | 45.5 |
| 2000 | 2,003,001 | **6.354802×10⁻³³** | 115.1 |

**Twenty-five orders of magnitude below anything this project had measured, and through
the founding claim's own regime.** Validated three ways:

1. **Against the established route** — identical to `1 − split` to 7–8 digits across the
   entire overlap (1e-2 down to 1e-11).
2. **Componentwise, not by a norm.** A norm residual is dominated by the large
   components and would not notice a garbage small one. Each row of the transient
   generator has ≤ 7 nonzeros, so the true residual is computable by exact summation per
   row; one refinement step gives the componentwise relative correction. At Ω = 2000,
   h = 6.35×10⁻³³ with **|δ/h| = 1.0×10⁻¹³ at the start state and 1.2×10⁻¹³ as the
   maximum over every positive component**.
3. **There is a reason.** The transient generator is an M-matrix, so its LU solve
   carries no subtractive cancellation and relative accuracy survives to arbitrarily
   small values. The floor was never in the physics or the linear algebra — only in the
   subtraction.

#### 35.1 Every collapse slope published here is a finite-Ω effective slope

**P2 confirmed.** The local slope is not constant. At γ = 0.20, ε-controlled, it drifts
monotonically −0.051603 → −0.049895 over 29.21 decades. So `P ~ A(Ω)·exp(−c·Ω)` with an
algebraic prefactor, and a straight-line fit returns a `c` contaminated by `A`.

Three ansätze, all fitted, all reported (rule 15):

| γ | decades | pure `c` | rms | prefactor `c` | exponent | rms | inverse `c` | rms | spread |
|---|---|---|---|---|---|---|---|---|---|
| 0.20 | **29.21** | −0.049730 | 0.0805 | −0.049064 | −0.4484 | **0.0016** | −0.049468 | 0.0133 | 1.35% |
| 0.25 | 21.20 | −0.035964 | 0.0730 | −0.035156 | −0.4394 | **0.0021** | −0.035588 | 0.0172 | 2.27% |
| 0.30 | 14.30 | −0.024311 | 0.0865 | −0.023643 | −0.4089 | **0.0091** | −0.023952 | 0.0186 | 2.79% |
| 0.35 | 8.79 | −0.014920 | 0.0769 | −0.014195 | −0.3964 | **0.0033** | −0.014523 | 0.0171 | **4.99%** |

The prefactor form beats the pure exponential by **10–45× in rms** for one extra
parameter. §5.1's lesson was that a three-parameter fit beating a two-parameter one by
3% is evidence *against* the extra parameter; a factor of 10–45 is the opposite case and
is read that way deliberately.

**The exponent sits near −0.4 and drifts toward −1/2 as the lever arm lengthens** (29
decades → −0.4484; 8.8 decades → −0.3964), which is WKB's predicted Ω^(−1/2). But
**fixing it at −1/2 costs a factor 2–6 in rms**, so it is not exactly the WKB value on
this data, and **no mechanism is asserted for the difference** (rule 17).

> ⚠ **§35.3 WITHDRAWS this reading.** The powers of 1/Ω are 90–99% correlated over any
> bounded Ω range, so b is an ill-conditioned projection, not a measured exponent: it
> swings −0.36 → −0.64 → −0.34 as the model order goes linear → quadratic → cubic. The
> numbers above are what that fit returned; the "drift toward −1/2" was conditioning
> noise read as convergence. **Nothing here constrains b.** What matters
downstream is that `c` shifts by less than **1.3%** between the free and fixed exponent,
so the rate is robust to the prefactor's exact form.

#### 35.2 The closed form's disagreement is larger than §28 measured, and much flatter in γ

> ⚠ **§36 REINTERPRETS EVERYTHING BELOW.** The 7.5–15.5% is not a property of §15's
> closed form. `V = ∫μ/D` integrates along the slaved manifold, while every exact run
> here starts with the pool at the attractor's `γ/(1+γ)` — off the manifold by 20–46%.
> Started on the manifold, the same parameter-free prediction agrees to **0.5–1.8%**
> (mean 0.9942) at the same γ, ε, threshold and Ω grid. **The rates below are correct
> measurements; they measure the cost of an off-manifold initial condition, not an
> approximation error.**

**P3 confirmed, in the uncomfortable direction predicted before the run.**

| γ | §28.3 ratio (P = 10⁻²→10⁻⁶) | asymptotic ratio | §28.3 excess | asymptotic excess | factor |
|---|---|---|---|---|---|
| 0.20 | 1.110 | **1.1553** | 0.110 | 0.1553 | 1.4× |
| 0.25 | 1.063 | **1.1331** | 0.063 | 0.1331 | 2.1× |
| 0.30 | 1.034 | **1.1053** | 0.034 | 0.1053 | 3.1× |
| 0.35 | 1.004 | **1.0747** | 0.004 | 0.0747 | **18.7×** |

> ⚠ **§28.3's zero crossing is WITHDRAWN as a statement about the asymptotic rate.**
> §28.3 fitted `excess = 0.2240 − 0.6276·γ` crossing zero near γ ≈ 0.357 and reported
> the parameter-free prediction as within **0.4%** at γ = 0.35 — "the closest agreement
> in the project". Asymptotically that cell is **7.5%** off, and the excess declines by
> only a factor of 2 across the γ range rather than a factor of 27. **The near-perfect
> agreement at γ = 0.35 was almost entirely finite-Ω contamination**, and it was
> flattering precisely where the window was shallowest.

> **What §28.3 got right, and what it was.** Its numbers stand as measured — they are
> correct *effective slopes on the P = 10⁻² → 10⁻⁶ window*, and §28.3's P2 (the
> ε-independence that attributed the residual to the slaved reduction) is untouched,
> since that was a comparison at fixed window. What is withdrawn is reading those
> effective slopes as the asymptotic rate. **This also answers T14-c-iv without needing
> γ = 0.38–0.44 at all:** the excess does not cross zero, because it was never as small
> as the shallow window made it look.

> **γ = 0.35 is UNRESOLVED by the criterion fixed in advance.** Its ansatz spread is
> 4.99% against the 3% threshold, because 8.79 decades is the shortest lever here. Three
> of four γ resolve; that one is reported as unresolved rather than quoted.

**What §35 settles.** T14's probability floor is gone — the founding claim's 1e-15
regime is now directly and exactly computable, limited only by the Ω²/2 state space and
double-precision underflow near 1e-308. And the first thing the unlock showed is that a
published headline was an artifact of the window it was measured in. Both halves of that
are the point.

### 35.2 The prefactor exponent is UNRESOLVED — and the rate does not care — `prefactor_exponent.py` — T14-d

§35 measured b = −0.4484 / −0.4394 / −0.4089 / −0.3964 at γ = 0.20…0.35 and flagged
that b's γ-dependence and its lever-arm dependence were confounded. **The confound was
tighter than "two things move together":** §35 used the same Ω grid at every γ, so
`decades = |slope|·ΔΩ/ln10` is a deterministic function of γ. They are one variable
with two names, and no fit to that data could separate them. Breaking it needs a
*different* Ω window per γ, sized to match decade counts — §30.2's method used in
advance rather than in hindsight.

**Sweep A — decades held fixed, γ varies.** Only one target achieved genuine matching:

| target | achieved decade span | spread in b |
|---|---|---|
| 9 | **38.6%** | 13.03% |
| 14 | **25.7%** | 5.20% |
| **21** | **7.5%** | **7.99%** |
| 29 | **36.9%** | 7.46% |

> ⚠ **My own matching tolerance was too loose and three of four targets are not
> interpretable.** `sub_window` takes the longest prefix under the target and the Ω grid
> has 14 points, so achieved decades jump coarsely; the 0.35×target tolerance then let
> a 19.55-decade cell count as "29 decades". Only the 21-decade target (7.5% span) is a
> real matched comparison. There b runs **−0.4519 / −0.4352 / −0.4172** — still
> monotone in γ, on three points.

**Sweep B — γ held fixed, decades vary.** b moves **10.76% / 6.88% / 11.95%** across
window length at γ = 0.20 / 0.25 / 0.30, lengthening toward more negative values.

> **So BOTH dependencies survive and P1 is not confirmed.** b tracks γ at (the one)
> matched decade count, and b tracks window at fixed γ. Neither is eliminated, so
> T14-d's question — is the γ-dependence physics or lever arm? — **is not answered.**

> ⚠ **P4 failed, and informatively.** Adding a `d/Ω` term was supposed to make b *more*
> window-independent if the three-term ansatz was simply incomplete. It does the
> opposite: b₄ ranges from **+0.0624 to −0.7362** across windows where b₃ spans only
> −0.37 to −0.45. `1/Ω` and `ln Ω` are strongly collinear over these ranges, so the
> four-term fit is ill-conditioned and the missing-term reading is **not** confirmed.
> By P4's own stated terms, that also means P1's interpretation is not safe.

**P3, reported for completeness and not leaned on:** extrapolating b against 1/decades
gives limits −0.5049 / −0.4654 / −0.4525 (R² = 0.558 / 0.928 / 0.921), mean **−0.4743**,
i.e. 94.9% of WKB's −1/2. Consistent with −1/2; not a demonstration of it, and the
γ = 0.20 fit rests on three points with R² = 0.56.

#### The result that matters: §35's conclusion does not depend on b

The reason to run this was not curiosity about the prefactor — it was whether §35's
headline survives b being unknown. It does, decisively:

| γ | windows | c spread | b spread | §35's `pred/c` | range over all windows |
|---|---|---|---|---|---|
| 0.20 | 6 | **0.19%** | 10.76% | 1.1553 | 1.1540 – 1.1561 |
| 0.25 | 8 | **0.17%** | 6.88% | 1.1331 | 1.1300 – 1.1320 |
| 0.30 | 8 | **0.17%** | 11.95% | 1.1053 | 1.1034 – 1.1053 |
| 0.35 | 2 | **0.12%** | 1.69% | 1.0747 | 1.0744 – 1.0757 |

**The rate `c` is stable to 0.12–0.19% across every window and every ansatz, while the
exponent it shares a fit with moves by up to 12%.** Propagated into §35's headline, the
closed-form disagreement moves in the fourth decimal. **§35's 7.5–15.5% asymptotic
disagreement, and the withdrawal of §28.3's zero crossing, stand independently of the
prefactor being unresolved.** That is the check rule 14 asks for — a withdrawal verified
as carefully as an assertion.

**T14-d stays open, with what it would take.** Matching decades to a few percent needs a
finer Ω grid than 14 points, and separating a γ-effect from a lever-arm effect needs
longer arms at *matched* arm length — which at fixed cost means larger Ω, where the
Ω²/2 state space bites (Ω = 2000 already costs ~500 s per γ for 14 points). The honest
present statement is **b ≈ −0.45 ± 0.05, consistent with −1/2, with neither dependence
eliminated.**

### 35.3 T14-d is ill-posed, and that can be proved — `local_slope_law.py`

§35.2 reported a null and blamed my own design: shared window edges and a collinear
basis. Both criticisms were correct and neither was the reason. Redone properly — the
prefactor read from the **local slope**, where `s(Ω) = d(lnP)/dΩ = −c + b/Ω` makes b the
slope of a two-parameter line with the constant differentiated away, ε-corrected before
differencing (§27) — the answer is that **the question cannot be asked of this data at
all.**

**P1 failed: the local slope is not linear in 1/Ω.** Adding a 1/Ω² term lifts R² from
0.942/0.950/0.947 to 0.993/0.988/0.991 and cuts the rms ~2.5–3×, at every γ. The
curvature is real, and by P4's terms it is an upward 1/Ω² term — a further algebraic
correction, `ln P = −cΩ + b·lnΩ − q/Ω + a`.

**But b does not converge with model order:**

| γ | b (linear) | b (quadratic) | b (cubic) | q (quad) | q (cubic) |
|---|---|---|---|---|---|
| 0.25 | −0.3560 | −0.6434 | −0.3404 | +44.2 | −77.7 |
| 0.30 | −0.3428 | −0.5826 | −0.3351 | +36.9 | −62.7 |
| 0.35 | −0.3261 | −0.5726 | −0.3068 | +37.9 | −69.0 |

b swings by a factor of ~2 and q changes sign. **The diagnosis is proved rather than
asserted:** over Ω = 150…1950 the reciprocal spans a factor of only 13, and on that
range the asymptotic basis functions are near-parallel —

| pair | correlation |
|---|---|
| 1/Ω vs 1/Ω² | **+0.961** |
| 1/Ω² vs 1/Ω³ | **+0.986** |
| 1/Ω vs 1/Ω³ | +0.905 |

with design condition numbers **6.7×10² → 3.6×10⁵ → 2.8×10⁸** at orders 1, 2, 3.

> **The function is determined; the decomposition is not.** Extrapolated to 1/Ω → 0 the
> three orders agree to **0.69% / 0.85% / 1.46%**, and extrapolated past the data to
> Ω = 4000 they agree to 0.49% / 0.60% / 1.03% — while the coefficients they are built
> from disagree by a factor of two. That is the exact signature of an ill-posed
> decomposition sitting inside a well-posed limit.

> ⚠ **So T14-d is not open — it is ILL-POSED for any bounded-Ω instrument**, and the
> decades of P are irrelevant to it. What matters is the range in 1/Ω, and no amount of
> depth in P widens that. This retroactively explains the whole arc: §35's b values,
> §35.2's inability to separate γ from lever arm, §35.2's P4 scatter, and the 40–139%
> spreads across sliding windows are all **one ill-conditioned projection reported four
> different ways.** §35.2 said "the basis was wrong"; the truer statement is that *every*
> basis is wrong for this question over a bounded range.

> ⚠ **§35.1's reading of b is WITHDRAWN.** It reported b = −0.4484 / −0.4394 / −0.4089 /
> −0.3964 and said the exponent "sits near −0.4 and drifts toward −1/2 as the lever arm
> lengthens". Those numbers stand as what that particular fit returned; **they are not a
> measurement of an exponent**, and the drift toward −1/2 was reading conditioning noise
> as convergence. Nothing in this project currently constrains b, including whether it
> is −1/2.

**What this buys, and it is the point.** §35's rate `c` is now verified across **four
bases and three model orders**: the ln P fit, and the local-slope fits at linear,
quadratic and cubic order, agreeing to **0.03–0.62%** against §35's published values.

| γ | c (§35, ln P) | c (slope, lin) | c (slope, quad) | c (slope, cubic) | shift |
|---|---|---|---|---|---|
| 0.25 | 0.035156 | 0.035286 | 0.035045 | 0.035213 | −0.32% |
| 0.30 | 0.023643 | 0.023837 | 0.023635 | 0.023773 | −0.03% |
| 0.35 | 0.014195 | 0.014314 | 0.014107 | 0.014255 | −0.62% |

**§35's 7.5–15.5% asymptotic disagreement with §15's closed form, and the withdrawal of
§28.3's zero crossing, stand — and stand more firmly than when they were published**,
because the one quantity they depend on is the one quantity this whole arc has proved
robust while everything around it was not.

**T14-e, open: derive b analytically instead of measuring it.** The splitting
probability is a ratio of scale-function integrals, `S(x) = ∫exp(2Ω∫μ/D)`, and Laplace
asymptotics on that ratio yields the algebraic prefactor in closed form — the Gaussian
widths at the dominant endpoints partially cancel, which is why the exponent need not be
−1/2. **How to kill:** derive it, then test it in absolute terms against the *function*
`s(Ω)`, which §35.3 shows is determined to ~1% even though its coefficients are not.
That is the right target: predict the curve, not the coefficients. Measuring b would
need 1/Ω decorrelated across ~100×, i.e. Ω ≈ 15,000 and ~10⁸ states — out of reach for
the exact solver, and not worth reaching for when the analysis is available.

### 36 The 7.5–15.5% discrepancy was an initial condition, not an approximation error

§22.4, §28, §28.1, §28.2, §28.3 and §35 all measured §15's parameter-free closed form
running too steep against the exact collapse — 4.1% at best, 7.5–15.5% asymptotically —
and every attempt to explain it named an approximation: the 1-D slaved reduction
(§28.3), the finite-Ω window (§35), the Gaussian truncation (§35.4, eliminated at
<0.62%). **None of them was the cause. The prediction and the measurement were about
different initial conditions.**

`V = ∫μ/D` integrates along the **slaved manifold** — it describes δ evolving with the
pool on its own nullcline. Every exact run since §12 has been started by `_setup`, which
places the pool at the **attractor's** value `b = γ/(1+γ)`. Those are not the same point:

| γ | pool at attractor | pool on nullcline at x₀ | gap | c (start OFF manifold) | c (start ON manifold) | pred/meas OFF | pred/meas ON |
|---|---|---|---|---|---|---|---|
| 0.20 | 0.16667 | 0.31045 | **46.3%** | −0.048838 | −0.057072 | 1.1678 | **0.9993** |
| 0.25 | 0.20000 | 0.31560 | **36.6%** | −0.035191 | −0.040347 | 1.1358 | **0.9907** |
| 0.30 | 0.23077 | 0.32004 | **27.9%** | −0.023880 | −0.026651 | 1.0962 | **0.9822** |
| 0.35 | 0.25926 | 0.32392 | **20.0%** | −0.014128 | −0.015194 | 1.0805 | **1.0047** |

**Same network, same γ, same ε, same threshold, same Ω grid, same exact solver, same
parameter-free prediction. The only change is where the pool starts.** Off the manifold:
pred/meas = 1.0805–1.1678, mean 1.1201. On it: **0.9822–1.0047, mean 0.9942.** The
discrepancy shrinks by a factor of **14.8**.

> **And the mechanism is quantitative, not just directional.** The pool gap orders
> exactly with the excess across γ — 46.3% → 16.8%, 36.6% → 13.6%, 27.9% → 9.6%,
> 20.0% → 8.1%. **That is the γ-dependence** which §28.3 could only describe with a
> straight line (`excess = 0.2240 − 0.6276γ`) and which §35 re-measured asymptotically
> without explaining. It was never a property of the closed form; it is the
> γ-dependence of how far `γ/(1+γ)` sits from the nullcline.

**The honest framing is not "§15 was wrong" and not "§15 was right".** Both starts are
legitimate physical preparations — `_setup`'s choice is a modelling convention, not a
law — and they pose different questions. The error was **comparing a prediction about
one against a measurement of the other**, for fourteen sections, while attributing the
gap to successively more refined approximations.

> **What this reinterprets.** §22.4's "κδ² is stiffer than the exact barrier", §28's
> 4.1%/24–36%, §28.1's ±15% scatter, §28.2's `sep⁻²·⁰³`, §28.3's zero crossing, and
> §35's 7.5–15.5% asymptotic disagreement are all measurements of **the cost of starting
> off the slow manifold**, reported as properties of the closed form. **Every number in
> those sections stands; what they were measurements *of* changes.**

> **What is untouched.** §35's instrument unlock (the direct rare-event solve, P = 6.35×10⁻³³,
> validated componentwise) — that is how this was found at all. §35.1's drift and §35.3's
> proof that the prefactor decomposition is ill-posed. §35.4's elimination of the
> Gaussian truncation. §28.3's ε-independence, which is consistent with this reading
> since the effect is a manifold property. And the §29–§31 identity arc, which never
> touched this.

> ⚠ **The residual is now 0.5–1.8% and is NOT one-signed** (0.9822, 0.9907, 0.9993,
> 1.0047 — two below 1, two above). At this size it is within reach of the genuine 2-D
> path correction, of the ε-lattice, and of the fit itself. **No mechanism is claimed
> for it** (rule 17), and in particular the two-sided sign means it is not obviously the
> one-sided minimum-action effect a 2-D correction would give.

**How this was found, because the route matters.** Not by looking for it. The
`slaving_axis` experiment was built to test §28.3's attribution on a separation axis
independent of γ, and its **P0 gate failed** — the T = 1 cell disagreed with §35's
published γ = 0.25 number. The gate was there to catch a broken instrument; what it
caught was a fourteen-section misattribution. The sweep it was gating had R² = 0.29 and
a non-monotone ratio, and would have been the weakest result of the session.

**T14-f, open: is the residual 0.5–1.8% the 2-D minimum-action correction?** It is the
right size and the only named candidate left, but its sign flips across γ where a
path-minimisation correction should be one-signed (a minimum over paths cannot exceed
the value along the slaved one, so the prediction should be uniformly ≥ the truth).
**How to kill:** compute the 2-D geometric minimum action with the full WKB Hamiltonian
and compare in absolute terms against the on-manifold rates above. If it is one-signed
and ~1%, it closes. If the measured residual keeps flipping sign, it is numerical and
the closed form is exact to the precision of this test.

### 37 An optimal drive exists and is protocol-robust; the cost per nat is not — T-COST

The founding question asks what restoration costs. This project has a minimum affinity
(§9.1), an error exponent (§1) and a fuel lifetime (§20), and no relation joining them.
Both sides of the trade are exact linear solves on the same generator:

    L(Ω)  = −ln P(error)                                  [nats of reliability]
    Σ(Ω)  = E[entropy produced before absorption]          [k_B]

`first_passage` solves `Q_tt T = −1`; the expected entropy production solves the same
system with the local entropy rate as its source, `Q_tt Σ = −σ_local` with
`σ_local(n) = Σ_j a_j(n)·ln[a_j(n)/a_rev(n+S_j)]`. So the cost of a decision is exact,
not a rate times a time. **R = Σ/L is the k_B spent per nat of reliability bought**, and
it needs §35's direct solve, §36's on-manifold start, and θ scaled with δ*.

**Gate: R is Ω-independent** — both sides are extensive, and R flattens to 0.25% between
Ω = 60 and 100. It is a property of the chemistry, not of the molecule count. The
asymptotic value is `R∞ = s/c` from *separate* linear fits of Σ(Ω) and L(Ω), since both
carry Ω-independent offsets (worst R² = 0.9992 over the reference grid).

**R diverges at both ends of the drive, and the components say why:**

| γ | 0.0025 | 0.01 | 0.04 | 0.16 | 0.32 |
|---|---|---|---|---|---|
| R∞ (k_B/nat) | 26.01 | 21.79 | **16.90** | 20.83 | 73.99 |
| s/ln(1/γ) | 0.791 | 0.796 | 0.767 | 0.857 | 1.480 |
| c | 0.1822 | 0.1682 | 0.1461 | 0.0754 | 0.0228 |

`s/ln(1/γ)` is **essentially constant at 0.77–0.80** below γ ≈ 0.08 — the entropy per
molecule tracks the cycle affinity A/3 — while `c` saturates near 0.19 as γ → 0 and
collapses at γ_c. So **R ≈ 0.79·ln(1/γ)/c(γ)**: drive too hard and every cycle dissipates
ln(1/γ) while the barrier has already saturated; drive too softly and the landscape
shallows faster than the saving. **The minimum is forced by two divergences with an
explicit mechanism, not fitted as a feature.**

> **This project has claimed an optimal drive before and withdrawn it.** THEORIES §4:
> *"Dissipation has a minimum near γ ≈ 0.3 — a clean U-shaped curve"*, killed because the
> threshold was held fixed while δ*(γ) shrank. **That withdrawal stands** — the minimum
> here is at γ ≈ 0.07, and §9.2's curve was an artifact. The coincidence of shape is
> exactly why the harder test was run before believing it.

**P2, the deciding test, SPLITS — and the split is the finding:**

| | γ* | R* |
|---|---|---|
| across θ ∈ {0.70, 0.80, 0.90} | 0.07400 / 0.07283 / 0.07290 — **1.6%** | — |
| across ε ∈ {0.25, 0.35, 0.50} | 0.08252 / 0.06327 / 0.07290 — 26.4%, non-monotone | — |
| whole 3×3 grid | 0.0624–0.0826, **factor 1.32** | 3.57–47.77, **203%** |

> **γ\* survives, and survives the specific way that matters.** §9.2 died because its
> optimum tracked θ. Here θ moves γ* by **1.6%** — that exact failure mode is ruled out.
> The ε variation is 26% and **non-monotone** (0.0825 → 0.0633 → 0.0729), consistent with
> parabolic refinement noise on a shallow minimum rather than a trend.

> ⚠ **R\* does NOT survive, and P4 is WITHDRAWN.** R* varies by 203% across the grid,
> collapsing with ε — obviously so in hindsight, since ε sets how hard the decision is
> and a decision from a wide margin is cheap per nat. **There is no universal cost per
> nat of reliability.** The arithmetic consequence "transistor-grade reliability costs
> R*·ln(10¹⁵) ≈ 564 k_BT per decision" is therefore **ε-specific and must not be quoted
> as a constant** — at ε = 0.50 the same arithmetic gives 123 k_BT and at ε = 0.25 it
> gives 1548. Reported here only to make the withdrawal concrete.

> **The optimum is BROAD, not sharp.** At the reference protocol R runs 16.97 / 16.73 /
> 16.40 / 17.83 across γ = 0.035 / 0.05 / 0.07 / 0.10 — within 5% of the minimum over
> **γ ∈ [0.03, 0.08]**, a factor of ~2.7 in drive (A ∈ [7.6, 10.5]). Quoting γ* to three
> digits would overstate what a minimum this flat can locate.

> ⚠ **§38 CORRECTS the optimum's location.** R = Σ/L is not a quantity at all — Σ falls
> with the input margin while L rises, so R varies by a factor of 36 for a trivial
> reason. Priced per **e-fold of gain** instead, the optimum sits at **γ ≈ 0.20**, not
> 0.07. The θ-robustness measured below is real; what it located was the optimum of a
> construction that conflates gain with margin.

**What §37 establishes.** A restoring chemical switch has an **optimal operating drive**,
located at γ ≈ 0.07 (A ≈ 8) with a broad basin, robust against the protocol axis that
destroyed the previous attempt. That is a design principle and it is new. **What it does
not establish is a universal price for reliability** — the cost per nat depends on the
margin the decision starts from, so the founding question's number remains
preparation-dependent rather than fundamental.

**T-COST-a, open: is there a cost that IS margin-independent?** R depends on ε because it
divides by the reliability bought from a particular start. The candidate invariant is the
cost per nat *at fixed margin-to-threshold ratio*, or the total Σ to traverse the whole
landscape (saddle to attractor), which has no free start point. **How to kill:** compute
Σ for the full traverse and check whether Σ/L is ε-free by construction; if it is, that
is the founding question's number and §37's R is a projection of it.

### 38 Restoration is priced per e-fold of GAIN, and that corrects §37 — T-COST-a

§37 measured R = Σ/L, entropy per nat of reliability, and found it varying 203% with
protocol. **That is structural, not a sensitivity.** At fixed γ and Ω, as the input
margin ε rises:

| ε | 0.20 | 0.30 | 0.40 | 0.50 | 0.60 |
|---|---|---|---|---|---|
| L (nats) | 15.0 | 29.7 | 50.5 | 74.7 | **101.7** |
| Σ (k_B) | 1094 | 759 | 514 | 339 | **205** |
| R = Σ/L | 72.8 | 25.5 | 10.2 | 4.5 | **2.0** |

**Σ falls while L rises**, so R collapses by a factor of **36** for a trivial reason: a
start nearer the threshold needs fewer reactions *and* is more reliable. **"Cost per nat
of reliability" is not a quantity.** Reliability is bought with input margin, which is
free. Dissipation buys something else.

**It buys GAIN.** A restoring switch takes ε·δ* and delivers θ·δ* — it amplifies — and
the entropy tracks the *logarithm* of that amplification, as an exponential amplifier
should, since it traverses margin multiplicatively:

> **G = Σ / (Ω · ln(θ/ε))   [k_B per molecule per e-fold of gain]**

preparation-free by construction, since ε enters only through the gain it defines.
Dividing rather than fitting leaves Σ's Ω-independent offset drifting through ln(gain),
so G is estimated as the **slope** of Σ against `Ω·ln(θ/ε)` fitted jointly across ε and Ω.

| γ | 0.02 | 0.04 | 0.07 | 0.12 | **0.20** | 0.30 | 0.40 |
|---|---|---|---|---|---|---|---|
| G | 3.788 | 3.125 | 2.645 | 2.264 | **1.986** | 2.113 | 3.271 |
| R² | 0.9995 | 0.9988 | 0.9979 | 0.9963 | 0.9942 | 0.9940 | 0.9979 |

**The improvement over R is the headline.** Naive spread across the ε×Ω grid falls from
R's **3600%** to **7–27%**.

> **An interior minimum, and it is θ-INVARIANT** — the test that killed §9.2 and that
> §37 passed:
>
> | θ | 0.70 | 0.80 | 0.90 |
> |---|---|---|---|
> | γ* | **0.20** | **0.20** | **0.20** |
> | G* | 2.0358 | 1.9927 | 1.9395 |
>
> γ* does not move at all (0.0%, i.e. within grid resolution) and G* spans **4.8%**,
> against R*'s 203%. **G* ≈ 1.94–2.04 k_B per molecule per e-fold of gain.**

> ⚠ **THIS CORRECTS §37, published the same session.** §37 located a design principle at
> γ* ≈ 0.07 by minimising R — and R is not a quantity. **The optimal drive for the cost
> of restoration is γ ≈ 0.20 (A ≈ 4.83), not 0.07.** §37's θ-robustness result stands as
> measured; what it was robustly locating was the optimum of a construction that
> conflates gain with margin.

> ⚠ **P1's gate MARGINALLY FAILS and G is not an exact invariant.** The joint linear fit
> gives R² = 0.9940–0.9995 against the 0.995 threshold fixed in advance, failing at
> γ = 0.20 and 0.30 — and the fitted intercepts are *negative* (−13 to −28), which is
> unphysical as Ω·ln(gain) → 0 and marks where the linear form gives out. Residual
> margin-dependence is 7–27%. **G is a very good description of how restoration is
> priced; it is not a universal constant, and the difference matters.**

> **G* ≈ 1.99 sits close to 2, and that is left as an observation, not a claim** (rule 15,
> and the P3 written before the run). The θ-trend is monotone — 2.036 → 1.993 → 1.939 —
> so it is drifting through 2 rather than converging on it. §28.2's power law and §35.1's
> −1/2 were both structure read into fitted quantities, and both were withdrawn.

**What §38 establishes.** The founding question has a well-posed form at last: not *what
does a bit cost* but **what does an e-fold of restoration gain cost**, and the answer is
≈ 2 k_B per molecule, minimised at γ ≈ 0.20 with a θ-invariant optimum. What it does not
establish is a universal constant — the residual 7–27% preparation-dependence is real and
unexplained.

**T-COST-b, open: what is the remaining 7–27%?** The negative intercepts say the linear
form fails at small `Ω·ln(gain)`, so the leading candidate is a finite-size offset that a
second term would absorb. **How to kill:** fit `Σ = G·Ω·ln(gain) + A·Ω + B·ln(gain) + C`
and check whether G stabilises and the intercept turns physical. §35.3 is the standing
warning: if the added terms are collinear over the available range, G will swing without
converging and the decomposition will be ill-posed rather than merely incomplete.

### 39 The cost has a closed form good to ~6%, and the optimum is predictable — T-COST-b

§38 left a 7–27% preparation-dependence with unphysical negative intercepts, and proposed
a four-term fit. **Physics forbids it:** as ε → θ the start *is* the threshold, absorption
is immediate and Σ → 0, so the constant and bare-Ω terms cannot exist. That also exposes
`ln(θ/ε)` as an approximation — it is the traversal time only under pure exponential
growth `δ̇ = λδ`, while the real drift saturates near the attractor. The unfitted
prediction is the ep rate integrated along the actual path:

> **Σ_pred = Ω · ∫ σ(δ)/μ(δ) dδ**, σ = Σ_r f_r ln(f_r/f_rev), both from the network's own
> fluxes on the slaved manifold. **No fitted parameter of any kind.**

**P1 FAILS, marginally.** Over 36 cells, pred/exact = 0.9050–1.1674, mean **1.0583**,
sd 0.070 — against the 5% gate fixed in advance. **P2 FAILS**: the ε-spread is 7.4–10.0%,
no better than §38's, so the traversal integral did *not* capture what `ln(θ/ε)` missed.

> ⚠ **P3 fires, and it is the interesting failure.** The residual does **not** shrink
> with Ω — mean |ratio−1| = 0.0784 / 0.0773 / 0.0771 at Ω = 150 / 300 / 450, flat over a
> 3× range. So it is **not** a finite-size effect. And it is not an error in σ either: a
> gate run afterwards against the CME's own `σ_local` shows my concentration-flux σ high
> by 8.0% at Ω = 150 falling to **1.1% at Ω = 1200** — that error converges properly,
> which the total residual does not. **There is a ~6% gap between the deterministic path
> cost and the exact stochastic cost that survives Ω → ∞, and I do not have it.**

> ⚠ **The gate against `σ_local` should have been run BEFORE the comparison, not after.**
> §36 was found by exactly such a gate failing, and the lesson written there —
> *gate every new instrument against the established one at a cell where they must
> agree* — was available and not applied here. It happened to exonerate σ; it might not
> have.

**P4 HOLDS, and it is what survives.** Minimising the closed-form integral over γ, with
**no CME solve at all**, puts the optimum at **γ = 0.240**, against §38's CME-measured
**γ* ≈ 0.20**. Two computationally independent routes — a deterministic path integral and
an exact master-equation solve with an entropy source — agree on the location of the
optimal drive.

| γ | 0.040 | 0.120 | 0.200 | 0.280 | 0.360 | 0.440 |
|---|---|---|---|---|---|---|
| ∫σ/μ (k_B/molecule) | 2.843 | 1.924 | 1.633 | 1.628 | 2.020 | 3.999 |

**What §39 establishes.** The thermodynamic cost of a restoring stage is predictable in
closed form to ~6% with no fitting, and **the optimal drive is predictable analytically**
— confirming §38's γ* ≈ 0.20 from an independent direction. What it does not establish is
an exact cost: the ~6% residual is flat in Ω, not attributable to σ, and unexplained.

**T-COST-c, open: what is the flat ~6%?** It is not finite-size and not σ. Remaining
candidates, none tested and none preferred (rule 17): (i) the mean first-passage time
differs systematically from the deterministic traversal even as Ω → ∞, because absorption
selects early-fluctuating trajectories; (ii) the exact path drifts off the slaved manifold
under noise, so ∫ along the manifold is the wrong contour; (iii) E[σ(state)] ≠ σ(E[state])
by a Jensen gap that does not close. **How to kill (i):** compare the exact MFPT to the
threshold against `∫dδ/μ` directly — that isolates the time from the entropy, and it is
one linear solve.

### 39.1 The cost residual is entirely a TIME residual — T-COST-c closed

§39 left a ~6% gap between the closed-form cost `Σ_pred = Ω∫σ/μ dδ` and the exact CME
cost, flat in Ω, not attributable to σ, with three untested candidates. Candidate (i) —
that the mean first-passage time differs from the deterministic traversal — was named as
the one to kill first, because it isolates the time from the entropy in a single solve.

**It is (i), and the test is decisive.** Comparing `T_det = ∫dδ/μ` against the exact MFPT
to `|δ| ≥ thr`, and both ratios against each other:

| | range | mean |
|---|---|---|
| `T_det / MFPT` | 1.0209 – 1.1764 | 1.1069 |
| `Σ_pred / Σ_exact` | 0.9367 – 1.1683 | 1.0796 |

**correlation +0.9513**, and their difference **shrinks monotonically with Ω in every
cell**:

| γ, ε | Ω=150 | Ω=300 | Ω=450 | Ω=700 |
|---|---|---|---|---|
| 0.07, 0.30 | 0.0346 | 0.0217 | 0.0148 | **0.0119** |
| 0.20, 0.50 | 0.0421 | 0.0207 | 0.0131 | **0.0073** |
| 0.30, 0.30 | 0.0843 | 0.0461 | 0.0344 | **0.0236** |
| 0.30, 0.50 | 0.0897 | 0.0422 | 0.0321 | **0.0207** |

> **The entire thermodynamic residual is kinematic.** The entropy *rate* along the path
> is right; the *clock* is wrong. `∫σ/μ dδ` overestimates the cost by exactly the factor
> by which `∫dδ/μ` overestimates the first-passage time. **Candidates (ii) — the path
> leaving the slaved manifold — and (iii) — a Jensen gap in σ — are not needed** and are
> withdrawn as explanations of the cost. Either may still explain the *time* gap, which
> is now a separate and cleaner question.

> **The time gap does not vanish with Ω**, converging to ~1.16 at γ = 0.07, ~1.11 at
> γ = 0.20, ~1.08 at γ = 0.30 rather than to 1. So the exact MFPT to an absorbing
> threshold is persistently **shorter** than the deterministic arrival along the slaved
> manifold — a first-passage effect, not a finite-count one, and larger where the
> landscape is deeper.

**What this buys.** §39's closed form is correct in its thermodynamics and wrong only in
its kinematics: `Σ = Ω · σ̄ · T`, with σ̄ right and `T` the deterministic traversal instead
of the MFPT. So **any improvement to the first-passage time carries straight through to
the cost** — the two are now one problem rather than two, which is why the optimum
survived (§39's P4: γ = 0.240 analytic against §38's γ* ≈ 0.20 measured) even though the
magnitude did not: a smooth multiplicative factor varying 1.16 → 1.08 across γ moves the
minimum's location far less than its value.

**T-COST-d, open: why is the MFPT persistently below the deterministic traversal?** The
gap converges to a nonzero limit and grows with landscape depth. §39's candidates (ii)
and (iii) remain live *for the time*, joined by (iv): absorption at a threshold selects
the leading edge of the packet, so the first-passage time sits below the mean arrival by
an amount set by the packet width relative to the drift — which need not vanish when the
threshold is crossed on the steep part of the drift. **How to kill (iv):** compare the
MFPT against `∫dδ/μ` for thresholds placed at different θ. If the gap tracks the local
drift steepness at the threshold rather than the path as a whole, it is an absorption
effect and the deterministic traversal is exact away from the boundary.

### 39.2 The closed form is EXACT in the slaved limit — T-COST-d closed

§39.1 reduced the cost residual to a time residual: `∫dδ/μ` overestimates the MFPT by
2–18%, flat in Ω. Four candidates remained. **The cause is (ii), the slow-manifold lag,
and the reasoning was available before the run:** the manifold is defined by `ds/dt = 0`
*at fixed δ*, but as δ evolves the manifold moves and the pool lags behind it. That lag
is O(1/sep), **not O(1/Ω)** — which is exactly why the gap survived Ω → ∞.

Tested on §36's independent separation axis (scaling the pool pair X+Y↔2B, which has
`U·S = 0`, at fixed γ = 0.20):

| sep | 7.00 | 12.48 | 22.54 | 63.64 | 208.04 | 620.73 |
|---|---|---|---|---|---|---|
| T_det/MFPT − 1 | +0.0914 | +0.0553 | +0.0285 | +0.0096 | +0.0063 | **+0.0002** |

> **The gap closes: intercept at 1/sep → 0 is 1.00089, R² = 0.9977.** And it is not
> finite-count — mean |ratio(Ω=400) − ratio(Ω=800)| = **0.0048** across the whole sweep.
> **The deterministic traversal along the slaved manifold is exact in the slaved limit.**

**So the cost of restoration has a closed form that is exact where the reduction is:**

> **Σ = Ω · ∫ σ(δ)/μ(δ) dδ**, with the entropy rate exact (§39.1) and the traversal exact
> as sep → ∞ (here), the leading correction being the slow-manifold lag.

**The correction scales as 1/sep, and its coefficient does NOT transfer between axes.**
Along the T axis `(T_det/MFPT − 1)·sep = 0.6465 ± 0.0285`, constant to **12.3%** over a
9× range in sep (the two largest sep points are excluded and reported: their gaps are
below 0.008 and at the numerical resolution). Carried to the γ axis it predicts 16.2% at
γ = 0.07 against 16% measured — but **5.4% at γ = 0.30 against 8% measured, 33% off**.

> ⚠ **Rule 9, and I checked rather than assumed.** The T axis and the γ axis both raise
> sep, and a law calibrated on one need not carry to the other, because scaling the pool
> pair deforms the network as well as the separation. **The 1/sep *scaling* is
> established on both; the coefficient 0.6465 is a T-axis value and is not universal.**
> Quoting `T_det/MFPT = 1 + 0.65/sep` as a general law would be exactly §28.2's error.

**What the cost arc now establishes, end to end.** §37 asked what reliability costs and
found the question malformed. §38 reframed it as cost per e-fold of gain and located an
optimal drive. §39 predicted the cost in closed form to ~6% and predicted the optimum
analytically. §39.1 showed the whole residual was the clock, not the entropy. §39.2 shows
the clock is exact in the slaved limit. **The result is that a restoring chemical stage
has a closed-form thermodynamic cost, exact as sep → ∞, with a measured O(1/sep)
correction — and an optimal drive at γ ≈ 0.20–0.24 confirmed independently by an exact
CME solve and by minimising the closed form.**

**Taken with §36, both halves of the founding question are now closed-form:**

| | closed form | accuracy |
|---|---|---|
| reliability | `−ln P = 2Ω·V(x₀)`, `V = ∫μ/D` | 0.5–1.8% on-manifold (§36) |
| cost | `Σ = Ω·∫σ/μ dδ` | exact as sep → ∞, ~9% at γ = 0.20 (§39.2) |

Neither was true at the start of this session: the first was thought to be 7.5–15.5% off
and the second did not exist.

### 36.1 T14-f refuted, and its premise was my own error — `manifold_residual.py`

§36 left a 0.5–1.8% residual and named T14-f — the 2-D minimum action — as the only
candidate, arguing the residual should be **one-signed and ≥ 1** because *a minimum over
all paths cannot exceed the value along the slaved one*. Two candidates were eliminated
first:

* **The pool wobble.** Realised `b` varies on the integer lattice exactly as ε does, and
  §36 showed the rate is far more sensitive to `b` than to ε — the obvious suspect.
  Measured: the wobble is **0.00–0.53%**, and adding a `b` regressor moves the rates by
  **under 0.03%**. Dead.
* **The fitting window.** Over six Ω windows spanning 200–1200 at 7 points, pred/meas
  ranges by **1.30–2.29%** per γ — *larger than the residual being explained*. At short
  windows the residual is entirely instrument.

Redone at §35 grade (Ω to 1800, 11 points, on-manifold start, ε **and** b controlled,
with a half-split precision check):

| γ | sep | decades | measured | WKB pred | pred/meas | half-split |
|---|---|---|---|---|---|---|
| 0.20 | 7.00 | 41.5 | −0.057243 | −0.057031 | 0.9963 | 0.42% |
| 0.25 | 9.00 | 29.7 | −0.040072 | −0.039970 | 0.9975 | 1.39% |
| 0.30 | 12.00 | 19.7 | −0.026387 | −0.026177 | 0.9921 | 1.27% |
| 0.35 | 17.00 | 11.6 | −0.015545 | −0.015266 | **0.9820** | **0.06%** |

> ⚠ **T14-f is REFUTED twice over.** The residual is one-signed **below** 1 — 0/4 above —
> which the minimum-action argument forbids outright. And the P4 ordering is **backwards**:
> the residual *grows* with separation (−0.37% at sep 7 → −1.80% at sep 17) where a
> deviation from the slaved manifold must *shrink* as the manifold becomes more strongly
> attracting.

> ⚠⚠ **The premise was wrong, and it was mine.** I argued the 1-D slaved result bounds
> the 2-D action from above because "a minimum over paths cannot exceed the value along
> the slaved one". **That is a category error.** The 1-D reduction is not the action along
> a particular 2-D path — it is the exact WKB action of a *different process*, a
> birth–death chain built from rates projected onto the manifold, whose momentum is
> conjugate to δ alone rather than to both coordinates. **No variational inequality
> relates them**, so a residual below 1 was never forbidden and the sign test I built the
> experiment around had no content. The measurement is fine; the reasoning that gave it
> meaning was not.

> **Precision, stated against the claim.** Half-split disagreement runs 0.06–1.39%.
> At γ = 0.25 and 0.30 the residuals (−0.25%, −0.79%) sit **at or below** their own
> precision and are not resolved. **Only γ = 0.35 is clearly resolved** (−1.80% against
> 0.06%) — and it is the cell with the fewest decades (11.6), which is exactly where a
> systematic would be most suspect.

**What stands.** §36's central finding is untouched and is if anything strengthened: the
on-manifold comparison gives **0.9820–0.9975, mean 0.9920**, against the off-manifold
1.0805–1.1678. §15's closed form agrees with the exact collapse to **within 2%, one-signed,
like-for-like** — the 14.8× improvement is not in question. What is withdrawn is the
*explanation* offered for the last percent, and the argument that made it seem necessary.

**T14-f is closed as refuted. T14-g, open: what is the ~1% one-signed deficit?** It grows
with sep and with γ, i.e. it is largest where the landscape is shallowest and the fit
shortest — so the leading candidate is now **instrumental, not physical**: the rate at
γ = 0.35 is extracted from 11.6 decades against 41.5 at γ = 0.20. **How to kill:** hold
the decade count fixed across γ by choosing per-γ Ω ranges (§35.2's matched-decade
construction, which is already written) and re-measure. If the deficit flattens, it is
lever-arm; if it survives at matched decades, it is physics and needs a mechanism that is
not the 2-D action.

### 40 The first external standard: AM sits ~5× from the thermodynamic bound

Every cost number in §37–§39 is a measurement with no external reference — nothing said
whether AM is a *good* decision element or merely a measurable one. A thermodynamic
uncertainty relation supplies one, because it is a **bound rather than a fit**. The
first-passage TUR bounds a decision's timing precision by its dissipation,
`Var(T)/⟨T⟩² ≥ 2/⟨Σ⟩`, so

> **Q = (Var(T)/⟨T⟩²) · ⟨Σ⟩ / 2 ≥ 1**, with Q = 1 at saturation.

Both sides were already exact and already built: `first_passage_moments` solves
`Q_tt m₂ = −2T` for Var(T) (added this session, cross-checked against the SSA), and §37's
`Q_tt Σ = −σ_local` gives ⟨Σ⟩. Same generator, same absorbing set, on-manifold start
(§36), θ scaled with δ*.

**P3: the bound HOLDS at all 32 cells**, Q ∈ [5.39, 180.9]. That is a genuine external
validation of the whole §37–§39 apparatus — an independent inequality that the exact
entropy solve and the exact first-passage moments had to satisfy together, and did.

> **P5: AM is not near the bound.** The closest approach is **Q = 5.39**, i.e. AM
> dissipates roughly **5.4× more than the thermodynamic minimum** for the timing
> precision it achieves. It is not a thermodynamically optimal decision element, and its
> ubiquity is presumably about something else — robustness, simplicity, or speed.

| γ | 0.05 | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.35 | 0.40 |
|---|---|---|---|---|---|---|---|---|
| Q | **5.43** | 5.74 | 6.76 | 8.79 | 12.48 | 20.79 | 44.01 | 152.6 |

> ⚠ **P4: the TUR optimum and the cost optimum are DIFFERENT, and that is informative.**
> Q falls monotonically toward small γ with its minimum **at the grid edge** (γ = 0.05,
> not bounded below), where §38's cost per e-fold has an **interior** minimum at
> γ ≈ 0.20. **They are two objectives with two optima**: driving harder makes the
> decision more deterministic (timing variance falls faster than Σ rises), which the TUR
> rewards without limit, while gain-per-dissipation pays a ln(1/γ) penalty that
> eventually dominates. §38's optimum is real and is *not* the thermodynamic one.

> ⚠ **P2 FAILS at large γ and the cells are flagged.** Q's spread across Ω is 1.7–4.7%
> for γ ≤ 0.20 but rises to 16.6 / 30.8 / 35.6 / 22.7% at γ = 0.25–0.40. Q should be
> Ω-independent (Var/⟨T⟩² ~ 1/Ω against Σ ~ Ω), so **the shallow-landscape cells are not
> converged** and their Q values are indicative only. The small-γ end, where the
> conclusion lives, is converged.

> **What the margin does and does not establish.** The bound holding at Q ≥ 5.4 is robust
> — applicability concerns could only matter if it were *violated*. But reading 5.4 as
> "5.4× from optimal" assumes this TUR form is **tight** for our setting, and ours has a
> **two-sided** absorbing set (|δ| ≥ thr) where the standard statement is one-sided.
> **If the correct bound for two-sided absorption is tighter, AM is closer to optimal
> than 5.4× suggests.** That is not established here and the factor should be read as an
> upper bound on the gap, not a measurement of it.

**What §40 settles.** The cost framework passes its first external check. AM is within an
order of magnitude of the thermodynamic limit on decision precision but does not
approach it, and the drive that minimises cost per e-fold of gain is **not** the drive
that comes closest to the bound — so §38's design principle is about amplification
economics, not thermodynamic optimality.

**T-TUR-a, open: is the two-sided bound tighter?** The gap of 5.4× rests on a one-sided
inequality applied to a two-sided absorbing set. **How to kill:** re-run with a one-sided
absorbing condition (`δ ≤ −thr` only, the error boundary) so the standard TUR applies
verbatim, and compare Q. If Q drops toward 1, the gap was the boundary convention and AM
is far closer to optimal than §40 reports.

### 40.1 The bound choice is robust, and AM is ribosome-grade

§40 used `CV² ≥ 2/⟨Σ⟩` and flagged the form as an assumption. Two literature results settle it.

**The bound form.** A first-passage TUR valid for **arbitrary initial conditions and
systems with absorbing states** (Pal, Reuveni & Rahav, arXiv:2103.16578) gives
`CV² ≥ 1/(Σ/2 + 1) = 2/(Σ+2)` — precisely our setting, where the standard steady-state
derivation does not apply. §40's form is its large-Σ limit, and with Σ ~ 900–1600 the
refinement moves Q by **0.13–0.22%**:

| γ | 0.05 | 0.20 | 0.35 |
|---|---|---|---|
| Q (§40's form) | 5.4251 | 8.7880 | 44.0074 |
| Q (absorbing-state form) | 5.4329 | 8.8072 | 44.0839 |

**§40's conclusion is unaffected: Q_min = 5.39 → 5.40.**

> **The calibration, which changes how §40 should be read.** Song & Hyeon
> (Phys. Rev. E **101**, 022415) measured the TUR distance for two real enzymes: **T7 DNA
> polymerase operates close to the bound, while the E. coli ribosome operates ~5× from
> it.** AM at **Q = 5.39** therefore sits **where the ribosome sits** — not "far from
> optimal", as §40 phrased it, but at the same distance from the thermodynamic limit as
> one of the most heavily optimised machines in biology.

> ⚠ **The comparison is suggestive, not exact.** Their Q is a steady-state KPR cycle
> quantity; ours is a first-passage decision. Both are dimensionless TUR distances, but
> they are not the same observable, and no claim of numerical equivalence is made. What
> the anchor establishes is the **scale**: a factor of ~5 is where real, selected
> biochemistry lives, so §40's number is unremarkable rather than damning.

**§40's wording is corrected**: "AM is not a thermodynamically optimal decision element"
stands; "its ubiquity is presumably about something else" does not follow, since the
comparison class it is being measured against also sits at ~5×.

### 41 An exact identity validates the entropy machinery — T-TUR-b

§37–§40 all rest on one object: `σ_local(n) = Σ_j a_j ln[a_j/a_rev]`, the reverse pairing
that defines it, and its sign convention. **Nothing had ever checked that object against
anything external.** §40 checked an *inequality*, which an entropy wrong by a constant
factor would still satisfy.

THEORIES §5 recorded Neri's martingale result and suggested testing `p₋ = exp(−ℓ₋)`, which
needs an augmented chain. **The same physics has a cheaper form:** the integral
fluctuation theorem `⟨e^(−S_tot)⟩ = 1` holds at *any* stopping time, including the
absorption this project already uses. And it is a linear solve, because tilting each
transition by its own entropy weight collapses:

> **a_j(n)·e^(−Δs_j) = a_j(n)·a_rev(n')/a_j(n) = a_rev(n')**

so the tilted generator is built from **reverse propensities**, at the same cost as every
other solve here. With `e^(−S_tot) = e^(−S_med)·π(n_T)/π(n₀)`, two solves differing only
in their absorbing boundary give the medium and total forms.

**P1 confirmed — the medium term alone comes nowhere near 1**, spanning
**3.7×10⁻¹¹ to 0.92** across the grid. Quoting it as though it satisfied the IFT is the
error this experiment exists to exclude.

**P2: the identity holds.** Over 36 cells (γ × Ω × ε × θ):

| | value |
|---|---|
| best cell (γ=0.30, Ω=40, ε=0.50) | **\|dev\| = 5.5×10⁻¹⁴** |
| median over all 36 cells | **1.33×10⁻⁹** |
| solve residual throughout | ~10⁻¹⁵ |

**σ_local, the reverse pairing and the sign convention are validated together against an
exact identity** — the first equality, as opposed to inequality, that this project's
entropy machinery has faced.

> ⚠ **It fails at γ = 0.10, and the failure is conditioning, not physics.** Deviation by
> γ runs **1.3×10⁻¹² (0.30) → 1.5×10⁻⁹ (0.20) → 4.2×10⁻³ (0.10)**, and one γ = 0.10 cell
> overflows to 5.9×10²⁸⁰. That tracks the dynamic range exactly: `E[e^(−S_med)]` runs
> 0.92 → 2.3×10⁻⁴ → 3.7×10⁻¹¹ over the same γ, and π spans correspondingly more decades
> as the barrier deepens. **The tilted solve carries the full range of π, so a deep
> landscape destroys it.** Reported, not dropped; the identity is verified where the
> arithmetic can represent it.

> **The test has teeth, demonstrated by catching me.** The first pass coded the boundary
> as `π(n₀)/π(n)` — the reciprocal of the correct `π(n)/π(n₀)`. The identity returned
> **0.0008 instead of 1**, a 1000× deviation, with a solve residual of 10⁻¹⁵ proving the
> linear algebra was exact and the *convention* was wrong. A sign error in the entropy
> convention is precisely the failure mode §37–§40 were exposed to, and this test finds
> it immediately.

**T-TUR-b closes.** The tight martingale route is available in this cheaper form and the
machinery passes it. **T-TUR-c, open: does `p₋ = exp(−ℓ₋)` hold?** That still requires
thresholding on accumulated entropy production rather than on δ — an augmented chain
(state × accumulated S) — and is now better motivated, since the identity above shows the
entropy bookkeeping is correct and the remaining question is genuinely about the
first-passage geometry rather than about conventions.

### 36.2 The 1-D reduction is ~1% shallow, and it is real — T14-g closed

§36.1 left a ~1% one-signed deficit and concluded "the leading candidate is now
instrumental", naming the matched-decade construction as the kill test. **Three sweeps
were needed because each one held a different thing fixed**, and the first two agreed only
by accident.

| sweep | held fixed | varies | γ=0.20 | 0.25 | 0.30 | 0.35 | spread |
|---|---|---|---|---|---|---|---|
| **B** (§36.1) | Ω window 150–1800 | decades 41.5→11.6 | 0.9963 | 0.9975 | 0.9921 | 0.9820 | 1.55 |
| **A** | decades ≈ 12 | Ω-ratio 4.1→12.4 | 1.0003 | 0.9999 | 0.9953 | 0.9839 | 1.64 |
| **C** | start depth **and** span | — | 0.9937 | 0.9949 | 0.9888 | 0.9871 | **0.78** |

> ⚠ **§36.1's "instrumental" reading is WITHDRAWN.** The deficit survives at matched
> decades (sweep A, spread 1.64 against B's 1.55), so it is not lever-arm. B and A
> confound oppositely — B fixes the Ω-window and varies decades, A fixes decades and
> varies the Ω-ratio — and the deficit tracks **γ** in both. Neither explains it.

> **But both silently shared the window's lower edge at Ω = 150**, where −ln P is 8.6 at
> γ = 0.20 and only **2.3** at γ = 0.35 — P ≈ 0.10, not a tail at all. Every large-γ fit
> was reaching into a region where the Laplace/WKB asymptotic has not taken hold, and
> matching *decades* constrains the span, not the starting depth. **Sweep C starts every
> window at P ≈ 10⁻⁴ and halves the spread, 1.64 → 0.78 points.** So half the apparent
> γ-dependence was asymptotic validity, and it took a third sweep to see it.

**What survives is a genuine result:** with all four γ compared in the regime where the
prediction is supposed to apply, **§15's 1-D closed form runs 0.9911 ± 0.0039 against the
exact 2-D collapse — one-signed, ~1% shallow.** Half-split precisions are 0.23 / 0.78 /
0.44 / 0.65%, so the deficit is resolved at **three of four** cells (γ = 0.25 is not).

> **A negative deficit is permitted, and §36.1 is why.** The 1-D reduction is not the
> action along a constrained 2-D path — it is the exact WKB action of a *different*
> process, a birth–death chain built from projected rates, so no variational inequality
> forces pred/meas ≥ 1. That correction was made when the sign test failed; here it is
> what licenses the result.

**T14-g closes.** The 1-D slaved reduction is shallow by ~1%, consistently and for real,
with a residual γ-dependence of 0.78 points now at the edge of resolution. **§36's
headline is strengthened, not weakened:** the on-manifold comparison is now 0.9911 ± 0.0039
on properly matched windows against 1.0805–1.1678 off-manifold — a like-for-like agreement
of **~1%** where fourteen sections read 7.5–15.5%.

> **Methodological, and it cost three runs.** "Lever arm" turned out to be *two* hidden
> variables — the span and the starting depth — and matching the first left the second
> free. Rule 9 says measure along an axis you did not choose; this says check how many
> axes you did not choose, because a matched-looking comparison can still share an
> unexamined edge.

---

### 42 The identity does not care about conservation laws — T15-c

§30 proved the pairwise identity for every n, γ and pair (i,j) in the AM family:

    d(n_i − n_j)/dt = (n_i − n_j)·(k/Ω)·[ n_B − Σ_{l≠i,j} n_l − γ(n_i + n_j − 1) ]

§31 then broke it with `am_asymmetric`'s tilt β. **Neither section could say which feature
the identity actually needs**, because in the whole AM family every network has exactly one
conservation law (X+Y+B = Ω). Physical networks carry several at once — charge, baryon
number, lepton number — so if the identity is a fact about conservation structure, §30's
restoration reading is far narrower than it sounds.

Reading the cancellations suggests otherwise. The identity survives because the pair's own
disagreement reaction consumes both equally, its reverse reaches both through the same
channels, and each recruitment is linear in the species it recruits **with the same
coefficient for i and j**. Only the last mentions the rate law, and it is exactly what β
breaks. Nothing in the derivation counts conservation laws.

`experiments/conservation_identity.py` crosses the two features. The cofactor family adds
conserved partner pairs *consumed and regenerated by the signal-moving reactions
themselves* — B + X + D → 2X + E and its mirror — so the extra laws are not decorative.
The probe is §30's own: hold everything fixed, hold s = n_X + n_Y fixed, vary the split,
and see whether d(δ)/dt ÷ δ moves. γ = 0.25, β = 0.20, Ω = 90, 28 random states each.

| | network | laws | exchange | worst spread | §30 identity |
|---|---|---|---|---|---|
| A | `am_reversible` | 1 | symmetric | 1.7×10⁻¹³ | **HOLDS** (control) |
| B | cofactor D/E on both | **2** | symmetric | 1.6×10⁻¹³ | **HOLDS** |
| C | cofactor D/E **+ F/G** | **4** | symmetric | 4.5×10⁻¹⁴ | **HOLDS** |
| D | cofactor, β = 0.2 on Y | 2 | rate broken | 1.9×10¹ | FAILS |
| E | cofactor on X only | 2 | struct. broken | 2.9×10² | FAILS |
| F | `am_asymmetric` | 1 | rate broken | 2.4×10² | FAILS (control) |

**P1–P4 all confirmed. The identity survives two and four conservation laws and dies the
moment exchange symmetry breaks.** Row **D is the load-bearing one**: it differs from B in
a *single rate constant* — same species, same reactions, same orders, same two conservation
laws — so nothing structural is available as an alternative explanation. E shows the same
break reached structurally instead.

D's tilt scales Y's recruitment **and its reverse** by the same (1+β), so the
forward/reverse ratio stays γ exactly as in B. The failure is therefore not a
broken-reversibility or thermodynamic-consistency artifact — the only thing that changed is
that X and Y stopped being each other's mirror image.

Per P6, the verdict is not "nothing changed": the implied bracket does pick up the cofactor
counts, spanning −0.184…0.694 (B) and −0.179…0.474 (C) against −0.164…0.517 (A). The
*bracket* moves; the *proportionality* does not.

---

### 43 Exchange symmetry forces divisibility — and §30's constant ratio is a degree accident

Reading §42 back exposes a distinction §30 never had to make, because in the AM family the
two properties coincide. **An antisymmetric polynomial in two variables is divisible by
their difference.** If a mass-action network is invariant under swapping i and j then
b_i(σn) = b_j(n), so b_i − b_j is antisymmetric under n_i ↔ n_j; propensities are
polynomials in the counts (falling factorials are polynomials); hence

> **b_i − b_j = (n_i − n_j)·P(n),  with P symmetric.**  (\*)

That is a **theorem, not a measurement**, and it holds at any reaction order and for any
number of conservation laws — including none.

> ⚠ **§65's literature check qualifies the novelty claim, not the statement.** The invariance
> half — that the diagonal n_i = n_j is flow-invariant — is a **standard result in equivariant
> dynamics**, the flow-invariance of the fixed-point subspace Fix(Z₂), called a *folk theorem*
> in that literature (Golubitsky, Stewart & Schaeffer 1988). Given it, divisibility follows in
> one algebraic step, because b_i − b_j is a polynomial vanishing on the irreducible variety
> {n_i = n_j}. **This section presented the invariance as the discovery; it is not.** What is
> not standard is the explicit *form* of P proved in §54 and its consequences in §56/§62 — see
> §65 for the theorem as it should be stated, with its prior art. Its consequence is the entire restoration
claim: δ = 0 is an invariant manifold of the deterministic flow, so **the sign of a lead is
a deterministic invariant and every reversal is a fluctuation.**

§30's identity is strictly stronger. P is symmetric, so P = P(s, δ²); independence of the
split needs P to have degree 0 in δ², i.e. the pair to enter every reaction at total degree
≤ 2. AM's highest pair term is 2X → B + X, giving P ~ (s−1) — which is why §30 could not
tell the two apart. **A cubic pair term separates them.** For the symmetric pair
X + 2Y → 3B and 2X + Y → 3B,

    b_X − b_Y ~ c·n_X n_Y[(n_Y−1) − (n_X−1)] = −c·n_X n_Y·δ,   P = −c(s² − δ²)/4

still divisible, no longer constant. `experiments/exchange_theorem.py`, γ = 0.25, Ω = 90:

| network | divisibility residual | divisible? | ratio spread | §30 ratio? |
|---|---|---|---|---|
| `am_reversible` (sym, degree 2) | **0.0** | YES | 1.2×10⁻¹³ | YES |
| `am_cubic` (sym, degree 3) | **2.8×10⁻¹⁷** | **YES** | 1.7×10¹ | **NO** |
| `am_asymmetric` | 1.9×10⁻¹ | NO | 8.4×10¹ | NO |

**P2: the properties separate.** Divisibility is the weaker one, it is what restoration
actually needs, and §30's constant ratio is the special case of pair-degree ≤ 2.

**P3 is an absolute check, not a fit (rule 16).** The cubic ratio is affine in δ² to a max
relative residual of **1.5×10⁻¹⁵**, and the fitted δ² coefficient is

| | value |
|---|---|
| fitted | 3.086420×10⁻⁵ |
| derived by hand, k/(4Ω²) = 1/(4·90²) | 3.0864198×10⁻⁵ |
| ratio | **1.00000008** |

**P4/P5 — the theorem on 200 random networks each way.** Networks with random species
counts, random orders up to 4 and random rate constants, symmetrised by construction versus
not:

| | divisible | worst residual | median residual |
|---|---|---|---|
| symmetrised | **200 / 200** | 3.0×10⁻¹⁷ | 0.0 |
| not symmetrised | 2 / 200 | 1.0 | 0.68 |

**No counterexample to (\*), and the probe demonstrably has power.** Per P6, broken out by
conservation-law count: **0 laws 157/157, 1 law 36/36, 2 laws 7/7.** The largest group
conserves *nothing at all* — so conservation is irrelevant to the identity in both
directions: extra laws do not break it and zero laws do not either.

> **Two limits, stated so the theorem is not read as more than it is.** (i) It is a
> statement about the **deterministic drift**. In the CME the same divisibility says only
> that the *mean* jump vanishes at δ = 0; the chain crosses δ = 0 constantly, which is
> precisely the "every reversal is a fluctuation" framing and nothing beyond it.
> (ii) Divisibility gives no-reversal, **not** amplification — P may be negative, and in
> the tables above it is (the bracket ranges start at −0.16, the sub-separatrix region).
> Restoration is divisibility *plus* sign(P) > 0, and only the first half is a theorem.

**What this answers.** The question was whether re-basing the model on particle-like
species, which carry several conserved charges simultaneously, would cost the restoration
result. It would not — and the reason is now sharper than "it survives": conservation
structure was never what the result rested on. What it rests on is that the two competing
species are each other's mirror image, which is a symmetry a particle model can keep or
break deliberately.

---

### 44 Real Arrhenius kinetics, and a free lever nobody had turned — T-COST-e

`crnl/cooling.py` maps the drive to a temperature, γ = exp(−ΔE/T), and its own docstring
flags what it leaves out: *"the forward rates have no temperature dependence, so this is
the minimal change that lets the balance move, not a thermochemistry."* Real forward rates
are Arrhenius too. If T moves them as well as the balance, then the optimal *temperature*
is not §38's optimal *drive* re-labelled, and the "which substrate should a restoring
element be built from" reading of §38 is unearned.

Writing k_r = A·exp(−E_r/T) splits T's job in two:

* A **uniform** activation energy is a pure clock change — Q → λQ and σ → λσ together, so
  Σ = Q⁻¹σ is exactly invariant. Temperature acting on all forward rates equally *cannot*
  move the cost optimum. §5.1's uniform-order argument, one level up.
* The channels need not share a barrier, and their ratio **ρ = k_dis/k_rec = exp(−ΔEa/T)**
  is a landscape change, not a clock change. **Every result in this document is at ρ = 1**,
  and nothing had ever checked that.

**The geometry, derived by hand and checked before use.** b\* = γ/(1+γ) is untouched (§30's
identity survives: ρ multiplies the disagreement channel, which cancels in the difference),
giving

    δ*(γ,ρ)² = (ρ − γ − 4ργ³) / [(1+γ)²(ρ − γ)]

which reduces to `delta_star(γ)` at ρ = 1. **It has a trap and the trap was hit.** For
ρ < γ both numerator and denominator are negative, so it returns a *plausible positive*
value — 0.806 at (γ,ρ) = (0.25, 0.05), and 6.37 and 11.4 near ρ = γ where δ\* ≤ 1 is a hard
bound — in a region where the ODE nullcline says the landscape is simply absent. The
landscape exists iff ρ > ρ_c = γ/(1−4γ³), and ρ_c > γ always, so one guard covers it.
Caught by checking the closed form against the nullcline at 30 (γ,ρ) pairs *before* using
it anywhere. Rule 10, in the cheapest possible place.

**P1a — clock invariance is exact.** Rescaling every rate by λ ∈ {0.1, 1, 7, 100} leaves
Σ = 311.092663 unchanged to a worst relative deviation of **4.6×10⁻¹⁵**. A uniform Arrhenius
factor therefore cannot move anything, as the algebra requires.

**P1b — the instrument reproduces §38.** At ρ = 1, sweeping γ (Ω = 200 / 300):

| γ | 0.08 | 0.12 | 0.16 | 0.20 | 0.24 | 0.28 | 0.32 | 0.36 | 0.40 |
|---|---|---|---|---|---|---|---|---|---|
| G (Ω=200) | 2.382 | 2.101 | 1.963 | 1.897 | **1.894** | 1.975 | 2.112 | 2.601 | 3.566 |
| δ\* | 0.925 | 0.889 | 0.854 | 0.817 | 0.777 | 0.732 | 0.681 | 0.619 | 0.541 |

γ\* = **0.2216** (Ω=200) and **0.2319** (Ω=300), sitting between §38's CME 0.20 and §39's
closed-form 0.240.

### 44.1 P3 was refuted, and in the opposite direction

**Predicted: cost has an interior minimum in ρ at ρ\* < 1**, on the reasoning that the
disagreement reaction moves δ by exactly zero — §30's first cancellation, it consumes both
species equally — so it produces entropy and no signal, and should be suppressed.

**Measured: G falls monotonically as ρ *rises*, at all three γ and both Ω, with no interior
minimum anywhere on ρ ∈ [0.35, 100].** The prediction was not merely wrong, it had the sign
backwards.

**The reasoning conflated flux with dissipation.** For the disagreement pair,
a_f/a_r = n_X n_Y/(γ n_B²) is *independent of ρ* — raising ρ scales both directions
equally. Driving a reaction fast in both directions pushes it onto its own local
equilibrium, where n_X n_Y → γ n_B² and ln(a_f/a_r) → 0, so the channel carries unbounded
flux at *bounded* net entropy production. Making the non-signal channel fast makes it
cheap, not expensive.

### 44.2 ρ is a free lever — cost halves, reliability doubles, time halves

At γ = 0.16, Ω = 200, sweeping ρ with everything else fixed:

| ρ | 0.35 | 1.0 | 1.8 | 6.0 | 20 | 100 |
|---|---|---|---|---|---|---|
| G (k_B/molecule/e-fold) | 4.135 | 1.963 | 1.535 | 1.247 | 1.167 | **1.117** |
| L (nats of reliability) | 7.67 | 17.18 | 22.65 | 29.03 | 32.14 | **34.38** |
| mean time | 17.68 | 5.41 | 3.84 | 2.98 | 2.80 | **2.70** |
| δ\* | 0.8490 | 0.8536 | 0.8543 | 0.8548 | 0.8549 | 0.8550 |

**All three move the right way at once, and P5's control shows why that is not a geometry
effect: δ\* moves by 0.16% while cost falls 43%.** Going from AM's ρ = 1 to the asymptote:

| γ | ΔG | L ratio | Δtime | Δδ\* |
|---|---|---|---|---|
| 0.16 | **−43.1%** | **2.00×** | −50.1% | +0.16% |
| 0.24 | **−46.9%** | **2.07×** | −54.0% | +0.93% |
| 0.32 | **−50.5%** | **2.05×** | −57.4% | +3.7% |

Cost halves, reliability doubles, time halves — consistently across γ.

**And this makes ρ a different *kind* of knob from γ.** Over the γ-sweep L falls
monotonically (32.1 → 13.2 → 1.7 as γ goes 0.036 → 0.209 → 0.449): **γ trades cost against
reliability, which is what makes §38's optimum an optimum.** ρ does not trade — it is
strictly dominating. The γ lever is worth ~21% in G from its worst reasonable setting;
**the ρ lever is worth 43–50% and costs nothing.** Nothing in this project had turned it.

**P9 confirmed — G asymptotes**, which is the signature the local-equilibration account
requires. The last step (ρ = 50 → 100) moves G by 0.5–0.6% at γ = 0.24 and 0.32, against
2.8% for the step before at γ = 0.16. The limit is the reduced model in which X+Y ⇌ 2B is
equilibrated and recruitment carries the whole net current.

> ⚠ **The last decade is at the resolution floor.** At ρ = 100 the Ω=200 and Ω=300 columns
> disagree by 1.2% at γ = 0.16 *in the wrong direction*, and the γ=0.32 Ω=200 column is
> non-monotone at the 1% level over ρ = 20…100. The generator's conditioning degrades as ρ
> grows. The asymptote is ~1.0–1.1; its third digit is not resolved, and the two "interior"
> ρ\* values the fitter reports (33.3 and 63.2) are **noise picking a point on a flat tail,
> not optima** — reported as such rather than as a located optimum.

> **Suspect, not result (rule 17).** ρ raises reliability while δ\* is frozen, which points
> at the timescale separation: faster disagreement means a faster pool, which is deeper
> slaving, which is §36's on-manifold condition and §39.2's 1/sep law. **Kill test:**
> compute sep(γ,ρ) and check whether §39.2's law predicts the L improvement quantitatively.
> If it does not, the mechanism is wrong and only the measurement stands.

### 44.3 The pre-registered P4 failed; its argmin was inadmissible

Sweeping T with γ(T) = exp(−ΔE/T) and ρ(T) = exp(−ΔEa/T), ΔE = 1, T_c = ΔE/ln2 = 1.4427:

**P2 holds.** At ΔEa = 0, T\* = 0.6367 / 0.6458 against the predicted ΔE/ln(1/γ\*) =
0.6636 / 0.6842 — ratios 0.96 and 0.94, the shortfall being parabolic interpolation on a
nonlinearly transformed grid rather than physics. **T\*/T_c = 0.441 / 0.448: the optimum
sits at ~44% of the temperature at which the landscape dies.**

**P4, as pre-registered, FAILS** — T\* did not split by sign(ΔEa): +75% and +0.7% at Ω=200.
**But its ΔEa = +0.6 argmin is not admissible.** It sits at T = 1.114, where δ\* = 0.227
and ρ/ρ_c = 1.04 — the landscape is 4% from death and L has fallen to 1.31 nats. **§9.2 was
withdrawn for exactly this**: G falls there not because the decision is cheap but because
it is *small*, the threshold θ·δ\*·Ω having collapsed with δ\*.

Under the admissibility floor δ\* ≥ 0.40 — fixed *before* the re-run, and chosen because the
ρ=1 γ-sweep that P1b validates against §38 itself spans δ\* = 0.40…0.97, so that is the
range where this instrument is known to agree with the established result:

| criterion | Ω | T\*(ΔEa=0) | ΔEa=+0.6 | ΔEa=−0.6 | verdict |
|---|---|---|---|---|---|
| P4 pre-registered | 200 | 0.6367 | 1.1143 (**+75.0%**) | 0.6415 (+0.7%) | FAILS |
| P4 pre-registered | 300 | 0.6458 | 1.1143 (**+72.5%**) | 0.6855 (+6.1%) | FAILS |
| P8, δ\* ≥ 0.40 | 200 | 0.6367 | 0.5839 (**−8.3%**) | 0.6415 (+0.7%) | HOLDS |
| P8, δ\* ≥ 0.40 | 300 | 0.6458 | 0.5910 (**−8.5%**) | 0.6855 (+6.1%) | HOLDS |

**Both verdicts stand in the record.** P4's is what the pre-registered criterion returned;
P8's is what the same data say once cells outside the instrument's validated range are
excluded. The substance — **temperature does not act through γ alone** — rests on P8.

> ⚠ **The split is carried by one arm, and read per cell it is asymmetric.** The +0.6 arm
> moves T\* down by 8.3% and 8.5%, consistent across Ω. The −0.6 arm moves it up by 0.7%
> and 6.1% — *not* resolved between Ω. So "opposite directions" is a one-sided result.
> **It is also the expected one**, and that is a consistency check rather than an excuse:
> §44.2 shows G falls steeply in ρ below ~5 and is flat above it. At ΔEa = −0.6, ρ ranges
> 7.4 → 1.6 across the sweep, entirely in the flat region, so the ρ channel has almost
> nothing to contribute; at ΔEa = +0.6, ρ ranges 0.135 → 0.58, entirely in the steep one.
> The asymmetry is predicted by the ρ-sweep, not an artifact of it.

**P6, arithmetic and labelled as such.** §16 pins the cycle affinity at A = 3 ln(1/γ) (it is
in `verify_base`), so γ\* = 0.2216 / 0.2319 gives **A\* = 4.52 / 4.39 k_BT**, an optimal fuel
drop of **0.117 / 0.113 eV at 300 K**, against ATP hydrolysis at ~20 k_BT. This is
arithmetic on §38 plus §16, not a new measurement. **It is not merged with §40's
Q_min = 5.39** — two numbers near 5 arriving along axes chosen for other reasons is rule 9's
trap, which has already sprung three times here.

---

### 45 ρ works partly through slaving — but sep is not the governing variable, and §39's 6% closes anyway

§44.2 named a suspect for ρ's free lever: a faster disagreement channel is a faster pool,
which is deeper slaving. This is its kill test, and it doubles as an attack on **T-COST-c**,
the project's outstanding unexplained number — §39's closed form sits ~6% off the exact
cost, flat in Ω over a 3× range and not attributable to σ.

**The design scan supplied the control before any prediction was written.** At γ = 0.20,
`sep` is **non-monotone in ρ**, with a minimum near ρ ≈ 1.5, while §44.2's cost G falls
monotonically across the same range. So the residual tracking *sep* and the residual
tracking *ρ* make opposite predictions, and the sweep can tell them apart.

**P1a — the instrument is anchored.** `sep_of` reproduces the closed form 3(1+2γ)/(1−2γ)
to **1.15×10⁻¹⁴** at ρ = 1 (3.976744186, 4.894736842, 7.000000000, 10.636363636,
17.000000000). **P1b** — pred/exact at ρ = 1 gives mean **1.0487**, range 0.9158–1.1357,
against §39's 1.0583 and 0.905–1.167. The closed form transplanted correctly.

**P2 HOLDS, and it is the real evidence.** The ρ knob at γ = 0.20:

| ρ | 0.5 | 0.75 | 1.0 | **1.5** | 2.0 | 3.0 | 4.0 | 6.0 | 8.0 | 16 | 32 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| sep | 11.70 | 8.04 | 7.00 | **6.73** | 7.19 | 8.73 | 10.56 | 14.45 | 18.47 | 34.82 | 67.76 |
| \|residual\| | 0.025 | 0.052 | 0.053 | **0.091** | 0.079 | 0.060 | 0.072 | 0.047 | 0.044 | 0.021 | **0.007** |

**The residual peaks at ρ = 1.5, exactly where sep bottoms out.** It is non-monotone in ρ
while cost is monotone, so it is tracking sep and not ρ. A monotone residual would have
looked like a confirmation and meant nothing; this is why the non-monotone design was worth
having.

**P6 survives: the residual closes.** |residual| averages **0.0826 at sep < 8** against
**0.0138 at sep > 20**, corr(|residual|, 1/sep) = **+0.78** over 32 cells. At ρ = 32 the
§39 residual is **0.7%**, down from 5.3% at ρ = 1 — a **7.9× reduction, and the first time
anything has moved T-COST-c.** This independently confirms §39.2 ("the closed form is exact
in the slaved limit") by reaching that limit with a different knob.

### 45.1 But the quantitative form is refuted, three ways

**P3 — matched pairs split, and the discriminating one fails.** The whole point was to move
ρ far at nearly fixed sep:

| pair | ρ ratio | sep difference | \|residual\| | verdict |
|---|---|---|---|---|
| ρ=1.0 (sep 7.00) vs ρ=2.0 (sep 7.19) | 2× | 2.7% | 0.041 vs 0.058 | agree |
| **ρ=0.5 (sep 11.70) vs ρ=4.0 (sep 10.56)** | **8×** | **10.3%** | **0.022 vs 0.067** | **DIFFER 3×** |

The 2× pair is a weak test and passes; the **8× pair is the one with power and it fails**.

**P4 — the two knobs do not collapse.** Over the overlap sep ∈ [6.73, 17.0]:

| knob | fit |
|---|---|
| ρ (γ = 0.20 fixed) | residual = **+0.441**/sep + 0.009 |
| γ (ρ = 1 fixed) | residual = **+1.147**/sep − 0.106 |

Slopes differ by **2.6×** and the intercepts have opposite signs. Two physically different
levers producing the same sep produce different residuals. **sep is not the governing
variable.**

**P5 — residual×sep is not constant.** It spans +0.200…+0.892 on the ρ knob (4.5× over a
10× range of sep). The 1/sep *form* is wrong even though the correlation is real.

> ⚠ **What the γ knob does at high γ is NOT resolved, and I am not claiming a sign change.**
> The γ-knob resid×sep range quoted by the run (−1.432…+0.540) is driven by one cell:
> γ = 0.35, where Ω=200 and Ω=300 give 0.9158 and 0.9951 — an **8% disagreement between Ω**,
> comparable to the residual being measured. §39's own range already straddled 1. So "the
> residual changes sign along γ" is a reading the data does not support at this resolution,
> and it is recorded here as unresolved rather than as a finding.

> ⚠ **§46 corrects the framing of this section in two places. The numbers stand; two
> readings do not.** (i) The closing paragraphs below hand the residual back to "T-COST-c's
> remaining candidates" — but **T-COST-c was closed by §39.1**, which showed the whole cost
> residual is a *time* residual and explicitly withdrew the off-manifold path and the Jensen
> gap as explanations of the cost. There were no remaining candidates to hand it to.
> (ii) The two-knob split reported as evidence *against* sep is **§39.2's published finding
> under rule 9** — that the 1/sep coefficient does not transfer between axes — arriving on a
> third axis. It is a confirmation read as a refutation. And this section tests
> `Σ_pred/Σ_exact` where §39.1/§39.2 test `T_det/MFPT`; they correlate at 0.95, not 1.
> **§46 re-runs the right quantity and finds something sharper.**

**Verdict. T-COST-f survives in its weak form and is refuted in its quantitative form.**
ρ does act partly through the timescale separation — the non-monotone signature is
unambiguous and the residual closes to 0.7% at large sep. But sep alone does not determine
the residual, so "ρ works by deepening slaving" is *incomplete*, not established.

> **The next suspect, with its kill test (rule 17).** `sep_of` measures the eigenvalue ratio
> at the **symmetric point** (x = y, b = b\*), but the traversal happens at δ ∈ [εδ\*, θδ\*],
> away from it. A separation evaluated *at one point* need not represent the separation
> *along the path*, and the ρ and γ knobs deform the manifold differently — which would
> explain both the matched-pair failure and the two-curve split without abandoning slaving.
> **How to kill:** compute a path-averaged separation over the actual traversal and re-run
> P3 and P4 against it. If the 8× matched pair then agrees and the two knobs collapse, the
> governing variable was the path separation all along. If they still split, slaving is not
> the mechanism and §44.2's lever keeps its measurement and loses its account.

---

### 46 The 1/sep scaling is axis-dependent too — a third axis breaks §39.2's law

§45 named the path separation as T-COST-f's successor: `sep_of` measures the eigenvalue
ratio at the *symmetric point*, while the traversal happens away from it. This tests it —
and on the way it corrects §45's framing, which was wrong about what was being asked.

**The averaging convention is forced, not chosen.** The slow eigenvalue crosses zero near
δ/δ\* ≈ 0.57, where the drift μ peaks, so sep(δ) *diverges mid-path*:

| δ/δ\* | 0.35 | 0.45 | 0.55 | 0.65 | 0.75 | 0.80 |
|---|---|---|---|---|---|---|
| sep (ρ=0.5) | 18.8 | 31.3 | **188.2** | 37.0 | 15.3 | 11.5 |

The arithmetic mean of sep is therefore meaningless and only a **harmonic** mean is
defined — which is also what the physics wants, since 1/sep corrections add along the path.
Three weightings computed, all reported (rule 15): uniform in δ, time-weighted (dδ/μ), and
cost-weighted (σ/μ dδ), with **cost nominated as primary in advance**.

**P1 gate passes — the distinction is real.** At γ = 0.20, ρ = 1 the path separations are
16.59 / 15.76 / 17.56 against a point value of 7.000, differing by **125–151%**. So the
experiment has content.

**And all three conventions fail, exactly as predicted (P5).**

| convention | ρ=0.5 vs ρ=4 sep | matched? | ρ slope | γ slope | collapse? |
|---|---|---|---|---|---|
| point | 11.70 vs 10.56 (1.11×) | yes | 0.441 | 1.147 | 2.60× split |
| path | 28.60 vs 26.46 (1.08×) | yes | 0.952 | 2.657 | 2.79× split |
| time | 27.27 vs 25.21 (1.08×) | yes | 0.901 | 2.515 | 2.79× split |
| cost | 29.61 vs 27.82 (1.06×) | yes | 1.007 | 2.909 | 2.89× split |

The 8× matched pair has residuals **3.11× apart** and no convention separates their
separations by more than 1.11×. The path conventions are marginally *worse* than the point
one. **The path separation is not the missing variable.**

**The reason is structural and was visible before the run.** Over ρ = 0.5 → 4 the fast
eigenvalue goes −1.13 → −3.83 (3.4×) and the slow one 0.060 → 0.221 (3.7×) — **both scale
together**, which is exactly why the point ratio matched. §44's P1a proved Σ is *exactly*
invariant under a uniform rate rescale, so under a perfectly uniform one the residual could
not move at all. No ratio, at a point or along a path, can manufacture a 3× from an 8%
departure from uniformity.

### 46.1 The correction — and the sharper result underneath it

**Reading §39.1 and §39.2 before writing §45 up showed the question was mis-posed.**

* **T-COST-c is closed.** §39.1 established that the entire cost residual *is* a time
  residual (correlation +0.9513) and **withdrew** the off-manifold path and the Jensen gap
  as explanations of the cost. §45's closing text hands the residual back to "remaining
  candidates" that do not exist.
* **§39.2 already recorded coefficient non-transfer, under rule 9 and in as many words:**
  the T-axis value 0.6465 predicts 16.2% at γ = 0.07 against 16% measured, but 5.4% at
  γ = 0.30 against 8% measured, 33% off. **So "the two knobs give different slopes" is that
  published finding on a third axis — a confirmation read as a refutation.**
* **The quantity was wrong.** §39.1/§39.2 test `T_det/MFPT`. §45 tested `Σ_pred/Σ_exact`.

So the corrected question is not whether the axes share a coefficient — §39.2 says they do
not — but **whether the 1/sep scaling itself holds along ρ, with its own coefficient.**

**P7 gate.** §39.2 publishes `T_det/MFPT − 1 = +0.0914` at γ = 0.20, sep = 7.00. This run
gives **+0.0805 (Ω=200) and +0.0938 (Ω=300)** — bracketing it. Reproduced.

**P8/P10 — and this is the result.** `(T_det/MFPT − 1)·sep`, the quantity §39.2 found
constant to 12.3% on the T axis, with §39.2's own resolution cut (gap < 0.008 excluded and
reported):

| ρ | 0.5 | 0.75 | 1.0 | 1.5 | 2.0 | 3.0 | 4.0 | 6.0 | 8.0 | 16 | 32 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| gap | .0808 | .0886 | .0805 | .1029 | .0771 | .0623 | .0674 | .0519 | .0423 | .0251 | .0160 |
| gap·sep | **0.946** | 0.712 | 0.563 | 0.692 | 0.555 | **0.544** | 0.712 | 0.750 | 0.781 | 0.874 | **1.082** |

| axis | coefficient | spread | verdict |
|---|---|---|---|
| T (§39.2) | 0.6465 | **12.3%** | scaling holds |
| γ (here, 4 cells + 1 excluded) | 0.5963 | **22.2%** | scaling holds |
| **ρ (here, 11 cells)** | 0.7465 | **72.1%** | **scaling FAILS** |

**§39.2's 1/sep law holds on the T axis, holds on the γ axis, and breaks on the ρ axis.**
And it does not break by drifting — `gap·sep` is **U-shaped in ρ**, bottoming near 0.545 at
ρ = 2–3 and rising to 0.946 and 1.082 at the ends. That is a different *functional form*,
not a different constant.

> ⚠⚠ **§46.1's headline is WITHDRAWN by §47. The instrument was not converged.** Every
> number in this subsection was computed at **Ω = 200**, while §39.2's T-axis result was
> established over Ω = 400–800 with convergence explicitly checked. §46 therefore compared
> an unconverged measurement against a converged published one — rule 13, applied to the
> model in §46 and not to the measurement. Re-run at Ω = 300/500/700, the ρ-axis spread
> falls **46.7% → 44.4% → 22.9%**, stops being anomalous against the γ axis (22.2%) and the
> T axis (12.3%), and heads toward the 7% band that §47's parameter-free prediction demands.
> **The claim "the 1/sep scaling fails on ρ" does not survive.** It is not replaced by "the
> scaling holds": several cells are still moving at Ω = 700 (ρ = 6 shifts 26% from Ω = 500),
> so the ρ axis is *unresolved*, not settled. The γ-axis 22.2% quoted below is under the
> same suspicion, having been measured at the same Ω = 200.
>
> **What survives §47 is §46's other half**: the path separation is not the missing
> variable. That rests on the matched pair differing 3× consistently at *both* Ω, and on the
> algebraic fact that ρ scales both eigenvalues together — neither of which is a convergence
> question.

> **This qualifies §39.2's central claim.** §39.2 established the 1/sep *scaling* on two
> axes and was careful to say the *coefficient* is not universal. **The scaling is not
> universal either.** §39.2's headline — that the closed form is exact in the slaved limit
> with an O(1/sep) correction — survives as a statement about the limit; the *rate* of
> approach is axis-dependent, and on ρ it is not 1/sep at all.

> **Suspect, with its kill test (rule 17).** The three axes move the two eigenvalues in
> *opposite proportions*: ρ takes fast ×26 and slow ×4, while γ takes fast ×1.4 and slow
> ÷3.3. A law in their ratio alone cannot capture both. The slow-manifold lag should
> physically go as *manifold velocity over pool relaxation rate*,
> `|db*/dδ · μ| / |λ_fast|` — dimensionally a different object from `|λ_slow|/|λ_fast|`,
> and equal to it only when the manifold's motion is set by the slow eigenvalue.
> **How to kill:** compute it on all three axes and check constancy. A partial check
> discourages the simplest version — `gap·|λ_fast|` alone drifts **5×** across the ρ axis
> (0.091 → 0.472) — so the `db*/dδ · μ` factor must carry the whole difference or the idea
> is wrong.

---

### 47 The lag, predicted absolutely — and the prediction caught the instrument, not the physics

Three sections have now **fitted** the `T_det/MFPT` gap: §39.2 as C/sep with C = 0.6465 on
the T axis, §46 as 0.5963 on the γ axis and a claimed failure on ρ. Rule 16 says a model
that is only ever fitted is never tested.

**Singular perturbation gives it in closed form and the fitted constant disappears.** In
(δ, s) with s = x+y, writing μ = dδ/dt and ν = ds/dt, the manifold is ν(δ, s\*(δ)) = 0. As δ
advances the manifold moves and the pool lags; in quasi-steady state
(ds\*/dδ)·μ = (∂ν/∂s)·Δ, so Δ = (ds\*/dδ)μ/(∂ν/∂s). The lag shifts the signal drift by
(∂μ/∂s)·Δ, and **μ cancels**:

> **ε(δ) = (∂μ/∂s)·(ds\*/dδ) / (∂ν/∂s)**, and **T_det/MFPT − 1 ≈ ⟨ε⟩_time**

a dimensionless property of the vector field alone, with **no constant and no fit**. Note
this is *not* |λ_slow|/|λ_fast|: the reduced slow eigenvalue is
λ_red = ∂μ/∂δ + (∂μ/∂s)(ds\*/dδ), so ε is λ_red − ∂μ/∂δ over the fast rate, and coincides
with the eigenvalue ratio only when the manifold's motion happens to be set by the slow
eigenvalue.

**P1 gate.** ε is built from three finite differences, so its step h is a second axis
(rule 13). ⟨ε⟩_time = 0.090945 at every h from 4×10⁻⁴ to 2.5×10⁻⁵, last relative change
**5.9×10⁻¹⁰**. Converged to machine precision.

**P2 — the absolute test.** `⟨ε⟩_time / (T_det/MFPT − 1)`:

| axis | cells | mean | range |
|---|---|---|---|
| ρ (γ=0.20, ρ = 0.5…32) | 11 | **0.921** | 0.777 – 1.219 |
| γ (ρ=1, γ = 0.07…0.35) | 5 | **0.931** | 0.780 – 1.025 |

**A parameter-free prediction of a quantity three sections fitted, right to 8% on average.**

**P5 — it explains the fitted constant.** On the γ axis ⟨ε⟩·sep spans 0.557–0.682, mean
**0.627**, against §46's fitted **0.5963** and §39.2's T-axis **0.6465**. The constant that
was measured on two axes and would not transfer is *computed* here, and it lands between
the two fitted values.

### 47.1 The prediction disagreed with §46, and §46 was wrong

P4 failed in a way that pointed at the instrument rather than the model:

| | across ρ = 0.5…32 |
|---|---|
| predicted ε·sep | 0.664 0.646 0.637 0.631 0.633 0.641 0.648 0.657 0.662 0.672 0.677 — **7% spread** |
| measured gap·sep | 0.845 0.788 0.657 0.737 0.745 0.571 0.685 0.539 0.720 0.773 0.871 — **62% spread** |

**The model says the 1/sep law should hold on ρ with C ≈ 0.65. §46 said it fails there.**
Rule 13 says look at the instrument first — and §46 computed its spread at **Ω = 200**,
never checking Ω-convergence of the *gap*, while §39.2's T-axis result was established over
Ω = 400–800 with convergence checked. Re-running at Ω = 300/500/700:

| | Ω=300 | Ω=500 | Ω=700 | predicted |
|---|---|---|---|---|
| gap·sep spread | 46.7% | 44.4% | **22.9%** | **7.0%** |
| mean | 0.7126 | 0.6839 | 0.6857 | 0.6541 |

**§46's "the 1/sep scaling fails on ρ" is withdrawn.** The spread halves, the mean settles
near 0.686 against the predicted 0.654, and the ρ axis stops being anomalous against the γ
axis (22.2%) and the T axis (12.3%).

> **The withdrawal is not "the scaling holds", and rule 14 forbids treating it as one.**
> Per cell, several are still moving at Ω = 700 — ρ = 6 shifts **26%** from Ω = 500 (0.0392
> → 0.0496), and ρ = 32's gap has fallen to 0.0095, approaching §39.2's own 0.008 resolution
> floor. So the ρ axis is **unresolved**, not settled. What is established is that §46's
> claimed *failure* rested on an unconverged instrument. **§46's γ-axis 22.2% is under the
> same suspicion**, having been measured at the same Ω = 200.

**The methodological point, and it cost a commit.** §46 applied rule 13 to the model —
checking three averaging conventions — and never applied it to the measurement. §47 applied
it to the model again (h-convergence to 10⁻¹⁰) and only found the problem because a
*parameter-free* prediction disagreed with the data. **A fitted model would have absorbed
the scatter into its constant and reported agreement.** That is the whole content of rule
16: the fit cannot tell you your instrument is broken, because it has a free parameter to
hide the breakage in.

> ⚠ **§48 resolves this: the shortfall is finite-Ω absorption bias and it extrapolates
> away.** Every one of five cells, with endpoints matched to the lattice, converges
> monotonically upward in Ω and extrapolates to pred/gap = **1.00 ± 0.02**. §47's refusal to
> attribute the shortfall to higher-order QSS was correct, and the "wrong direction" puzzle
> below is explained: the shortfall is not proportional to the gap (absolute 0.002–0.010
> while the gap spans 0.012–0.143), so its *fraction* grows as the gap shrinks. **The lag
> model is exact as Ω → ∞.**

> **Suspect, and it does NOT fit (rule 17).** The 8% systematic shortfall was predicted in
> advance as the expected size of the next-order quasi-steady-state correction, which is
> O(gap) in the same small parameter. **But it does not behave like one.** A higher-order
> correction must shrink where the gap is small; measured, pred/gap is **0.936 at γ = 0.07
> where the gap is 0.150, and 0.780 at γ = 0.35 where the gap is 0.051** — larger shortfall
> at smaller gap, the wrong direction. So the shortfall is *not* explained by higher-order
> QSS, and attributing it there would be exactly the move this project keeps having to undo.
> **How to kill:** the γ = 0.35 cell is also the one nearest the resolution floor, so
> re-measure the shortfall against the gap at Ω = 700+ before treating its sign as real.

---

### 48 The shortfall extrapolates away — the lag model is exact as Ω → ∞

§47 predicted the `T_det/MFPT` gap with no fitted constant and landed 8% short, then
declined to attribute the shortfall to next-order quasi-steady state because it ran the
wrong way: pred/gap was 0.936 where the gap was 0.150 and 0.780 where it was 0.051. Two
instrument candidates, both measured before predictions were written: **(a)** `T_det`
integrates the *nominal* [εδ\*, θδ\*] while the CME first-passage runs between the *lattice*
endpoints d₀/Ω and thr/Ω — rule 11, literally — and **(b)** the absolute shortfall clustered
near +0.003 across cells whose gaps differed 7×.

**P2 is refuted. The endpoints are not the cause.** Matching them moves `T_det` by
**−0.09% to +1.02%** against a shortfall of 10–25%, and the per-cell prediction failed
outright: γ = 0.35, which had the largest overshoot and the worst ratio, got **worse** by
0.043, while γ = 0.07 improved by only 0.006.

**But matching them turned noise into a signal, and that is what mattered.** Nominal
endpoints give erratic series; matched endpoints give monotone ones — in **all five cells**:

| cell | Ω=300 | Ω=500 | Ω=700 | Ω=1000 |
|---|---|---|---|---|
| γ=0.07 nominal | 0.936 | 1.052 | 1.000 | 0.991 |
| γ=0.07 **matched** | 0.9415 | 0.9655 | 0.9741 | **0.9745** |
| γ=0.20 nominal | 0.969 | 0.937 | 0.983 | 0.941 |
| γ=0.20 **matched** | 0.8836 | 0.9172 | 0.9421 | **0.9558** |
| γ=0.35 **matched** | 0.7369 | 0.7397 | 0.7541 | **0.7959** |
| ρ=0.5 **matched** | 0.7226 | 0.7950 | 0.8214 | **0.8594** |
| ρ=32 **matched** | 0.6065 | 0.7159 | 0.7776 | **0.8281** |

The lattice rounding varies erratically with Ω, so the nominal ratios bounce and no
convergence can be read off them at all. **The endpoint fix bought monotonicity, not
accuracy** — and monotonicity is what makes P3 and P4 testable.

**P4 holds, and its sign was forced in advance.** Absorption at a threshold selects the
leading edge of the packet, so the MFPT sits *below* the mean arrival, the measured gap
*exceeds* the pure lag, and pred/gap must approach 1 **from below**. Every cell is below 1
and every cell rises with Ω: means 0.778 → 0.827 → 0.854 → 0.883.

> ⚠ **§52 qualifies this.** These fits use four points over Ω = 300–1000 with a free
> exponent, on data now known to be non-asymptotic in six of eight cells — cme_coeff is still
> rising at Ω = 1400 everywhere except the two fast-pool cells. **The Ω → ∞ conclusion is
> properly established only at ρ = 32**, where the coefficient is flat and (gap − pred) ∝ 1/Ω
> exactly as claimed. Elsewhere it is consistent with the data, not demonstrated by it.

**P3/P6 — it extrapolates to 1.** Per cell, fitting `1 − pred/gap = c·Ω^(−p)`:

| γ | ρ | Ω=300 | Ω=1000 | exponent | intercept (1/√Ω) | intercept (free p) |
|---|---|---|---|---|---|---|
| 0.07 | 1 | 0.0585 | 0.0255 | 0.72 | 1.0212 | **1.0038** |
| 0.20 | 0.5 | 0.2774 | 0.1406 | 0.55 | 1.0206 | **0.9995** |
| 0.20 | 1 | 0.1164 | 0.0442 | 0.82 | 1.0460 | **0.9982** |
| 0.20 | 32 | 0.3935 | 0.1719 | 0.69 | 1.0980 | **0.9968** |
| 0.35 | 1 | 0.2631 | 0.2041 | 0.20 | 0.8427 | **0.9838** |

**All five extrapolate to pred/gap = 1.00, spanning 0.984–1.004, mean 0.996.** The lag
model — `T_det/MFPT − 1 = ⟨ε⟩_time`, no fitted parameter — is **exact in the Ω → ∞ limit**,
and §47's 8% shortfall is finite-Ω absorption bias: §39.1's candidate (iv), which §39.2 left
live for the time and which nothing had measured until now.

> **Rule 15: the two ansätze disagree, so the decay LAW is unresolved.** Free-p gives
> 0.984–1.004; a fixed 1/√Ω gives **0.843–1.098**, and at γ = 0.35 it extrapolates to 0.843,
> not 1. The exponents span **0.20–0.82** (median 0.69) and are fitted from four points with
> two parameters, so they are weakly determined. **The intercept is robust; the exponent is
> not.** Convergence to exactly 1 is established only under a free exponent.

> ⚠ **Rule 18, caught in the act.** Averaged across cells, the shortfall looks like a clean
> law: `(1 − mean ratio)·√Ω` = 3.842, 3.876, 3.866, 3.708 — **constant to 4% over a 3.3×
> range in Ω**, which would have been reported as a discovered 1/√Ω absorption law. Per
> cell, the same coefficient spans **0.7 to 6.5**. The mean is an artifact of averaging cells
> with different coefficients, and the clean law is not claimed. This is precisely §23.5's
> failure — a fit across an axis averaging a localised effect into something tidy.

**§47's puzzle is resolved.** The shortfall is *not* proportional to the gap — in absolute
terms it spans 0.002–0.010 at Ω = 1000 while the gap spans 0.012–0.143 — so its ratio to the
gap necessarily grows as the gap shrinks. That is the "wrong direction" §47 observed, and
§47 was right to refuse to attribute it to higher-order QSS.

---

### 49 The absorption exponent is 1, derived — the coefficient is not the bulk term

§48 closed T-COST-i but left the decay law unresolved, fitting per-cell exponents of
0.20–0.82 and declining to claim the tidy 1/√Ω that appeared after averaging. Fitting an
exponent is what rule 16 says never settles anything, so this derives it.

On the slaved manifold δ jumps by ±1 with total fluxes `up` and `dn`, both already computed
by `updown`, so in concentration units μ = up − dn and **D₀ = (up+dn)/2**. Writing
g = μ/D₀ so that ψ′ = Ω·g, the Laplace expansion of the first-passage integral gives

> **T_det/MFPT − 1 = ⟨ε⟩_time + K/Ω**,  **K = −(1/T_det)∫ g′/(μg²) dδ**

**The exponent is 1, not 1/2, and K is fully computable with no fit.**

> ⚠ **§52 extends this to Ω = 1400 and the picture splits.** ρ=32 holds superbly (0.1% last
> step, 6.7% total drift over 4.7× in Ω); **γ=0.20 does not** — it rises to 4.899, a 14.1%
> last step. And γ=0.07's Ω=1000 point, excused below as "a numerical floor, not physics",
> **is not a floor** — the trend keeps rising to 4.870 at Ω = 1400. That diagnosis is
> withdrawn. The exponent-1 confirmation stands in the fast-pool cells and is premature
> elsewhere.

**The exponent is confirmed where it is testable.** `(gap − pred)·Ω` over Ω = 300…1000:

| cell | T_det | 300 | 500 | 700 | 1000 | |
|---|---|---|---|---|---|---|
| γ=0.20, ρ=32 | 3.06 | 1.95 | 1.98 | 2.00 | 2.08 | **constant to 6%** |
| γ=0.20, ρ=1 | 6.65 | 3.58 | 4.11 | 3.90 | 4.21 | **constant to 17%** |
| γ=0.07 | 4.89 | 2.61 | 2.49 | 2.60 | 3.65 | constant to 4% for Ω ≤ 700 |
| γ=0.20, ρ=0.5 | 13.57 | 6.55 | 7.31 | 8.62 | 9.28 | rising 42% |
| γ=0.35 | 12.91 | 4.33 | 7.00 | 9.13 | 10.26 | rising 137% |

**§48's fitted 0.20–0.82 was averaging asymptotic cells together with pre-asymptotic ones.**
The exponent is 1 where the expansion applies; the cells where it does not are still
approaching it.

**But P2 fails: K is not the coefficient.** In the two cleanly asymptotic cells, K/measured
is **0.333 and 0.293** — mean 0.313 with a spread of only **12.9%**, so the bulk term has
the right shape and scale but is short by a factor of about 3. That crosses the "off by more
than 3×" bar fixed in advance, so **the bulk diffusion term is refuted as a complete
account** even though it is clearly a real component.

**And at γ = 0.07 the sign is wrong — K = −0.263 against a measured +3.65.** No
multiplicative factor repairs a sign.

> ⚠ **§50 withdraws the 31%.** An *exact* tridiagonal solve of the same 1-D process — bulk,
> boundary layer and discreteness together — gives 0.217 at γ = 0.20, ρ = 32, against the
> 0.293 the bulk term alone gives here. A complete account cannot be smaller than one of its
> parts, so K is wrong, for the reason recorded just below: it is a near-cancellation. **The
> exponent, which does not depend on K, stands.**

> **The reason K is fragile, found by inspecting it rather than by fitting around it.**
> g rises then falls across the traversal — sign pattern `+ + + + + − − − −` at γ = 0.07,
> `+ + + + − − − − −` at γ = 0.20 and 0.35 — so **K is a near-cancellation between two
> comparable halves.** A quantity that is the small difference of larger pieces is the
> *least* robust member of an O(1/Ω) budget, not the most, and at γ = 0.07 the rising half
> simply wins. This is why the 31% deficit should not be read as "the bulk term times a
> constant".

**Two same-order terms are missing by construction, one of which was not named in advance.**
The derivation extends the inner integral to −∞, so it is the bulk term only and omits the
absorbing boundary layer of width D/μ — that was P4. The one I did not account for: **the
CME is a jump process, and its diffusion approximation is itself wrong at O(1/Ω)**, so the
Kramers–Moyal truncation contributes at exactly the same order as everything else here.

> **The tempting constant is declined (rule 9).** 1/0.313 = 3.19, and π = 3.14159. The two
> cells individually give 3.00 and 3.41. **Two cells cannot distinguish π from 3 from 3.2**,
> and naming one would be precisely the move this project keeps having to undo.

> ⚠ **The T_det ordering proposed in the design scan does NOT survive.** It predicted the
> longest-traversal cells would be the pre-asymptotic ones, which holds for ρ=0.5 (13.57) and
> γ=0.35 (12.91) — but γ=0.07 has T_det = 4.89, *shorter* than the asymptotic γ=0.20 (6.65),
> and is still flagged. Its flag rests entirely on the Ω=1000 point, where the residual has
> fallen to 0.0037 and stopped decaying between Ω = 700 and 1000 — the signature of a
> numerical floor, not of physics. **Reported as a failed ordering, not quietly dropped.**

**What §49 settles and what it leaves.** Settled: the absorption correction is **O(1/Ω)**,
derived rather than fitted, and confirmed to 6% and 17% in the two cells where the expansion
applies. Open: its coefficient. The bulk diffusion term supplies ~31% of it with a
consistent offset, and the remainder must come from the boundary layer and the
Kramers–Moyal truncation, neither of which is computed here.

**T-COST-k, open: what supplies the other ~69%?** Both missing terms are O(1/Ω) and both are
computable. **How to kill:** the Kramers–Moyal term is the cheaper of the two — it is the
difference between the exact jump-process MFPT and its diffusion approximation on the *same*
1-D slaved chain, which is one tridiagonal solve against one quadrature and needs no CME at
all. If that difference alone closes the gap to the measured coefficient, the boundary layer
is negligible and the budget is complete.

---

### 50 A complete 1-D account supplies only ~21% — the absorption correction is mostly 2-D

§49 left T-COST-k: the bulk diffusion term K supplies ~31% of the absorption coefficient
and has the wrong sign at γ = 0.07, with two same-order terms missing — the absorbing
boundary layer and the Kramers–Moyal truncation. **Computing them separately is
unnecessary: the 1-D slaved birth–death chain contains all three at once.** δ hops ±1 with
rates Ω·up and Ω·dn, both already returned by `updown`, so its exact MFPT is a tridiagonal
solve — no expansion, no CME, no fitted parameter.

**P1 fails as pre-registered, and the first pass failed for a reason §48 had already
taught me.** `T_det` integrated the unrounded limits while the chain ran between the
*rounded* lattice points m₀/Ω and thr/Ω. That mismatch is O(1/Ω) in δ — which, multiplied
by Ω, is **exactly the size of the effect being measured**. It produced sign-flipping
garbage:

| γ, ρ | Ω=300 | 500 | 700 | 1000 | 2000 | spread |
|---|---|---|---|---|---|---|
| 0.35, 1 — **before** | −3.587 | 0.920 | −0.583 | 0.363 | 4.167 | **3029%** |
| 0.35, 1 — after | −3.587 | −0.403 | 0.738 | 1.686 | 2.842 | 2518% |
| 0.07, 1 — after | 0.483 | 0.500 | 0.538 | **0.545** | **0.546** | 12.2% |
| 0.20, 32 — after | 0.354 | 0.381 | 0.436 | **0.447** | **0.450** | 23.0% |

Matching the endpoints removed every sign oscillation and made all five series monotone.
**Rule 11, twice in three sections, in the same experiment family.**

**The gate still fails**: three of five cells have not converged by Ω = 2000. Restricting
to the two whose last two Ω agree within 1%:

| cell | bd_coeff | cme_coeff | **bd/cme** |
|---|---|---|---|
| γ = 0.07, ρ = 1 | 0.546 | 2.568 | **0.213** |
| γ = 0.20, ρ = 32 | 0.450 | 2.076 | **0.217** |

**The two agree to 2%.** A *complete* one-dimensional account — bulk, boundary layer and
jump discreteness together, computed exactly — supplies **~21%** of the measured
correction.

> **That selection is post-hoc and it weakens the claim.** The convergence criterion was
> chosen after seeing which cells converged, unlike §47's δ\* ≥ 0.40 which was fixed in
> advance. The three unconverged cells give 0.200, 0.269 and 0.277, and their bd_coeff is
> still rising, so those are lower bounds. The result is provisional pending a run at larger
> Ω.

**P4 holds and is diagnostic.** At γ = 0.07 the bulk term K = −0.263 against a measured
+3.65 — a sign error no factor could repair. The complete chain gives **+0.546**. So §49's
sign failure was the near-cancellation in K, exactly as diagnosed, and not something
structural in the slaved reduction.

**P2, P3, P5 all fail.** The 1-D chain is not better than K (0.243 mean vs 0.313), does not
supply the whole correction, and does not exceed the still-rising measurements in the
pre-asymptotic cells.

> **~79% is genuinely two-dimensional, and the suspect was named nine sections ago
> (rule 17).** The pool `b` fluctuates about its slaved value, so the effective drift felt
> by δ is E[μ(δ, b)] rather than μ(δ, b\*) — a Jensen term, O(1/Ω), that no one-dimensional
> reduction can carry. **That is §39.1's candidate (iii)**, which §39.1 withdrew as an
> explanation of the *cost* while explicitly recording that it "may still explain the *time*
> gap, which is now a separate and cleaner question." It is now the leading explanation of
> precisely that gap. **How to kill:** compute E[μ(δ,b)] − μ(δ,b\*) with the pool variance
> from the linear-noise approximation about the manifold — one 2×2 Lyapunov solve per δ —
> and test the resulting coefficient absolutely against the missing 79%. It must also
> reproduce the *sign*, which is fixed: absorption makes the MFPT shorter, so the term must
> be positive.

> ⚠ **§49's 31% and §50's 21% are inconsistent, and §50's supersedes it.** At the *same*
> cell (γ = 0.20, ρ = 32) the bulk term alone gives 0.293 while the *complete* 1-D account
> gives 0.217 — **a whole that is smaller than one of its parts.** Since the boundary-layer
> and discreteness terms have no reason to be negative, one of the two is wrong, and it is
> K: K is an asymptotic expansion resting on a near-cancelling integral (§49 showed g′
> changes sign mid-path and that the cancellation flips K negative at γ = 0.07), while
> bd_coeff is an *exact* tridiagonal solve of the actual one-dimensional process. **§49's
> "the bulk term supplies 31%" is withdrawn; the exact 1-D chain supplies ~21% and
> supersedes it.** §49's surviving contribution is the *exponent*, which is independent of K.

**Where the cost budget stands.** `T_det/MFPT − 1 = ⟨ε⟩_time + C/Ω`, with ⟨ε⟩ the
deterministic lag (§47, exact as Ω → ∞ by §48) and the exponent on C derived and confirmed
(§49). Of C, an exact one-dimensional account supplies **~21%**; **the remaining ~79% is
2-D pool noise, unmeasured, and is the only piece of the cost arc still without a number.**

---

### 51 The Jensen term is exactly zero — killed by §30's identity

§50 left ~79% of the absorption coefficient unaccounted and named the suspect: the pool
fluctuates about the manifold, so δ feels E[μ(δ,s)] rather than μ(δ,s\*), a Jensen term that
no 1-D reduction can carry. **This is §39.1's candidate (iii)**, withdrawn there as an
explanation of the *cost* while explicitly left live for the *time*. It is now dead for both.

Expanding at fixed δ, the first-order term E[s] − s\* is the deterministic lag — already
§47's ε. The new piece is second order: with the pool fast, s holds a quasi-stationary
Gaussian of variance `Var(s) = D_s/|∂ν/∂s|`, giving

    ε_J(δ) = (1/2μ)·(∂²μ/∂s²)·Var(s),   J = Ω·⟨ε_J⟩_time

**Measured: J = 0 to machine precision in all five cells.** Not small — zero. And the
reason is exact:

> **μ = k·δ·(1 − (1+γ)s)**, verified against the network's own fluxes to **9.8×10⁻¹⁶** over
> 81 states × 3 γ × 3 ρ.

**That bracket is §30's identity at n = 2**, `k(b − γs)`, written in concentrations. The
drift is **exactly linear in the pool coordinate**, so ∂²μ/∂s² ≡ 0 and the Jensen term
vanishes identically. `(bd+J)/cme` is unchanged at **0.213 and 0.217** — the same 21% §50
reported, with the Jensen term contributing nothing at all.

**Note also that ρ does not appear in μ.** The disagreement channel moves δ by exactly
zero — §30's first cancellation — which is why the whole ρ family shares one drift law. The
theorem that opened this session is what closes this question.

> **P1's "failure" is a gate artifact and is reported as one.** The h-convergence gate
> measured the *relative* change of a quantity that is exactly zero, so it divided by zero
> and reported 8×10⁻¹ to 1.04. The absolute values were 0.00000 at every h. The gate was
> written for a nonzero quantity and does not apply; it is not evidence against anything.

> **P2's sign test was the cheap way to be wrong, and it was not the way this failed.** The
> prediction was that ∂²μ/∂s² > 0 was *required* for the Jensen term to have the right sign.
> It is neither positive nor negative — it is zero — so the sign test never got to fire.

**What survives as a suspect (rule 17).** The *noise* is curved even though the drift is not:
`up + dn = k·b·s + γk(s²+δ²)/2`, so **∂²(up+dn)/∂s² = k(γ−2) ≠ 0**. But a Jensen term in the
*diffusion* enters the MFPT only through the correction itself, making it O(1/Ω²) — second
order, and so almost certainly not the missing 79% either.

**The surviving candidate is noise-induced drift from the δ–s cross-correlation.** Adiabatic
elimination of a fast variable produces an effective slow drift with a term set by
⟨δ-fluctuation · s-fluctuation⟩, distinct from both the deterministic lag and the Jensen
curvature. It survives here precisely because **∂μ/∂s = −k(1+γ)δ ≠ 0 even though
∂²μ/∂s² = 0** — the drift is linear in s, not independent of it. **How to kill:** compute the
stationary 2×2 covariance from the Lyapunov equation about the manifold, take the δs entry,
and test `(∂μ/∂s)·⟨δs⟩/μ` absolutely against the missing 79%, with the sign again forced
positive.

**The budget after §51.** `T_det/MFPT − 1 = ⟨ε⟩_time + C/Ω`, lag exact as Ω → ∞ (§47/§48),
exponent derived (§49), and of C: **~21% exact 1-D discreteness (§50), 0% Jensen (§51), ~79%
still unaccounted.** Two named candidates are now eliminated rather than merely unfavoured.

---

### 52 The 21% survives on new cells — but §49's convergence claims do not

§50's "a complete 1-D account supplies ~21%" rested on **two cells selected post-hoc**, one
of them (γ = 0.07) with a top point §49 had excused as a numerical floor. §51 then
eliminated two candidate mechanisms for the remaining 79% without touching that measurement.
Rule 8 says verify before building further, so this adds Ω = 1400 and three new cells rather
than a third mechanism.

**P2 holds, and on a cell §50 never saw.** Applying §50's convergence criterion unchanged
(last two Ω within 1%):

| γ | ρ | bd | cme | **bd/cme** | cme converged? |
|---|---|---|---|---|---|
| 0.20 | 32 | 0.446 | 2.079 | **0.214** | **YES** (0.1%) |
| **0.28** | **8** | 0.658 | 3.264 | **0.202** | **YES** (1.0%) |
| 0.20 | 1 | 1.173 | 4.899 | 0.239 | no |
| 0.35 | 1 | 2.450 | 10.930 | 0.224 | no |
| 0.12 | 1 | 0.782 | 3.757 | 0.208 | no |
| 0.20 | 4 | 0.658 | 3.191 | 0.206 | no |
| 0.20 | 0.5 | 1.766 | 9.093 | 0.194 | no |
| 0.07 | 1 | 0.557 | 4.870 | 0.114 | no |

**γ = 0.28, ρ = 8 is a new cell on a new γ and a new ρ, it converged, and it gives 0.202**
against ρ=32's 0.214 — mean 0.208, spread 6.1%. And **seven of eight cells give 0.194–0.239
regardless of convergence status.** §50's headline is strengthened, not weakened: the
one-dimensional share really is about a fifth.

### 52.1 But only 2 of 8 cells are Ω-converged, and §49 was wrong about why

**P3 resolves against §49.** γ = 0.07's series is **2.608, 2.494, 2.601, 3.653, 4.870** — the
Ω = 1000 point that §49 attributed to "a numerical floor, not physics" is not a floor at all.
**The trend genuinely rises**, and keeps rising at Ω = 1400. That diagnosis is withdrawn.

> ⚠ **§49's exponent confirmation holds only in the fast-pool cells.** §49 reported
> `(gap − pred)·Ω` "constant to 6% (ρ=32) and 17% (γ=0.20)" over Ω = 300–1000. Extended to
> 1400: **ρ=32 holds superbly** — 1.948, 1.980, 1.996, 2.076, 2.079, a 0.1% last step and
> only 6.7% total drift over a 4.7× range in Ω. **γ=0.20 does not** — it rises to 4.899, a
> 14.1% last step. The claim was right for one cell and premature for the other.

**The pattern is coherent and it is about sep, not Ω alone.** The two cells that converge are
the two with the fastest pools — ρ = 32 (sep 67.8) and γ = 0.28 with ρ = 8. Every slow-pool
cell is still climbing at Ω = 1400. **The absorption coefficient is asymptotic in both Ω and
sep, and at moderate sep, Ω ~ 10³ is not yet in the asymptotic regime.**

> ⚠ **This qualifies §48 as well.** §48's per-cell extrapolations to pred/gap = 1.00
> (0.984–1.004) were fitted from four points over Ω = 300–1000 with a free exponent — on data
> now known to be non-asymptotic in six of eight cells. **The "lag model is exact as Ω → ∞"
> conclusion is properly established only where the data is asymptotic**, which is ρ = 32,
> where cme_coeff is flat and therefore (gap − pred) ∝ 1/Ω → 0 exactly as claimed. Elsewhere
> it is consistent with the data but not demonstrated by it.

**What this leaves.** The 1-D share of the absorption coefficient is **~0.21, now measured on
two independently converged cells and clustered across seven of eight** — so the "~79% is
two-dimensional" statement is on much firmer ground than it was in §50. What is *not*
established is the Ω-asymptotics of the coefficient itself outside the fast-pool limit.
**Running this before proposing T-COST-m was the right order**: a third mechanism fitted
against six non-asymptotic cells would have been fitted against a moving target.

---

### 53 The amplification half is NOT structural — and P has a closed form — T15-e

§43 proved `b_i − b_j = (n_i − n_j)·P(n)` for every exchange-symmetric mass-action network,
making δ = 0 an invariant manifold, and was explicit that this is only half of restoration:
**divisibility gives no-reversal; amplification additionally needs P > 0, and only the first
half is a theorem.** This closes the second half — negatively, with a closed form as
consolation.

**P is the symmetry-breaking eigenvalue, and §43 did not name it.** Differentiating
b_i − b_j = δ·P at δ = 0 gives P(symmetric state) = d(b_i − b_j)/dδ, which at a symmetric
*fixed point* is the Jacobian eigenvalue along (1,−1,0) — the quantity THEORIES **T7** and
**§14** built the n-winner barrier on. Verified to **5.8×10⁻¹⁰** over 52 states on AM and on
random symmetrised networks, and exactly at AM's symmetric fixed point.

**And it has a closed form.** Measured at 11 values of γ from 0 to 0.9:

> **P(symmetric fixed point) = (1 − 2γ)/3**, worst deviation **8.9×10⁻¹²**

| γ | 0 | 0.20 | 0.35 | 0.49 | 0.55 | 0.90 |
|---|---|---|---|---|---|---|
| P | 0.333333 | 0.200000 | 0.100000 | 0.006667 | −0.033333 | −0.266667 |

It vanishes at **γ_c = 1/2**, which is the pitchfork, and at γ = 0 it is **1/3 = 1/(2n−1)
at n = 2 — exactly T7/§14's symmetry-breaking eigenvalue λ(n).** An absolute check (rule 16)
linking two sections that had never been connected. It also sits in fixed ratio to §12's
wall coefficient: **P/κ = (1+γ)/4.5** exactly, both carrying the same (1−2γ) zero.

**P2 — the attractors lie on P's zero set.** |P| at δ\* is **7.5×10⁻¹⁶, 1.7×10⁻¹⁷,
1.5×10⁻¹⁶, 2.4×10⁻¹⁶** across γ. Since dδ/dt = δ·P, an off-symmetric fixed point *must* have
P = 0, and `delta_star(γ)` reaches the same points by an independent route. Note this is not
the separatrix — §43's theorem already makes that δ = 0.

**P5 — no interior sign flip.** P falls monotonically from the symmetric state to the
attractor at every γ, so AM's amplifying region is exactly (0, δ\*) and "restoring" needs no
region qualifier here.

### 53.1 sign(P) is not combinatorial, and the two halves have different status

**P3, the kill test: hold the topology fixed, vary only the rate constants.** Over 188
random symmetrised topologies, **17 flip the sign of P** at a fixed state — examples
spanning **−11.46…+17.20**, **−4.82…+3.11**, and **−38.26…+2.04**. AM itself flips at
γ_c. **One counterexample suffices**, so:

> **sign(P) is not determined by the stoichiometry.** §43's theorem covers the
> cannot-reverse half structurally, and that is the whole of what is structural.
> Amplification is a linear-stability condition on the rate constants.

> ⚠ **The pre-registered verdict rule was badly designed and is reported as such.** It
> declared "not combinatorial" only if *more than half* of topologies flipped; 17/188 = 9%
> flip, so the script printed "sign(P) looks topology-determined". **That is the wrong
> logic for a universality claim** — refuting "topology determines sign" needs one
> counterexample, not a majority. The verdict above is the correct reading of the same data,
> and the criterion, not the data, is what changed.

The 91% that do *not* flip are informative in the other direction: **topology carries
substantial information about sign(P) without determining it.**

**P4 — and this is the sharpest thing here.** Over 200 random symmetrised networks, how many
amplify at all?

| property | fraction |
|---|---|
| divisibility (§43) — no-reversal | **200 / 200** |
| P > 0 somewhere — amplification | **21 / 200 = 10.5%** |

**Divisibility is generic; amplification is rare.** The two halves of restoration are not
merely different in logical status — they are different in *prevalence*. Almost every
exchange-symmetric network preserves the sign of a lead, and almost none of them grows it.

**What T15-e settles.** Restoration = divisibility + sign(P) > 0. The first is a theorem
about stoichiometry and symmetry, holds universally, and needs no conservation law (§43).
The second is a rate-constant condition, equals the symmetry-breaking eigenvalue, has the
closed form (1−2γ)/3 in AM, and is **rare** among symmetric networks. **A restoring element
is not a structural accident — it is a tuned one**, and the project's founding claim about
the transistor being a near-ideal restoring switch is a statement about tuning, not topology.

---

### 54 Where the combinatorics runs out — P is a rate-weighted sum of integers — T15-f

§53 closed T15-e negatively (sign(P) is not determined by stoichiometry) and left the
obvious question unasked: **what do the amplifying 10.5% have in common?** They have a
decomposition.

Group an exchange-symmetric network's reactions into mirror pairs {r, r̄}, and let r be the
X-heavy member (X-power p > Y-power q in its reactants). Symmetry forces
S_X(r̄) = S_Y(r), so with **d_r = S_X(r) − S_Y(r)** the pair contributes
d_r c_r O_r (x^p y^q − x^q y^p), and factoring x^p y^q − x^q y^p = (xy)^q (x−y) Σ_m x^m y^…
gives

> **P = Σ_pairs d_r · c_r · [ O_r (xy)^q Σ_m x^m y^(p−q−1−m) ]**, every bracket **≥ 0** on x,y ≥ 0

A self-mirror reaction (p = q) forces S_X = S_Y and contributes **nothing**. **Verified to
1.9×10⁻¹⁴** on AM, on `am_cubic`, and on 120 random symmetrised networks.

**sign(P) is therefore a rate-weighted sum of integers with non-negative weights**, and
three regimes follow — two of them purely combinatorial:

| class | prediction | networks | amplify | violations |
|---|---|---|---|---|
| all d_r ≤ 0 | **P ≤ 0 everywhere**, whatever the rates | 113 | **0** | **0** |
| all d_r ≥ 0, some > 0 | **P > 0 everywhere** | 3 | 1 | **0** |
| mixed | rates decide, and only here | 183 | 33 | — |

> ⚠⚠ **§56 corrects this table: `classify` was counting p = q pairs.** Those have an empty
> bracket and contribute identically zero, but `mirror_pairs` records their d as an arbitrary
> ±1 (§55 found this in the n-winner counter and the fix was never propagated here).
> **148 of 232 networks were misclassified "mixed" when they are genuinely all-d ≤ 0.**
> Corrected, the table reads **all≤0 211 (0 amplify), all≥0 7 (7 amplify), mixed 63 (27
> amplify), trivial 19** — and the "two of three all≥0 networks read P ≈ 0" caveat below
> **was not a threshold artifact, it was this bug**: corrected, all≥0 amplifies **7/7**.
> The flip cross-check strengthens too: 24/37 mixed flip, 0 of 263 unanimous.

**P2 and P3 hold with zero violations.** Topology *can* rule restoration out, and *can*
rule it in — it is silent only in the mixed case. (Two of the three all-≥0 networks read
P ≈ 0 rather than P > 0 at the sampled states; their positive bracket falls below the 10⁻⁹
threshold at small concentrations, which is a threshold artifact, not a sign violation.)

**P4 — the cross-check, and it is the sharp one.** §53 measured, for a different purpose,
that 17/188 topologies flip sign(P) under rate changes alone. The classification here
predicts that **only mixed topologies can flip**:

| class | networks | flipped |
|---|---|---|
| all ≤ 0 | 159 | **0** |
| all ≥ 0 | 8 | **0** |
| mixed | 132 | **24 (18%)** |

**Zero of 168 unanimous topologies flipped.** And the rates reconcile: mixed are 44% of the
sample and 18% of those flip, giving 8% overall against §53's independently measured 9%.

**d_r > 0 is NOT autocatalysis, which is the obvious guess and is wrong.** `B + X → 2X`
gives d = +1 (X makes more X), but so does **`2X + Y → 2X + B`, where S_X = 0 and
S_Y = −1** — X catalysing *Y's destruction*. The governing notion is positive feedback on
the **difference**, not on either species.

### 54.1 AM decomposed, and what γ_c actually is

| d_r | k | reaction |
|---|---|---|
| **+1** | 1.0 | `f2: B + X → 2X` (recruitment) |
| **−1** | γ | `r2: 2X → B + X` (its reverse) |
| — | — | `f1: X + Y → 2B` is **self-mirror → contributes 0** |

Evaluating the brackets gives **P = k(b − γs)** — **§30's identity, recovered term by
term.** So:

* **AM is a *mixed* network**, which is why γ can flip its sign at all.
* **γ is literally the weight on the contracting term**, and **γ_c = 1/2 is the point where
  the mixed sum changes sign.** §53's closed form P = (1−2γ)/3 is that sum evaluated at the
  symmetric fixed point.
* **The disagreement channel is self-mirror and contributes exactly zero** — which is the
  same fact as §30's first cancellation, and the same fact as §51's discovery that ρ does
  not appear in μ at all. Three sections, one structural cause.

**P5 — §53's 10.5% decomposed.** This draw gives 34/300 = 11.3% amplifying, of which
**11.0 points come from mixed networks and 0.3 from guaranteed ones**. The
topology-guaranteed class is ~1% of networks. **So nearly all restoration in this family is
rate-tuned rather than topology-forced**, which is §53's "a restoring element is a tuned
object" with a mechanism attached rather than an observation.

**What T15-f settles.** §53 said sign(P) is not combinatorial. §54 says *precisely where*
the combinatorics runs out: it decides the unanimous cases completely and the mixed case not
at all. Restoration is impossible for 38% of these networks on stoichiometric grounds alone,
guaranteed for 1%, and a matter of tuning for the remaining 61%.

---

### 55 The decomposition predicts γ_c(n) and λ(n) with no fit — T15-g

§54 decomposed P into a rate-weighted sum of integers, but only at n = 2. §30 proved the
pairwise identity for **every** n, so the decomposition should follow — and at n > 2 it makes
predictions against two numbers this project published long ago, with nothing fitted.

**Sorting `n_winner_reversible`'s reactions by the (i,j) swap:**

| reaction | p, q | d_r | weight | bracket | count |
|---|---|---|---|---|---|
| `B + X_i → 2X_i` | 1, 0 | **+1** | k | b | 1 |
| `2X_i → B + X_i` | 2, 0 | **−1** | γk | x_i + x_j | 1 |
| `X_i + X_k → 2B`, k ∉ {i,j} | 1, 0 | **−1** | k | x_k | **n−2** |
| `X_i + X_j → 2B` | 1, 1 | — | — | **self-mirror → 0** | 1 |
| `2B → X_i + X_k` | 0, 0 | — | — | **p = q → cancels** | n−2 |

Summing at a symmetric state (all x_l = x, b = 1 − nx) gives **P = k[b − γ·2x − (n−2)x]**,
and §30's published bracket is `(k/Ω)[n_B − Σ_{l≠i,j} n_l − γ(n_i + n_j − 1)]` — **the same
three terms, same order, same signs.**

**P1** — the identity and the decomposition both survive n > 2, worst **3.8×10⁻¹⁰** over
n = 2…6 at three γ each.

**P2 — absolute, against T7/§14.** At γ = 0 the symmetric fixed point has x = 1/(2n−1), so
b = (n−1)/(2n−1) and P = [(n−1) − (n−2)]/(2n−1) = **1/(2n−1)**:

| n | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|
| P measured | 0.3333333333 | 0.2000000000 | 0.1428571429 | 0.1111111111 | 0.0909090909 |
| 1/(2n−1) | 0.3333333333 | 0.2000000000 | 0.1428571429 | 0.1111111111 | 0.0909090909 |

Worst deviation **8.9×10⁻¹²**. **T7/§14's symmetry-breaking eigenvalue λ(n), recovered from
a stoichiometric decomposition that knows nothing about it.**

**P3 — absolute, and the sharper one.** P = 0 is the loss of amplification, so it must vanish
at γ_c(n), which `gamma_critical(n)` finds by an entirely independent bracketed root-find:

| n | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|
| γ_c | 0.50000000 | 0.20222577 | 0.06807612 | 0.02613332 | 0.01219827 |
| **P at γ_c** | **0.0** | 3.5×10⁻¹² | **0.0** | 8.7×10⁻¹³ | **0.0** |

**The decomposition predicts γ_c(n) with no fitted parameter**, and γ_c(3) = 0.202226 against
the published 0.2023.

**P4 — and this explains why γ_c(n) falls.** Every n-winner pair is **mixed**, with the
contributing terms counting as:

| n | +1 | −1 | zero-bracket |
|---|---|---|---|
| 2 | 1 | 1 | 0 |
| 3 | 1 | 2 | 1 |
| 4 | 1 | 3 | 2 |
| 5 | 1 | 4 | 3 |
| 6 | 1 | 5 | 4 |

**Exactly one amplifying term against n−1 contracting ones**, plus n−2 that vanish. γ_c(n)
falls with n because the single d = +1 must outweigh one more contracting term at every step
— which is the radix penalty of §3/§6.1 read off the stoichiometry.

> ⚠ **The first pass of P4 reported the wrong counts, and the bug was mine.** It read
> `1 / 2n−3 / 0` — three negatives at n = 3 where the derivation says two. The cause: a
> p = q pair has an **empty** bracket sum and contributes exactly zero, but `mirror_pairs`
> records its d as an arbitrary ±1, because "the X-heavy member" is undefined when p = q.
> Counting raw d therefore reported n−2 phantom contracting terms. **P1–P3 were unaffected**,
> since they use the decomposition itself, where the empty sum correctly gives zero — the
> bug was in the *reporting*, not the *algebra*, and it was caught by the counts disagreeing
> with a derivation written before the run.

**P5 — two different routes to zero, and §54 saw only one.** At n = 4 the self-mirror
reactions are `X1+X2→2B`, `X3+X4→2B`, `B+X3→2X3`, `B+X4→2X4` and their reverses; separately,
the four `2B → X_i+X_k` reactions form genuine mirror *pairs* whose contributions cancel
through p = q. **Self-mirroring and pairwise cancellation are distinct mechanisms**, and
conflating them is what produced the P4 bug.

**What §55 adds.** §54's decomposition is not an n = 2 curiosity: it reproduces §30's
n-winner bracket term by term, recovers T7/§14's λ(n) = 1/(2n−1) exactly, and predicts
γ_c(n) at five values of n with no free parameter. The radix penalty, the symmetry-breaking
eigenvalue, and the critical drive are all **the same integer count** — one amplifying
reaction against n−1 contracting ones — evaluated at different places.

---

### 56 The restoration trichotomy — capability is combinatorial, realisation is tuned

§53 concluded "restoration is tuned, not topological". That is half right, and the
decomposition already contained the other half.

**P is LINEAR in the rate constants.** Writing v_r(x) = d_r B_r(x), §54's decomposition is
P(x) = ⟨c, v(x)⟩. So a network fails to restore exactly when ⟨c, v(x)⟩ ≤ 0 for every
accessible x — and that condition is closed under addition and positive scaling:

> **THEOREM. The non-restoring rate vectors form a CONVEX CONE** — the polar of the cone
> generated by {v(x) : x accessible}.

Verified: over networks admitting two independent non-restoring rate vectors, **0 violations**
of c₁, c₂ non-restoring ⟹ c₁ + c₂ non-restoring. Linearity itself holds to **2.7×10⁻¹⁶**.

**Two corollaries follow, and they complete the classification.**

> **(a)** all d_r ≤ 0 ⟹ v(x) ≤ 0 componentwise ⟹ the *whole* positive orthant of rate
> constants lies in the cone. Never restores, for any rates.
>
> **(b)** some d_r > 0 ⟹ load c onto that reaction and choose x where its bracket
> dominates ⟹ **the network restores for SOME c.**

**So capability is decidable from stoichiometry alone:**

| class | prediction | networks | capable | flips under rate change |
|---|---|---|---|---|
| all d_r ≤ 0 | never restores | 157 | **0 (0%)** | 0 of 240 |
| all d_r ≥ 0, some > 0 | always restores | 2 | **2 (100%)** | 0 of 4 |
| **mixed** | **capable, realisation tuned** | 34 | **34 (100%)** | 24 of 37 |

> **A network can restore for some rate constants ⟺ some d_r > 0.**
> **It restores for every rate constant ⟺ all d_r ≥ 0.**
> **Otherwise the failing rates are a convex cone, and only there does tuning decide.**

**P6 — compactness does not bite.** The domination argument in (b) needs room to scale x, and
a conservation law confines the state to a compact simplex where that room might vanish. It
does not: **35/35 mixed networks are capable on the simplex**, the same 100% as on the open
orthant. So a conservation law does not forbid restoration a network is otherwise capable of
— which is the third time (with §42 and §43) that conservation structure has turned out not
to matter to this question.

**AM as the worked case.** d = {+1, −1} → **mixed → capable**, and γ decides whether it
realises it. γ_c = 1/2 is where the rate vector crosses out of the non-restoring cone, and
§55 showed the same crossing at every n with the cone boundary at γ_c(n).

### 56.1 This corrects §54, and the correction sharpens it

**`classify` was counting p = q pairs.** §55 found that those have an *empty* bracket sum —
contributing identically zero — while `mirror_pairs` records their d as an arbitrary ±1,
because "the X-heavy member" is undefined when p = q. §55 fixed the n-winner counter and
**the fix was never propagated to `classify`**, which §54's whole table rests on.

| class | §54 as published | corrected |
|---|---|---|
| all ≤ 0 | 113 | **211** |
| all ≥ 0 | 3 (1 amplifies) | **7 (7 amplify)** |
| mixed | 183 (33 amplify) | **63 (27 amplify)** |
| trivial | 1 | 19 |

**148 of 232 networks were misclassified.** And §54's one soft spot — "two of the three
all≥0 networks read P ≈ 0 rather than P > 0 … a threshold artifact, not a sign violation" —
**was not a threshold artifact. It was this bug.** Corrected, all≥0 amplifies **7 of 7**, and
P3 holds without an excuse attached.

The cross-check strengthens rather than weakens: **0 of 263 unanimous topologies flip, 24 of
37 mixed ones do (65%)**, against §54's 0/168 and 24/132.

> **What did not change:** the decomposition itself (1.9×10⁻¹⁴), the AM term-by-term
> recovery of §30's bracket, §55's γ_c(n) and λ(n) predictions. Those use the *brackets*,
> where the empty sum correctly gives zero. **The bug was in the classifier, not the
> algebra** — the same distinction §55 drew about its own counting bug, in the same place,
> which is why it should have been propagated then.

**Where this leaves the founding question.** §53 said a restoring element is a tuned object
rather than a topological one. §56 splits that:

* **Capability is topological.** 70 of 300 random symmetric networks (23%) can restore at
  some rates; the other 77% are forbidden by stoichiometry, no tuning possible.
* **Realisation is tuned.** Of the capable ones, 34/70 = 49% restore at randomly drawn rates,
  and the failing rates form a convex cone.

So "the transistor is a near-ideal restoring switch" decomposes into a structural claim that
is decidable and a tuning claim that is not — and 77% of the symmetric networks one might
write down cannot be tuned into a switch at all.

---

### 57 An exhaustive search rediscovers AM — and the optimum is what §54 predicts

The founding question, with every piece finally in place. §40 measured AM at **Q_min = 5.39**
against the thermodynamic floor Q ≥ 1. §56 made the search *tractable*: a network can restore
for some rates iff some d_r > 0, so the capable ones can be enumerated rather than stumbled
on.

**The family is finite.** On {X, Y, B} with X+Y+B conserved and bimolecular reactions there
are exactly **30 conservative reactions**, in **16 exchange-symmetric classes** (2
self-mirror, 14 mirror pairs). Networks from 1–3 classes number **696**, and **AM is one of
them**. Every network is parameterised exactly as AM is: forward rate 1, reverse rate γ.

**P1 GATE — and it earned its place.** AM's best in this pipeline is **Q = 5.4750 at
γ = 0.05**, against §40's published **5.39** — 1.6% apart, on a quantity computed through
entirely separate code.

> ⚠ **The gate caught a screening bug that had silently discarded the answer.**
> `slaving_axis.delta_star_of` brackets sign changes between fixed grid points, but `slaved`
> returns None as δ → 1, so a large δ\* falls in a dead zone: at γ = 0.03 and 0.05 the closed
> form gives 0.971 and 0.952 and `delta_star_of` returns **None**. Since AM's Q minimum sits
> at γ ≈ 0.05, **the first pass threw away exactly the best cells** and reported AM's optimum
> at γ = 0.10 with Q = 5.79. It was caught only because §40's published number disagreed.
> Fixed with a local bracket over consecutive *finite* samples, gated against
> `delta_star(γ)` to **2.1×10⁻¹³**. `delta_star_of` itself is left untouched — §36 and §39.2
> rest on it.
>
> **Scope limit, stated rather than fought:** below γ ≈ 0.05 `slaved` itself dies, because
> the attractor sits where the losing species is 4×10⁻⁵ of the tank. The γ grid starts at
> 0.05, which is also §40's own grid.

**P3 — the answer, and it is a small number.** Over 227 capable networks (32.6% of 696) and
39 cells that produce both a landscape and a valid Q:

| rank | Q | γ | reactions |
|---|---|---|---|
| **1** | **5.0045** | 0.05 | `X+Y→2B`; `B+X→2X`; **`2B→X+Y`** |
| 2 | 5.4750 | 0.05 | `X+Y→2B`; `B+X→2X` — **this is AM** |
| 3 | 5.4803 | 0.08 | `X+Y→2B`; `B+X→2X` |
| 4 | 5.7703 | 0.08 | `X+Y→2B`; `B+X→2X`; `2B→X+Y` |
| 5 | 6.1120 | 0.12 | `X+Y→2B`; `B+X→2X` |

> **Nothing beats AM by more than 9%.** And **every one of the top ten contains AM's two
> classes** — the search does not find an alternative motif, it finds AM with a decoration.

**P4 holds: 0 of 39 cells read Q < 1**, so §40's pre-registered instrument warning never
had to fire.

### 57.1 The winner is what §54 predicted, not a surprise

Written out, rank 1 is:

| reaction | rate |
|---|---|
| `X + Y → 2B` | 1.0 |
| `2B → X + Y` | **1.0** |
| `B + X → 2X`, `B + Y → 2Y` | 1.0 |
| `2X → B + X`, `2Y → B + Y` | γ = 0.05 |

**The disagreement channel is at detailed balance; only recruitment is driven.**

And §54 predicts exactly this. `X + Y → 2B` is **self-mirror**, so d = 0 and it contributes
**identically zero** to P — the same fact as §30's first cancellation and §51's discovery
that ρ never appears in μ. **A channel that contributes nothing to the signal but carries a
thermodynamic drive is spending entropy for nothing.** Setting it to equilibrium removes that
cost while leaving the drift untouched, and buys the 9%.

This also converges with **§44** from a different direction: §44 found ρ — speeding the
disagreement channel — a free lever worth 43–50%, for the same structural reason. Here an
*exhaustive* search, told nothing about §44, arrives at the same channel and neutralises its
drive instead of its speed.

> **Three limits on the headline, all real.** (i) **39 valid cells from 696 networks.** Most
> capable networks never form a landscape on this rate slice, which §56 explains rather than
> excuses: capability is combinatorial but realisation needs rates inside the complement of
> the non-restoring cone, and **a single γ traces a one-dimensional curve through a
> high-dimensional rate space.** (ii) The search is over 1–3 classes, bimolecular, three
> species. (iii) Ω = 200, where Q is converged to ~1.6% on AM (5.826 → 5.735 over
> Ω = 150…600) but to only ~25% at γ = 0.30.

> ⚠⚠ **§58 corrects this.** "Nothing beats AM by more than 9%" is true **on the slice
> only**. With free per-class rates, Q reaches **1.25–1.33** against AM's 5.475 — a 4.3×
> improvement §57's one-parameter curve could not see. But every network that beats AM does
> so **by being slower**: the Q-optimum rides the bifurcation at γ_rec = 0.49 (γ_c = 0.50)
> with a mean decision time of **783 against AM's 4.09**. On the (time, Q) Pareto frontier
> **AM is the fastest point and sits on the frontier**, so §57's conclusion survives in a
> stronger form: *AM is Pareto-optimal at its own speed*, not merely near-optimal on a slice.

**So the claim is bounded and it is still the answer to the founding question.** Within the
696-network family and a single-γ rate slice, **AM is within 9% of the best chemical decision
element there is, and every network that matches or beats it contains AM.** The canonical
consensus motif is not merely one workable design among many — at this size it is
essentially *the* design, and the only improvement available is to stop wasting drive on the
channel that carries no signal.

**T-OPT-a, open: does the full rate space change the answer?** §56 characterises the
restoring rates as the complement of a convex cone, so per-network rate optimisation is a
well-posed problem this search did not solve — it sampled one curve. **How to kill:** for the
top networks, optimise Q over independent per-class rates rather than a shared γ. If AM's
margin survives free rates, the optimality claim strengthens from "on a slice" to "in the
family"; if some network overtakes it, the slice was hiding the answer, exactly as the
screening bug above did.

---

### 58 AM is on the speed–optimality frontier — and §40's ruler had a factor of 10⁵⁷ in it

§57 concluded "nothing beats AM by more than 9%" on a one-parameter rate slice, and opened
T-OPT-a because §56 says that slice is a curve through a high-dimensional space. Optimising
freely changes the answer, and checking the ruler changes it again.

### 58.1 The one-sided setting is inapplicable, by a factor of 10⁵⁷

§40 named its own leading suspect in advance: *"our absorbing set is TWO-SIDED, |δ| ≥ thr,
where the standard statement is one-sided."* Every "how far is AM from the bound" number
rests on that. **Nothing had measured it.**

Measured now: one-sided absorption (at +thr only) makes the first passage **exponentially
long**, because reaching the *correct* threshold alone requires waiting out excursions into
the wrong basin — which need a barrier crossing back.

| Ω | 10 | 14 | 18 | 22 |
|---|---|---|---|---|
| ⟨T⟩ two-sided | 4.47 | 4.77 | 4.92 | 5.00 |
| ⟨T⟩ **one-sided** | 3.5×10⁵ | 4.9×10⁶ | 7.1×10⁷ | 1.1×10⁹ |

> **ln(T_one/T_two) = 0.630·Ω + 5.06**, so at Ω = 200 the ratio is **e¹³¹ ≈ 10⁵⁷**.

**The two-sided set is not a sloppy substitute for the one-sided one — the one-sided
quantity diverges exponentially in Ω.** The TUR's standard form is a statement about a
current reaching *one* threshold; a bistable decision element that has committed to the
wrong answer does not wait 10⁵⁷ time units to spontaneously fix itself. **§40's Q = 5.39 is
measured against a bound whose derivation assumes a setting this system exponentially cannot
occupy**, and that is now quantified rather than flagged.

### 58.2 My own optimiser failed, and a direct scan caught it

Nelder-Mead over free log-rates (4 restarts, 120 iterations) reported AM's optimum at
**ρ = 0.997, Q = 4.84** — i.e. that §44's ρ lever does not help Q. A direct scan says
otherwise:

| ρ | 1 | 2 | 4 | 8 | 20 |
|---|---|---|---|---|---|
| Q | 4.856 | 3.697 | 3.205 | 2.994 | **2.919** |

**The optimiser missed a 40% improvement lying along a single axis.** Replaced with grid
search. This is the third instrument failure in two sections — after §57's screening bug and
§55's counting bug — and like both, it was caught only by checking against something
independent.

### 58.3 The frontier, which is the actual answer

Grid search over (ρ, γ_dis, γ_rec), 274 cells across three refinements, drives Q down to
**1.2531** — against the thermodynamic floor of 1, and against AM's 5.4750. But the optimum
sits at **γ_rec = 0.49 with γ_c = 0.50**, riding the bifurcation, where δ\* has collapsed to
0.67 and the mean decision time is **783**.

**Q is dimensionless and therefore time-blind.** Plotting every cell as (time, Q) gives a
clean Pareto frontier:

| mean time | 4.04 | 5.91 | 9.88 | 18.80 | 45.51 | 112.7 | 241.0 | 783.3 |
|---|---|---|---|---|---|---|---|---|
| **Q** | **5.40** | 3.38 | 2.43 | 2.04 | 1.89 | 1.84 | 1.37 | **1.25** |

> **AM sits ON that frontier, at its fast end.** AM at γ = 0.05 is (t = 4.09, Q = 5.475);
> the frontier's fastest point is (4.04, 5.400). **AM is Pareto-optimal at its own operating
> speed.**
>
> **Buying 4.3× closer to the thermodynamic bound costs 194× in decision time.**

> ⚠ **§59 corrects the frontier's far end.** §58's parameterisation had too little rate
> freedom: pooling every rate assignment of AM's own reaction set reaches **Q = 1.115 at
> t = 3747**, not 1.253 at 783. The shape is unchanged and the trade-off sharpens —
> **5.5× closer to the bound costs 950× in time.** §59 also shows no *different* topology
> reaches AM's frontier at any speed.

**This reframes §40 and corrects §57.** §40's "Q_min = 5.39, about 5× from the bound" is not
a measure of AM's inefficiency — **it is the price of deciding in 4 time units instead of
800.** And §57's "nothing beats AM by more than 9%" was true only on the slice: with free
rates plenty of things beat AM *on Q*, but every one of them does it by being slower, and
none of them beats AM at AM's speed.

> **What Q is not.** §38 already found that "cost per nat of reliability" is not a quantity,
> because reliability is bought with free input margin. §58 is the same lesson one level up:
> **Q is not a design objective either, because it is blind to time.** A figure of merit that
> a system can improve arbitrarily by slowing down is measuring a trade-off, not a quality.
> The frontier is the quantity; Q alone is a coordinate on it.

> **Limits.** Ω = 200 throughout, with the best cell checked at Ω = 150/200/300/450 giving
> Q = 1.321 / 1.253 / 1.319 / 1.326 — non-monotone at the 5% level, so the best Q is
> **1.25–1.33, not 1.25 exactly**. Both γ_dis → 0 and γ_rec → γ_c are at their grid edges,
> so the frontier's far end is a bound, not a located optimum. The family is AM's two
> classes with free rates; §57's 696-network enumeration was not re-run freely.

**T-OPT-b, open: is the frontier's shape universal?** Q ~ 1 + a/t^b would say the approach
to the bound has a rate; the measured points span Q = 5.40 at t = 4 to Q = 1.25 at t = 783,
which is close to Q − 1 ~ t^(−0.5). **How to kill:** fit the frontier per network family and
check whether the exponent transfers — §39.2 and §46 both found coefficients that did not
transfer between axes, so the prior is that it will not, and a shared exponent would be the
surprise.

---

### 59 No other topology beats AM — the winners are AM re-parameterised — T-OPT-b

§58 named its own limit: §57's 696-network enumeration was never re-run with free rates, so
"AM is Pareto-optimal" was a statement about AM's own two classes. This tests it.

**§54 buys back the search dimensions.** A network with m classes has 2m rates. But §54
classifies every class by d_r: those with d_r ≠ 0 carry signal, those with d_r = 0 contribute
**identically zero** to the drift (§51, §54, §55 — one fact in three places) and act only
through the pool and the noise. §58's optimum has exactly that shape, so the search collapses
to (ρ_ns, γ_ns, γ_s) whatever m is. **Of the 16 classes, 14 are signal-carrying and 2 are
not.**

**P1 gate holds.** AM's frontier reproduces §58's slow end exactly — (783.30, 1.253) — and
its fast end at (3.93, 6.19) against §58's (4.04, 5.40), the reduced parameterisation being
coarser where §58's grid was fine.

### 59.1 The networks that beat AM are AM

Raw, five candidates beat AM somewhere. **Two survive an honest comparison, and both are AM.**

> ⚠ **My P2 test had an extrapolation flaw and it inflated the result.** It scored each
> network against AM by `np.interp` on AM's frontier, which returns the **endpoint value**
> outside the traced range. AM's frontier stops at t = 783, so *any* network slower than that
> automatically "beat" it. Restricting the comparison to AM's own traced range removes three
> of the five.

The two that survive are `AM+cls1` and `{cls1, dis}`, where

> **cls1 = {2X → B+X, 2Y → B+Y}** — which is AM's recruitment pair with the *forward
> direction relabelled*. Its reverse is `B+X → 2X`, the recruitment itself.

Checked directly: `{cls1, dis}` and `{dis, rec}` generate **identical reaction sets**. So they
are not different networks — they are AM with the recruitment pair's two rates decoupled from
the shared γ, which is rate freedom §58's parameterisation did not have.

> **P2's real answer: no structurally different topology beats AM at any speed.** Every
> candidate with genuinely different chemistry — `AM+revdis`, `AM+cls0`, `AM+cls2` — is worse
> at every overlapping time, by factors from 2× to 90×. **§58's headline strengthens from
> "AM is on its family's frontier" to "no other topology in the enumeration reaches it."**

**And it corrects §58's frontier**, which was traced with too little rate freedom. Pooling
every parameterisation of AM's reaction set:

| time | 3.9 | 4.4 | 39.5 | 80.6 | 189 | 782 | 1594 | 3747 |
|---|---|---|---|---|---|---|---|---|
| **Q** | 6.19 | 5.02 | 1.52 | 1.51 | 1.36 | 1.196 | 1.139 | **1.115** |

**AM's reaction set reaches Q = 1.115**, not §58's 1.253 — closer to the thermodynamic floor
of 1, and still by going slower: t = 3747 against 4.1 at the fast end. The trade-off §58
identified is unchanged in shape and sharpened in extent: **buying 5.5× closer to the bound
costs 950× in time.**

### 59.2 The frontier exponent does not transfer — a third failure

Fitting `Q − 1 = a·t^(−b)` per network:

| network | AM | AM+revdis | AM+cls0 | {cls1,dis} | AM+cls2 |
|---|---|---|---|---|---|
| **b** | 0.584 | 0.078 | 0.300 | 0.402 | 0.493 |

**0.078 to 0.584, a 136% spread.** There is no shared exponent, so **there is no universal
time-cost law** for approaching the thermodynamic bound in this family — the approach rate is
a property of the network, not of the bound.

> **P4 predicted this, and said so to avoid the flattering reading.** §39.2 found a 1/sep
> coefficient that did not transfer between axes; §46 found the *scaling* did not transfer
> either. **This is the third attempt at a transferable exponent in this project and the
> third failure.** At this point the prior should be explicit: in this system, exponents
> fitted on one axis do not carry to another, and any future claim of one needs a
> cross-axis test before it is written down, not after.

**P6: 0 cells read Q < 1**, across every network and every rate combination tried. §40's
instrument rule never had to fire, and §58 has since measured the suspect it was written for.

**Where the founding question stands.** A chemical decision element can approach the
thermodynamic bound to within **1.12×**, and AM's own chemistry is what does it. Nothing else
in the enumerated family comes close at any speed. The 4× gap §40 measured is not
inefficiency and not a better motif waiting to be found — **it is the price of deciding
quickly**, and AM occupies the fast end of the only frontier there is.

---

### 60 The fluctuation theorem factorises over outcomes — so it cannot bound the error

§41 verified `⟨e^(−S_tot)⟩ = 1` at the absorption time to 5.5×10⁻¹⁴ and **never split that 1
by outcome**. The split is the founding question in its sharpest form: if error paths carry
exponentially little entropy, the identity might pin the error probability to the
dissipation. Writing Φ_o = p_o·⟨e^(−S_tot)|o⟩, the hypothesis was **Φ_e = p_c**, i.e.
⟨e^(−S)|error⟩ = p_c/p_e — the odds of being right *are* the exponentiated error entropy.

Computed with §41's tilted generator (built from reverse propensities) under §35's
outcome-selective boundary. Two extra solves, no new machinery.

**P1 gate, after one failure.** On the first grid it failed spectacularly — γ = 0.10, Ω = 90
returned Φ_c ≈ 10²⁷⁹. That is **the conditioning failure §41 documented at exactly γ = 0.10**
(*"it fails at γ = 0.10, and the failure is conditioning, not physics"*), reproduced here.
Restricted to γ ≥ 0.20 and Ω ≤ 60, the gate holds: **|Φ_c + Φ_e − 1| ≤ 5.7×10⁻⁹, median
6.6×10⁻¹³** over 32 cells — §41's identity rebuilt from two independent solves.

**P2 is refuted by three orders of magnitude.** Φ_e/p_c spans **3.75×10⁻⁴ to 0.205**, never
near 1, and drifts systematically with every axis swept (γ: 0.0057→0.105; Ω: 0.071→0.015;
ε: 0.045→0.014), so P4's fallback — a constant ratio — fails too.

**The actual structure was in a column the prediction did not name.**

| | range | max deviation from 1 |
|---|---|---|
| ⟨e^(−S_tot) \| correct⟩ | 0.998567 – 1.002775 | **2.8×10⁻³** |
| ⟨e^(−S_tot) \| error⟩ | 0.949877 – 1.034502 | **5.0×10⁻²** |

> **⟨e^(−S_tot)⟩ = 1 holds for each outcome SEPARATELY**, not merely in aggregate — and the
> deviations carry **random sign**, flipping between adjacent cells, which is the signature
> of numerical error rather than a physical effect.

### 60.1 Why that is a negative result, and a sharp one

If ⟨e^(−S)|o⟩ = 1 for each outcome, then **Φ_o = p_o identically** and the aggregate identity
is p_c + p_e = 1 wearing a disguise.

> **The integral fluctuation theorem carries no information about which outcome occurred.**
> It cannot bound the error probability, because it factorises over outcomes. **Reliability
> is not dissipation in this exact sense**, and the hope that motivated this experiment is
> dead rather than merely unconfirmed.

> **The aggregate's exactness does not validate the split, and this is the trap worth
> naming.** Φ_c + Φ_e = 1 and p_c + p_e = 1 together force
> p_c(⟨e^(−S)|c⟩ − 1) + p_e(⟨e^(−S)|e⟩ − 1) = 0 **by construction**. The two deviations are
> exactly anti-correlated whatever they are, so §41's 5.5×10⁻¹⁴ says nothing about whether
> the outcome-wise identity holds. It is established here only to ~5% on the error branch —
> eleven orders worse than the aggregate — because Φ_e is a small number extracted from a
> tilted generator with an enormous dynamic range.

> ⚠ **P6's verdict was a false alarm from my own criterion, and that is now a pattern.** I
> wrote that ⟨e^(−S)|error⟩ *must* exceed 1 "for error paths", so anything below 1 would
> "suspect a sign error in sigma_local". The premise was wrong: error paths do **not** carry
> strongly negative total entropy production, because the two absorbing boundaries are
> exchange images of one another and so carry equal stationary weight — the system term
> cancels. Measured, ⟨e^(−S)|e⟩ ≈ 1. **The criterion was wrong, not the instrument.**
>
> **That is the fourth badly-designed verdict rule in this session**: §53's demand for a
> *majority* of topologies to flip before declaring non-combinatoriality; §55's counting of
> p = q pairs; §59's `np.interp` extrapolating flat beyond the traced range; and this. In
> every case the measurement was sound and the *summary rule* was not. **Pre-registering a
> prediction does not protect against pre-registering the wrong test of it**, and this
> project's rules — which are all about the measurement — do not currently say so.

**What survives.** §41's aggregate identity, reconfirmed through an independent decomposition
(5.7×10⁻⁹). The outcome-wise identity, to ~5%, consistent with exact and not established as
exact. And a closed door: **entropy production and decision outcome do not couple through the
fluctuation theorem**, so the exact reliability–dissipation relation this project has been
circling since §37 is not to be found there.

---

### 61 The reduction's exponent is exact and its prefactor is a START effect — T14-e

§35.3 proved the prefactor cannot be extracted **numerically** from the 2-D problem: the
candidate basis functions are collinear over any accessible Ω range (correlations 0.961 and
0.986, condition number 2.8×10⁸), so a slowly-varying prefactor and a slightly wrong exponent
cannot be told apart by fitting. That stands. **But §50's exact 1-D slaved chain is an
instrument §35.3 did not have**, because a birth–death chain's splitting probability is
closed-form:

    P(hit a before b | i) = Σ_{k=i}^{b−1} π_k / Σ_{k=a}^{b−1} π_k,   ln π_k = Σ [ln μ_j − ln λ_j]

evaluated in logs — **exact at any Ω, no solve, no fitting**, and naming the wrong outcome
directly as §35 requires. Verified against a sparse solve of the same chain to <10⁻⁸.

**And §44's ρ is the knob that separates the two.** An exponent error shows as a *slope* in
ln(P_1D/P_2D) versus Ω; a prefactor shows as an *intercept*. §39.2 says the reduction becomes
exact as sep → ∞, so driving ρ up kills the slope and leaves the intercept clean. §35.3 had
no such knob — ρ did not exist until §44.

| ρ | exponent error | ln(P_1D/P_2D) span | intercept |
|---|---|---|---|
| 64 | 0.107% | 0.069 | 0.3209 |
| 256 | — | 0.022 | 0.3201 |
| **1024** | **0.006%** | **0.0138** | **0.3188** |

At ρ = 1024 the ratio reads **0.3143, 0.3249, 0.3112, 0.3146, 0.3140** across Ω = 150…700 —
a span of **0.0138 nats while ln P itself spans ~70**. **The exponent is exact in the slaved
limit**, confirming §39.2 on a new axis and for the error probability rather than the time,
and **a constant prefactor ratio survives.**

### 61.1 The prefactor is a property of the START, not the boundary

The surviving constant is **not universal** — checked before writing it down rather than
after:

| γ | ε=0.35, θ=0.80 | ε=0.50, θ=0.80 | ε=0.35, **θ=0.70** |
|---|---|---|---|
| 0.10 | 1.6002 | 2.3790 | **1.6002** |
| 0.20 | 1.3720 | 1.8299 | **1.3720** |
| 0.30 | 1.2556 | 1.5484 | **1.2556** |

It varies with γ (1.26 → 1.60) and strongly with ε (1.37 → 1.83 at γ = 0.20). **So there is
no constant to name** — and per §49's precedent that declined π on two cells, the check came
first.

> **But the θ column is identical to four decimals, and the thresholds genuinely differ** —
> thr = 98 vs 86 at Ω = 150, 328 vs 287 at Ω = 500, a 14% change that moves nothing.
>
> **The prefactor discrepancy between the 1-D reduction and the exact 2-D CME depends on the
> START and not at all on where the threshold is read.**

That localises T14-e. Both descriptions relax into the same quasi-stationary escape mode; the
threshold only reads that mode out, so it cannot affect the ratio. What differs is **how much
of the initial condition feeds into the mode** — set by the start ε and the landscape shape
γ, and by nothing downstream. **The 1-D reduction gets the escape mode exactly right and its
excitation wrong.**

> ⚠ **Rule 19's first live catch, on the run that motivated it.** P3's verdict criterion
> thresholded the *absolute* drift at 0.2 nats. The four values were 0.138, 0.315, 0.185,
> 0.053 — straddling the threshold, so the verdict flip-flopped non-monotonically in ρ
> (FLAT, LINEAR, FLAT, FLAT) on data that is in fact monotone. **The criterion was measuring
> noise around its own threshold.** Normalising the drift by the range of ln P turns it into
> an exponent error and the picture becomes monotone and interpretable. Written one commit
> after the rule, and caught by it.

**What T14-e now has.** Not the prefactor — that remains uncomputed, and §35.3's proof that
fitting cannot supply it is untouched. What is new is its **structure**: the slaved
reduction's error is a pure start-side factor, exactly θ-independent, on top of an exponent
that is correct to 0.006% in the slaved limit. **A prefactor calculation now has a target
with a known form and a known set of arguments — (γ, ε) and not θ** — which is what an
analytic route (Assaf–Meerson, or the cat-qubit path integral of arXiv:2507.18714 noted in
THEORIES §5) would have to reproduce.

---

### 62 The restoration boundary is a closed form — T15-i, and §56's cone does not survive it

§56 proved the non-restoring rate vectors form a convex cone and that **capability** is
combinatorial: some d_r > 0 suffices. It said nothing about **realisation**, and left the
trichotomy half-finished — "tuned" is the absence of a criterion, not a criterion.

**§55 already contained the missing half.** It measured P at the symmetric fixed point and
got (1−2γ)/3 for AM, vanishing exactly at γ_c = 1/2. So the state that decides realisation
is not an arbitrary accessible x — **it is the symmetric steady state, the decision point the
dynamics actually occupies** — and §54 gives P there explicitly:

> **RESTORES ⟺ Σ_r c_r d_r B_r(x\*) > 0**,  d_r = S_X(r) − S_Y(r),  B_r(x\*) ≥ 0
>
> a **single linear inequality in the rate constants at fixed x\***, for any
> exchange-symmetric mass-action network — §43 makes P exist, §54 makes it explicit, and no
> simulation is involved.

**Gates.** §54's closed form equals §53's `P_at` on the symmetric line to **7.8e-11** over 40
random networks. For AM it reproduces §55's published (1−2γ)/3 to **8.3e-15**, and
root-finding the criterion returns **γ_c = 0.50000000000000**, |diff| **1.1e-15** — an
absolute check against a stored number, not a fit (rule 16).

### 62.1 Tested against the dynamics, not against another rule

sign(P at x\*) versus the full mass-action ODE integrated from x\* with an antisymmetric
kick. **120 decidable networks drawn from 434, 60 predicted restoring and 60 decaying, and
0 disagreements.** Excluded and counted rather than dropped: 32 with more than one symmetric
fixed point, 7 solver, 275 refused because their branch was already full. None marginal,
none ambiguous, none with a species going negative (rule 10's guard: the integrator never
clips, it discards).

> ⚠ **The first version of this test could not fail, and that is rule 19 again.** It drew 20
> networks, the criterion predicted **20 decaying and 0 restoring**, and it printed "P2 HOLDS:
> the criterion IS the boundary" — off a branch that never ran. Random class-sets are almost
> never capable, so the restoring half of a two-sided claim was simply absent. The sampling is
> now stratified and **the verdict is gated on both branches having ≥ 30 cases**, printing
> UNDER-TESTED otherwise. The 60/60 above is what that gate is for.

### 62.2 The realising set is NOT convex — and the counterexample is AM itself

§56's cone theorem needs P linear in c at **fixed** x. Here x\*(c) moves with c, so the
argument cannot transport. Searching 6000 pairs, of which **210 had both members restoring**
(only those can test convexity), found one counterexample — in **AM's own two classes**,
{X+Y→2B} and {B+X→2X, B+Y→2Y}, with free forward and reverse rates:

| | k_dis^f | k_dis^r | k_rec^f | k_rec^r | u\* | P | ODE |
|---|---|---|---|---|---|---|---|
| c₁ | 0.8273 | 7.9309 | 5.2147 | 0.4279 | 0.45718 | **+5.535e−2** | restores, spread ×6236 |
| c₂ | 1.6947 | 0.0782 | 2.8352 | 1.3938 | 0.32696 | **+6.975e−2** | restores, spread ×5659 |
| **c₁+c₂** | 2.5220 | 8.0091 | 8.0499 | 1.8217 | 0.41858 | **−2.142e−1** | **decays, ×1.1e−12** |

Verified three independent ways per rule 14 — the closed form, §53's `P_at` (agreeing to all
printed digits), and the ODE — and none of the three values is marginal.

> **Two restoring elements can sum to a non-restoring one.** §56's convexity is a property of
> **capability** and does not extend to **realisation**; the mechanism is visible in the table,
> since x\* moved from 0.457 and 0.327 to 0.419 and P is not monotone along that path.

### 62.3 Critical points this project did not have

The criterion is constructive, so γ_c follows for any class-set by one root-find:

| network | γ_c | u\* at γ_c |
|---|---|---|
| AM {X+Y→2B, B+X→2X} | **0.50000000** (gate: published) | 0.33333333 |
| AM + {X+Y→2X, X+Y→2Y} | **0.16666667** (= 1/6) | 0.33333333 |
| AM + {X+Y→B+X, X+Y→B+Y} | **0.38829144** | **0.30585428** |
| AM + five other classes | — | restore at **no** γ (P from −8.0e−5 to −3.3e0) |

The third is the interesting one: adding that class **moves the symmetric fixed point off
1/3**, which is why its γ_c is not a simple fraction. The first two keep x\* at 1/3 and give
exact rationals.

### 62.4 Scope, and a second criterion caught by rule 19

The criterion is a linearisation, so rule 9 demands an axis I did not choose: the kick size.
Over d₀ = 10⁻⁶ … 10⁻¹, **40/40 verdicts agree at every decade** — d₀-independent, so on this
family the criterion is global and not merely local.

> ⚠ **It did not read that way at first, and the fault was again the verdict rule.** The
> original P4 declared restoration at "spread ratio ≥ 10". At d₀ = 0.1 that printed **20
> ambiguous**, which looks like the criterion failing at large kicks. Reading the cells
> (rule 18) instead of the summary: all 20 were restoring networks whose spread had
> **saturated at s₁ = 0.70–1.00** — they restored completely — while s₀ was already 0.11–0.43.
> The normalised spread cannot exceed 1, so the largest attainable ratio is 1/s₀, and that was
> **2.3–9.5 in every one of the 20**. *The threshold was unreachable by construction.* A
> conservative geometric estimate caught only 5 of the 20; the cells had to be read one by
> one. Replaced by a ceiling-aware rule, s₁ ≥ min(10 s₀, 0.5), which is satisfiable at every
> d₀ and identical at the small d₀ where §62.1 runs — §62.1's 60/60 is unchanged by the fix.

**Stated scope.** Conservative bimolecular networks on {X, Y, B} with a **unique** symmetric
steady state; 7.4% of draws had more than one and are outside the claim, since which fixed
point decides is then undetermined. The dynamical test is mass-action, not CME — this is a
statement about the deterministic restoring element, which is what §56's trichotomy was about.

---

### 63 How sharp is a restoring threshold made of Ω molecules? — T15-k

§9.1 settled the deterministic picture: the AM landscape dies at γ_c = 1/2 and above it "no
population size Ω can restore, because there is nothing to restore toward". §62 then made
that boundary exact and free. **So the deterministic side has zero uncertainty, and anything
measured against it is entirely the CME's** — the setup this project has wanted since §9.

The device question §9.1 never asked: *below* γ_c, how many molecules does a switch need
before its threshold is sharp? A transistor restores well because its threshold is abrupt.
One made of Ω molecules cannot be.

**The instrument is exact and was implied by §53 all along.** §53 defined P as the
antisymmetric eigenvalue of the *Jacobian*. The CME generator has an antisymmetric sector
too: exchange X↔Y commutes with Q for any exchange-symmetric network (§43's premise), so Q
block-diagonalises, and on antisymmetric observables

    Q_A[p,q] = Q[s_p, s_q] − Q[s_p, σ(s_q)]        over representatives with n_X > n_Y

with the n_X = n_Y states dropping out identically. **λ_A, its leading eigenvalue, is the
exact stochastic counterpart of §53's P** — no simulation, no threshold, no first-passage
definition. Gate: spec(Q_A) ⊂ spec(Q) to **5.7e-13**; λ_A < 0 at every finite Ω, as ergodicity
requires — **there is no sign change at finite size, and that is the point.**

**The two sides meet (absolute, no fit).** Above γ_c the CME's antisymmetric relaxation must
converge to §62's exact (1−2γ)/3:

| γ | Ω=20 | 40 | 80 | 120 | approach |
|---|---|---|---|---|---|
| 0.60 | 1.1647 | 1.1289 | 1.0857 | 1.0643 | ~Ω^−0.52 |
| 0.75 | 0.9290 | 0.9664 | 0.9838 | 0.9893 | ~Ω^−1.06 |

Converging to 1 from both sides — and **the approach is twice as slow near γ_c**, which is
the critical slowing this section goes on to measure.

### 63.1 The width of the threshold goes as Ω^(−1/2)

Metastability is measured as e-folds of excess, and the width as the **γ-interval between two
fixed metastability levels** — a level *difference*, so it needs no "E ≈ 0" endpoint. Two
instruments, one anchored to §62's exact deterministic rate and one mentioning no reference
at all (D = ln|λ_A(2Ω)| − ln|λ_A(Ω)|, finite everywhere):

| levels | with the exact reference | **with no reference** |
|---|---|---|
| 1→2 | −0.4744 | **−0.4934** |
| 1→4 | −0.4856 | **−0.5069** |
| 2→4 | −0.4915 | **−0.5133** |

Stable across levels (3.5% spread, rule 13 satisfied *before* comparing across Ω), and
**ρ-independent** across §44's lever (−0.5069, −0.4789, −0.4748; 6.6%), so the width belongs
to the threshold and not to the timescale separation.

> **w(Ω) ~ Ω^(−1/2).** The blur on a chemical switch's threshold falls as the square root of
> its molecule count: **to halve the blur, use four times the molecules.**

Per rule 15, the extrapolation is *not* resolved even though the exponent is: at Ω = 1000 the
power law gives w = 0.0174, c/ln Ω + b gives 0.0141, and a/Ω + b gives 0.0336 — a 90% spread.
The exponent over the measured decade stands; the value beyond it does not.

### 63.2 The exponent is *nearly* the pitchfork — and 2 is excluded

The exponent is not free. If the escape action vanishes as A ∝ (γ_c−γ)^ν then w ~ Ω^(−1/ν),
so −1/2 says ν = 2 — and **ν = 2 is what §9.1 predicts**, since it identified the bifurcation
as a pitchfork (δ\* ∝ √(γ_c−γ), 3 fixed points below and 1 above from 1830 starts) and a
pitchfork's barrier goes as (γ_c−γ)². Fitting ln A against ln(1/2 − γ) with **γ_c held fixed
at the value §62 proves exactly** — so this is a test, not a two-parameter fit:

| window | n | ν |
|---|---|---|
| γ ∈ [0.20, 0.45] | 8 | 1.9507 |
| γ ∈ [0.30, 0.45] | 6 | 1.9465 |
| γ ∈ [0.38, 0.45] | 4 | 1.9498 |
| γ ∈ [0.41, 0.45] | 3 | 1.9517 |

> ~~**ν = 1.9496 ± 0.0026, flat and non-monotone across nested windows, with 2 lying 9.7×
> the scatter away. The pitchfork value is excluded over the accessible window.**~~
>
> **WITHDRAWN by §64.** The gate admitted points whose remaining error was 4× the drift it
> measured — ~6% at γ = 0.45 against ~0.1% at γ = 0.38 — and that bias runs in exactly the
> direction that manufactures a deficit below 2. Left as first printed per rule 7; the
> corrected reading is ν ≈ 2 ± 0.1, and no instrument here excludes the pitchfork.

**Scope, and it is the whole caveat.** The narrowest window still stops at γ = 0.45, and the
γ = 0.46 point was excluded because its action had not converged in Ω (4.33% drift). So
**whether ν → 2 in the last 10% of the approach to γ_c is not resolved here** — that needs
Ω well beyond 500. What the data exclude is ν = 2 *over γ ∈ [0.20, 0.45]*, and they show no
drift toward it.

**The two measurements agree, and that agreement is NOT a second confirmation.** 1/ν = 0.5129
against the reference-free width exponent 0.5133 at levels 2→4 — 0.1%. But both are readings
of the same measured surface λ_A(γ, Ω), one along Ω and one along γ, related *by the scaling
ansatz itself*. **So their agreement confirms the ansatz — that the width is set by the action
alone — and does not measure ν twice.** This experiment's own P4 predicted that trap and it
is honoured here rather than quietly enjoyed.

### 63.3 Three verdict rules wrong in one section

Rule 19 fired three times, on measurements that were all sound:

> ⚠ **P1(c) printed FAILS on converging data.** It demanded |ratio−1| < 0.05 at Ω = 120 and
> saw 1.0643. The series was 1.165 → 1.129 → 1.086 → 1.064 and 0.929 → 0.966 → 0.984 → 0.989 —
> converging cleanly from both sides. *A fixed-Ω tolerance tests the size of a finite-Ω
> correction, not whether it vanishes.* Replaced by monotone convergence with its exponent.
>
> ⚠ **P2's crossover statistic could not find a crossover.** γ\* = argmax|d ln|λ_A|/dγ|, chosen
> because it is parameter-free "so it cannot be an artifact of a threshold". In the metastable
> branch the slope is −ΩA′(γ), *largest deep inside that branch and smallest at the crossover* —
> so the argmax runs to the low-γ edge of whatever window is swept, and it duly printed "max
> slope at the EDGE" for every Ω. **Being free of a threshold does not make a statistic measure
> the thing you named it after.** No data would have made it print a crossover.
>
> ⚠ **P3 read a 0.001 wobble as a trend.** It compared |ν−2| at the widest window against the
> narrowest and printed **"P2 HOLDS: ν → 2"** because 0.0483 < 0.05. The values were 1.9507,
> 1.9465, 1.9498, 1.9517 — flat, non-monotone, scatter 0.005. A rule that cannot tell a drift
> from a constant is the only thing it was there to do. Replaced by comparing the trend against
> the scatter, which prints the opposite verdict on the same numbers.

And rule 13 caught the estimator before any of this: **A = −ln|λ_A|/Ω drifted 8–39% across
Ω**, because dividing by Ω leaves the whole WKB prefactor as an O(ln Ω/Ω) contamination. The
gate refused all 7 γ and reported no ν. The repair is the local slope in Ω, which cancels the
constant exactly — and deliberately **not** a three-parameter fit, since §35.3 proved Ω and
ln Ω collinear over any bounded window here and §35.1's b was withdrawn for exactly that.

---

### 64 §63.2's ν = 1.95 was an unconverged estimator — the exclusion of 2 is withdrawn

**§63.2 gated the wrong quantity.** It extracted A as the local slope of −ln|λ_A| in Ω and
required the *drift* between the last two slopes to be under 2%. If A_eff(Ω̄) = A + c/Ω̄, that
drift is c(1/Ω̄₁ − 1/Ω̄₂) while the **remaining error is c/Ω̄₂** — and for those Ω ladders the
second is 4× the first. Measured, not assumed:

| γ | 0.38 | 0.41 | 0.43 | 0.45 | 0.46 |
|---|---|---|---|---|---|
| drift (what §63.2 gated) | 0.03% | 0.12% | 0.20% | 1.49% | 4.33% |
| **remaining error** | 0.14% | 0.48% | 0.79% | **5.95%** | 17.31% |

So §63.2 admitted a point carrying ~6% error at γ = 0.45 against ~0.1% at γ = 0.38. **That
bias is systematic in γ and signed** — A_eff approaches A from above, worse toward γ_c — which
flattens ln A against ln(γ_c−γ) and *lowers* ν. Predicted correction, written before the run:
ν should move **up** by ≈ 0.067, to 2.02 ± 0.03.

**The direction was right, the value was not.** After extrapolating A_eff → A:

| window | §63.2 (unextrapolated) | **extrapolated** | ladder B (disjoint Ω) |
|---|---|---|---|
| γ ≥ 0.24 | — | 2.0425 | 2.0569 |
| γ ≥ 0.35 | 1.9498 | 2.0836 | 2.1030 |
| γ ≥ 0.41 | 1.9517 | **2.1941** | **2.2141** |

ν moved **+0.15**, up as predicted (P2 hit); the predicted 2.02 ± 0.03 was **missed** (2.10).
A disjoint Ω ladder agrees to 2.5% in A, so the extrapolation is not ladder-dependent.

### 64.1 Three extractions, and rule 14 refuses the correction

A correction is a claim, so it gets an independent instrument. The action is also the
quasipotential barrier, and −ln π/Ω → V for the CME whether or not detailed balance holds —
computable from the stationary distribution alone, with no eigenvalue, no antisymmetric block
and no local slope. And §63.1's **width** exponent gives 1/ν with *no action extrapolation at
all*. Three routes:

| route | extrapolation done | ν |
|---|---|---|
| width of the transition (§63.1) | **none** | 2.027, 1.973, **1.948** |
| stationary distribution π | mild (Ω ≤ 220) | 1.852 → 1.986 → **1.991**, drifting toward 2 |
| λ_A action, extrapolated | heavy (10.9% and 25.4% shifts at γ = 0.45, 0.46) | 2.04 → **2.19**, rising |

> **They disagree, so the corrected value is NOT established.** The route doing the most
> extrapolation is the outlier. What all three agree on is that ν is **near 2**, and none of
> them supports 1.95 as a distinct value.

**~~§63.2: "ν = 1.9496 ± 0.0026 … 2 is excluded at 9.7× the scatter"~~ — WITHDRAWN.** The
number was produced by an estimator whose remaining error was 4× the gate that admitted it,
and the bias ran the right way to manufacture exactly that deficit. **Present statement:
ν ≈ 2 ± 0.1, consistent with §9.1's pitchfork, and not determined more precisely than that by
any instrument here.**

**What survives untouched: §63.1.** The width exponent −0.4934 / −0.5069 / −0.5133 was measured
directly from λ_A's Ω-scaling with no action extrapolation anywhere, so nothing in this section
touches it. **w(Ω) ~ Ω^(−1/2) stands**, as does the device reading — four times the molecules
for half the blur.

> ⚠ **And §63.1's "agreement confirms the ansatz" is weakened.** §63.1 noted 1/ν = 0.5129
> against the width's 0.5133 and read it as confirming that the width is set by the action
> alone. Both numbers came from the same biased ν. With the bias removed the action route
> gives 2.10–2.19 against the width's 1.95–2.03, so **the ansatz is no longer confirmed — it
> is untested**, and the tension is itself the open question.

### 64.2 Two more verdict rules, and one withdrawn kill test

> ⚠ **P1's residual was normalised by the range of the points being fitted** — a range that
> shrinks toward zero as convergence improves, so it punished exactly the γ it should pass:
> γ = 0.43 shifts by 0.62% and scored 0.41, γ = 0.46 shifts by 25% and scored 0.019. It
> rejected all 8 γ and printed no ν. A residual must be normalised by the quantity being
> estimated, not by the spread of the estimates.
>
> ⚠ **The three-way verdict helper had an unreachable branch.** It tested `trend > scatter`
> with trend = |ν_first − 2| − |ν_last − 2|; for any sequence *monotone toward 2* those are
> identically equal, so a drifting sequence could never print "(c) still drifting" and fell
> through to "(a) ν = 2 within the scatter" — which it duly did on 1.8522, 1.9857, 1.9909.
> Replaced, and **checked against §63.2's own numbers, where the new rule still returns that
> section's verdict** — so it is not tuned toward the answer now expected.

> **T15-m's kill test is withdrawn as unreachable.** §63 proposed settling ν by pushing the
> action measurement to γ = 0.47–0.49, calling it "reachable rather than hypothetical" on the
> strength of Ω = 640 costing 7s. At γ = 0.47 the slope still drifts **4.2% at Ω = 1000**
> (17.8s, 5×10⁵ states), so 1% needs Ω ≈ 4000 and ~8×10⁶ states. **The estimate was made from
> wall-clock without asking what Ω the convergence required.** A kill test costed by runtime
> instead of by the convergence it must reach is not a kill test.

---

### 65 The Symmetric Restoration Theorem — stated once, with its prior art

§43, §54, §55, §56, §62 and §63 have been reported as six results. They are six clauses of
one statement, and it was never written down. It is written here, with the scope it actually
has and with the literature it sits next to — **which had not been checked for any of §42–§64
before now**, in violation of this project's own rule that the search happens when a question
is *named*, not when it is finished.

> **THE SYMMETRIC RESTORATION THEOREM.**
> Let `N` be a mass-action reaction network on species including X and Y, invariant under the
> exchange σ: X ↔ Y, with rate vector `c > 0`. Write `b = S·v(n)` for the mass-action field,
> `d_r = S_X(r) − S_Y(r)` for the signed stoichiometric asymmetry of reaction `r`, and
> `x*` for a symmetric steady state (`x*_X = x*_Y`). Then:
>
> **(1) Invariance.** `b_X − b_Y = (n_X − n_Y)·P(n)` with `P` symmetric — so `{n_X = n_Y}` is
> invariant, **the sign of a lead is a deterministic invariant, and every reversal is a
> fluctuation.**
>
> **(2) Decomposition.** `P(x) = Σ_r c_r · d_r · B_r(x)`, where each `B_r(x) ≥ 0` on the
> positive orthant and is an explicit monomial expression fixed by the reaction's
> stoichiometry alone. Reactions with `d_r = 0` — self-mirror, or entering X and Y at equal
> order — contribute **identically zero**.
>
> **(3) Capability is combinatorial.** `N` restores for *some* rate vector ⟺ **some `d_r > 0`**;
> it fails to restore for *every* rate vector ⟺ all `d_r ≤ 0`. The non-restoring rate vectors
> form a **convex cone** (the polar of the cone generated by {v(x)}), because `P` is linear
> in `c` at fixed `x`.
>
> **(4) Realisation is a single linear inequality.** `N` restores ⟺ `Σ_r c_r d_r B_r(x*) > 0`.
> The realising set is a **cone but is *not* convex**, because `x*(c)` moves with `c` — two
> restoring rate vectors can sum to a non-restoring one.
>
> **(5) Finite Ω.** The CME generator commutes with σ, so it block-diagonalises; on
> antisymmetric observables the block is `Q_A[p,q] = Q[s_p,s_q] − Q[s_p, σ(s_q)]`. Its leading
> eigenvalue `λ_A` is the exact stochastic counterpart of `P`, is **negative at every finite
> Ω**, and the deterministic sign change of (4) becomes a metastability transition of width
> **~ Ω^(−1/2)**.

**Scope, which is part of the statement.** Mass-action; exchange symmetry (clause 1 fails
without it — §42 measured 1.9e1 and 2.9e2 residuals when the symmetry is broken); a *unique*
symmetric steady state for clause 4 (7.4% of sampled networks had more than one and are
outside it, T15-j); clause 4 is **linear stability of the symmetric branch**, and should be
stated as such rather than as "the restoration boundary"; clause 5's Ω^(−1/2) is measured on
AM, and the barrier exponent behind it is only `ν ≈ 2 ± 0.1` (§64).

### 65.1 What is prior art, and what is left

**Clause 1 is standard.** Flow-invariance of the fixed-point subspace `Fix(Z₂)` is a *folk
theorem* of equivariant dynamics (Golubitsky, Stewart & Schaeffer, *Singularities and Groups
in Bifurcation Theory* II, 1988). Given it, the polynomial divisibility follows immediately,
since `b_X − b_Y` vanishes on the irreducible variety `{n_X = n_Y}`. §43 framed the invariance
as the discovery and that framing is withdrawn.

**The nearest neighbour is not CRNT — it is spontaneous mirror-symmetry breaking.** An
exchange-symmetric network whose two symmetric species are enantiomers L and D, asked when it
leaves the racemic state, is *the* central problem of the homochirality literature, from Frank
(1953) through Ribó & Hochberg to the algorithmic work of Montoya, Cruz & Ágreda (*Life* 9:74,
2019) and the `Listanalchem` tool for searching mechanisms by linear stability. **Clause 4 is,
in substance, linear stability of the racemic steady state, which that literature does
routinely.** Reported there, the resulting conditions are *semialgebraic and hard to sample*,
linearised where possible via Clarke's stoichiometric network analysis. Adjacent CRNT anchors:
Craciun–Feinberg's species-reaction graph, Joshi & Shiu on multistationarity, and
"Switching in mass action networks based on linear inequalities" (arXiv:1002.1054).

**So the contribution is not that a criterion exists — it is that the criterion has a fixed
sign structure.** Clause 2 says every reaction enters `P` as `c_r · d_r · (something ≥ 0)`,
which is what makes clauses 3 and 4 possible: it **separates topology from rates**, reduces
capability to the signs of a list of integers, and attributes amplification to *named
reactions*. A generic linear-stability calculation gives a number; this gives the number's
decomposition.

> ⚠ **NOVELTY UNVERIFIED, and labelled as such per rule 17.** Whether clause 2's decomposition
> is already in the mirror-symmetry-breaking literature has *not* been settled — that needs
> the Montoya et al. and MATCH algebraic-analysis papers read properly, not their abstracts.
> **How to kill:** read them and check whether the racemic-stability condition is ever written
> as a nonneg-weighted sum over reactions with combinatorial signs. Until then this section
> claims priority for nothing; it claims only that the statement above is true, tested, and
> was worth assembling into one object.

**Clause 5's exponent is also expected, not discovered.** Ω^(−1/2) rounding of a Z₂ pitchfork
is standard mean-field finite-size scaling (Curie–Weiss). §63's contribution is the *exact
antisymmetric-sector instrument*, which gives the stochastic order parameter with no
simulation and no first-passage definition — not the exponent it produced.

---

### 66 §60's closure was tested only where symmetry guarantees it — T-TUR-e

§60 split §41's fluctuation-theorem identity by outcome, found `⟨e^(−S_tot)|o⟩ = 1` separately
for each outcome, and concluded the identity cannot bound the error rate. That closed the
founding question's sharpest form. **But §60's own explanation of *why* is a fact about
β = 0:** it reasoned that the two absorbing boundaries are exchange images carrying equal
stationary weight, so the system term cancels outcome by outcome. **For any tilted network —
that is, for every real device, since an inverter drives toward one rail and not the other —
the boundaries are not exchange images and the cancellation has no reason to occur.** §60
stated a general closure and tested the one case whose symmetry guarantees the answer.

Repeated on `am_asymmetric` at γ = 0.25, every cell strictly below β_c = 0.6365 (P4), with the
start taken from the **separatrix** (§31's matched rule) rather than a symmetric point, and
`b` on the slow manifold solved generically so no `am_reversible` closed form leaks in:

| β | β/β_c | mean \|r_c − 1\| | mean \|r_e − 1\| | cells |
|---|---|---|---|---|
| 0.00 | 0.00 | 0.000468 | **0.001071** | 6 |
| 0.10 | 0.16 | 0.002088 | 0.005204 | 5 |
| 0.20 | 0.31 | 0.002030 | 0.006522 | 4 |
| 0.30 | 0.47 | 0.010045 | 0.028819 | 3 |
| 0.40 | 0.63 | 0.003916 | **0.042831** | 2 |

**`|r_e − 1|` rises 40× with the tilt.** Whether that is the physics is *not* settled, and
§66.1 is why.

### 66.1 The instrument fails exactly where the effect should be largest

**The confound here is the dangerous kind: the nuisance is caused by the variable under test.**
Tilt is what makes the two boundaries' stationary weights differ — and that weight ratio `w` is
precisely what breaks the solve, because §41's convention weights the boundary by π(n)/π(n₀).
At β = 0.40, Ω = 90 the error boundary carries **ln w = −38.9**, i.e. e^−38.9 ≈ 10⁻¹⁷ of the
correct boundary's weight — below double precision. So `Φ_all` itself comes back as 1.161 or
0.813 instead of 1. **A deviation rising with β is predicted by both hypotheses.**

A per-cell precision gate (`|Φ_all − 1| < 10⁻⁶`) was added and 37 cells survive with a maximum
of 9.4e−7. Breaking the confound the way this project has broken confounds before — sweeping
the nuisance *against* the cause, since |ln w| also rises with Ω at fixed β, so matched-|ln w|
cells exist at different β:

| \|ln w\| band | cells | β range | corr(β, \|r−1\|) |
|---|---|---|---|
| 0.0 – 4.4 | 7 | 0.05–0.10 | +0.440 |
| 4.4 – 11.3 | 13 | 0.10–0.30 | +0.660 |
| 11.3 – 20.4 | 13 | 0.15–0.40 | +0.637 |

> **SUGGESTIVE ONLY, and stated as such.** Mean within-band correlation +0.579 on ~10 cells
> per band is too weak to separate tilt from solve error. §60's mechanism fares no better:
> corr(|ln w|, |r−1|) = +0.546, so the boundary-weight account remains a **suspect**, not a
> result (rule 17).

**The precision wall is itself part of the finding.** The accessible β range is bounded by
arithmetic, and bounded *exactly where the predicted effect is largest*, because the tilt that
would break the factorisation is the same tilt that makes the two boundary weights differ by
more than 10¹⁷. This method cannot answer the question it was built for beyond β/β_c ≈ 0.6.

> **What does NOT depend on resolving any of that: §60's scope.** It claimed a general closure
> from a case whose symmetry guarantees the result, and offered a mechanism that cannot apply
> once the symmetry is gone. **§60 is rescoped to exchange-symmetric elements.** For tilted
> ones the factorisation is *untested*, not established — and the founding question's sharpest
> form is **open where real devices live**, not closed.

### 66.2 The verdict rule was unit-tested before the experiment ran, and it was broken

§63 and §64 lost four verdict rules to defects a five-line test would have caught. Per that
lesson this section's rule was written into `tests/test_outcome_split_tilted.py` **first**,
with data engineered to trigger each branch. It failed two of them immediately:

> ⚠ It printed **DEVIATION** on `[0.30, 0.31, 0.32]` with a β = 0 residual of 0.30 — a 6% rise
> sitting entirely on its own baseline. And it could not distinguish a **physical baseline**
> (the β = 0 residual: a signal inside it means *no effect*) from **instrument noise** (the Ω
> scatter: a signal inside it means *not measurable*), so it would have printed "no deviation"
> on data that was merely too noisy to read. Both were fixed before a single cell was solved.

**This is the first section where no verdict rule had to be corrected after seeing results** —
five sections in a row had one. The cost was about ten minutes.

> ⚠ One thing the tests did not cover, caught only by reading the output: **P5's offset sweep
> applied no precision gate.** Its first pass reported |r_e − 1| ≈ 0.988 at β = 0.40 across all
> three offsets and concluded the effect "survives the offset sweep" — and all three of those
> cells are the ln w = −38.9 precision failure. Gated, **all three are excluded and P5 has
> nothing to say at β = 0.40.** A gate applied in one place and not another is the same class
> of defect as a gate on the wrong quantity (§64), and it took reading the cells to see it.

---

### 67 The rulers on an element with no symmetry to break — T-THM-b

Every quantitative result in this project is measured on Approximate Majority or its
exchange-symmetric relatives, §65's theorem is conditioned on that symmetry outright, and the
founding object — an inverter driving toward one rail — is asymmetric. So the founding
question's quantitative answer has **n = 1**. Schlögl's model is the textbook bistable CRN,
chemostatted, with **one dynamic species and therefore no exchange symmetry to break**:

    2X ⇌ 3X   (k1a = k1[A], k1r)        ∅ ⇌ X   (k2b = k2[B], k2r)

It is 1-D, so the chain is **exact** — splitting probability in closed form, MFPT, its
variance, and the mean entropy production all tridiagonal solves. A failure to transfer
cannot be blamed on numerics here.

Placing the fixed points at x₀−m, x₀, x₀+m gives `k1a/k1r = 3x₀`, `k2r/k1r = 3x₀²−m²`,
`k2b/k1r = x₀(x₀²−m²)`; m → 0 is where bistability dies, the analogue of γ → γ_c.

**A structural difference, before any number.** In Schlögl the affinity is *determined* by the
landscape — `A = ln[3(3x₀²−m²)/(x₀²−m²)]`. There is **no free drive knob at fixed landscape**,
unlike AM's γ. That is a fact about the substrate.

### 67.1 Both substrates have an affinity floor, and they are 3 ln 2 against 2 ln 3

Derived analytically, then confirmed against the engine's `cycle_affinity` — which takes the
null space of the per-pair stoichiometry and knows nothing of the formula — agreeing to
**8.9e−16** across seven cells. As m → 0:

> **Schlögl:  A_c = ln 9 = 2 ln 3 = 2.1972245773, independent of x₀** (checked at x₀ = 0.4,
> 1.0, 2.5).
> **AM (§9.1): A_c = 3 ln 2 = 2.0794415417.**
>
> **Different — so the floor is not universal — but 5.66% apart, from chemistry with nothing
> in common.** Both are ln(small integer): ln 8 against ln 9.

§9.1 read its own floor as "3 reactions × ln 2" and called the Landauer resemblance
arithmetic coincidence. This substrate has **2 reversible pairs and gives 2 × ln 3**, so the
pattern is `(pairs) × ln(pairs+1)` in both cases — which is now a testable statement rather
than a decorative one, and a third substrate would settle it.

### 67.2 §38's cost per e-fold has no counterpart here — two definitions, both refuted

The raw numbers are G = 8.3–27.1 against AM's ~2, and Q = 388–7102 against AM's 5.475.
**Reporting those as "does not transfer" would have been a rule-11 error**: AM is closed and
conservative, so all its dissipation belongs to the decision, while Schlögl is chemostatted
and burns entropy sitting still. Σ = 14293 at Ω = 800 against a cycle affinity of 2.26 is
housekeeping, not restoration.

| repair attempted | result |
|---|---|
| subtract housekeeping, `Σ − σ_local(x₀)·⟨T⟩` | **negative in 5/5 cells** (−1316 … −5878). Not a small number — an invalid definition. |
| dimensionless cycles per molecule, `Σ/(ΩA)` | 2.56–9.67 against AM's 0.342 — but it **grows with Ω** (4.52, 6.46, 7.44 at Ω = 200, 400, 800), so it is not intensive either. |

The reason no subtraction works is structural: a 1-D birth–death chain has **zero stationary
probability current**, so essentially all of σ_local is adiabatic and the non-adiabatic part
is the system term `ln[π(n₀)/π(n_f)]`, which is *negative* along a trajectory running from the
unstable point down to a rail.

> **NO dissipation-based cost measure ported.** The chemostatted element's dissipation is
> dominated by a term proportional to the decision *time*, for which a closed element has no
> analogue. **§38's G is not merely a different number on this substrate — it is not the same
> quantity.** The founding question's "price of restoration" is defined relative to closed
> conservative bookkeeping, and that bookkeeping does not survive the move to a driven device.

**P5 is therefore NOT DECIDED**: Q inherits the same non-comparability and is reported, not
compared. (No cell fell below §40's floor of 1, so the instrument is at least not broken.)

### 67.3 What this does to the programme's answer

The four-currency reading — existence bought with affinity, gain with entropy, reliability
with molecules, speed with time — survives only in its **first** clause across substrates.
The affinity floor is a real structural feature of both elements and the two values are within
6%. The **gain** clause is AM-specific bookkeeping. Nothing here touches the reliability
clause, which was always about molecule count and is substrate-neutral by construction.

> ⚠ **Rule 19's convention caught a third defect before any cell ran.** P1's gate demanded the
> chain's drift equal the ODE field to 1e−12. **That is false and could never have passed** —
> the chain uses falling factorials n(n−1) where mass action uses x², so the two differ by
> O(1/Ω) *by construction*, and that gap is the discreteness this whole project measures. The
> gate is now that Ω·|drift − field| is Ω-independent, which holds to 3e−3. A gate demanding
> something false is the same defect class as a verdict that cannot fail, and the pre-run test
> cost ten minutes.

---

### 68 The affinity floor is not a law — T-THM-c killed one commit after it was opened

§67 found affinity floors of **3 ln 2** (AM) and **2 ln 3** (Schlögl), noted they are 5.66%
apart and both ln(small integer), and opened T-THM-c as a *pattern* rather than a law. Two
readings fit those two points equally well:

    (a)  A_c = (pairs) × ln(pairs + 1)        3 ln 2  and  2 ln 3
    (b)  A_c = (pairs) × ln(max order)        AM's max order 2, Schlögl's 3

**Both are wrong, and a derivation settles it without any fitting.** Generalise Schlögl's
autocatalysis to order p at fixed pair count — `pX ⇌ (p+1)X`, `∅ ⇌ X` — and impose a triple
root (f = f′ = f″ = 0):

    k1a/k1r = (p+1)x₀/(p−1),   k2r/k1r = (p+1)x₀^p/(p−1),   k2b/k1r = x₀^(p+1)

so that `A = ln[k1a·k2r/(k1r·k2b)]` gives, with **x₀ cancelling exactly**,

> **A_c(p) = 2 ln[(p+1)/(p−1)]**

Confirmed against the engine's `cycle_affinity` — which takes the null space of the per-pair
stoichiometry and has no access to the formula — to **4.4e−16**, and x₀-independent to
**8.9e−16** across x₀ = 0.6, 1.0, 2.5:

| p | 2 | 3 | 4 | 5 | 6 | 8 |
|---|---|---|---|---|---|---|
| **A_c** | 2.1972 | **1.3863** | 1.0217 | 0.8109 | 0.6729 | 0.5026 |
| reading (a) predicts | 2.1972 | 2.1972 | 2.1972 | 2.1972 | 2.1972 | 2.1972 |
| reading (b) predicts | 2.1972 | 2.7726 | 3.2189 | 3.5835 | 3.8918 | 4.3944 |

**p = 3 kills both on the spot** — one counterexample, no tolerance written (rule 19).

> **The floor spans 2.1972 down to 0.5026 within ONE family at FIXED pair count, and tends to
> 0 as p → ∞.** So §67's near-agreement between 3 ln 2 and 2 ln 3 is a coincidence of two
> points that happen to sit 5.66% apart, and **there is no substrate-independent affinity
> floor.** The floor is a property of the specific element, not of restoration.

### 68.1 Two things §67 asserted without checking

> **§67 called A_c a "floor" without verifying it is a minimum.** It is the affinity at which
> bistability dies — a *boundary* value. If A(m) dipped below it, a restoring element could
> run cheaper than "the floor" and the word would be wrong. Swept here: A(m) rises
> monotonically away from m = 0 for every p (e.g. p = 2: 2.1972, 2.1989, 2.2039, 2.2246,
> 2.2611). **It is a genuine minimum, so the word survives** — but it was luck, not evidence.

> **P5 found the pair count never enters, and cannot.** Adding a third reversible pair opens a
> *second independent cycle*, so `cycle_affinity` refuses with a `ValueError` — by design; it
> checks the cycle space is one-dimensional rather than summing an arbitrary combination.
> **So "A_c = (pairs) × something" was malformed as well as wrong**: at more than two pairs a
> single affinity is not defined at all, and both readings were extrapolating a quantity that
> stops existing.

> ⚠ **The P1 gate failed twice, and both times the probe was wrong, not the claim.** First it
> perturbed the *constant* term k2b — but near a triple root f ≈ −k(x−x₀)³ + δ has one real
> root for any δ. Then it perturbed k2r alone, which also leaves the family. §67's actual knob
> adds m²(x−x₀), moving **k2r and k2b together**; with that, every p gives three positive roots
> straddling x₀ and the gate holds. **A gate that fails because the probe is wrong looks
> exactly like a gate that fails because the claim is wrong**, and the only thing that
> distinguishes them is checking the probe against a case whose answer is already known — here
> §67's p = 2, which the first two probes also failed.

### 68.2 What survives of §67

The measurements stand: Schlögl's floor *is* 2 ln 3, AM's *is* 3 ln 2, both exact. What is
withdrawn is the invitation to read them as instances of a law. **T-THM-c is closed as
refuted**, and the honest cross-substrate statement is now:

> Every restoring element examined has an affinity floor — a nonzero drive below which no
> bistable landscape exists — but **its value is set by the element's own stoichiometry and
> ranges over at least a factor of four**. There is no universal price of admission.

---

### 69 The outcome-wise factorisation is deeper than symmetry — T-TUR-f

§66 could not settle whether §60's closure survives asymmetry, because the π(n)/π(n₀) boundary
convention exhausts double precision exactly where the predicted effect is largest. **Rather
than fight that wall in 2-D, ask the question where the instrument is exact.** §67's Schlögl
element has **one dynamic species — no exchange symmetry exists to appeal to** — and is 1-D,
so every quantity is a tridiagonal solve.

**And the tilted generator collapses.** For an up-jump the two channels contribute
`a₁(n)·[a₂(n+1)/a₁(n)] + a₃(n)·[a₄(n+1)/a₃(n)] = μ(n+1)`, and for a down-jump `λ(n−1)`. So the
reverse-weighted generator is **just the chain with λ and μ swapped and shifted** — explicit,
tridiagonal, and free of the enormous factors that broke §66. ln π comes from the exact
birth–death product formula in logs.

Asymmetry is built in by placing the fixed points at arbitrary r₁ < r₂ < r₃, with skew
s = (r₃−r₂)/(r₂−r₁):

| skew | mean \|Φ_lo/p_lo − 1\| | mean \|Φ_hi/p_hi − 1\| | max \|ln w\| |
|---|---|---|---|
| 1.0 | 9.07e−14 | 9.96e−14 | 3.7 |
| 1.7 | 1.32e−14 | 1.82e−13 | 1.6 |
| 2.2 | 3.54e−14 | 5.27e−14 | 9.4 |
| 3.0 | 1.01e−13 | 5.58e−14 | 29.2 |
| 4.0 | 4.56e−14 | 4.81e−14 | 32.3 |

P1 holds on all 17 surviving cells (max |Φ_lo+Φ_hi−1| = 2.0e−13) — **the gate §66 could not
hold**. One cell was excluded on precision, at ln w = 67.

> **Φ_o = p_o to 2e−13 on an element with no symmetry whatsoever, while the two boundaries'
> stationary weights differ by a factor of e³². §60's RESULT generalises. §60's MECHANISM for
> it is wrong.**

§60 attributed the cancellation to the boundaries being exchange images of equal stationary
weight. That account is simply unavailable here — there is one species, the rails are images
of nothing, and ln w runs over tens — yet the factorisation is exact. **So the outcome-wise
identity does not come from symmetry**, and what does produce it is now an open question with
no candidate mechanism attached (rule 17: it gets none until something independent supplies
one).

### 69.1 This resolves §66, against my own prediction

§66 measured |r_e − 1| rising 40× with the tilt and could not separate that from the solve,
because the nuisance grew with the cause. **§69 breaks that confound by construction**: ln w
spans 0 → 32.3 here while |r−1| stays at 10⁻¹³. A boundary-weight ratio of e³² does not
disturb the factorisation when the arithmetic can carry it.

> **So §66's rise was the conditioning, not the tilt** — and §66's own within-band correlation
> of +0.579 was reading solve error. §66 declined to call it either way and recorded
> "suggestive only"; that caution was correct, and the effect it was suggestive *of* is now
> refuted.
>
> **My prediction was wrong in the informative direction.** §66 and §69 both pre-registered
> that the factorisation would fail once symmetry was gone. It does not. §60's closure of the
> founding question's sharpest form therefore **stands, and stands more broadly than §60 could
> claim** — but for a reason nobody has yet given.

### 69.2 Where this leaves the founding question

The reliability–dissipation link is now closed on much firmer ground than in §60. The
fluctuation theorem factorises over outcomes on:

* a symmetric two-species element at every γ (§60),
* an asymmetric one-species element at every skew up to 4:1 (§69),

with the boundary weights differing by up to e³² and the identity holding to 2e−13. **The
error rate is not obtainable from the entropy production by this route, and the obstruction is
not a symmetry accident.**

That closes the last exact route this project had to a thermodynamic price for reliability,
and it agrees with the four-currency reading arrived at independently: reliability is bought
with **molecules**, not free energy. §1, §35 and §38's exp(−cΩ) is the whole story, and §68
has already removed the affinity floor's claim to universality.

> ⚠ **Scope.** Both tests are single restoring elements with two absorbing outcomes and a
> unique separatrix. Nothing here addresses cascades, where §12's depth ceiling lives, nor
> n > 2 outcomes. The identity is an identity; what is shown is that *splitting it by outcome
> yields no information*, on the two substrate families examined.

---

### 70 T-THM-a settled: clause 4 is prior art, and in more generality than §65 claimed

§65 stated the Symmetric Restoration Theorem and recorded its novelty as **UNVERIFIED**,
noting that whether clause 2's decomposition is already in the mirror-symmetry-breaking
literature "has *not* been settled — that needs the Montoya et al. and MATCH
algebraic-analysis papers read properly, not their abstracts." Read properly:

**Montoya, Cruz & Ágreda, *Life* 9:74 (2019), Theorem 1 — the "MM-condition":**

> *"Let s be a racemic steady state of the pseudochiral network Ω, we have that s is
> symmetry-breaking, if and only if, the characteristic polynomial of A_Ωs − B_Ωs is
> unstable."*

with `A_Ω` the k×k block of the Jacobian on the first enantiomer set and `B_Ω` the block
coupling it to the mirror set. **For a Z₂-symmetric Jacobian `[[A, B], [B, A]]` the spectrum
splits into A+B (symmetric) and A−B (antisymmetric), so `A_Ω − B_Ω` IS the antisymmetric
block.** Their Definition 7 makes the same distinction CRNL's P does, requiring an unstable
eigenvector with `v_i ≠ v_{i+k}` — an antisymmetric mode, not just any instability.

> **So clause 4 is prior art, and §65 understated how much.** §65 said it was "in substance
> linear stability of the racemic state, which that literature does routinely". It is more
> specific than that: **it is their Theorem 1**, stated for k enantiomeric pairs where CRNL's
> §62 has k = 1. The generality runs the other way from what §65 implied.

**What that literature does *not* contain**, on the reading of both papers:

* **No per-reaction signed decomposition.** The criterion is Jacobian-eigenvalue-based and
  tested through Hurwitz–Routh inequalities on symbolic determinants. Nothing of the form
  `P = Σ_r c_r·d_r·B_r(x)` with `d_r` integral and `B_r ≥ 0` appears.
* **No purely combinatorial capability result.** No statement of the form "the network can
  break symmetry for *some* rate constants ⟺ [condition on stoichiometry alone]". Their one
  structural theorem is negative and case-by-case (the Calvin model cannot break symmetry,
  proved by ruling out Hurwitz inequalities). Clarke's extreme currents enter only to make
  sampling cheaper, not to decide capability.

### 70.1 The corrected novelty statement

| clause | status |
|---|---|
| **1** invariance / divisibility | **prior art** — flow-invariance of Fix(Z₂), a folk theorem of equivariant dynamics (§65.1) |
| **2** decomposition `P = Σ_r c_r d_r B_r`, `B_r ≥ 0` | **not found** in either paper — the candidate contribution |
| **3** capability ⟺ some `d_r > 0`; non-restoring rates a convex cone | **not found** — follows from clause 2 |
| **4** realisation as one inequality at x\* | **prior art** — Montoya et al. Theorem 1, at k ≥ 1 |
| **5** exact stochastic antisymmetric sector, λ_A | **unchecked** against literature |

**So the theorem's contribution narrows to clauses 2 and 3, exactly as §65 hedged** — not that
a criterion exists, but that it has a fixed sign structure separating topology from rates. That
hedge was written before the papers were read and it survives them.

> **And the agreement is itself worth recording.** Two independent routes — this project's
> §53/§54 decomposition and their Jacobian block analysis — arrive at the same object. Their
> `A_Ω − B_Ω` is the deterministic antisymmetric block; §63's λ_A is its **exact CME
> counterpart**, which is where clause 5 sits and where nothing comparable was found.

> ⚠ **Scope of this literature check, stated because it is not exhaustive.** Two papers read:
> Montoya/Cruz/Ágreda (clear) and the MATCH algebraic-analysis paper (PDF extraction poor, its
> findings weakly supported and not relied on here). The homochirality field is seventy years
> old and this is not a survey. **Clause 2 is "not found in two papers", which is weaker than
> "novel"**, and §65's refusal to claim priority stands until someone reads more widely.

---

### 71 The depth ceiling on an asymmetric element — NOT DECIDED, and the gate says so

Everything in §62–§70 is a single element. **The founding claim is about composition**, and
§12.1 priced that for AM: the ceiling is set by the *inter-stage channel*, not by molecule
count, and at large Ω the per-stage error saturates so that

    D_max ≈ exp(δ*²/2σ²) / 4

AM ran about 3× that (9 against 3.0, 50 against 14.8, 489 against 147).

**The prior here was the opposite of §67's and §68's, deliberately.** Those asked whether
*thermodynamic* quantities transfer, and neither did — the affinity floor is element-specific
(§68) and the cost per e-fold has no counterpart on a driven element (§67). But the depth
ceiling is **information-theoretic**, and its formula mentions only the geometry of the two
rails against the channel noise. Nothing in it refers to bookkeeping, affinity, or how the
element is driven. So this one *should* transfer, and writing the prediction that way is the
test: three prior failures make "does not transfer" the cheap answer and it must not become
the automatic one.

Run on §67/§69's Schlögl element (one species, no symmetry, 1-D so the stage kernel is an exact
matrix exponential), rails at 0.1 and 1.9, channel-then-chemistry as §7 orders it:

| σ/Δ | predicted | Ω=400 | Ω=600 | Ω=900 |
|---|---|---|---|---|
| 0.45 | 2.95 | 4.39 | 4.76 | 5.03 |
| 0.35 | 14.81 | 14.85 | 18.58 | 21.77 |
| 0.28 | 147.12 | 58.52 | 97.84 | 143.20 |

Both gates hold: the kernel is stochastic to 7e−13 and holds its rail (0.99973 after 20
noiseless stages), and with the channel **off** the mutual information stays at 0.9637 through
600 stages — so the decay being measured is the channel's, not the harness's.

> **But P4 fails: D_max keeps rising with Ω** (13.6%, 37.6%, 84.8% spread), and pushed further
> it is still climbing at Ω = 1300 (f = 0.45: 4.39, 4.76, 5.03, 5.22). **§12's ceiling is an
> Ω-*saturated* quantity; this one has not saturated, so it is not yet the same object and the
> ratios are not comparable.**

**Verdict: NOT DECIDED.** The comparison is withheld rather than reported.

### 71.1 What can honestly be said, and what cannot

Cannot: that the depth ceiling transfers. The measurement never reached the regime where §12's
formula is a statement.

Can, and it is worth recording as a *qualitative* contrast: **the ceiling is of the right
order.** Schlögl's measured/predicted runs 0.97–1.70 where AM's is 3.0–3.4 — the same
quantity within a factor of ~2, against a formula containing no free parameter and fitted to
nothing. Set beside the thermodynamic comparisons, where §67's cycles-per-molecule missed AM
by **20×** and had no intensive definition at all, and §68's affinity floor varies by **4×
inside one family**, that is a visibly different kind of disagreement.

> That is a suggestive contrast, **not a result** (rule 17). It is one substrate pair, at
> unsaturated Ω, on a formula whose own AM ratio is 3 rather than 1. The right reading is that
> the information-theoretic quantity is the *most promising* transfer candidate found so far,
> and that saying more requires the saturated regime.

**What it would take.** f = 0.45 is decelerating toward roughly 5.7 (increments 0.37, 0.27,
0.19), so it is within reach at Ω of a few thousand — cap = 2·r₃·Ω states and a dense matrix
exponential, so ~10⁴ states is the practical wall for `expm`. The narrower channels need
much more. **Opened as T-CASC-a with that cost stated**, rather than left as an implied
success.

> ⚠ **Two of my own criteria failed here, both caught.** P2's first version tested
> `d_max(I) is None` and printed **HOLDS on a run whose I was 0.000000** — `d_max` returns
> None both when I never falls through 0.5 and when it started *below* 0.5 and so never
> crossed it, which are the two opposite cases the gate exists to separate. It is now
> `min(I) > 0.5`. And P3 printed **"THE DEPTH CEILING TRANSFERS"** for f = 0.45 in the first
> full run **while P4 had already reported the ceiling still moving with Ω** — a verdict
> announced from a comparison its own precondition had just invalidated. P3 is now gated on
> P4 and prints WITHHELD.

### 71.2 The wall is structural, and the obvious way around it makes things worse

§71 estimated the cost of reaching saturation as "dense `expm` walls out near 10⁴ states". The
obvious fix is to stop forming a dense matrix: the chemistry step only needs `exp(Qt)·v`
(sparse `expm_multiply`) and the channel is a Gaussian convolution rather than a (cap+1)²
matrix. Built as `run_fast`, and **validated against the exact instrument to 3e−6 in I, with
D_max identical to four decimals**.

**It is 50–90× slower.**

| Ω, σ/Δ | max \|ΔI\| | D_max dense vs sparse | dense | sparse |
|---|---|---|---|---|
| 150, 0.45 | 5.2e−6 | 3.0451 / 3.0450 | 2.1 s | 95.5 s |
| 200, 0.35 | 7.4e−8 | 8.5255 / 8.5255 | 1.6 s | 139.4 s |

The reason is structural and settles the cost question properly: **propensities are extensive,
so `‖Qt‖` grows linearly with Ω** — 8.5e4 at Ω = 200, 3.8e5 at Ω = 900, 1.5e6 at Ω = 3600.
Every Krylov or Taylor exponential costs O(‖Qt‖) work, while dense scaling-and-squaring absorbs
the norm in log time but costs O(n³) and O(n²) memory. **Both methods wall out around 10⁴
states, for opposite reasons**, so §71's estimate was right about the wall and wrong about
which method sets it.

> **T-CASC-a therefore needs a different instrument, not a bigger machine.** The stage map is a
> fixed stochastic kernel applied repeatedly, so the natural route is its dominant eigenvalues
> rather than its repeated action — D_max is set by the second eigenvalue of `C·K`, and a
> sparse partial eigensolve is O(nnz) per iteration with no dependence on ‖Qt‖ at all. That is
> the reachable version of this measurement and it is what the open question should have asked
> for.

> ⚠ **Two bugs in the fast instrument, both caught by the gate against the exact one.** The
> channel was normalised globally where the reference normalises **per source row** (the dense
> map truncates the Gaussian asymmetrically at the lattice edges), which disagreed by |ΔI| up
> to **0.30**; and the two steps were applied in the wrong order — `run` computes
> `p @ (C @ K)`, channel then chemistry — which left a half-stage offset worth **~1%** in
> D_max. A faster instrument that disagreed with the exact one would have been worse than no
> instrument, and neither bug is visible in the output of the fast path alone.

---

### 72 The depth ceiling TRANSFERS — the first quantity in this project that does

§71 left T-CASC-a undecided and §71.2 showed the wall was structural: `‖Qt‖ ~ Ω` because
propensities are extensive, so every exponential method — dense or sparse — walls out near 10⁴
states. §71.2 proposed the stage map's second eigenvalue, **which inherits the same wall**,
since applying `C·K` still needs `exp(Qt)v`.

**The right reduction avoids the exponential entirely, and it is exactly the regime in
question.** §12's ceiling describes the limit where each stage *fully restores*. In that limit
the chemistry contributes one number per rail — the probability of committing to the wrong one
— so the per-stage error is

    ε = Σ_n q(n)·p_cross(n)

with `q` the channel's output from a rail and **`p_cross` the exact birth–death splitting
probability**, which §61/§69 already compute in closed form from the scale function, in logs,
with no solve. **O(cap) per cell, no exponential anywhere.** The element is asymmetric, so
ε_hi ≠ ε_lo and the cascade is a binary *asymmetric* channel — handled exactly by a two-state
chain, where `(1−2ε)^D` would be wrong.

| σ/Δ | Ω=400 | 1800 | 7200 | 28800 | **saturated** | predicted | **ratio** | AM's ratio |
|---|---|---|---|---|---|---|---|---|
| 0.45 | 6 | 7 | 7 | 8 | ✓ | 2.95 | **2.71** | 3.00 |
| 0.35 | 30 | 42 | 44 | 45 | ✓ | 14.81 | **3.04** | 3.38 |
| 0.28 | 200 | 447 | 528 | 551 | ✓ | 147.12 | **3.75** | 3.33 |

> ~~**All three channels saturate in Ω, and all three ratios land on AM's — 10%, 10%, 13%.**
> §12's ceiling `D_max ≈ exp(δ*²/2σ²)/4`, which AM overshoots by ≈3×, is overshot by ≈3× on an
> element with one species, no symmetry, chemostats, and a different reaction order.**~~
>
> **DEFLATED by §73.** The measurement stands; the reading does not. A **step function with no
> dynamics at all** gives 2.71/3.11/3.80 on the same geometry, so the agreement is not evidence
> about substrates — the ceiling does not depend on the element beyond where it puts its rails.
> Left as first printed per rule 7.

**This is the first quantity in the project to transfer across substrates**, and it is the one
that was *predicted* to — §71 wrote that prior deliberately, against three prior failures, on
the grounds that the depth ceiling is information-theoretic rather than thermodynamic. Set
against §67 (cost per e-fold: no counterpart at all) and §68 (affinity floor: varies 4× inside
one family), the contrast is now measured rather than suggested.

### 72.1 The caveat that keeps this at "tens of percent"

**The two sides are not perfectly like-for-like, and the mismatch runs one way.** This
measures the t_stage → ∞ idealisation; AM's published 9 / 50 / 489 were measured at *finite*
stage time. P1 quantifies the gap: at Ω = 400, σ/Δ = 0.45 the finite-time cascade plateaus at
**5.10** (t = 6 and t = 18 agree) against the saturated **6–7** — the saturated value is an
upper bound, ~25% high.

So a like-for-like correction would raise AM's effective ratios by roughly that much, moving
the disagreements from 10/10/13% to something nearer 25/25/10%. **The honest claim is that the
ceiling transfers to within tens of percent, not to 10%** — which is still a different
category from a quantity that has no counterpart or varies fourfold.

> ⚠ **Two more criteria of mine, both caught.** The saturation test demanded spread < 5%, but
> D_max is an **integer** depth, so at D_max ≈ 7 one unit of quantisation is 14% and the test
> could never pass there — the f = 0.45 column read 6,7,7,7,7,7,8 and was declared "still
> moving". Now `max(5%, 1.5 units/mean)`. And `p_cross` was written with a per-state
> log-sum-exp loop, O(cap²), which could not reach the Ω this experiment exists to test;
> replaced by a suffix cumulative sum in O(cap), **gated against the O(cap²) version to
> 9.8e−15**.

### 72.2 What the founding claim now says

| you buy | you pay | transfers? |
|---|---|---|
| existence of a landscape | a nonzero drive | **no** — element-specific (§68) |
| gain (e-folds) | entropy per molecule | **no counterpart** off a closed element (§67) |
| reliability | molecules, exp(−cΩ) | yes, by construction |
| **composition depth** | **channel noise against rail separation** | **YES — §72** |

**The thermodynamic accounting of restoration is substrate-specific; the information-theoretic
accounting is not.** That is the sharpest form the founding claim has taken here: what makes a
restoring element good at composing is *how far apart its rails are relative to the noise* —
and that quantity, unlike every price attached to it, is the same across chemistry with
nothing in common.

---

### 73 §72's transfer was never about the element — T-CASC-b

§72 read the saturated depth ceiling landing on AM's ratios as **"the first quantity in this
project to transfer across substrates."** There is a nastier reading it did not test, and §72's
own reduction makes it obvious in hindsight: **the chemistry enters only through `p_cross`**,
the probability of committing to the wrong rail — a monotone function near 1 at the low rail
and near 0 at the high one. If *any* such function with the same geometry gives the same ratio,
the ceiling was never a property of the element.

Holding the rails, the channel and the depth criterion fixed, and swapping **only** the
commitment function:

| commitment function | f=0.45 | f=0.35 | f=0.28 |
|---|---|---|---|
| **(a)** Schlögl, exact birth–death splitting | 2.37 | 3.04 | 3.69 |
| **(b)** Langevin double well, exact — different physics, no counting noise | 2.71 | 3.04 | 3.79 |
| **(c)** bare sigmoid, slope matched at the saddle, **no dynamics** | 2.37 | 3.04 | 3.66 |
| **(d)** step function — **no element at all** | 2.71 | 3.11 | 3.80 |
| AM (published) | 3.00 | 3.38 | 3.33 |

> **A step function reproduces the ceiling.** Commit to whichever rail is nearer, with no
> chemistry, no noise, no dynamics of any kind, and the ratio to `exp(Δ²/2σ²)/4` is
> 2.71/3.11/3.80 — inside the band Schlögl's exact solution gives. The sigmoid matches Schlögl
> to **three significant figures**.

**~~§72: "the depth ceiling TRANSFERS — the first quantity in this project that does"~~ —
DEFLATED.** Left as first printed per rule 7. The agreement between Schlögl and AM is real, and
it is not evidence about substrates: **both agree with a function that has no substrate.**

### 73.1 The corrected claim, which is sharper than the one it replaces

The right reading is not that the ceiling transfers but that **the ceiling does not see the
element**:

> **Composition depth is fixed by the readout geometry — rail separation against channel noise
> — and the restoring element's entire contribution is where it puts its rails.** Reaction
> order, symmetry, drive, chemostatting, the shape of the commitment function: none of it
> survives into D_max once Δ is fixed.

That is a stronger statement than §72's and a more useful one. It says what a restoring element
is *for* in a cascade: it is a rail-placer. Everything this project measured about how the
element gets to its rails — §67's dissipation, §68's affinity floor, §63's threshold blur — is
invisible to the depth it can support.

It also explains the four-currency table's shape rather than adding a row to it. §72 recorded
composition depth as a fifth currency that transfers; it is not a currency at all, because
**nothing is being bought**. `SYNTHESIS.md` is corrected accordingly.

### 73.2 What this does and does not license

**Does not license:** any claim that the chemistry is irrelevant to restoration. The element
determines Δ, whether rails exist at all (§68's floor), how fast it decides (§58), and how
sharply it thresholds at finite Ω (§63). All of that is upstream of the geometry and none of it
is touched here.

**Does license:** dropping the search for a substrate-independent *price*. Three currencies
were tested across substrates and none transferred (§67, §68); the fourth appeared to and
turned out not to be a currency. **The honest summary is that no cost of restoration measured
here is substrate-independent, and the one quantity that is substrate-independent is not a
cost.**

> ⚠ **A broken instrument nearly became the flattering evidence.** Row (b) first read
> 0.68/0.14/0.01 — comfortably "OUTSIDE", which would have *saved* §72 by suggesting the
> commitment function's shape matters after all. It was wrong twice over: the scale density for
> `dX = −U′dt + √(2D)dW` is `exp(U/D)` and I wrote `exp(2U/D)`, and a plain `cumsum` over an
> integrand spanning e⁹⁴⁵⁰ produced 0/0 which `np.clip` silently turned into a plausible
> column. Redone in logs with a suffix log-sum-exp, (b) joins the others at 2.71/3.04/3.79.
> **The bug pointed toward the conclusion I had just published, and a RuntimeWarning was the
> only thing that flagged it.**

---

### 74 What Δ costs: closed elements pay affinity, open ones pay material — T-DEPTH-a

§73 collapsed the founding question to one line. If composition depth is fixed by the readout
geometry and the element's entire contribution is where it puts its rails, then **the cost of
depth, if it exists, lives in Δ**.

**§12's convention had hidden the question.** §12, §71 and §72 all set the channel noise as a
*fraction* of the element's own rails, σ = f·Δ. Then Δ/σ = 1/f identically, D_max depends only
on f, and every element gives the same answer — which is exactly why the predictions matched
across substrates in §72 and why §73's step function reproduced them. **To ask what Δ buys, σ
must be held fixed in physical units.** That change of convention is what this section rests
on, and everything below uses σ = 0.15 in concentration units.

### 74.1 A conservative element has a maximum composition depth

AM conserves X + Y + B = Ω, so concentrations are normalised and δ\* ≤ 1. Sweeping the drive:

| γ | A = −3 ln γ | δ\* | D_max |
|---|---|---|---|
| 0.45 | 2.40 | 0.4005 | 32.7 |
| 0.20 | 4.83 | 0.8165 | 4.75e6 |
| 0.05 | 8.99 | 0.9521 | 1.14e9 |
| 0.002 | 18.64 | 0.9980 | 8.68e9 |
| 1e−4 | 27.63 | 0.9999 | 9.46e9 |

> **δ\* → 1 as γ → 0, so D_max saturates: at σ = 0.15 no amount of drive can push a
> conservative element past D_max = 9.50e9.** Spending 27.6 k_BT of affinity per cycle instead
> of 8.99 buys a factor of 8 in depth; spending infinitely more buys another factor of 1.04.
> **Drive buys rail separation, rail separation is capped by conservation, and depth is capped
> with it.**

### 74.2 An open element buys depth with material, and affinity is free

Schlögl's cycle affinity is `ln[e₁e₂/e₃]`, and under `r → λr` the elementary symmetric
polynomials go as `λe₁`, `λ²e₂`, `λ³e₃` — so **A is exactly invariant** while Δ scales
linearly. Measured over 2.5 decades of λ, the affinity moves by **0.00e+00**:

| λ | rails | A | Δ | D_max |
|---|---|---|---|---|
| 0.25 | 0.125 / 0.375 | 2.3978952728 | 0.125 | 1 |
| 1 | 0.5 / 1.5 | 2.3978952728 | 0.5 | 289.5 |
| 4 | 2 / 6 | 2.3978952728 | 2.0 | **>1e18** |
| 64 | 32 / 96 | 2.3978952728 | 32.0 | **>1e18** |

**At matched affinity — Schlögl's A = 2.3979, i.e. AM at γ = 0.4496 — AM reaches D_max = 33.5
and cannot exceed 9.50e9 at any drive, while Schlögl passes 10¹⁸ by scaling its rails.** The
affinity is identical. What Schlögl spends instead is *material*: its rails sit at r₁Ω and r₃Ω
molecules, so λ = 64 is 64× the molecules for the same thermodynamic force.

> ~~**Closed and open elements pay for depth in different currencies — affinity and material —
> and that is why §67 and §68 found no substrate-independent price. They were pricing the wrong
> thing.** The founding question's cost is not in the dissipation, the gain, or the affinity
> floor; it is in Δ, and what Δ costs depends on whether the element is conserved.~~
>
> **DEFLATED by §75.** True at fixed σ, which this section flagged as load-bearing. With σ
> intrinsic to the chemistry, the rail fluctuation is Poissonian and Δ/σ depends on the
> molecule count λΩ **alone**, to 0.000% — so scaling rails and scaling Ω are the same act,
> there is no separate "material" currency, and AM's ceiling disappears too. One currency:
> molecules at the rail.

### 74.3 Scope, and the convention that decides everything

**This is geometry, not dynamics** — §73 licenses that, since the commitment function's shape
does not enter D_max. Every number above uses a step commitment deliberately, to isolate the
one thing that matters.

**σ fixed in physical units is doing all the work, and it should be argued rather than
assumed.** Under σ = f·Δ the question cannot even be posed. The physical case for fixed σ is
that the inter-stage channel is a property of the wiring, not of the gate — but a real chemical
cascade might well have channel noise scaling with the signal, in which case §12's convention
is the right one and *no* element beats any other. **Which convention is physical is not
settled here, and the entire content of §74 hangs on it.** Stated as the load-bearing
assumption rather than buried.

> ⚠ **Two instrument faults, both caught.** `depth_at` clamped its grid at zero — correct for
> Schlögl's concentrations, **wrong for AM**, whose coordinate is the signed lead δ = x − y
> with its low rail at −δ\*. That cut off the entire low rail and returned `None` in all ten AM
> cells. And `d_max_saturated` iterates O(depth) per evaluation, which does not return at the
> depths here (~10¹⁰); replaced by a closed form, `T^D = π ⊕ λ^D(I − π)`, **gated against the
> iterative version to 1e−12 in I(D) with matching ceilings**, so §72's published numbers are
> untouched.

---

### 75 There is one currency after all: molecules at the rail — T-DEPTH-b

§74 rested entirely on holding the inter-stage noise σ fixed in physical units, and said so.
**In a cascade whose coupling is chemical there is no such freedom**: stage 1's output species
*is* stage 2's input, so the corruption between stages is stage 1's own fluctuation about its
rail, which the CME fixes exactly.

**My prediction for how that noise scales was wrong.** I argued that under Schlögl's rescaling
`r → λr` every propensity scales as λ³, so the quasipotential stretches by λ, its curvature
falls as 1/λ², and σ ∝ λ — which would have made Δ/σ λ-invariant and deflated §74 outright.
Measured: **the exponent is 0.486, not 1.**

The right reading is simpler and is the ordinary one: **the rail fluctuation is Poissonian in
the molecule count.** σ_n ~ √n̄, so in concentration units

    σ_x ~ √(λ/Ω)     measured  d ln σ / d ln λ = +0.486,   d ln σ / d ln Ω = −0.52

and therefore **Δ/σ ~ √(λΩ) — the square root of the number of molecules at the rail.**

### 75.1 The collapse is exact

| λ | Ω | λΩ | Δ/σ |
|---|---|---|---|
| 0.50 | 1600 | 800 | 3.0037 |
| 2.00 | 400 | 800 | 3.0037 |
| 1.00 | 800 | 800 | 3.0037 |
| 0.25 | 12800 | 3200 | 6.8718 |
| 0.50 | 6400 | 3200 | 6.8718 |
| 2.00 | 1600 | 3200 | 6.8718 |
| 8.00 | 400 | 3200 | 6.8718 |

> **Δ/σ depends on the product λΩ alone, to 0.000% spread**, across a 32× range in each. Not a
> fit — an exact scaling symmetry: **scaling an element's rails and increasing its molecule
> count are the same act.**

And AM behaves identically. Its rail width scales as Ω^(−0.5035) and Ω^(−0.5000) at γ = 0.20
and 0.05 — Poisson again — so δ\*/σ = 14.74, 20.92, 29.62 at Ω = 60, 120, 240: **√Ω, without
bound.**

### 75.2 §74 is deflated, and the answer gets simpler

> **~~§74: "closed elements pay affinity, open ones pay material" / "a conservative element has
> a maximum composition depth of 9.5e9"~~ — DEFLATED.** Both statements were true *at fixed
> σ*, which §74 flagged as load-bearing and untested. With σ intrinsic:
>
> * AM's ceiling disappears — δ\*/σ ~ √Ω grows without bound, so the cap was a fact about the
>   convention, **not about conservation**;
> * Schlögl's "material" currency is not a separate currency — scaling rails by λ is
>   *identical* to scaling Ω by λ.

Left as first printed per rule 7. What replaces it is shorter:

> **Depth is bought with molecules at the rail, in both substrates, and with nothing else.**
> `D_max ~ exp(Δ²/2σ²)` with `Δ/σ ~ √N_rail`, so `D_max ~ exp(c·N_rail)`.

That is the *same currency and the same functional form* as reliability — §1, §35 and §38's
exp(−cΩ) error probability. **The founding question's two halves collapse into one statement:
a restoring element buys both its error rate and its composition depth with molecules,
exponentially, and buys neither with free energy.**

### 75.3 What this does to the affinity

Affinity has not become irrelevant — it has become a *gate* rather than a *price*. Drive is
required for rails to exist at all (§9.1's floor, §68's `2 ln[(p+1)/(p−1)]`), and it sets where
the rails sit. But **once the rails exist, buying more depth is a matter of counting, not of
dissipation.** §74's own table already showed the shape of this without my seeing it: going
from A = 8.99 to A = 27.6 bought a factor of 8 in depth, and going to infinity bought 1.04.

> ⚠ **Three criteria, and the pattern is now familiar.** P1 first demanded |exact/LNA − 1| <
> 0.10 at every cell and failed on the coarsest one while the series ran 0.1538 → 0.0200 →
> 0.0044 — a fixed tolerance on a converging quantity, **the same error as §63's P1(c)**.
> Replaced by a convergence test; that then failed too, because Ω = 400 sits outside the
> asymptotic regime the LNA is about, so the tail is asserted and the exclusion is printed. And
> `ln_pi` used an interpreted O(cap) loop that did not return at cap ≈ 10⁶; vectorised.

---

### 76 Depth *is* the error rate: `D_max · ε = c*` — T-DEPTH-c

§75 found depth and reliability bought in the same currency with the same functional form. That
invites a sharper question, and it can be answered before measuring anything. For a binary
symmetric channel the D-step bias is (1−2ε)^D and the mutual information falls through ½ at
b\*, defined by H((1+b\*)/2) = ½. So

> **D_max = ln b\* / ln(1−2ε) → c\*/ε,   c\* = −ln(b\*)/2 = 0.124266404564**

**D_max · ε is a pure number.** Verified against the closed-form chain:

| ε | 1e−1 | 1e−2 | 1e−3 | 1e−4 | 1e−6 | 1e−8 |
|---|---|---|---|---|---|---|
| \|D·ε/c\*−1\| | 1.0e−1 | 1.0e−2 | 1.0e−3 | 1.0e−4 | 1.0e−6 | **8.9e−9** |

And it is the **same constant on every substrate** — Schlögl at λ = 1 and 4, Ω = 3600 and
14400, and §73's step function at two channel widths — all giving D·ε/c\* within **0.8%**, the
residual being the asymmetry correction below rather than anything about the substrate.

> **So §75's unification is an identity, not a coincidence.** The depth ceiling contains *no
> information beyond the per-stage error rate*. §12's entire depth apparatus — cascade kernels,
> mutual-information decay, the ceiling formula — reduces to one number the single element
> already determines.

**The asymmetric case.** Real elements have ε_hi ≠ ε_lo. The decay rate is λ = 1 − ε_hi − ε_lo,
which depends on the sum alone, so the **arithmetic** mean is the right ε — confirmed against
geometric and harmonic. It is exact only in the symmetric limit, and the deviation grows with
the asymmetry ratio: 3× → 1.050, 9× → 1.165, 20× → 1.247, 1000× → 1.385. Every element measured
in this project is mildly asymmetric and sits inside 1%.

### 76.1 What §12's "factor of 3" actually was

§12.1 measured its ceiling running ≈3× the prediction and concluded: *"The ceiling's 'factor
≈3' is not a prefactor — it is a 7% error in the exponent."* **That reading is wrong.**

`exp(Δ²/2σ²)/4` is the Gaussian tail with its algebraic prefactor discarded. The exact ceiling
is c\*/Φ(−Δ/σ), and since Φ(−z) ≈ exp(−z²/2)/(z√2π),

    exact / §12's formula  ≈  4c*√(2π) · (Δ/σ)

— **an algebraic prefactor that grows linearly in Δ/σ, not a constant.**

| σ/Δ | Δ/σ | §12's formula | exact c\*/ε | ratio | predicted 4c\*√2π·(Δ/σ) |
|---|---|---|---|---|---|
| 0.45 | 2.222 | 2.95 | 9.46 | **3.204** | 2.769 |
| 0.35 | 2.857 | 14.81 | 58.14 | **3.926** | 3.560 |
| 0.28 | 3.571 | 147.12 | 700.01 | **4.758** | 4.450 |

The ratio rises exactly as the dropped prefactor requires, and §73's clean step-function
measurements rise with it (2.71, 3.11, 3.80). §12.1's own AM numbers (3.00, 3.38, 3.33) do not
rise monotonically, but they were measured at finite Ω in a different setup and were read as a
*constant* 3 — which is the reading this section overturns.

> **~~§12.1: "the factor ≈3 is not a prefactor — it is a 7% error in the exponent"~~ —
> WITHDRAWN.** It is precisely a prefactor: the Laplace factor that `exp(Δ²/2σ²)/4` drops. The
> exponent was never wrong. Left as first printed per rule 7.

### 76.2 What the founding claim reduces to

Chaining §73 → §76: composition depth does not see the element's dynamics (§73), the element
enters only through Δ (§73), Δ/σ is fixed by molecules at the rail (§75), and depth is exactly
c\*/ε (§76). So the whole cascade apparatus collapses to a single element quantity:

> **A restoring element is characterised, for every purpose this project set out to measure, by
> one number: its per-stage error probability ε.** Reliability is ε. Depth is c\*/ε. There is
> no second quantity, no trade-off between them, and no thermodynamic price on either — ε is
> bought with molecules, exponentially.

That is the founding question answered in one line, and it is a *deflationary* answer: the
depth advantage of a restoring switch is not an additional property beyond its error rate but
the same property counted twice.

> ⚠ **Two more gates too strict, and the class is now stable enough to name.** P1 demanded
> |D·ε/c\*−1| < 1e−3 for every ε ≤ 1e−3 and rejected ε = 1e−3 at *exactly* 1e−3; P3 demanded
> 0.1% across substrates when the physical residual is the 0.8% asymmetry correction. **Both
> are fixed tolerances applied to converging quantities — the same error as §63's P1(c) and
> §75's P1, now three sections running.** The fix each time is to test convergence rather than
> a level, and the reason it keeps recurring is that writing a tolerance is easier than asking
> what the residual is made of.

---

### 77 The last free number: nats per molecule, and it does not transfer — T-DEPTH-d

§76 reduced everything this project set out to measure to one element quantity: the per-stage
error ε. Reliability is ε, depth is c\*/ε, there is no third thing. **So the whole founding
question rests on one coefficient** — how fast ln(1/ε) grows with molecule count:

    η  =  d ln(1/ε) / dΩ        [nats of reliability per molecule]

**This is the transfer test that matters, and the earlier ones did not.** §67 priced
dissipation per e-fold and §68 the affinity floor; §73 later showed both are invisible to
composition. η is the only quantity left.

**It exists.** ln(1/ε) becomes linear in Ω on both substrates, gated on convergence rather
than a tolerance (rule 20):

| element | landscape | ln(1/ε) at successive Ω | **η** | last-step drift |
|---|---|---|---|---|
| AM | γ = 0.20 | 112.2, 222.7, 443.0, 736.6 | **1.8346** | 0.09% |
| AM | γ = 0.05 | 597.2, 1190.3, 2376.2, 3957.2 | **9.8813** | 0.01% |
| Schlögl | λ = 1 | 24.7, 97.2, 383.7, 1527.6 | **0.059574** | 0.20% |
| Schlögl | λ = 4 | 97.2, 383.7, 1527.6, 6100.7 | **0.238185** | 0.05% |

**P3's prediction held.** Schlögl's ratios looked superlinear in §75's data (2.288 against
√4 = 2), and I predicted that was small-Ω contamination of the rail width rather than different
physics, *because §75's exact λΩ collapse cannot coexist with two exponents*. It straightens:
η settles to 0.20% by Ω = 25600. And λ = 4 gives **0.238185 / 0.059574 = 3.998** — exactly the
λ-scaling §75's collapse requires, so the two sections agree where they must.

### 77.1 η does not transfer, and that is the answer

| element | landscape | η |
|---|---|---|
| AM | γ = 0.30 | 0.6813 |
| AM | γ = 0.20 | 1.8346 |
| AM | γ = 0.05 | **9.8813** |
| Schlögl | spread 0.6 | **0.015617** |
| Schlögl | spread 0.9 (λ=1) | 0.059574 |
| Schlögl | λ = 4 | 0.238185 |

> **η spans a factor of 633**, and it moves with the landscape *inside* each substrate — γ from
> 0.30 to 0.05 buys 14×, and Schlögl's root spacing from 0.6 to 0.9 buys 3.8×. **It does not
> transfer, and it was never going to: η is the one place the chemistry survives.**

So the programme's closing arithmetic is:

> Everything about a restoring element washes out of composition — its reaction orders, its
> symmetry, its drive, its dissipation, the shape of its commitment function, even whether it
> is conserved — **except one number, η, and that number is the element's own.** A restoring
> element is a device for converting molecules into nats at an exchange rate set by its
> landscape, and *nothing else about it matters* to what it can compute.

That is the founding question answered, and the answer is deflationary in the same direction
§73 and §76 pushed: the transistor's advantage is not a special thermodynamic property. It is a
good exchange rate.

### 77.2 What is unresolved

**η's dependence on the landscape is measured, not derived.** γ → η runs 0.6813, 1.8346, 9.8813
at γ = 0.30, 0.20, 0.05 — plainly not linear, and WKB should give it in closed form as the
barrier action per molecule. That is the natural next quantity and this section does not claim
it. **T-DEPTH-e.**

**And the definition of ε is a convention.** It is the Gaussian readout of the element's own
intrinsic noise, which is what §75 established the chemically-coupled cascade actually uses —
but a different readout would give a different η. The *ratios* above would survive; the
absolute numbers might not.

> ⚠ **A bare `except` was hiding an exclusion.** The first version wrapped the AM loop in
> `except Exception: continue`, and one cell (γ = 0.30, Ω = 240) was silently dropped when the
> engine's own guard refused an untrustworthy stationary solve. The guard was working; my
> handler discarded its message. Failures are now counted and printed — **1 cell excluded** —
> which is what rule 10 asks for and what a bare `except` is structurally incapable of doing.

---

### 78 η is derived, and the founding question lands in the ODE — T-DEPTH-e

§77 left η undetermined and proposed deriving it from WKB as the barrier action per molecule.
**That named the wrong theory.** WKB gives the *escape* probability — the chance the element
spontaneously crosses its own saddle — but §75/§77's ε is not an escape. It is the **Gaussian
readout of the rail's own fluctuation**, ε = Φ(−Δ/σ), which is what a chemically-coupled
cascade actually applies. The relevant object is therefore the **linear-noise variance**:

    σ² = V/Ω  (Lyapunov at the fixed point)  ⟹  ln(1/ε) ≈ Δ²Ω/2V  ⟹  **η = Δ²/(2V)**

**The LNA is the exact rail width in the large-Ω limit**, converging on every element (rule 20
— convergence, not a tolerance): AM at γ=0.20 gives |ratio−1| running 0.0062 → 0.0028 → 0.0013
→ 0.0008; at γ=0.05 it is already 1.0000 to four places; Schlögl runs 0.0084 → 0.0020 → 0.0005.

Against §77's stored numbers, with nothing fitted:

| element | Δ | V | η predicted | η measured | ratio |
|---|---|---|---|---|---|
| AM γ=0.05 | 0.95213 | 0.045880 | 9.8797 | 9.8813 | **0.9998** |
| AM γ=0.20 | 0.81650 | 0.181858 | 1.8329 | 1.8346 | 0.9991 |
| AM γ=0.30 | 0.70741 | 0.369014 | 0.6781 | 0.6813 | 0.9952 |
| Schlögl s=0.6 | 0.60000 | 11.555556 | 0.015577 | 0.015617 | 0.9974 |
| Schlögl s=0.9 | 0.90000 | 6.802469 | 0.059537 | 0.059574 | 0.9994 |

**Worst 0.48%.** P4's prediction held: the error is largest at γ = 0.30, nearest γ_c where the
rail is shallowest and the LNA's harmonic assumption is weakest, and vanishes at γ = 0.05.
P5's holds too — η(λ=4)/η(λ=1) = **4.0000**, exactly what §75's λΩ collapse requires.

### 78.1 How much of this is content and how much is algebra

**Stated plainly, because P2 is mostly P1 restated.** η_measured was extracted from ε = Φ(−Δ/σ)
with σ from the exact CME, so η_measured ≈ Δ²/(2V_exact) up to the algebraic corrections in Φ.
The test `η_pred/η_meas ≈ 1` is therefore, at bottom, **V_LNA = V_exact** — which is P1 — plus
algebra. What P2 adds is that the Φ prefactor corrections, which enter the *slope* in Ω, are
small enough not to disturb it at the 0.5% level.

The content is P1's, and it is worth having: **the linear-noise variance is the whole story for
this ε.**

### 78.2 The founding question is in the ODE

`η = Δ²/(2V)` needs the rail separation, the Jacobian, and the diffusion matrix `S·diag(a)·Sᵀ`
at a fixed point. **All three are deterministic-plus-LNA quantities.** No master equation, no
stationary solve, no simulation, no entropy production.

> Chaining §73 → §78: depth does not see the element's dynamics; the element enters only
> through Δ; Δ/σ is molecules at the rail; depth is exactly c\*/ε; ε is fixed by η; and
> **η = Δ²/2V is computable from the deterministic field and its linearisation.**
>
> **The entire founding question — how deeply can you compose a noisy restoring element, and
> what does it cost — is answered by the ODE and its linear-noise correction.** The exact CME
> was needed to *establish* that, and is not needed to *use* it.

That is the strongest form of the deflationary answer this programme has been converging on
since §73. It also explains, retrospectively, why every thermodynamic price failed to transfer:
entropy production is not in the formula. Affinity enters only by setting where the rails sit
and how deep the wells are — that is, through Δ and V — and never appears in its own right.

> ⚠ **§77's kill test named the wrong theory, and I wrote it.** "WKB should give η in closed
> form as the barrier action per molecule" is a well-formed proposal about a *different* ε than
> the one §77 measured. Both quantities are real and both are exponential in Ω; they are simply
> not the same exponent. The lesson is narrow and worth keeping: **when a section reduces to
> "one number", check which definition of that number the reduction actually used before
> proposing how to derive it.**

---

### 79 Out of sample: the ODE predicts systems it never saw — and §78 needs one qualification

§78 was checked on exactly the systems it was built from. **A formula that fits the data it was
derived on is not a formula that works** — rule 16 exists here because §22 fitted a convolution
for three subsections and was out by 3688× against an exactly-computable quantity. So: predict
first, from the ODE and the Lyapunov equation alone, then check against the exact CME.

| system | never used before? | \|σ_pred/σ_exact − 1\| across Ω |
|---|---|---|
| AM γ = 0.10 | γ never used | 0.04% → 0.02% → **0.01%** |
| AM γ = 0.35 | γ never used | 13.96% → 5.04% → 1.98% → **1.12%** |
| **Schlögl QUARTIC** (3X⇌4X) | **different reaction order** | 2.63% → 0.58% → **0.14%** |
| Schlögl asymmetric rails | unequal basins | 0.70% → 0.17% → **0.04%** |
| *AM γ = 0.45* | *predicted to fail* | *10.76% → 8.88% → 23.72% → 23.33%* |

**All three out-of-sample systems converge.** The quartic is the one that matters most: it
changes the *reaction order*, not a parameter inside a family — a quartic deterministic field,
`3X ⇌ 4X`, never touched by §73–§78 — and the ODE predicts its rail width to 0.14%.

**And the case predicted to fail, failed.** AM at γ = 0.45 does not converge at all
(0.8924 → 1.0888 → 1.2372 → 1.2333), landing **20.9× worse** than the worst passing case.
That is §78's P4 diagnosis confirmed on a system chosen in advance to break it: the rail is
shallowest near γ_c, so the LNA's harmonic assumption fails there and nowhere else. **Including
a case the theory should get wrong is what makes the other three mean anything.**

### 79.1 The qualification §78 needs

`D_max ~ exp(Δ²Ω/2V)`, so a relative error δ in V becomes a factor `D_max^δ` in the depth. **A
1% error in σ is not a 1% error in D_max.**

| system | ln D predicted | ln D exact | ratio | σ error |
|---|---|---|---|---|
| AM γ = 0.10 | 1135.55 | 1135.32 | **1.26×** | 0.01% |
| Schlögl asymmetric | 1226.23 | 1225.20 | **2.8×** | 0.04% |
| Schlögl quartic | 918.10 | 915.52 | **13.2×** | 0.14% |
| AM γ = 0.35 | 143.53 | 140.40 | **22.8×** | 1.12% |

> **§78's "the exact CME is not needed to use it" holds for η and not for D_max.** The ODE route
> predicts the *exponent* to well under a percent on systems it never saw — that is the real
> claim, and it survives — but the depth itself only to a factor of 1.3–23×, because depth is
> exponentially sensitive to V. Anyone wanting D_max to better than an order of magnitude needs
> the master equation after all.

That is a genuine limit rather than a failure. `ln D_max` is predicted to 0.02–2%; `D_max` is
not. The distinction matters because §72–§77 reported depths as numbers, and those numbers
carry the exponential amplification with them.

### 79.2 Rule 20 did not stop me writing the same gate

> ⚠ **I added rule 20 two sections ago — never gate a converging quantity with a fixed
> tolerance — and then wrote exactly that gate here.** P1's first version demanded
> |ratio − 1| < 1% at the largest Ω and called AM γ = 0.35 a **MISS at 1.12%** while its series
> ran 1.1396 → 1.0504 → 1.0198 → 1.0112. Replaced by a convergence test, under which it passes
> and γ = 0.45 correctly does not. **Writing the rule down did not stop me writing the gate**,
> which is worth recording because it says something about how these rules actually work: they
> catch the error on re-reading, not on first drafting.

> ⚠ **And one instrument failure that looked like a physics failure.** The quartic system first
> read 0.966 → 1.253 → 1.047, non-monotone, which would have been evidence that the formula
> breaks under a change of reaction order. It was my landscape: at m = 0.35 the rails sit only
> Δ/σ ≈ 2.4 apart at Ω = 1600, so the basin-restricted second moment was contaminated by the
> *other* rail's tail. Separating them (m = 0.8, roots 0.298/1.0/1.505) gives the clean
> 2.63% → 0.58% → 0.14%. **Non-monotone convergence is a signature worth trusting — it meant the
> instrument, not the theory.**

---

### 80 §75–§79 priced the wrong failure mode — T-DEPTH-g

§75 argued that a chemically-coupled cascade has no abstract channel: stage 1's output species
*is* stage 2's input, so the inter-stage noise is stage 1's own rail fluctuation. §76–§79 built
everything on the resulting per-stage error, the **Gaussian readout** of that fluctuation.

**An element in a cascade has a second way to fail, and the arc never compared them.** During
the stage time it can spontaneously escape its rail. Both are exponential in Ω:

    ε_read ~ exp(−ηΩ),   η = Δ²/(2V)        [LNA, §78]
    ε_esc  ~ t·exp(−AΩ), A = the escape action [quasipotential]

**Whichever exponent is smaller dominates absolutely.** Measured on §75's own element, with the
escape action from the exact mean first-passage time (ln T linear in Ω, local A settling
0.025923 → 0.026046):

| Ω | ln ε_read | ln ε_esc | difference | binds |
|---|---|---|---|---|
| 400 | −24.73 | −10.53 | 14.20 | **ESCAPE** |
| 1600 | −97.22 | −41.70 | 55.52 | **ESCAPE** |
| 6400 | −383.75 | −166.70 | 217.05 | **ESCAPE** |
| 12800 | −765.14 | −333.40 | **431.74** | **ESCAPE** |

> **A = 0.026046 against η = 0.059537** — the escape exponent is 2.29× smaller, so at Ω = 12800
> spontaneous escape is **e⁴³² ≈ 3×10¹⁸⁷ times more likely** than misreading. Readout would bind
> only for stage times below **6×10⁻¹⁸⁸**.
>
> **§75's own premise — that the physical cascade is chemically coupled — selects the regime in
> which §75–§79's ε is not the physical one.**

### 80.1 What survives, and one thing that survives in better shape than expected

**§76 stands as mathematics.** `D_max = c\*/ε` is a statement about a binary channel with
per-stage error ε, verified to 8.9e−9. It does not care where ε comes from. What changes is
which ε to put in it.

**§77's η is the wrong coefficient for this regime.** "The last free number" is the right *kind*
of statement — one exponent governs everything — but for a chemically-coupled cascade the number
is A, not η. §77's measurement of η is correct and its transfer conclusion is unaffected (A is
just as substrate-specific), but η is not what sets the depth.

**§78's headline survives, with a different formula.** Its claim was that the founding question
is answered by deterministic-side quantities with no master equation. That is still true — the
escape action is an integral of the *propensity densities*:

    A = −∫ ln(μ(x)/λ(x)) dx    from the saddle to the rail

Computed that way, from the rate functions alone with no chain and no lattice: **0.026047,
against the exact MFPT's 0.026046 — agreement to 1×10⁻⁴.** So the answer is still on the
deterministic side; it is the quasipotential integral rather than the linear-noise variance.

> **Scope on that.** The integral form is a 1-D fact. For a multi-dimensional element like AM
> the quasipotential requires solving a Hamilton–Jacobi equation rather than an integral —
> harder than the LNA, still deterministic-side, and not done here.

**§71/§72 are unaffected, and this was checked rather than assumed.** They used an *external*
channel with σ = f·Δ, 4.8× wider than the intrinsic width at Ω = 1600. There ln ε_read = −6.15
against escape at −41.70, so the readout term binds by 35 nats. **Engineered wiring and chemical
coupling are physically different cascades, and the readout analysis is right for the first.**

### 80.2 Why I missed it for five sections

The arc from §73 was a chain of reductions, each correct given the last, and **the error entered
at the one step that changed the physical setup rather than simplifying the mathematics.** §75
replaced an external channel with the intrinsic rail width — a change of *model*, not of
algebra — and that is precisely the step at which a second failure mode became available.
Every later section then inherited the assumption without restating it.

> The generalisable form: **when a reduction replaces one physical mechanism with another,
> re-enumerate the failure modes.** Simplifying algebra cannot introduce new ones; changing what
> the model represents can.

> ⚠ **The first escape instrument returned a negative time.** A banded MFPT solve gave
> −4.98×10¹³ at Ω = 3200 and a value that did not grow with Ω at all — the reflecting-boundary
> row was written at the wrong lattice site, since the band spans lo+1…hi−1 and the reflecting
> site is hi. A mean first-passage time cannot be negative, and nothing was read off it.
> Replaced by the exact log-space sum, which gives ln T linear in Ω to 0.01%.

---

### 81 §77–§79 redone against the escape action — T-DEPTH-h

§80 showed the coefficient governing depth in a chemically-coupled cascade is the escape action
A, not the linear-noise η. Three sections were built on η, and each needed redoing.

**§79 redone (out of sample), and it is the strongest form yet.** A is an integral of the
propensity densities, so it can be *predicted* with no chain, no lattice and no master equation
— and checked against the exact mean first-passage time:

| landscape | A from ∫ | A from exact MFPT | ratio | used in §80? |
|---|---|---|---|---|
| 0.1 / 1.9 | 0.026047 | 0.026044 | 1.0001 | yes |
| 0.4 / 2.2 | 0.034764 | 0.034761 | 1.0001 | **no** |
| 0.3 / 1.7 | 0.010138 | 0.010132 | 1.0006 | **no** |
| 0.05 / 2.5 | 0.093421 | 0.093419 | 1.0000 | **no** |
| quartic 3X⇌4X | 0.009985 | *(1-D integral)* | — | **no** |

**Worst 0.06%, on four landscapes never used**, including a different reaction order. So §78's
"deterministic-side, no master equation" claim holds for A at least as well as it did for η —
and A is the exponent that actually governs the depth.

**§77 redone: A does not transfer either.** It spans 0.0101 to 0.1358, a factor of 13, moving
with landscape shape inside each substrate. The same conclusion as η's, reached independently.

**And §80 is not a relabelling.** A/η spans **0.2002 to 0.7560 — a factor of 3.78** across
shapes and substrates. The two coefficients respond differently to the landscape, exactly as
predicted: η is a *local harmonic* quantity at the rail, A is a *global integral* over the
barrier. Had A/η been constant, §77's conclusions would have carried over unchanged and §80
would have been a change of name.

**§72's depths restated.** With escape dominant, `ln D_max = AΩ − ln(t/c*)`:

| element | Ω | ln D (escape) | ln D (readout) |
|---|---|---|---|
| Schlögl 0.1/1.9 | 6400 | **163.9** | 381.0 |
| Schlögl 0.4/2.2 | 6400 | **219.7** | 305.9 |

The escape route gives a far shallower ceiling. §72's published depths were computed with an
external channel and are a different regime (§80 P5), not these numbers.

### 81.1 An underflow was driving both headline spreads

> ⚠ **AM at γ = 0.20 first returned A = 5.7137, 3.8145, 2.8591 at Ω = 120, 180, 240.** Those
> products with Ω are 685, 686, 686 — the barrier was not growing at all. `π(saddle)` had reached
> the 1e−300 floor, so `ln π(rail) − ln π(saddle)` saturated at ln(10³⁰⁰) ≈ 691 and "A" became
> exactly 691/Ω. **That single cell was inflating A's reported span from a factor of 13 to a
> factor of 282, and A/η's from 3.78 to 7.79.**
>
> The tell was that the value scaled as 1/Ω — a quasipotential barrier per molecule must
> *converge*, not decay. Cells at the floor are now excluded rather than fitted, and γ = 0.20 is
> reported as having no usable stationary solve.

Both conclusions survive the exclusion, with smaller numbers. **That is the point of catching
it: the verdicts were right and the magnitudes were not**, and a factor of 282 quoted from an
underflow would have been the most quotable wrong number in the record.

---

### 82 The escape action has no thermodynamic price — T-COST-n

The founding question was pursued four ways and closed each time (§3 of `SYNTHESIS.md`). **Every
one of those closures priced a quantity §80 has since shown governs nothing.** §80/§81
established that the exponent setting reliability, and therefore depth, is the escape action A.
A had never been asked the thermodynamic question.

**The test is exact, because the drive of these networks is one number.** AM's reversible pairs
span a one-dimensional cycle space, so its entire non-equilibrium force is the cycle affinity
(§16, derived generically by `cycle_affinity`, not assumed); Schlögl's two pairs likewise span
one cycle, with affinity `ln(k1a·k2r / k1r·k2b)`. **So the force can be pinned to a constant by
construction while the kinetics move** — rule 9's opposite sweep in its strongest form.

**Why the saddle is pinned at r₂ = 1.** The Schlögl affinity `ln(e₁e₂/e₃)` is invariant under
rescaling all three roots, while A scales *linearly* with concentration (§75, §78). A uniform
rescaling therefore moves A at fixed affinity for free — and would be correctly dismissed as a
change of units. Pinning r₂ = 1 removes that freedom.

**P2, 1-D and exact to quadrature (errors 1e−17 to 1e−10):**

| r₁ | r₃ | affinity | A | Δ |
|---|---|---|---|---|
| 0.08 | 1.3480 | 3.543245 | 0.002240 | 0.634 |
| 0.10 | 1.9000 | 3.543245 | 0.026047 | 0.900 |
| 0.20 | 4.3505 | 3.543245 | 0.439982 | 2.075 |
| 0.50 | 9.6411 | 3.543245 | 2.073935 | 4.571 |

The affinity varies by **4.9e−14 nats** across the family. **A spans a factor of 926.**

**P3, a second substrate.** AM with per-reaction reverse ratios (γ₁, γ₂, γ₂) on the level set
γ₁γ₂² = const — which holds the affinity fixed to 4.4e−16 nats *and* keeps X↔Y exchange symmetry
exact, so the saddle stays on the diagonal. A spans **0.01770 to 0.08808, a factor of 5.0**, with
local `d ln T/dΩ` converging to four figures in every retained cell. One cell (γ₂ = 0.25) failed
the convergence gate and was excluded, not fitted.

**P4, and the dissociation runs both ways.** Elements tuned to *identical* A = 0.026047 (to
0.0000%) carry affinities from **3.226 to 4.143 nats**. So A does not bound the drive either.

> **Neither implication holds. The thermodynamic force does not set the escape action, and the
> escape action does not report the force.** This is the fifth closure of the founding question,
> and the first one aimed at the exponent that actually governs anything.

### 82.1 P5's verdict was right for a range and a half, and then wrong

P5 asked what *does* set A, predicting that Δ — one scalar against a whole field — would not be
enough. **Swept over r₁ ∈ [0.05, 0.35], A spanned 3% and the verdict printed "Δ alone reproduces
A".** That sweep moves the saddle's relative position by only 1.4×. Widened to the full
admissible range 0.02 < r₁ < 0.97 — a **32.7×** move — at three separate Δ:

| Δ | saddle position moved | A spans | affinity spans |
|---|---|---|---|
| 0.9 | 32.7× | **1.754×** | 2.558 nats |
| 1.5164 | 32.7× | 1.189× | 2.731 nats |
| 3.0 | 32.7× | 1.426× | 2.919 nats |

> **This is rule 9 caught in flight, in the section that invokes rule 9 as its own method.** The
> narrow sweep was not wrong about its own data — A really does vary by 3% there. It was wrong
> that the range it covered was the range that mattered, and nothing in the output said so.

Both numbers are kept, and both are in the test. **The corrected reading is sharper than either
prediction: Δ and the affinity both fail to determine A, by margins three orders of magnitude
apart — 1.75× against 926×.** The drive has essentially no grip on the escape action; the
landscape's headline number has most of one. The suspect (rule 17) remains *"A is a functional
of the whole field, not of any single number extracted from it"*, and it is not yet confirmed —
what would confirm it is a second field with identical Δ **and** identical saddle placement
giving a different A.

---

### 83 Two networks with identical mass-action ODEs have different reliability — T-COST-o

§82 left one suspect standing: A is a functional of the whole field and of no scalar summary of
it. The kill test as opened was to match Δ and the saddle placement and vary the curvature.
**A sharper version was available, and it kills more than the suspect.**

`A = −∫ln(μ/λ)dx` depends on λ and μ **separately**. The mass-action ODE depends only on their
*difference*. So any reaction pair adding the same function to both leaves the deterministic
dynamics exactly unchanged while moving A — and mass action supplies such pairs directly, since
a birth and a death with identical **reactant complexes** have identical propensities:

| substrate | the neutral pair | shared propensity |
|---|---|---|
| Schlögl | `X→2X` and `X→∅` at equal rate c | `c·n` |
| AM | `X+Y→2X` and `X+Y→2Y` at equal rate c | `c·n_X·n_Y/Ω` |

**P1(a), the gate the section means nothing without.** Schlögl: `max|f_c(x) − f_0(x)| = 1.07e−14`
over 41 points × 6 values of c, against a field scale of 2.16 — a relative perturbation of
**4.9e−15**. AM: `1.11e−16` over 200 random simplex points × 3 values of c. **The deterministic
dynamics is identical on both substrates**, every fixed point, Δ, and the saddle included.

**P1(b), rule 16, absolute.** The quadrature is checked against the *exact* first-passage time of
the modified chain at every c — worst **0.215%**. Without this the whole result could have been a
quadrature drifting with c while the true barrier did not.

**P2/P3/P4 — the measurement, at a fixed deterministic field:**

| c | A | A/A(0) | η | η/η(0) | ln D_max at Ω=6400 |
|---|---|---|---|---|---|
| 0.0 | 0.026047 | 1.0000 | 0.059537 | 1.0000 | **163.9** |
| 1.0 | 0.021108 | 0.8104 | 0.050782 | 0.8529 | 132.3 |
| 5.0 | 0.012044 | 0.4624 | 0.031974 | 0.5370 | 74.3 |
| 20.0 | 0.004629 | 0.1777 | 0.013384 | 0.2248 | **26.8** |

Both exponents fall, monotonically and in the same direction — A by 5.6×, η by 4.4×. They are
computed by completely different routes (quadrature over the barrier; Lyapunov at the rail), and
a disagreement in *direction* would have meant one of them was wrong.

> **At Ω = 6400 the depth ceiling spans e¹³⁷ between elements that no ODE measurement can tell
> apart.**

**P5, rule 9, on a conservative two-species element.** AM's neutral pair collapses it the same
way: at γ = 0.35, A falls 0.07158 → 0.02034, a factor of **3.52**, monotone, with local
`d ln T/dΩ` converging to four figures in every retained cell; at γ = 0.40, a factor of 2.27 with
the c = 1.0 cell excluded on the convergence gate. AM's c = 0 value reproduces §82's 0.07158
exactly, which pins the two sections' instruments together.

**What this costs and what it does not.** §78/§81's operational claim survives untouched — A is
still obtained by quadrature with no master equation. **Its gloss did not: "a property of the
deterministic field" was too loose.** The Lyapunov route needs `D = S·diag(a)·Sᵀ`, which is the
propensities and not the drift, and so does A. The correct statement is that A is a functional of
the **propensity pair**, and the ODE is a strictly coarser object than the thing that sets
reliability.

**And it sharpens §82 rather than qualifying it.** Both neutral pairs are irreversible as
written, so their dissipation is formally infinite. The element therefore pays an *unbounded*
thermodynamic cost and gets **worse** reliability for it — drive and reliability moving in
opposite directions inside a single network, which is a stronger statement than §82's
independence.

---

### 84 ν = 2 exactly — T15-n closed, and every route was measuring a window

§64 extracted the barrier exponent in `A ~ (γ_c − γ)^ν` three ways and got three answers — width
1.95–2.03, stationary 1.99, extrapolated action 2.10–2.19 — and left T15-n open on *why*. §64's
own present statement was **"ν ≈ 2 ± 0.1, and not determined more precisely than that by any
instrument here."**

**Every one of those is a fit, and rule 16 is exactly about this: nobody computed what the normal
form predicts.** So this section derives A(γ) in closed form and checks it in *absolute* terms.

**The derivation.** Near the pitchfork b is fast (symmetric eigenvalue −(1+2γ)) and the lead
u = x−y is slow. The lead's drift factors exactly — this is §43's invariance — as
`du/dt = k·u·[b(1+γ) − γ]`, so the saddle sits at `b* = γ/(1+γ)`. Eliminating b from `ds/dt = 0`
gives `b(u) = b₀ + b₂u²`, and **b₀ = 1/3 exactly for every γ** (the symmetric fixed point never
moves), so

    ε = b₀ − b* = (1−2γ)/(3(1+γ)),     b₂ = −(1−γ)/(2(1+2γ))

Per **§83**, the noise must come from the propensities and not the drift — a normal form carrying
only the drift would be precisely the error §83 identified. Summing (Δu)² over the four
lead-changing reactions gives `D_u(0) = 2k(1+γ)/9`, and `A = 2∫₀^{u*} F/D du` collapses to

> **A = 9ε²/(4|b₂|) = (1+2γ)(1−2γ)² / (2(1−γ)(1+γ)²)**
>
> manifestly ∝ (γ_c−γ)², so **ν = 2 exactly.** k, D_u and the prefactor all cancel: the formula
> has no free parameter of any kind.

**P1, the gate.** `k(1+γ)ε` reproduces the module's independently-computed `lambda_antisym(γ)` to
**5.6e−17** at eight γ, and the normal form's δ* converges to the exact `delta_star`:
0.1835 → 0.1136 → 0.0646 → 0.0320 → 0.0117 → 0.0050 → 0.0009 → **0.0001**. (Rule 20: it is a
leading-order form, 18% off at γ=0, and that is not a failure — convergence is the criterion.)

**P2, the test, absolute — closed form against the exact first-passage action:**

| γ | A measured | A closed form | ratio |
|---|---|---|---|
| 0.30 | 0.125373 | 0.108199 | 0.8630 |
| 0.38 | 0.046444 | 0.042929 | 0.9243 |
| 0.42 | 0.021153 | 0.020138 | 0.9520 |
| 0.44 | 0.012080 | 0.011657 | 0.9650 |
| 0.45 | 0.008458 | 0.008215 | 0.9713 |
| 0.46 | 0.005460 | 0.005338 | **0.9776** |

Monotone, from below, with residual/(γ_c−γ) = 0.685 → 0.559 — an O(γ_c−γ) correction, as an
omitted next-order term must be. **Zero cells excluded; zero fitted parameters.**

**P3/P4, why the three routes disagreed. ν is a limit, and every finite window reads below it.**
The closed form's own effective exponent runs **1.8300** on §63.2's window [0.20, 0.45] and
climbs to **1.9898** on [0.48, 0.499]. The *measured* A gives 1.9424 / 1.9416 / 1.9430 over three
windows — **squarely on §64's width route (1.95–2.03) and its stationary route (1.99)**.

> **§63.2's "2 is excluded over [0.20, 0.45] with no drift toward it" was reading the window
> bias.** An effective exponent below 2 over a finite window is exactly what ν = 2 with a
> correction to scaling produces. The exclusion was already withdrawn by §64; this says what it
> was. §64's present statement "ν ≈ 2 ± 0.1, not determined more precisely by any instrument
> here" is now superseded — **ν = 2 exactly, by derivation, validated absolutely.**

### 84.1 Two criteria that were wrong, in a section about wrong criteria

> **P3 as first written demanded the closed form land within 0.05 of §64's width route on that
> window, and printed FAILS off 1.8300.** But §64's routes fit a width, a stationary distribution
> and an extrapolated action over their own γ grids — **none of them this quadrature**. A
> threshold between effective exponents of *different observables fitted by different protocols*
> cannot be satisfied only by the thing it claims to test (rule 19). What is comparable is the
> sign and direction of the window bias, and that is what the verdict now tests.
>
> **P4 gated the measured-vs-closed-form exponent gap at 0.10 and passed at 0.091** — a fixed
> tolerance on a quantity that converges, which is rule 20 verbatim, two sections after rule 20
> was last invoked.

The gap it was hiding is real and is reported instead of gated: the closed form's window bias is
**larger** than the true one (1.8513/1.8691/1.8814 against 1.9424/1.9416/1.9430), in the same
direction at every window. That is a leading-order artifact — the closed form's subleading term is
not the true one, which is also why P2's ratio sits below 1 away from γ_c. **ν = 2 rests on P1 and
P2, the absolute check, and on none of these fits.**

**Still open (T15-n.1):** the quantitative gap between the closed form's window exponent (1.83)
and §64's measured routes (1.95–2.19) on the same window is *not* explained here.
