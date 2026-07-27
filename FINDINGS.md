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

### 17.2 The coefficient rises with Ω, and this data cannot say where it stops

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

P3 itself shows the same incompleteness directly: `β* ∝ Ω^−x` with x measured at
**0.508 / 0.603 / 0.678 / 0.752** across consecutive Ω pairs — drifting toward the
predicted 1 and nowhere near it yet.

### 17.3 The careful refinement was worse than the crude argument

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
