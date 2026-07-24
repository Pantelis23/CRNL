# CRNL — Chemical Reaction Network Landscape

**Chemical reaction networks: from logic to landscape.**

## 0. What this is, and what it is deliberately not

CRNL is a small simulation rig whose purpose is epistemic, not performative. It
is not a chemical computer, not a fast solver, and not a library anyone needs.
It exists to make one property measurable: **signal restoration** — the ability
of a physical system to keep its states distinguishable against noise,
indefinitely, across a deep cascade.

The claim to be felt from the inside: binary did not win because 2 is a special
number. It won because the transistor is a near-ideal restoring switch.
Chemistry, given the right network motif, can restore too — but only by running
away from equilibrium and paying free energy for it. Alphabet size was never the
bottleneck; restorable, composable, addressable interactions were.

The name records the arc. The project starts from the question *can you build a
logic out of chemistry — a richer vocabulary than binary or ternary?* and the
expected finding is that the vocabulary framing dissolves: the universe does not
compute with symbols, it computes by relaxation — rolling downhill on a
landscape until it settles, with the answer being *which attractor*, not *which
gate output*. Logic in, landscape out.

**Design principle that governs everything below:** the value of this simulation
is proportional to how much noise and dissipation it deliberately keeps in. This
inverts normal engineering discipline. A clean, well-behaved simulation of
chemical logic is a simulation that lies — it abstracts away exactly the
nonidealities where the entire lesson lives. Resist every instinct to smooth,
clamp, or idealize.

## 1. Core method

Run the same reaction network two ways, and measure the gap.

- **Deterministic:** mass-action ODEs, infinite population, continuous
  concentrations.
- **Stochastic:** exact Gillespie SSA, finite molecule count Ω, integer counts.

The deterministic path gives perfect readout, a perfectly sharp decision
threshold, and reproducibility. The stochastic path puts the noise back.
Everything CRNL teaches lives in the difference between them.

## 2. The physics the engine encodes

### 2.1 Mass action is polynomial ODEs

A CRN is a set of species plus a set of reactions. Each reaction *j* has a rate
constant *k_j* and a flux given by the law of mass action: *k_j* times the
product of reactant concentrations, each raised to its stoichiometric
coefficient.

Collect fluxes into a vector **v(x)**, and stoichiometry (products − reactants,
per species per reaction) into a matrix **S**. The dynamics is entirely:

    dx/dt = S · v(x)

Every flux is a monomial in concentrations, so the RHS is a vector of
polynomials. The near-converse (Hárs–Tóth) means this class of dynamical system
and this class of chemistry are essentially the same object — which is why a
sandbox built on CRNs teaches something general rather than something about
beakers.

**Convention, fixed once:** flux *v* is *k* × product of concentrations to their
stoichiometric powers, and the stoichiometric coefficient lives in **S**. So for
2A → B, *v = k[A]²* and *S_A = −2*, giving *d[A]/dt = −2k[A]²*. Every rate
constant in the codebase means this and nothing else.

### 2.2 The first instance: Approximate Majority (AM)

Two committed species X, Y and a blank B. Three reactions, all *k = 1*:

    r1:  X + Y → 2B      disagreement — a collision cancels both to blank
    r2:  B + X → 2X      autocatalysis — the leader recruits blanks
    r3:  B + Y → 2Y      autocatalysis, mirror

Every reaction is 2→2, so total count is conserved. Normalize to 1, set
*b = 1 − x − y*:

    dx/dt = x(1 − x − 2y)
    dy/dt = y(1 − 2x − y)

Three reactions, and it computes consensus while restoring its own signal.

### 2.3 Why AM restores — the landscape

Jacobian: `J = [[1−2x−2y, −2x], [−2y, 1−2x−2y]]`

| Fixed point       | Eigenvalues | Type          | Meaning              |
|-------------------|-------------|---------------|----------------------|
| (1, 0, 0) all-X   | −1, −1      | stable node¹  | clean output "X"     |
| (0, 1, 0) all-Y   | −1, −1      | stable node¹  | clean output "Y"     |
| (⅓, ⅓, ⅓)         | −1, +⅓      | saddle        | the decision threshold |
| (0, 0, 1) all-blank | +1, +1    | repeller      | undecided; a whisper blows it up |

¹ Technically a degenerate (improper) node — the repeated eigenvalue −1 has a
single eigenvector, not a full eigenbasis. Still asymptotically stable.

The saddle is the whole point. A restoring device is a system with multiple
stable outputs separated by an unstable threshold that pushes inputs away from
the boundary. The unstable eigenvector at the saddle is (1, −1): the difference
X−Y is what grows, so any initial bias is exponentially amplified to a clean
rail. And the negative eigenvalues at the rails mean noise that knocks you off an
output decays: errors are corrected rather than accumulated — the one property a
deep cascade requires.

Generic chemistry cannot do this. A system at thermodynamic equilibrium sits in a
single free-energy minimum — one stable state, no threshold, no restoration. AM
escapes only because the autocatalytic steps are effectively irreversible; it
runs away from equilibrium. Restoration costs dissipation. (Same positive-
feedback motif Cardelli identified at the core of the biological cell-cycle
switch.)

### 2.4 Where the deterministic view lies, quantitatively

Deterministically the separatrix is a perfect wall. At finite Ω the saddle is a
hill of finite height, and fluctuations of order √Ω can push you over it — wrong
decisions, and spontaneous flips.

Erroneous-crossing probability scales as

    P(error) ~ exp(−c(ε)·Ω)

(Kramers / large-deviation escape), where *c(ε)* is a quasipotential barrier
height set by network geometry and the initial bias *ε*.

The barrier *c(ε)* is the **noise margin**; Ω is the **multiplier** that scales
your protection. Restoration is never deterministic — only exponentially
reliable. This is the radix-vs-margin tradeoff from the chemistry side: pack more
distinguishable states into the same population and the basins crowd, lowering
the barrier per symbol. Exactly why ternary keeps losing to binary in silicon.

Formal backstop: finite stochastic CRNs are Turing-universal (Soloveichik, Cook,
Winfree & Bruck, 2008), and AM's n-winner generalization is plurality consensus
in the population-protocol sense (Angluin).

## 3. Engine design

General, not AM-specific. The engine takes species, reactions, and rate
constants, and derives both dynamics from that same data. AM is the first
network loaded into the engine — it is not the engine.

### 3.1 Two paths

- **Deterministic.** Build S and v(x); integrate `dx/dt = S·v(x)` with
  `scipy.integrate.solve_ivp`. Default LSODA (adaptive, switches to implicit) so
  a future network with rate constants spanning orders of magnitude doesn't
  explode and masquerade as diverging physics.
- **Stochastic.** Hand-write the Gillespie loop. Compute propensities *a_j*;
  total *a₀ = Σa_j*; draw *τ = −ln(u₁)/a₀*; select reaction *j* with probability
  *a_j/a₀*; apply stoichiometry column *j*; advance *t += τ*; repeat until
  *a₀ = 0* (absorption) or a step budget is exhausted.

### 3.2 Units convention — decided once, at engine level

A single scalar Ω = molecule count at concentration 1. Then
concentration = count / Ω, and stochastic constants derive from deterministic
ones as:

| Reaction order        | Deterministic | Stochastic constant |
|-----------------------|---------------|---------------------|
| unimolecular A → …    | k             | c = k               |
| heterobimolecular A+B | k             | c = k / Ω           |
| homodimer A+A         | k             | c = 2k / Ω          |

Bimolecular constants carry the Ω scaling; unimolecular do not — and the two
kinds of bimolecular do not share the same constant.

The homodimer 2 is the same fact as the ½ in §3.3, from the other side.
Propensity for A+A is *c·n(n−1)/2*; each firing changes *n_A* by −2; so the mean
drift is *−2 · c n²/2 = −c n²*, i.e. *d[A]/dt = −cΩ[A]²*. Matching the
deterministic *−2k[A]²* forces *c = 2k/Ω*. Write *k/Ω* for a homodimer and you
have silently halved its noise relative to its own ODE.

Getting any of this wrong rescales the noise without touching the deterministic
trajectory — the most insidious bug class in the project. Enforced at engine
level; no individual network gets a vote.

The engine generalizes the table with one rule: for a reaction of total order
*m*, *c = k · s / Ω^(m−1)*, where *s = Π_i (coeff_i!)* is the product of the
factorials of the reactant coefficients. (unimolecular s=1 m=1; heterobimol. s=1
m=2; homodimer s=2!=2 m=2; trimolecular 3A gives s=6 m=3.)

### 3.3 Combinatorics, at engine level

Propensity = *c* × number of distinct reactant combinations:

    heterobimolecular A + B:  a = c · n_A · n_B          with c = k/Ω
    homodimer A + A:          a = c · n_A(n_A − 1)/2      with c = 2k/Ω

The ½ and the 2 must be owned together by the propensity builder. AM has no
homodimer and never triggers this — which is precisely why it must be handled at
engine level rather than discovered later.

### 3.4 Classifier — attractor-count-agnostic, two distinct criteria

The classifier returns however many end-states a network has, with no hardwired
count. But "end-state" is two different objects in the two engines.

- **Stochastic absorption:** a configuration where every propensity is zero
  (*a₀ = 0*). The chain physically halts. Cheap and exact. AM has exactly three
  such configurations: (Ω,0,0), (0,Ω,0), (0,0,Ω). Three is a theorem about AM,
  not a constant to hardwire.
- **Deterministic stable end-state:** a stable root of *S·v(x) = 0* (Jacobian
  eigenvalues with negative real part). A settling criterion, not a halting one.

These need two code paths. A general CRN can have a stable interior fixed point
at which the stochastic system never absorbs — it fluctuates around it forever.
For those the stochastic classifier needs a **dwelling test**: the trajectory has
remained within a small ball of a known fixed point for longer than some multiple
of the local relaxation time. AM is the lucky case where all three attractors are
simultaneously zero-propensity corners and stable roots, so the *a₀ = 0* test
suffices and the dwelling path goes unused — but the engine must carry both, or
it is AM-only wearing a general costume.

**The third bin matters.** All-blank (0,0,Ω) is stochastically absorbing and
genuinely reachable at finite Ω, when the last committed pair annihilates via r1
before autocatalysis amplifies a lead. The deterministic view calls (0,0,1) a
repeller and would never rest there; the stochastic view lands there with nonzero
probability precisely because it can hit the corner exactly at integer count.
Give it its own bin and report it.

> **Stability inside the stoichiometric subspace.** The full three-species
> Jacobian carries a spurious zero eigenvalue per conservation law (AM's 2→2
> reactions make (1,1,1) a left null vector of S). The classifier projects the
> Jacobian onto Im(S) before taking eigenvalues, recovering the reduced 2-D
> picture of §2.3 exactly — rails (−1, −1), saddle (−1, +⅓), all-blank (+1, +1).

## 4. The experiment: measuring the restoration wall

Goal — watch *exp(−c(ε)Ω)* fall out of data, and read off *c(ε)*.

1. Fix a small fractional bias *ε*. Hold the **fraction** constant across Ω —
   never the absolute count difference. Only the fractional convention keeps the
   starting distance from the separatrix Ω-independent. The measured exponent is
   *c(ε)*, not a universal constant; *ε* is part of the result.
2. Sweep Ω across the observable window — roughly 20 to a few hundred. Ω in the
   thousands drives error to ~e⁻²⁵, so every trial comes out correct and you
   measure zero.
3. Run many Gillespie trials per Ω, to absorption. Bin into {X-wins, Y-wins,
   all-blank}.
4. Report **two** curves: raw error fraction (Y / all) and conditional-on-
   decision (Y / (X + Y)). At low Ω the all-blank bin is non-negligible, so these
   differ; collapsing them hides the finite-count effect. The clean exponential
   lives in the conditional fraction.
5. Fit on the linear region only. The error-rich window overlaps the region where
   an algebraic prefactor still curves the log plot; treat a curved low-Ω tail as
   prefactor, not exponent.
6. Contrast against the ODE, which glides to the X rail every time, at every Ω.
   The deterministic curve is the lie; the stochastic error fraction is the
   truth; *c(ε)·Ω* is the restoration wall as a number.

> **A note on the bias, from actually running it.** At the literal 51/49
> (ε = 0.02) the quasipotential barrier *c(ε)* is intrinsically tiny, so *c·Ω*
> stays below ~1 across the entire observable window and the wall never clears
> the algebraic-prefactor crossover until Ω ~ thousands — where errors fall to
> ~e⁻²⁵ and read as exactly zero. This is the §4 squeeze made concrete. A
> somewhat larger (still small) bias — the default ε = 0.10 (55/45) — puts a
> clean *exp(−c·Ω)* squarely in the accessible window: measured
> **c(0.10) ≈ 0.018, R² ≈ 0.94**. The experiment stays parameterized by `--bias`
> so 51/49 remains reproducible; it just shows the crossover rather than the
> wall.

## 5. Traps (all bite regardless of network)

- **Homodimer factor.** Propensity *c·n(n−1)/2* (not *c·n²*) and constant
  *c = 2k/Ω* (not *k/Ω*). Same fact, two places.
- **Ω scaling lives only on bimolecular constants.** Misplace it and the noise
  rescales invisibly while the ODE still looks right.
- **Boundary drift — do not naively clamp.** The polynomial RHS can push a
  concentration slightly negative near *x = 0*. Clamping injects mass and breaks
  conservation at the most delicate point. Monitor the conserved total; use tight
  tolerances. (CRNL uses an *analytic* Jacobian precisely because a central
  difference gets corrupted at the boundary, halving one-sided derivatives.)
- **Overlay normalization.** ODE lives in [0,1], Gillespie in [0,Ω]. Divide
  counts by Ω before comparing.
- **Stiffness later.** Default to LSODA so a future wide-spread-rate network
  doesn't blow up and get mistaken for divergent physics.
- **RNG discipline.** Seed per-trial deterministically from (Ω, trial index) so
  any single anomalous trajectory can be replayed exactly.

## 6. Layout

    crnl/
      reactions.py     species + reaction data model; builds S; owns the units
                       convention and the propensity builder (§3.2, §3.3)
      deterministic.py scipy LSODA path; S·v(x); analytic Jacobian; conservation
      stochastic.py    hand-written Gillespie SSA (§3.1) — the lesson
      classify.py      absorption test + dwelling test + fixed-point classifier
      networks/am.py   AM as data: 3 species, 3 reactions, k=1
    experiments/
      restoration_wall.py   the §4 protocol
      phase_portrait.py     the §2.3 landscape, made visible
    tests/
      test_engine.py   the physics checks that verify each build stage
    docs/
      design.md        this file

**Build order:** reactions → deterministic (verify AM reproduces the §2.3 fixed
points) → stochastic (verify it converges to the ODE as Ω grows — the single best
test that the units convention is right) → classifier → experiment.

The large-Ω agreement check catches a §3.2 error, because a wrong stochastic
constant leaves the deterministic path untouched and shows up only as the two
engines failing to converge.

## 7. What success looks like

- the same three reactions produce a smooth deterministic decision and a noisy
  probabilistic one;
- the saddle's role as a restoring threshold is visible as basin structure, not
  as a slogan (`phase_portrait.py`);
- the restoration wall appears as a measured *exp(−c(ε)Ω)*, with the barrier
  *c(ε)* as the margin and Ω as the multiplier (`restoration_wall.py`);
- the all-blank outcome shows up exactly where the deterministic repeller said
  nothing could rest.

The universe does not compute with a vocabulary. Physical systems compute by
relaxation — the logic is emergent from dynamics, and the output is an attractor,
not a gate result.

## 8. Extensions, deliberately out of scope for v1

- **n-winner AM (the radix experiment).** *n* committed species plus blank,
  pairwise disagreement Xᵢ+Xⱼ→2B and per-species autocatalysis B+Xᵢ→2Xᵢ. As *n*
  grows at fixed Ω, the basins carve the simplex into smaller pieces,
  separatrices crowd, and the per-symbol barrier falls — the radix-vs-margin
  tradeoff made measurable. (The engine already handles arbitrary species/
  reaction counts, so this is a new `networks/*.py` plus a sweep.)
- **Analytic saddle height.** Predict *c(ε)* from the quasipotential before
  measuring it, then watch the data land on the line.
- **Free-energy accounting.** Dissipation cost per restoration event, and the
  Landauer-adjacent floor beneath it.
