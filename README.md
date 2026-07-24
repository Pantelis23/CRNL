# CRNL — Chemical Reaction Network Landscape

**Chemical reaction networks: from logic to landscape.**

A small simulation rig whose purpose is epistemic, not performative. It is not a
chemical computer, not a fast solver, and not a library anyone needs. It exists
to make one property measurable: **signal restoration** — the ability of a
physical system to keep its states distinguishable against noise, indefinitely,
across a deep cascade.

Binary did not win because 2 is a special number. It won because the transistor
is a near-ideal restoring switch. Chemistry, given the right network motif, can
restore too — but only by running away from equilibrium and paying free energy
for it. CRNL runs one such motif (Approximate Majority) two ways —
**deterministic** mass-action ODEs and **exact stochastic** Gillespie SSA — and
measures the gap. Everything it teaches lives in that difference.

> Full rationale and derivations: [`docs/design.md`](docs/design.md).

## The one-paragraph physics

Approximate Majority is three reactions on two committed species X, Y and a blank
B, all with rate 1:

    r1:  X + Y → 2B      disagreement cancels both to blank
    r2:  B + X → 2X      the leader recruits blanks (autocatalysis)
    r3:  B + Y → 2Y      mirror

Deterministically this has two stable rails (all-X, all-Y) separated by a saddle
at (⅓,⅓,⅓) — the decision threshold. Any bias is amplified to a clean rail; noise
that knocks you off a rail decays. That is restoration. But the deterministic
separatrix is a perfect wall only at infinite population. At finite molecule
count Ω, fluctuations of order √Ω can push a decision over the saddle, and the
error probability follows a **restoration wall**:

    P(error) ~ exp(−c(ε)·Ω)

The barrier `c(ε)` is the noise margin; Ω is the multiplier on it. Restoration is
never deterministic — only exponentially reliable.

## Results

Run from the repo root (see setup below).

### The restoration wall — `experiments/restoration_wall.py`

![restoration wall](experiments/restoration_wall.png)

Left: the conditional error fraction `Y/(X+Y)` falls log-linearly with Ω, fitting
`exp(−c(ε)·Ω)` with measured **c(0.10) ≈ 0.018, R² ≈ 0.94**. The deterministic
ODE glides to the X rail every time — error exactly 0 at every Ω, off the log
floor: that curve is the lie. Right: the **all-blank** outcome, which the
deterministic view calls an impossible repeller, genuinely occurs at low Ω (~5%
at Ω=6) and vanishes as Ω grows — a pure finite-count effect with its own Ω
scaling.

```bash
python -m experiments.restoration_wall                 # default: eps=0.10, clean wall
python -m experiments.restoration_wall --trials 20000  # tighter statistics
python -m experiments.restoration_wall --bias 0.02     # the literal 51/49 (shows the crossover, not the wall — see note below)
python -m experiments.restoration_wall --quick         # fast smoke run
```

### The landscape — `experiments/phase_portrait.py`

![phase portrait](experiments/phase_portrait.png)

The deterministic flow (blue streamlines) rolls downhill to a rail; the red
diagonal is the separatrix (the restoring threshold); grey lines are exact
Gillespie trajectories fluctuating around the smooth green ODE from the same
start. The saddle's role as a threshold is *basin structure*, not a slogan.

```bash
python -m experiments.phase_portrait
```

## Why the default bias is 55/45, not 51/49

The design doc's illustrative 51/49 (ε = 0.02) has an *intrinsically tiny*
barrier: `c·Ω` stays below ~1 across the whole observable window, so the wall
never clears the algebraic-prefactor crossover until Ω reaches the thousands —
where errors collapse to ~e⁻²⁵ and every trial reads as correct (you measure
zero, not a small number). This is the squeeze §4 of the design doc warns about,
made concrete by actually running it. A slightly larger — still small — bias puts
a clean `exp(−c·Ω)` squarely in the accessible window. The experiment stays fully
parameterized by `--bias`, so 51/49 is one flag away; it simply shows the
crossover rather than the wall.

## Radix experiment (n-winner AM)

Generalizes AM to n committed symbols to measure **radix vs. margin** — the
project's core claim — as data.

- `experiments/radix_wall.py` — champion-vs-field. One symbol leads each rival by
  a fixed pairwise margin δ (55/45 at n=2). Measures the barrier `c(n)` and the
  population `Ω_required(n)` to hold a fixed reliability as the alphabet grows.
- `experiments/radix_discovery.py` — symmetric start. Characterizes the outcome
  distribution (single winner / all-blank / coexistence / undecided) and consensus
  time vs n, under fixed-total-Ω and fixed-density conventions — built to reveal
  high-n failure modes (blank collapse, long coexistence) rather than confirm a
  prediction. At the tested (n, Ω) grid the system stays robust (single-winner
  ≈ 1.0); the cost of radix shows up instead as a falling barrier c(n) and rising
  consensus time, not as collapse.

```bash
python -m experiments.radix_wall --quick
python -m experiments.radix_discovery --quick
```

The engine reaches n≈100 via a NumPy-vectorized SSA path (`crnl/vectorized.py`)
validated to match the readable reference propensities to 1e-12 (rtol), including
the boundary states where naive fast paths diverge.

## Setup

Requires Python 3.10+.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Verify (the physics checks are the tests)

```bash
pytest -q
```

The suite verifies each build stage: AM's stoichiometry is conservative; the RHS
matches the reduced ODEs; the fixed-point eigenvalues reproduce §2.3 exactly
(−1,−1 rails; −1,+⅓ saddle; +1,+1 repeller); the homodimer and heterobimolecular
units conventions are correct; **the SSA converges to the ODE as Ω grows** (the
single best test that the units convention is right); every AM trajectory absorbs
into one of three bins including all-blank; and seeded trajectories replay
exactly.

## Layout

| Path | Role |
|------|------|
| `crnl/reactions.py` | species/reaction data model; builds S; **owns the units convention and propensity builder** (§3.2, §3.3) |
| `crnl/deterministic.py` | scipy LSODA path; `S·v(x)`; analytic Jacobian; conservation monitoring |
| `crnl/stochastic.py` | hand-written Gillespie SSA — the lesson |
| `crnl/classify.py` | absorption test + dwelling test + fixed-point classifier (stoichiometric-subspace aware) |
| `crnl/networks/am.py` | AM as data: 3 species, 3 reactions, k=1 |
| `crnl/networks/n_winner.py` | n-winner AM as data: n committed species + blank, pairwise disagreement + per-species autocatalysis |
| `crnl/vectorized.py` | NumPy-vectorized SSA path validated against the reference propensities, letting the radix experiments reach n≈100 |
| `experiments/restoration_wall.py` | the §4 protocol |
| `experiments/phase_portrait.py` | the §2.3 landscape, made visible |
| `experiments/radix_wall.py` | champion-vs-field barrier c(n) and population cost Ω_required(n) as the alphabet grows |
| `experiments/radix_discovery.py` | symmetric-start outcome distribution and consensus time vs alphabet size |
| `tests/test_engine.py` | the verification suite |
| `tests/test_n_winner.py` | n-winner network construction and stoichiometry checks |
| `tests/test_radix_experiments.py` | radix_wall / radix_discovery helper and fit checks |
| `docs/design.md` | full design rationale |

The engine is general: it takes species, reactions, and rate constants and
derives both dynamics from that same data. AM is the first network loaded into
the engine — it is not the engine. n-winner AM / the radix experiment is now
implemented (see the Radix experiment section above and
`experiments/radix_wall.py` / `radix_discovery.py`). The remaining
out-of-scope-for-v1 extensions — analytic saddle height and free-energy
accounting — are sketched at the end of `docs/design.md`.
