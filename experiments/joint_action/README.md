# joint_action — the least joint action of a never-dipped downstream fall (§124.7)

Freidlin–Wentzell / WKB minimum action for the two-stage cascade: the cheapest joint path from the
stable point (both stages high) to the joint deterministic separatrix Σ (the boundary between
initial states whose coupled rate equations end with stage 2 high vs low), with stage 1 kept above
its saddle s₁. Γ = V₁(s₁) − W* is the exponential growth rate of the never-dipped term relative to
the dip rate.

Written by an independent agent as a cross-check, reviewed before inclusion. Run from this directory.
Pair is selected by environment: `R2_STAGE1` (default 1.005, E-070), `R2_STAGE2` (default 0.606,
E-084); rails 0.15 / 3.0 throughout.

| script | what |
|---|---|
| `model.py` | rates, Hill coupling, the WKB Hamiltonian and its equations, V₁ |
| `shoot.py` | rays launched from the linear-noise ellipse (Lyapunov covariance) around (3, 3) |
| `sigma.py` | Σ by bisection on the coupled rate equations |
| `main_sigma.py` | ray family vs Σ; `p20.py` the free-x₂ (p₂ = 0) ray, which gives a bound |
| `graze.py EPS RTOL NSIG` | bisects for the limiting ray into the joint saddle; W* extrapolated to \|p\| → 0 |
| `mam.py N init.npy` | independent direct minimum-action check (converges from above) |

Converged values are stored in `results/joint_action.json`. `tests/test_joint_action.py` re-runs the
free-x₂ bound for X (7 s).
