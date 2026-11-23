"""§134.31 (T-CASC-bb): the finite-gap correction in closed form, E = E0 F(z), against the linear model and the chain.

Derivation (written before running). Source at stretched (lam g0, -lam L0), lam = sqrt(Omega lambda_u). After the
deterministic inflow reaches the inner scale (t* = ln(lam L0)/mu), the mean's normal coordinate is lam g0 e^t and
Sigma_11 = (e^{2t} - 1)/2, so the Gaussian at H carries exp(-lam^2 g0^2 / (e^{2t} - 1)) beyond the action S(H) =
lam^2 g0^2. With u = e^{-t}, v(H) e^{S(H)} ~ int_0^{u*} exp(-(lam g0)^2 u^2) du ~ erf(z) / (lam g0), u* = (lam L0)^{-1/mu},
z = lam g0 u* = g0 L0^{-1/mu} lam^{1 - 1/mu}. The reference point on the path gives lam^{-1}, so
Q = ln erf(z) + const, and
    E = dQ/dlnOmega = E0 F(z),   E0 = (mu - 1)/(2 mu),   F(z) = 2 z e^{-z^2} / (sqrt(pi) erf z).
Limits: z -> 0 gives E0 (§134.29); z -> inf gives 0 (the decaying control, §134.28).

Predictions:
  F1 (the derivation, against the model it came from) the window-averaged E from the formula matches the linear model's
     converged-grid E_lin (results/saddle_graze_grid4.json) within 0.01 in every window, for all four geometries
     (XK-like, X-like, S-like, control). Failure means the derivation misses a term.
  F2 (the chain) then the formula with each pair's own (mu, lambda_u, D22, gap, L0) gives the chain's E at H within
     0.02 in the largest windows: XK 0.374 (2800->4000), S 0.318 (2800->4000), X 0.359 (1200->1600).
"""
import json
import math
from pathlib import Path

from scipy.special import erf

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

GEOM = {   # mu, lambda_u, D11, D22 (from saddle_graze.params), gap (§134.18 at x1 = 2.7)
    "XK_like": (5.2083, 1.09166, None, None, 0.0107),
    "X_like": (5.2083, 1.09166, None, None, 0.0024),
    "S_like": (3.3333, 1.70573, None, None, 0.0122),
    "control_gap0.10": (5.2083, 1.09166, None, None, 0.10),
}


def Q(om, mu, lu, g0, L0):
    lam = math.sqrt(om * lu)
    z = g0 * L0 ** (-1.0 / mu) * lam ** (1.0 - 1.0 / mu)
    return math.log(erf(z)), z


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(HERE))
    import saddle_graze as sg
    grid = json.loads((ROOT / "results" / "saddle_graze_grid4.json").read_text())
    oms = (400, 800, 1600, 2800, 4000)
    out = {}
    for name, (mu, lu, _, _, gap) in GEOM.items():
        hk, r2 = {"XK_like": ("1.2", "0.606"), "X_like": ("1.0", "0.606"), "S_like": ("1.0", "1.005"),
                  "control_gap0.10": ("1.2", "0.606")}[name]
        mu_, c, lu_, D11, D22 = sg.params(hk, r2)
        g0, L0 = gap / math.sqrt(D22), 1.1 / math.sqrt(D11)
        E0 = (mu_ - 1) / (2 * mu_)
        Qs = [Q(o, mu_, lu_, g0, L0) for o in oms]
        Ef = [(Qs[i + 1][0] - Qs[i][0]) / math.log(oms[i + 1] / oms[i]) for i in range(len(oms) - 1)]
        El = grid[name]["E_grid4"]
        out[name] = {"E_formula": Ef, "E_linear": El, "z": [q[1] for q in Qs]}
        print(f"{name}: z {Qs[0][1]:.3f}..{Qs[-1][1]:.3f}  formula " + " ".join(f"{e:+.4f}" for e in Ef) +
              "  | linear " + " ".join(f"{e:+.4f}" for e in El), flush=True)
    (ROOT / "results" / "graze_formula.json").write_text(json.dumps(out, indent=1))

# OUTCOME (§134.31): F1 FAILED as registered. Formula vs linear model (converged grid), per window:
#   XK_like 0.400/0.398/0.393/0.389 vs 0.397/0.388/0.384/0.381 (max 0.009); X_like 0.404/0.404/0.404/0.403 vs
#   0.401/0.393/0.392/0.392 (max 0.011 > 0.01); S_like 0.348/0.347/0.346/0.345 vs 0.346/0.341/0.339/0.339 (max 0.007);
#   control 0.170/0.085/0.027/0.007 vs 0.204/0.132/0.077/0.047 (max 0.050).
# The structure is right -- the scaling variable z, both limits, the ordering -- but two subleading terms are missing:
# a ~0.01 finite-lam offset at small z (the g = 0 model itself reaches E0 only slowly, §134.29) and a slower large-z
# tail than erf's. F2 is therefore not credited to the derivation (formula XK 0.389 vs chain 0.374).

# REGISTERED BEFORE RUNNING (§134.32): keep the stable direction's smooth factor instead of a step at u*:
#   v(H) e^{S(H)} ~ I(a, b) = int_0^inf exp(-a u^2 - b u^{2 mu}) du, a = lam^2 g0^2, b = mu lam^2 L0^2;  Q = ln I + (1/2) ln Omega.
# At a = 0 this is E0 exactly at every lam, so it CANNOT produce the small-z finite-lam offset (that must come from early
# times u ~ 1 or from the reference point). Predictions: G-tail -- it closes the control's large-z residual to <= 0.01;
# G-offset -- it leaves the small-z offset (~0.01) unchanged. Either failing says the split of the residual is wrong.
# OUTCOME (§134.32): I-form vs linear, max |d|: XK 0.0129, X 0.0115, S 0.0081, control 0.0240 (from 0.050).
# G-offset HELD (small-z offset unchanged, ~0.011). G-tail FAILED as registered (control halved, not <= 0.01). The
# residual has an early-time / reference-point term (~0.01) and a tail term only partly fixed here. Refinement stopped.
