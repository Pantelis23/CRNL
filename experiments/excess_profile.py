"""§134.21 (registered in inner_linear_box.py before running): L = ln P + Omega W along the normal through H at fixed
STRETCHED distance r. W off H is not W*: use W(H + delta e_n) = W* + (1/2) Z_u,22 delta^2 on B's unstable manifold
(Z_u = diag(0.1502, 0) in x -- zero along x2), i.e. W = W* to second order along x2; the BVP's own Z22 = 0.66 is on the
noise-stable side (the incoming one, above H). Both are reported: L_flat (W = W*) and L_in (W* - (1/2) 0.66 delta^2,
the incoming-side quadratic above H)."""
import json, math, sys
from pathlib import Path
import numpy as np
import experiments.episode_committor as ec
from experiments.sigma_flux import ARMS, element
ec.SOLVER = sys.argv[2] if len(sys.argv) > 2 else ec.SOLVER
W_STAR, Z22 = 0.1484037485, 0.6597751589951498
A, B = element(ARMS["XK"][0]), element(ARMS["XK"][1])
out = {}
for om in [int(v) for v in sys.argv[1].split(",")]:
    m = A.cap(om) + 1
    sad1, sad2 = int(np.ceil(A.r2 * om)), int(np.ceil(B.r2 * om))
    ec.qprocess_J(A, B, om, m, sad1, sad2, lower="rail")
    L = ec.qprocess_J.last
    P, n1s, base2 = L["P"] / ec.PIN_SCALE, L["n1s"], L["base2"]
    i = int(round(3.0 * om)) - int(n1s[0])
    rec = {}
    for r in (-10, -3, -1, 0, 1, 3, 10):
        d = r / math.sqrt(om * Z22)
        x2 = (B.r2 + d) * om - base2
        j0 = int(math.floor(x2)); f = x2 - j0
        lnP = (1 - f) * math.log(P[i, j0]) + f * math.log(P[i, j0 + 1])
        rec[str(r)] = {"delta": d, "L_flat": lnP + om * W_STAR,
                       "L_in": lnP + om * (W_STAR - 0.5 * Z22 * d * d) if r > 0 else None}
    out[str(om)] = rec
    print(om, " ".join(f"r{k}:{v['L_flat']:+.3f}" for k, v in rec.items()), flush=True)
Path("results/excess_profile.json").write_text(json.dumps(out, indent=1))
ks = sorted(out, key=int)
for r in ("-10", "-3", "-1", "0", "1", "3", "10"):
    s = [(out[b][r]["L_flat"] - out[a][r]["L_flat"]) / math.log(int(b) / int(a)) for a, b in zip(ks, ks[1:])]
    print(f"r={r:>3}: dL_flat/dlnOmega per window " + " ".join(f"{v:+.3f}" for v in s))

# OUTCOME (§134.21): last-window slopes dL_flat/dlnOmega at fixed stretched r along x2 through H:
#   r = -3: +0.088 | -1: +0.279 | 0: +0.368 | +1: +0.454 | +3: +0.618 | +10: +1.128 (still rising); r = -10 is in the
#   low basin's steep side (values meaningless there). A W quadratic near H adds only an Omega-independent constant at
#   fixed r (Omega delta^2 = r^2 / Z), so these slopes measure the excess, not the W = W* approximation.
#   As registered (symmetric branches): UNRESOLVED. The pattern: the excess rises monotonically across H from ~0 below
#   to > 1 several tube widths up the INCOMING side -- not confined to the inner zone.
