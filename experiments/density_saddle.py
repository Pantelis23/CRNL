"""§134.12: which half of the exit carries the anomalous Omega-power -- the DENSITY at the joint saddle, or the flux /
committor part? A split of §134.11's conflict (measured core b = 0.72-0.80 against Maier-Stein's 1/2).

WKB for the stationary law pinned at the rails A = (r3, r3) (pin = 1 there, as qprocess_J does): P(x) ~ K(x) e^{-Omega
(W(x) - W(A))} with K(A) = 1, so at the joint saddle H = (3, 0.606) the combination
    L(Omega) = ln P(H) + Omega W*,   W* = 0.1484037485 (BVP, §134.8)
must tend to ln K(H), a constant -- no power of Omega -- if the density is WKB-regular at H (Maier-Stein, mu > 1:
K(H) finite). The Q-process law is psi * (original law); psi -> 1 near x1 = 3, so the read is of the original law.

Predictions, written before running (XK, cap 1.1, lower = "rail", Omega = 200, 400, 800, 1200, 1600, krsolve):
  D1  dL/dlnOmega per window tends to 0 (read the sequence, rule 20; rule 21: windows reported).
      If instead it sits near +0.25 (the core's excess) and does not fall, the anomaly is in the DENSITY at H;
      if it tends to 0, the anomaly is in the flux / committor half, and the density is regular.
  D2  (instrument) P at the pin equals PIN_SCALE exactly, and the Q-process residual is <= 1e-10.
  My expectation: D1 tends to 0 -- the anomaly sits in the committed-flux half (the boundary layer of q near H).

Usage: CRNL_CAP_MULT=1.1 PYTHONPATH=. python experiments/density_saddle.py 200,400,800,1200,1600 kr
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

import experiments.episode_committor as ec
from experiments.sigma_flux import ARMS, element

W_STAR = 0.1484037485

if __name__ == "__main__":
    if len(sys.argv) > 2:
        ec.SOLVER = sys.argv[2]
    A, B = element(ARMS["XK"][0]), element(ARMS["XK"][1])
    out_p = Path("results/density_saddle.json")
    out = json.loads(out_p.read_text()) if out_p.exists() else {}
    for om in [int(v) for v in sys.argv[1].split(",")]:
        m = A.cap(om) + 1
        sad1, sad2 = int(np.ceil(A.r2 * om)), int(np.ceil(B.r2 * om))
        _, res = ec.qprocess_J(A, B, om, m, sad1, sad2, lower="rail")
        L = ec.qprocess_J.last
        P, n1s, base2 = L["P"] / ec.PIN_SCALE, L["n1s"], L["base2"]
        i = int(round(3.0 * om)) - int(n1s[0])
        x2h = B.r2 * om                                            # the saddle's stage-2 coordinate, in molecules
        j0 = int(math.floor(x2h)) - base2
        f = x2h - math.floor(x2h)
        lnP = (1 - f) * math.log(P[i, j0]) + f * math.log(P[i, j0 + 1])     # log-linear between the two sites
        k1 = int(np.clip(round(A.r3 * om) - sad1, 0, P.shape[0] - 1))
        k2 = int(np.clip(round(B.r3 * om) - base2, 0, P.shape[1] - 1))
        rec = {"lnP_H": lnP, "L": lnP + om * W_STAR, "pin": float(P[k1, k2]), "residual": float(res),
               "lnP_sites": [math.log(P[i, j0]), math.log(P[i, j0 + 1])]}
        out[str(om)] = rec
        print(om, rec, flush=True)
        out_p.write_text(json.dumps(out, indent=1))
    ks = sorted(out, key=int)
    for a, b in zip(ks, ks[1:]):
        print(f"window {a}->{b}: dL/dlnOmega = {(out[b]['L'] - out[a]['L']) / math.log(int(b) / int(a)):+.4f}")

# OUTCOME (§134.12), Omega = 200 / 400 / 800 / 1200 / 1600:
#   L = -0.549 / -0.340 / -0.107 / +0.039 / +0.145; dL/dlnOmega = +0.302 / +0.336 / +0.360 / +0.368 per window.
#   D1 FAILED as the regular-density prediction: the slope does not fall toward 0, it rises toward ~0.37-0.4. The
#   anomaly is in the DENSITY at H. My expectation (flux half) was WRONG. Not a W* error (that would grow linearly in
#   Omega; these saturate).
#   D2 as written FAILED on the residual (2.9e-10 / 8.0e-10 at 800-1600 vs a fixed 1e-10 -- rule 20: a fixed tolerance on
#   a solve whose relative residual grows mildly with size; pin exact at 1.0 throughout).
#
# ADDED BEFORE THE HARNESS CHECK (D3): the same L at points ALONG the BVP orbit, L(x) = ln P(x) + Omega W(x) with W(x)
# the orbit's own action. If L's slope at a mid-path point (far from H) is also ~0.3, the Omega-power is a
# normalisation artefact (pin / Doob psi), not physics at H. Prediction: at the mid-path point the slope tends to 0,
# and it grows as the point approaches H -- the anomaly is local to the saddle.
# -> D3 HELD (harness clean, anomaly local to H): along the orbit L converges to a constant ln K(x) at every point to 3
#    decimals by Omega = 800 (last-window slopes +0.003 / +0.001 / +0.000 / -0.001 / -0.005 / -0.009 / -0.013 going
#    toward H); the point (2.68, 0.62) just short of H relaxes with slope -0.33 -> -0.055. Only AT H does the density
#    carry a growing positive power (+0.37). results/density_along_path.json; experiments/density_along_path.py.
#
# ADDED BEFORE THE LAPTOP CELLS (Omega = 1600 check, 2800, 4000; §134.14), the reading of T-CASC-ax kill test (a):
#   1600 must reproduce the local L = +0.145193 to <= 1e-6 (laptop route check).
#   E(1600->2800) or E(2800->4000) > 0.5             -> suspect 4 (ridge collection, excess <= 1/2) is DEAD.
#   E falling in both windows, below 0.35 by 2800->4000 -> supports suspect 4 (the turn-down it needs).
#   E within [0.35, 0.42] in both windows, no trend    -> supports suspect 5 (a settled local exponent).
#   Anything else -> UNRESOLVED, written as such. None of these is decisive alone: 4 needs an eventual decline it may
#   not show by 4000, and 5's "settled" is a trend read over two windows (rule 21: windows reported).
# -> §134.14 (laptop): 1600 reproduced bit-for-bit (L = +0.14519328219140). L(2800) = +0.353721, L(4000) = +0.487019.
#    E(1600->2800) = +0.3726, E(2800->4000) = +0.3737; increments +0.0046 then +0.0011 -> flattening near 0.374.
#    Registered reading: both windows in [0.35, 0.42] with a vanishing trend -> SUPPORTS suspect 5 (a settled local
#    exponent); suspect 4 not dead (<= 0.5) but shows no turn-down by 4000. Residuals 1.5e-9 / 3.7e-9 (D2's fixed 1e-10
#    again fails as written; rule 20).
