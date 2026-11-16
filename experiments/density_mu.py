"""§134.15 (T-CASC-ax, suspect 5's kill test): does the density excess at the joint saddle move with mu?

Suspect 5 (§134.13-14): ln P(H) + Omega W* carries Omega^beta with beta set by the saddle's linearisation; XK and X (both
mu = 5.208, differing only in c) converge onto E ~ 0.37 (XK: +0.360 / +0.368 at 800-1600, settling at 0.374 by 4000;
X: +0.393 / +0.359).

Pair S = E-070 -> E-070 (Hill K = 1): saddle H = (3, 1.005), lambda_s = 5.686, lambda_u = 1.706, so mu = 3.333. W* by the
heteroclinic BVP: 0.14457152566 (bvp_gamma.py, with A's projection from a sorted Schur form -- S's identical stages make
A's Jacobian defective; the change reproduces X's and XK's W* to every printed digit).

Criterion, written before running (Omega = 200, 400, 800, 1200, 1600; cap 1.1; krsolve; local):
  K1  E_S in the 1200 -> 1600 window within +-0.02 of XK's +0.368 AND |E_S(1200->1600) - E_S(800->1200)| < 0.01
      (settled)                                  -> beta does not move with mu (5.21 -> 3.33): suspect 5 is DEAD.
  K2  E_S(1200->1600) differs from +0.368 by >= 0.05 and the last increment is < 0.02 in size
      (settled elsewhere)                        -> beta moves with mu: suspect 5 SURVIVES.
  Otherwise UNRESOLVED (e.g. not yet settled by 1600), written as such.
  Caveat stated in advance: S's Sigma ridge also differs, and §134.13 showed the ridge changes how E approaches its
  limit; a not-yet-settled S is therefore expected to be possible.
  My expectation, stated so it can be wrong: K2 -- beta moves with mu.

Usage: CRNL_CAP_MULT=1.1 PYTHONPATH=. python experiments/density_mu.py S 200,400,800,1200,1600 kr
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

import experiments.episode_committor as ec
from experiments.sigma_flux import ARMS, element

W_STAR = {"X": 0.1540536705, "XK": 0.1484037485, "S": 0.14457152566}

if __name__ == "__main__":
    pair = sys.argv[1]
    if len(sys.argv) > 3:
        ec.SOLVER = sys.argv[3]
    A, B = element(ARMS[pair][0]), element(ARMS[pair][1])
    out_p = Path(f"results/density_mu_{pair}.json")
    out = json.loads(out_p.read_text()) if out_p.exists() else {}
    for om in [int(v) for v in sys.argv[2].split(",")]:
        m = A.cap(om) + 1
        sad1, sad2 = int(np.ceil(A.r2 * om)), int(np.ceil(B.r2 * om))
        _, res = ec.qprocess_J(A, B, om, m, sad1, sad2, lower="rail")
        L = ec.qprocess_J.last
        P, n1s, base2 = L["P"] / ec.PIN_SCALE, L["n1s"], L["base2"]
        i = int(round(3.0 * om)) - int(n1s[0])
        x2h = B.r2 * om
        j0 = int(math.floor(x2h)) - base2
        f = x2h - math.floor(x2h)
        lnP = (1 - f) * math.log(P[i, j0]) + f * math.log(P[i, j0 + 1])
        out[str(om)] = {"L": lnP + om * W_STAR[pair], "lnP_H": lnP, "residual": float(res)}
        print(om, out[str(om)], flush=True)
        out_p.write_text(json.dumps(out, indent=1))
    ks = sorted(out, key=int)
    for a, b in zip(ks, ks[1:]):
        print(f"window {a}->{b}: E = {(out[b]['L'] - out[a]['L']) / math.log(int(b) / int(a)):+.4f}")

# OUTCOME (§134.15), S: L = -0.227 / -0.813 / -0.778 / -0.669 / -0.584; E_S = -0.845 / +0.050 / +0.270 / +0.293.
# K1 FAILS (1200->1600 is 0.075 below XK's 0.368). K2 FAILS narrowly (difference >= 0.05, but the last increment is
# 0.023, not < 0.02: S has not settled). Verdict UNRESOLVED as written. S starts far lower and is still climbing at 1600
# -- the not-yet-settled case flagged in advance. Needs Omega = 2800 and 4000.
