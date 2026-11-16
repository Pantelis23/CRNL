"""§134.13 (T-CASC-ax kill test (b)): does the density excess at the joint saddle track the Sigma ridge, across pairs?

Suspect 4 (§134.12): P(H) collects the near-flat ridge of Sigma (W within ~1e-4 of W*) through Sigma's free inflow into
H, so ln P(H) + Omega W* gains an excess power that depends on how much of the exit is spread along the ridge.

Pair X (E-070 -> E-084, Hill K = 1) has the same stage-2 saddle H = (3, 0.606) but a different Sigma. W* for X by the
heteroclinic BVP: 0.1540536705 (§134.8, B1).

Measured per pair and Omega (200, 400, 800, 1200, 1600; cap 1.1; krsolve; local):
  E(Omega)  = dL/dlnOmega per window, L = ln P(H) + Omega W*           (the excess at H; XK: +0.30 / .34 / .36 / .37)
  R(Omega)  = 1 - (Sigma-cut flux within |x1 - 3| <= 0.05) / total     (ridge occupancy; XK: 0.956 / 0.927 / 0.904 / 0.880)

Criterion, written before running (a single counterexample decides a universality claim; rule 19):
  R1  If X's R at Omega = 1600 differs from XK's 0.880 by a factor >= 2 in (1 - R) (i.e. X's core share is <= 0.06 or
      >= 0.24), the pairs' ridges differ enough to test. Then: X's E within +-0.03 of XK's in EVERY window -> the
      excess does not track the ridge -> suspect 4 is DEAD.
  R2  If X's core share is within a factor 2 of XK's, the test cannot discriminate -> UNRESOLVED, written as such.
  R3  (supporting only) if R1's premise holds and E differs, the direction must be: larger ridge occupancy -> larger E.
      The opposite direction also kills suspect 4.
  My expectation, stated so it can be wrong: X's core share differs by more than 2x and its E differs in R3's direction.

Usage: CRNL_CAP_MULT=1.1 PYTHONPATH=. python experiments/ridge_x.py 200,400,800,1200,1600 kr
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

import experiments.episode_committor as ec
from experiments.sigma_flux import ARMS, cut_flux, element, sigma_table

W_STAR = {"X": 0.1540536705, "XK": 0.1484037485}

if __name__ == "__main__":
    if len(sys.argv) > 2:
        ec.SOLVER = sys.argv[2]
    pair = "X"
    A, B = element(ARMS[pair][0]), element(ARMS[pair][1])
    sig = sigma_table(A, B)
    out_p = Path("results/ridge_x.json")
    out = json.loads(out_p.read_text()) if out_p.exists() else {}
    for om in [int(v) for v in sys.argv[1].split(",")]:
        m = A.cap(om) + 1
        sad1, sad2 = int(np.ceil(A.r2 * om)), int(np.ceil(B.r2 * om))
        x1, f = cut_flux(A, B, om, "sigma", sig, raw=True)            # also leaves the Q-process law cached
        Lq = ec.qprocess_J.last
        P, n1s, base2 = Lq["P"] / ec.PIN_SCALE, Lq["n1s"], Lq["base2"]
        i = int(round(3.0 * om)) - int(n1s[0])
        x2h = B.r2 * om
        j0 = int(math.floor(x2h)) - base2
        fr = x2h - math.floor(x2h)
        lnP = (1 - fr) * math.log(P[i, j0]) + fr * math.log(P[i, j0 + 1])
        core = float(f[np.abs(x1 - 3.0) <= 0.05].sum() / f.sum())
        out[str(om)] = {"L": lnP + om * W_STAR[pair], "core_share": core, "ridge": 1 - core, "total": float(f.sum())}
        print(om, out[str(om)], flush=True)
        out_p.write_text(json.dumps(out, indent=1))
    ks = sorted(out, key=int)
    for a, b in zip(ks, ks[1:]):
        print(f"window {a}->{b}: E = {(out[b]['L'] - out[a]['L']) / math.log(int(b) / int(a)):+.4f}")

# OUTCOME (§134.13): X core share 0.040 / 0.032 / 0.029 / 0.028 / 0.032 (XK 0.044 / 0.073 / 0.096 / 0.120 at
# 400-1600) -> at 1600 XK/X = 3.7x: R1's premise HOLDS. E(X) = +0.888 / +0.580 / +0.393 / +0.359 against
# E(XK) = +0.302 / +0.336 / +0.360 / +0.368: not within +-0.03 in every window -> suspect 4 NOT killed (registered).
# R3: early windows differ in the predicted direction (more ridge -> larger E). My expectation held as worded.
# NOT PREDICTED: the two pairs' E converge from opposite sides onto ~0.36-0.37 by 1200-1600 despite a 3.7x different
# ridge -- a pattern pointing to a LOCAL exponent at H (both pairs share H's linearisation up to c), which suspect 4
# does not predict. Recorded as a pattern, not a verdict.
