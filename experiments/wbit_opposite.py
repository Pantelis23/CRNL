"""§132 -- T-WBIT-a⁗: §131's c-dependence, is it c or h? The opposite sweep.

§131 held b2 fixed and moved c, so the healthy drive h moved with it. Here c is held at 0.34 and h is moved
(0.14, 0.16, 0.18, 0.20). That moves b2 = h - c (2 HI1 - 1) and so Gamma = -S2(b2), which is exact (§130), and
which is removed by looking at the residual
    rho(h) = ln R + S2(b2) / sigma^2.

THE WINDOW SUSPECT'S PREDICTION. ln R = -S2/sigma^2 - k ln(2 c S2'(b2) kappa / sigma) + const, so at fixed c the
residual moves only through the window width's S2'(b2):
    rho(0.20) - rho(0.14) = -k [ln S2'(b2(0.20)) - ln S2'(b2(0.14))] = -k x 0.362,
with S2' = 1.234 / 1.437 / 1.613 / 1.772 at h = 0.14 / 0.16 / 0.18 / 0.20 (computed beforehand by finite
differences of the exact 1-D action). k is taken from §131, which is independent data: 7.06 (sigma = 0.15) and
7.69 (0.12). The prediction is a drop of 2.56 and 2.78 nats.

PREDICTIONS, WRITTEN BEFORE RUNNING.
  O0  CLEAN. Cell 2's own route (positive 1-D rate at drive h over the dip rate, §131's method) is below 0.1 R at
      every cell. Any cell failing this is excluded before the residuals are read.
  O1  THE RESIDUAL FALLS BY THE WINDOW'S AMOUNT. rho(0.20) - rho(0.14) is negative, with magnitude within a factor
      2 of k x 0.362: [1.28, 5.12] nats at sigma = 0.15 and [1.39, 5.57] at 0.12.
      Verdicts:
        inside the band              -> the window suspect SURVIVES the opposite sweep;
        |drop| < 0.5 nat             -> KILLED as stated: the prefactor does not follow S2', and §131's
                                        dependence is on c itself through something else;
        positive, or beyond 2x       -> KILLED: h acts through something other than the window width.
  O2  MONOTONE. rho decreases at each step in h at both sigma.

OUTCOMES (appended after running; predictions above left as written).
  O0 HELD: own/R <= 6.6e-4 (sigma = 0.15) and <= 7.9e-7 (0.12) at every cell.
  O1 THE WINDOW SUSPECT IS KILLED AS STATED. The drop is -1.117 (sigma = 0.15) and -1.205 (0.12) against predicted
     -2.554 / -2.781. Both are outside the factor-2 bands [1.28, 5.12] / [1.39, 5.57], on the low side: 0.437 and
     0.433 of the prediction. The direction is right and the size is not. In the window's form, c and S2' enter only
     as the product c S2' with a single k; the data give c an exponent of about 7 (§131) and S2' an effective
     exponent of about 3.1-3.3.
  O2 HELD: steps -0.449 / -0.386 / -0.282 (0.15), -0.471 / -0.415 / -0.320 (0.12). FINDINGS §132.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from experiments.wbit_action import escape_action
from experiments.wbit_cascade import DC, HI1
from experiments.wbit_qsd import exits_pos
from experiments.wbit_window import rate_pos

C = 0.34
HS = (0.14, 0.16, 0.18, 0.20)
K131 = {0.15: 7.058809534803738, 0.12: 7.685718595949311}
DLNS2P = 0.5721 - 0.2103


def b2_of(h):
    return h - C * (2 * HI1 - 1)


def main():
    out = Path("results/wbit_opposite.json")
    rep = json.loads(out.read_text()) if out.exists() else {}
    for s in (0.15, 0.12):
        for h in HS:
            key = f"{h}:{s}"
            if key in rep:
                continue
            b2 = b2_of(h)
            r = exits_pos(320, s, C, b2, h=h)
            S2 = escape_action(b2)
            own = rate_pos(3200, s, h) / rate_pos(3200, s, 0.0) / r["R"]
            rep[key] = {"h": h, "b2": b2, "edge": HI1 - (DC + h) / (2 * C), "R": r["R"], "S2_b2": S2,
                        "rho": float(np.log(r["R"]) + S2 / s ** 2), "own_over_R": own, "mass_check": r["mass_check"]}
            print(key, rep[key], flush=True)
            out.write_text(json.dumps(rep, indent=2))
        clean = [h for h in HS if rep[f"{h}:{s}"]["own_over_R"] < 0.1]
        rho = [rep[f"{h}:{s}"]["rho"] for h in HS]
        rep[f"verdict:{s}"] = {"clean": clean, "drop": rho[-1] - rho[0], "predicted": -K131[s] * DLNS2P,
                               "steps": [rho[i + 1] - rho[i] for i in range(3)]}
        print("verdict", s, rep[f"verdict:{s}"], flush=True)
        out.write_text(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main()
