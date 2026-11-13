"""§134 -- why the total's prefactor stays flat (a ~ 0) while the exit on Sigma widens: the naive degenerate model.

§133: on Sigma the committed exit spreads over x1 ~ 1.2 - 3 and widens with Omega, and the action along Sigma slopes
toward the saddle, eps(x1) ~ 1e-4 per unit x1 (ray family). §124.10 / §124.14 measured the total
R ~ C Omega^a e^{Gamma Omega} with a = -0.012 (800 -> 1600) and -0.010 (1600 -> 2000).

THE MODEL (no free parameter after one calibration). Each point of Sigma carries flux
    F(x1; Omega) = C(x1) Omega^{1/2} e^{(Gamma - eps(x1)) Omega}
  - the local Omega^{1/2} is the naive degenerate count: a normal crossing at each point, with the density along
    Sigma flat to leading order;
  - eps(x1) comes from the rays (§133), interpolated on the lower envelope;
  - C(x1) is calibrated from the Sigma-cut density at ONE Omega (800) by inverting the formula.
Then F_tot(Omega) = sum over bins is a prediction at every other Omega, and so is the total's local a.

PREDICTIONS, WRITTEN BEFORE RUNNING.
  M1  The model's a over 800 -> 1600 lies within 0.10 of the measured total's a over the same window (from §133's
      Sigma-cut totals, which equal the committed fall rate). If the model gives a >= +0.25, the local Omega^{1/2}
      on Sigma is wrong, and a ~ 0 is not explained by a spread exit with this eps.
  M2  The same over 1600 -> 2000 (Sigma cut at Omega = 2000, cap 1.1, laptop), within 0.10.
  M3  The model's a drifts DOWN with Omega: it is lower over 1600 -> 2000 than over 800 -> 1600, because e^{-Omega eps}
      progressively removes the low-x1 shoulder. The asymptote at Omega >> 1e4 depends on eps near the saddle, which
      this model does not resolve.

ADDED AFTER M1-M3 WERE JUDGED, BEFORE OMEGA = 2800 RAN. M1 held (model +0.061, measured -0.013, 0.074 apart).
M2 held (+0.026 against -0.012, 0.038 apart). M3 held (the model drifts down). The model's next windows:
+0.0136 (2000 -> 2800), -0.0008 (2800 -> 4000), -0.018 (4000 -> 8000), -0.038 (8000 -> 16000).
  M4  FORWARD. The measured a over 2000 -> 2800 (Sigma cut, cap 1.1, laptop) lies within 0.05 of the model's +0.0136,
      and the model-minus-measured offset is below the previous window's 0.038.

OUTCOMES (appended after running; predictions above left as written).
  M1 HELD: model +0.061, measured -0.0127 (800 -> 1600).
  M2 HELD: model +0.026, measured -0.0118 (1600 -> 2000). The cap 1.1 Sigma-cut totals reproduce cap 1.25 at 1600
     to 5.9e-8.
  M3 HELD: the model's a drifts down.
  M4 HELD: measured -0.0073 over 2000 -> 2800, against the model's +0.0136 (0.021 apart). The offsets run 0.074 /
     0.038 / 0.021: the model and the exact chain converge. Omega = 2800 is the first exact cell past 2000
     (PARDISO, 13 min per cut; the two cut totals agree to 2e-13). FINDINGS §134.

ADDED AFTER OMEGA = 4000 FAILED, BEFORE OMEGA = 3400 RAN. T-CASC-au's kill test asked for Omega = 4000. PARDISO
failed there with error -2 (insufficient memory), in-core and with out-of-core requested, on the 64 GB laptop.
Stand-in, weaker than the registered test because the window is shorter: Omega = 3400.
  M5  The measured a over 2800 -> 3400 lies within 0.03 of the model's +0.0024 (the model's value, computed before
      the run), and the offset is below 0.021. The registered asymptotic question (does the measured a turn
      down?) stays open unless the measured a here drops below -0.0073.
  -> M5 HELD: measured -0.0026 against +0.0024 (offset 0.005). The offsets run 0.074 / 0.038 / 0.021 / 0.005. The
     measured a is still rising (-0.0127 / -0.0118 / -0.0073 / -0.0026), so the turn-down stays open. On Sigma the
     spread narrows: 0.704 / 0.692 / 0.679 at 2000 / 2800 / 3400 (29 min per cut on the laptop).

ADDED 2026-11-13, WHILE OMEGA = 4000 WAS RUNNING (out-of-core PARDISO, 16 GB in-core cap), BEFORE IT RETURNED.
  M6  The unchanged model (the one M1-M5 tested) gives a = -0.0047 over 3400 -> 4000. The measured a lies within
      0.03 of that, and |offset| <= 0.01.
  An exact bound, computed today: at the joint saddle (3, 0.606), the Lyapunov form gives h = 0.1500 along Sigma's
  tangent (equal to stage 1's own well curvature 2|f1'|/D11). Climbing back along Sigma after reaching the saddle
  is a feasible path, so eps(3 - delta) <= 0.075 delta^2 + O(delta^3). The ray envelope EXCEEDS this bound for
  delta < ~0.012 (9.6e-6 against 6.1e-6 at delta = 0.009), so the envelope is not W there. The model's
  extrapolation to negative a (-0.018 at 4000 -> 8000, -0.046 at 16000 -> 64000) rests on that region, and does
  not hold asymptotically. A first attempt to substitute the bound for eps as an estimate (a -> -0.59) was wrong:
  a bound may only be used as a bound. It was discarded.
  -> M6 NOT RUN. Omega = 4000 failed a third time with PARDISO -2, with out-of-core mode confirmed requested
     (pypardiso's set_iparm is 1-based, so iparm(60) = 2 reached MKL) and a 16 GB in-core cap. The test stands,
     written, for a larger machine. FINDINGS §134.2.
  (M6 was later run with krsolve, out-of-core: HELD, +0.0009 against -0.0047. FINDINGS §134.5.)

ADDED 2026-11-13 AFTER OMEGA = 4000, BEFORE OMEGA = 5000 RAN.
  M7  The unchanged model gives a = -0.0114 over 4000 -> 5000. The verdict depends on where the measured a lands:
        a_measured >= 0                        -> the model's negative asymptotics are REFUTED on data;
        |a_measured - (-0.0114)| <= 0.01       -> the model SURVIVES this window;
        -0.0014 < a_measured < 0               -> UNRESOLVED (between the two), and written as such.
      (The bands cannot overlap: -0.0114 + 0.01 = -0.0014.)
  M8  On Sigma the exit keeps narrowing: sd(x1) at Omega = 5000 is below 4000's 0.665.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

GAMMA = 0.00734


def eps_fn(s):
    r = np.array(s["ray_profile_XK"]["rows"])
    d = json.loads(Path("results/joint_action.json").read_text())
    ws = d["V1_s1"]["E-070"] - d["XK"]["Gamma"]
    x, S = r[:, 1], r[:, 3] - ws
    xe = np.linspace(1.15, 2.995, 400)
    env = np.array([S[(x >= a - 0.01) & (x <= a + 0.01)].min() if ((x >= a - 0.01) & (x <= a + 0.01)).any()
                    else np.nan for a in xe])
    ok = np.isfinite(env)
    return lambda z: np.interp(z, xe[ok], env[ok], right=0.0)          # beyond the saddle: eps -> 0


def main():
    s = json.loads(Path("results/sigma_flux.json").read_text())
    eps = eps_fn(s)
    bins = np.array(s["XK:800:sigma"]["flux_bins_x1"])
    F800 = np.array(s["XK:800:sigma"]["flux_density"])
    C = F800 / (800 ** 0.5 * np.exp((GAMMA - eps(bins)) * 800))
    model = lambda om: float((C * om ** 0.5 * np.exp((GAMMA - eps(bins)) * om)).sum())
    out = {}
    for (o1, o2), suf in (((800, 1600), ""), ((1600, 2000), ":cap1.1")):      # each window at one cap
        k1, k2 = f"XK:{o1}:sigma{suf}", f"XK:{o2}:sigma{suf}"
        a_mod = (math.log(model(o2) / model(o1)) - GAMMA * (o2 - o1)) / math.log(o2 / o1)
        rec = {"a_model": a_mod}
        if k1 in s and k2 in s:
            rec["a_measured"] = (math.log(s[k2]["total"] / s[k1]["total"]) - GAMMA * (o2 - o1)) / math.log(o2 / o1)
        out[f"{o1}-{o2}"] = rec
        print(o1, o2, rec, flush=True)
    s["S134"] = out
    Path("results/sigma_flux.json").write_text(json.dumps(s, indent=2))


if __name__ == "__main__":
    main()
