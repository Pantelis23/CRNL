"""T10b-iii-c: is the hazard integral's hard stop at theta = 1/3 the remaining error?

§23.4 built a parameter-free hazard integral that reproduced the measured depths to
1.09-1.20x at the large budgets but gave exponent 0.8925 against a measured 0.6474.
§23.5 eliminated the obvious suspect: imposing quasi-staticity on the state by hand
moves the exponent only to 0.7077, 1.6 sigma from doing nothing.

THE SURVIVING SUSPECT. The integral kills the bit the instant `theta` crosses 1/3,
and not as a modelling choice -- as an artifact of how `q(theta)` was measured.
`hazard_at` seeds the state with `initial_counts`, which past `gamma_c` puts zero
separation between X and Y because `delta*` does not exist there. So the seeded state
has already lost the bit and `q = 1` by construction. Meanwhile §23.4 measured that
the real cascade holds a `delta ~ 0.19-0.23` bit through a median 3 of its 7 stages in
exactly that region at the smallest budget, and ~0 stages at the largest.

THE FIX AND WHAT IT COSTS. Past `gamma_c` the state is seeded at an IMPOSED
separation `delta_past` instead of the nonexistent rail, then the stage runs as
before. `delta_past` is not derivable -- it is an empirical input read off §23.4 --
so **this is a one-parameter model and not the absolute test §23.4 ran.** It is
therefore swept, not fitted, and every value is reported (rule 15). Seeding a state
the simulation was OBSERVED to occupy is not a free restoring element: no
separation is created that the chemistry did not already carry, and the conditional
hazard measured from it is a real conditional hazard.

PREDICTION, written before running, and it is quantitative rather than directional.
§23.4's stage counts say the past-`gamma_c` region contributes a median 3 extra
stages out of 7 at Phi/Omega = 25 and 0 out of 46 at Phi/Omega = 400. Restoring
those stages multiplies the small-tank depth by ~10/7 and the large-tank depth by
~1, so it adds `ln(10/7) = 0.357` to the small end of a log-log fit spanning
`ln 16 = 2.77`, i.e. it should LOWER the exponent by about `0.357/2.77 = 0.129`:

  P1  The exponent falls from 0.8925 to roughly **0.76**, at the delta_past that
      matches §23.4's observed 0.19-0.23.
  P2  So the hard stop accounts for about HALF the 0.25 gap and no more, leaving
      ~0.11 unexplained. If instead the exponent lands at or below 0.68, the hard
      stop was the whole error and T10b-iii-c closes.
  P3  The exponent falls monotonically as `delta_past` rises -- a bigger imposed
      bit survives more stages past `gamma_c`, and those stages are worth
      proportionally more to a small tank.
  P4  WHAT IT WOULD MEAN IF P1 FAILS UPWARD. If the exponent barely moves from
      0.8925, then the past-`gamma_c` stages are not what the integral is missing
      either, and the remaining suspect is the survival product itself: it treats
      stages as independent, and a trial that drifts to small delta is more likely
      to lose next stage. Correlated stages die sooner than an independent product
      predicts, which is the wrong sign to explain an integral that already
      OVER-predicts small-tank depth -- so that outcome would leave the gap with no
      suspect at all, which is worth knowing.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.networks.am_fueled import am_fueled, gamma_effective, initial_counts
from crnl.networks.am_reversible import GAMMA_C
from crnl.vectorized import compile_network, gillespie_fast
from experiments.cascade_fuel_vs_noise import _kick, stage_time
from experiments.fuel_hazard import GAMMA0, GAMMA_INF, predict_depth

FUEL_CONCS = [25.0, 50.0, 100.0, 200.0, 400.0]


def seed_past_gc(omega: int, phi: int, w: int, delta_past: float) -> np.ndarray:
    """State at burn fraction w/phi carrying an imposed separation delta_past.

    Mirrors `initial_counts` -- n_B at the gamma_eff attractor, X+Y+B = omega, tank
    exact -- but sets the separation by hand because past gamma_c there is no rail
    to sit on. The value comes from what the cascade was measured to carry (§23.4).
    """
    g = float(gamma_effective(phi - w, w, GAMMA_INF))
    nb = int(round(omega * g / (1.0 + g)))
    rest = omega - nb
    sep = int(round(delta_past * omega))
    if (rest - sep) % 2:
        sep -= 1
    sep = max(sep, 2)
    nx = (rest + sep) // 2
    return np.array([nx, rest - nx, nb, phi - w, w], dtype=np.int64)


def hazard_at(theta: float, omega: int, phi: int, sigma_counts: float,
              t_stage: float, trials: int, seed: int,
              delta_past: float) -> dict:
    w = int(round(theta * phi))
    g = float(gamma_effective(phi - w, w, GAMMA_INF))
    if g >= GAMMA_C:
        n_start = seed_past_gc(omega, phi, w, delta_past)
        imposed = True
    else:
        n_start = initial_counts(omega, phi, w, gamma_inf=GAMMA_INF)
        imposed = False
    if int(n_start[0]) <= int(n_start[1]):
        return {"theta": theta, "q": 1.0, "c": float("nan"),
                "imposed": imposed, "trials": 0}

    comp = compile_network(am_fueled(GAMMA_INF), float(omega))
    rng = np.random.default_rng(seed)
    lost, waste = 0, []
    for _ in range(trials):
        n = _kick(n_start.copy(), sigma_counts, rng)
        nf = gillespie_fast(comp, n, rng, t_max=t_stage).n_final
        if int(nf[0]) <= int(nf[1]):
            lost += 1
        waste.append(int(nf[4]) - int(n_start[4]))
    return {"theta": theta, "q": lost / trials, "c": float(np.mean(waste)),
            "imposed": imposed, "gamma_eff": g, "trials": trials}


def depth_continuous(curve: list[dict], phi: int, theta0: float,
                     max_depth: int = 200000) -> float:
    """Fractional depth where survival crosses 1/2 -- NOT the integer stage.

    `predict_depth` counts whole stages, and that discreteness is the
    integrator's own numerical parameter (rule 13). It has to be removed before
    comparing models: at Phi/Omega = 25 the predicted depth is ~5 stages, so one
    stage is 20% of it, and a swept `delta_past` that genuinely lowers the
    past-gamma_c hazard from 0.429 to 0.302 moved the integer depth not at all --
    four visibly different models returned byte-identical exponents. Interpolating
    the crossing in ln S makes that dependence visible.
    """
    ths = np.array([c["theta"] for c in curve])
    qs = np.array([c["q"] for c in curve])
    cs = np.array([c["c"] for c in curve])
    ok = np.isfinite(cs)
    theta, lnS, d = theta0, 0.0, 0
    hist = [(0, 0.0)]
    target = np.log(0.5)
    while lnS > target and d < max_depth:
        if not ok.any() or theta > ths[ok][-1]:
            break
        c = float(np.interp(theta, ths[ok], cs[ok]))
        if c <= 0:
            break
        q = min(float(np.interp(theta, ths, qs)), 1.0 - 1e-12)
        lnS += np.log1p(-q)
        theta += c / phi
        d += 1
        hist.append((d, lnS))
    for (d0, s0), (d1, s1) in zip(hist, hist[1:]):
        if s1 <= target <= s0:
            return d0 + (s0 - target) / (s0 - s1) * (d1 - d0)
    return float(d)


def fit(x, y) -> tuple[float, float]:
    lx, ly = np.log(x), np.log(y)
    p, c = np.polyfit(lx, ly, 1)
    r = ly - (p * lx + c)
    se = np.sqrt(max(r @ r, 1e-30) / (len(x) - 2) / np.sum((lx - lx.mean()) ** 2))
    return float(p), float(se)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=40)
    ap.add_argument("--sigma-frac", type=float, default=0.03)
    ap.add_argument("--deltas-past", type=float, nargs="+",
                    default=[0.10, 0.15, 0.20, 0.25])
    ap.add_argument("--n-theta", type=int, default=14)
    ap.add_argument("--trials", type=int, default=800)
    ap.add_argument("--seed", type=int, default=20260730)
    ap.add_argument("--measured", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_depth_scaling.json"))
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_hazard_pastgc.json"))
    args = ap.parse_args()

    om, sf = args.omega, args.sigma_frac
    from crnl.networks.am_reversible import delta_star
    sigma_counts = sf * delta_star(GAMMA0) * om
    theta0 = GAMMA0 / (1.0 + GAMMA0)
    thetas = np.linspace(theta0, 0.48, args.n_theta)
    phis = np.array([int(round(fc * om)) for fc in FUEL_CONCS], float)

    meas = json.load(open(args.measured)) if args.measured.exists() else []
    dmeas = []
    for phi in phis:
        g = [r for r in meas if r.get("arm") == "fueled" and r.get("phi") == int(phi)
             and abs(r.get("sigma_frac", -1) - sf) < 1e-9]
        dmeas.append(g[0]["depth_median"] if g else float("nan"))

    t0 = time.time()
    print(f"Omega={om} sigma/delta*={sf} {args.trials} trials/cell")
    print(f"  reference: §23.1 measured 0.6474 +- 0.022; "
          f"§23.4 hard-stop integral 0.8925 +- 0.017; §23.5 rail-reseeded 0.7077")
    print(f"  P1 predicts ~0.76 at delta_past matching §23.4's observed 0.19-0.23\n")

    rows = []
    print(f"{'d_past':>7} {'depths (continuous)':>30} {'ratios vs measured':>26} "
          f"{'exponent':>16} {'integer':>9}")
    for dp in args.deltas_past:
        preds, conts = [], []
        for fc, phi in zip(FUEL_CONCS, phis):
            ts = stage_time(GAMMA0, fc, 2.0)
            curve = [hazard_at(float(t), om, int(phi), sigma_counts, ts,
                               args.trials, args.seed + i * 31 + int(fc), dp)
                     for i, t in enumerate(thetas)]
            preds.append(predict_depth(curve, int(phi), theta0)["depth_pred"])
            conts.append(depth_continuous(curve, int(phi), theta0))
            rows.append({"delta_past": dp, "phi": int(phi), "sigma_frac": sf,
                         "depth_pred": preds[-1], "depth_continuous": conts[-1],
                         "curve": curve})
        p, se = fit(phis, np.array(conts))
        pi, _ = fit(phis, np.array(preds))
        ds = " ".join(f"{v:.2f}" for v in conts)
        rs = " ".join(f"{a/b:.2f}" if b == b else "-" for a, b in zip(conts, dmeas))
        print(f"{dp:>7.2f} {ds:>30} {rs:>26} {p:>9.4f}+-{se:<5.3f} {pi:>9.4f}")

    print(f"\n  measured depths for comparison: "
          + " ".join(f"{v:.0f}" for v in dmeas))
    print(f"  the 'integer' column is why this matters: it is IDENTICAL across a "
          f"2.5x sweep of delta_past")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
