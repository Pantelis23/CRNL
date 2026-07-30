"""T10b-iii-d: is the hazard integral's residue survivor selection?

§23.6 closed 39.6% of the integral's exponent error by removing its hard stop at
`theta = 1/3`, leaving +0.1445 (4.6 sigma) and a SHAPE error across the budget --
predicted/measured 0.78, 0.72, 0.88, 1.00, 1.10 -- so no normalisation fixes both
ends. The named suspect: the survival product is mean-field. It applies the
UNCONDITIONAL hazard `q(theta)` at every stage, while among trials that have already
survived several stages past `gamma_c` the separation is selected upward -- survivors
are the ones that happened to keep a big bit, and they face a lower hazard than `q`
says. The bias grows with the number of past-`gamma_c` stages, 3 of 7 at the smallest
tank and 0 of 46 at the largest, i.e. exactly where the integral under-predicts.

This measures the conditional hazard directly. For each trial the stages with
`gamma_eff >= gamma_c` are indexed k = 1, 2, 3, ...; the conditional hazard at k is
P(lost during the k-th such stage | survived k-1 of them), and the carried separation
at entry to each is recorded alongside.

THE CONFOUND, and the control for it. `gamma_eff` keeps RISING with k because the
tank keeps draining, so the landscape is getting worse at the same time the survivor
pool is getting fitter. A raw hazard that rises with k therefore proves nothing. The
control is to divide by the unconditional `q(theta)` measured at the same theta in
§23.6 (`delta_past = 0.20`, the imposed-separation curve the integral actually used).
That RATIO is the mean-field error, and it is what has to move with k.

PREDICTIONS, written before running:

  P1  The conditional/unconditional hazard ratio is BELOW 1 and FALLS with k. Below
      1 because the integral's imposed `delta_past = 0.20` is the median over all
      past-`gamma_c` stages while survivors sit above that median; falling because
      the selection compounds with every stage survived.
  P2  The carried separation conditional on survival RISES with k, which is the
      direct signature and does not depend on the unconditional curve being right.
  P3  The effect is large enough to matter: a ratio reaching ~0.6-0.7 by k = 3.
      Anything above ~0.9 throughout means selection exists but is too small to be
      the 0.145 residue, and the independence assumption itself becomes the suspect.
  P4  WHAT WOULD KILL IT. If the ratio is flat in k, or rises, survivor selection is
      not the residue. Then the honest position is that the 0.145 has no surviving
      suspect, and §23.6's closing note applies -- a non-Markovian repair is a
      different object from the integral §23.4 built, so the integral should be left
      where it is rather than patched further.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.networks.am_fueled import am_fueled, gamma_effective, initial_counts
from crnl.networks.am_reversible import GAMMA_C, delta_star
from crnl.vectorized import compile_network, gillespie_fast
from experiments.cascade_fuel_vs_noise import _kick, stage_time

GAMMA0 = 0.3
GAMMA_INF = 1.0


def run(omega: int, fuel_conc: float, sigma_frac: float, trials: int,
        seed: int, max_depth: int = 600) -> dict:
    sigma_counts = sigma_frac * delta_star(GAMMA0) * omega
    phi = int(round(fuel_conc * omega))
    w0 = int(round(phi * GAMMA0 / (1.0 + GAMMA0)))
    ts = stage_time(GAMMA0, fuel_conc, 2.0)
    comp = compile_network(am_fueled(GAMMA_INF), float(omega))
    n_start = initial_counts(omega, phi, w0, gamma_inf=GAMMA_INF)
    rng = np.random.default_rng(seed)

    # per past-gamma_c index k: entered, lost, and the (theta, delta) on entry
    entered: dict[int, int] = {}
    lost: dict[int, int] = {}
    thetas: dict[int, list] = {}
    deltas: dict[int, list] = {}
    for _ in range(trials):
        n = n_start.copy()
        k = 0
        for _d in range(1, max_depth + 1):
            w = int(n[4])
            g = float(gamma_effective(phi - w, w, GAMMA_INF))
            past = g >= GAMMA_C
            if past:
                k += 1
                entered[k] = entered.get(k, 0) + 1
                thetas.setdefault(k, []).append(w / phi)
                deltas.setdefault(k, []).append((int(n[0]) - int(n[1])) / omega)
            n = _kick(n, sigma_counts, rng)
            n = gillespie_fast(comp, n, rng, t_max=ts).n_final
            if int(n[0]) <= int(n[1]):
                if past:
                    lost[k] = lost.get(k, 0) + 1
                break
    ks = sorted(entered)
    return {"omega": omega, "fuel_conc": fuel_conc, "phi": phi, "t_stage": ts,
            "sigma_frac": sigma_frac, "trials": trials,
            "by_k": [{"k": k, "entered": entered[k], "lost": lost.get(k, 0),
                      "hazard": lost.get(k, 0) / entered[k],
                      "theta_median": float(np.median(thetas[k])),
                      "delta_median": float(np.median(deltas[k])),
                      "delta_mean": float(np.mean(deltas[k]))} for k in ks]}


def unconditional_q(pastgc_json: pathlib.Path, delta_past: float, phi: int):
    """q(theta) as §23.6's integral used it, for the same Phi."""
    if not pastgc_json.exists():
        return None
    d = json.load(open(pastgc_json))
    rs = [r for r in d if abs(r["delta_past"] - delta_past) < 1e-9
          and r["phi"] == phi]
    if not rs:
        return None
    cur = rs[0]["curve"]
    th = np.array([c["theta"] for c in cur])
    q = np.array([c["q"] for c in cur])
    return lambda t: float(np.interp(t, th, q))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=40)
    ap.add_argument("--sigma-frac", type=float, default=0.03)
    ap.add_argument("--fuel-concs", type=float, nargs="+", default=[25.0, 50.0, 100.0])
    ap.add_argument("--trials", type=int, default=3000)
    ap.add_argument("--delta-past", type=float, default=0.20)
    ap.add_argument("--seed", type=int, default=20260730)
    ap.add_argument("--pastgc", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_hazard_pastgc.json"))
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/fuel_survivor_bias.json"))
    args = ap.parse_args()

    t0 = time.time()
    print(f"Omega={args.omega} sigma/delta*={args.sigma_frac} "
          f"{args.trials} trials   gamma_c={GAMMA_C}")
    print(f"unconditional q from §23.6's delta_past={args.delta_past} curve\n")
    rows = []
    for fc in args.fuel_concs:
        r = run(args.omega, fc, args.sigma_frac, args.trials, args.seed + int(fc))
        uq = unconditional_q(args.pastgc, args.delta_past, r["phi"])
        print(f"=== Phi={r['phi']}  (Phi/Omega={fc:.0f})")
        print(f"{'k':>3} {'entered':>8} {'lost':>6} {'hazard':>8} {'+-':>7} "
              f"{'theta':>7} {'q_uncond':>9} {'ratio':>7} {'delta_med':>10}")
        for b in r["by_k"]:
            se = np.sqrt(max(b["hazard"] * (1 - b["hazard"]), 1e-12) / b["entered"])
            qu = uq(b["theta_median"]) if uq else float("nan")
            ratio = b["hazard"] / qu if qu and qu == qu else float("nan")
            b["q_uncond"] = qu
            b["ratio"] = ratio
            b["hazard_se"] = float(se)
            print(f"{b['k']:>3} {b['entered']:>8} {b['lost']:>6} "
                  f"{b['hazard']:>8.4f} {se:>7.4f} {b['theta_median']:>7.4f} "
                  f"{qu:>9.4f} {ratio:>7.3f} {b['delta_median']:>10.4f}")
        rows.append(r)
        print()

    print("=== P2: does the carried separation rise with k? (survivor selection's "
          "direct signature)")
    for r in rows:
        ds = [b["delta_mean"] for b in r["by_k"] if b["entered"] >= 30]
        print(f"  Phi={r['phi']}: " + " -> ".join(f"{v:.4f}" for v in ds))

    print("\n=== P1/P3: conditional / unconditional hazard by k")
    for r in rows:
        rs = [b["ratio"] for b in r["by_k"] if b["entered"] >= 30]
        print(f"  Phi={r['phi']}: " + " -> ".join(f"{v:.3f}" for v in rs))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
