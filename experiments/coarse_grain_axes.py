"""Two axes of coarse-graining: keeping the NOISE vs keeping the COORDINATE

THE BIGGER PICTURE THIS TESTS. §21 measured the project's central claim at the
trajectory level: against an exact CME reference, every approximation that keeps ANY
noise recovers the restoration error exponent to 2-12%, while the deterministic ODE
is categorically wrong. The headline was "noise is what matters."

§23 built a different approximation -- a stage-level hazard integral -- and spent six
subsections failing to repair it. The decomposition (§23.10) says its error is 56%
inter-stage MEMORY and 17% theta-dispersion. But that model is not merely
deterministic: it also **collapses the state**, throwing away `delta` and keeping only
the burn fraction `theta`. §21's levels (CLE, tau-leaping) keep every species, so they
keep both the noise AND the coordinates.

So §21's cliff may not be about noise at all -- it may be about what a level RETAINS,
with noise the thing those levels happened to retain. That predicts something sharp
and previously untested:

  **A model that keeps `delta` but is deterministic in `theta` should beat a model
  that keeps `theta`-noise but discards `delta`.** If it does, "keep the noise" is the
  wrong lesson and "keep the coordinate that carries the memory" is the right one --
  and §21's result is the special case where those coincide.

WHAT THIS IS NOT, established by a smoke run BEFORE the real one. The collapsed cells
here are NOT §23.4's integral or §23.9's ensemble: the collapsed cell gives 82.9
stages at Phi/Omega = 400 where §23.4's integral gave 48.2. The collapse here forces
`delta` onto the mean the FULL model occupies, a more favourable choice than §23.4's
rail seeding, so absolute depths do not compare across sections. **This experiment can
support the ORDERING between two collapses of one kernel, not the absolute recovery of
any of them.** An earlier draft of this docstring advertised the collapsed cells as a
built-in control on §23.4; they are not, and P2 below is correspondingly weak.

THE DESIGN is a 2x2 on ONE measured kernel, so nothing differs between cells except
what the propagation retains. Measure `K(theta, delta)` -- per-stage loss probability,
the distribution of outgoing `delta`, and the burn -- by seeding states at each grid
point. Then propagate four ways:

    delta kept    x  theta stochastic      <- retains everything the kernel has
    delta kept    x  theta deterministic
    delta collapsed x theta stochastic     <- this is §23.9's ensemble
    delta collapsed x theta deterministic  <- this is §23.4's integral

Collapsing `delta` means replacing it, at each `theta`, by the mean `delta` the FULL
model actually carries there -- a mean-field collapse of the same data, not a
different measurement.

PREDICTIONS, written before running:

  P1  THE ONE THE SYNTHESIS RIDES ON. Keeping `delta` while propagating `theta`
      deterministically recovers MORE of the gap than keeping `theta`-noise while
      collapsing `delta`. Quantitatively: §23.10 puts memory at 56% and dispersion at
      17%, and `delta` is the coordinate memory lives in, so I predict the
      delta-kept/theta-deterministic cell lands nearer the plain simulation
      (0.6325) than the delta-collapsed/theta-stochastic cell's 0.7641.
  P2  The full cell (both retained) reaches ~0.63-0.68, i.e. close to the plain
      simulation, because between them the two coordinates are most of the state
      that matters. It cannot be exact -- the kernel is measured on a lattice in
      `delta` and interpolated in `theta`.
  P3  The two axes are NOT additive. Every correction in this thread has partly
      cancelled against another (§22.4, §23.8, §23.10's anti-correlation), so I
      expect keeping both to recover less than the sum of keeping each.
  P4  WHAT KILLS THE SYNTHESIS. If the delta-kept/theta-deterministic cell sits at
      ~0.79 -- no better than the integral -- then `delta` is not where the memory
      lives, §23.10's 56% is not attributable to the collapsed coordinate, and the
      reframing of §21 has no support. The honest fallback would be that §21's
      "noise is what matters" stands unqualified and §23's integral simply fails for
      reasons specific to it.
  P5  A WEAKER FAILURE that must not be read as success: if the full cell reaches
      0.63 but the delta-kept/theta-deterministic cell does NOT beat
      delta-collapsed/theta-stochastic, then both axes matter but the ordering claim
      -- coordinate over noise -- is unsupported, and only the unsurprising half of
      the synthesis survives.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.networks.am_fueled import am_fueled, gamma_effective
from crnl.networks.am_reversible import delta_star
from crnl.vectorized import compile_network, gillespie_fast
from experiments.cascade_fuel_vs_noise import _kick, stage_time
from experiments.fuel_hazard import GAMMA0, GAMMA_INF
from experiments.fuel_hazard_pastgc import FUEL_CONCS, fit

DELTAS = np.array([0.05, 0.125, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70])


def seed(omega: int, phi: int, w: int, delta: float) -> np.ndarray:
    """State at burn fraction w/phi carrying separation `delta`. Tank exact."""
    g = float(gamma_effective(phi - w, w, GAMMA_INF))
    nb = int(round(omega * min(g, 0.999) / (1.0 + min(g, 0.999))))
    rest = omega - nb
    sep = int(round(delta * omega))
    if (rest - sep) % 2:
        sep -= 1
    sep = max(min(sep, rest), 2)
    nx = (rest + sep) // 2
    return np.array([nx, rest - nx, nb, phi - w, w], dtype=np.int64)


def kernel_cell(theta: float, delta: float, omega: int, phi: int,
                sigma_counts: float, t_stage: float, trials: int,
                seed_i: int) -> dict:
    """q, outgoing-delta distribution on DELTAS, and burn mean/sd."""
    w = int(round(theta * phi))
    n0 = seed(omega, phi, w, delta)
    comp = compile_network(am_fueled(GAMMA_INF), float(omega))
    rng = np.random.default_rng(seed_i)
    lost, burns, outs = 0, [], []
    for _ in range(trials):
        n = _kick(n0.copy(), sigma_counts, rng)
        nf = gillespie_fast(comp, n, rng, t_max=t_stage).n_final
        burns.append(int(nf[4]) - int(n0[4]))
        if int(nf[0]) <= int(nf[1]):
            lost += 1
        else:
            outs.append((int(nf[0]) - int(nf[1])) / omega)
    b = np.array(burns, float)
    if outs:
        idx = np.abs(np.asarray(outs)[:, None] - DELTAS[None, :]).argmin(axis=1)
        hist = np.bincount(idx, minlength=len(DELTAS)).astype(float)
        hist /= hist.sum()
    else:
        hist = np.zeros(len(DELTAS))
    return {"theta": theta, "delta": delta, "q": lost / trials,
            "c": float(b.mean()), "c_sd": float(b.std(ddof=1)),
            "delta_out": hist.tolist(), "n_surv": len(outs)}


def propagate(K, thetas, phi, theta0, keep_delta, stoch_theta, n_traj, rng,
              max_depth=4000, dbar_idx=None, record_occupancy=False):
    """K[i, j] over (theta_i, delta_j). Returns the continuous median depth.

    `keep_delta=False` is the MEAN-FIELD COLLAPSE of the same kernel: at each
    `theta` the trajectory is forced onto the mean `delta` that the FULL model
    carries there, supplied as `dbar_idx` from a first pass. Collapsing against a
    delta the full model never occupies would be a different (and unfair) model,
    which is why this needs two passes rather than a fixed constant.
    """
    nt, nd = len(thetas), len(DELTAS)
    q = np.array([[K[i][j]["q"] for j in range(nd)] for i in range(nt)])
    c = np.array([[K[i][j]["c"] for j in range(nd)] for i in range(nt)])
    csd = np.array([[K[i][j]["c_sd"] for j in range(nd)] for i in range(nt)])
    out = np.array([[K[i][j]["delta_out"] for j in range(nd)] for i in range(nt)])

    theta = np.full(n_traj, theta0)
    # start on the rail of the initial gamma, as the cascade does
    d_idx = np.full(n_traj, int(np.abs(DELTAS - delta_star(GAMMA0)).argmin()))
    alive = np.ones(n_traj, bool)
    surv = [1.0]
    occ = np.zeros((nt, nd))
    for _ in range(max_depth):
        if not alive.any():
            break
        idx = np.flatnonzero(alive)
        t = theta[idx]
        ti = np.clip(np.searchsorted(thetas, t) - 1, 0, nt - 1)
        if keep_delta:
            dj = d_idx[idx]
            if record_occupancy:
                np.add.at(occ, (ti, dj), 1.0)
        else:
            dj = dbar_idx[ti]          # mean-field collapse, from the full pass
        died = rng.random(t.size) < q[ti, dj]
        cc = c[ti, dj]
        step = (np.maximum(rng.normal(cc, csd[ti, dj]), 0.0) if stoch_theta else cc)
        theta[idx] = t + step / phi
        if keep_delta:
            p = out[ti, dj]
            tot = p.sum(axis=1)
            ok = tot > 0
            u = rng.random(t.size)
            nxt = np.zeros(t.size, int)
            cum = np.cumsum(np.where(ok[:, None], p / np.maximum(tot[:, None], 1e-12),
                                     1.0 / nd), axis=1)
            nxt = (u[:, None] > cum).sum(axis=1).clip(0, nd - 1)
            d_idx[idx] = nxt
        alive[idx[died]] = False
        alive[idx[theta[idx] > thetas[-1]]] = False
        surv.append(alive.sum() / n_traj)
    s = np.array(surv)
    depth = float(len(s) - 1)
    for i in range(len(s) - 1):
        if s[i] >= 0.5 >= s[i + 1]:
            depth = i + (s[i] - 0.5) / (s[i] - s[i + 1])
            break
    if not record_occupancy:
        return depth
    tot = occ.sum(axis=1, keepdims=True)
    mean_d = np.where(tot[:, 0] > 0,
                      (occ * DELTAS).sum(axis=1) / np.maximum(tot[:, 0], 1e-12),
                      delta_star(GAMMA0))
    return depth, np.abs(mean_d[:, None] - DELTAS[None, :]).argmin(axis=1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=40)
    ap.add_argument("--sigma-frac", type=float, default=0.03)
    ap.add_argument("--n-theta", type=int, default=12)
    ap.add_argument("--trials", type=int, default=600)
    ap.add_argument("--n-traj", type=int, default=100000)
    ap.add_argument("--seed", type=int, default=20260730)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/coarse_grain_axes.json"))
    args = ap.parse_args()

    om, sf = args.omega, args.sigma_frac
    sigma_counts = sf * delta_star(GAMMA0) * om
    theta0 = GAMMA0 / (1.0 + GAMMA0)
    thetas = np.linspace(theta0, 0.48, args.n_theta)
    phis = np.array([int(round(fc * om)) for fc in FUEL_CONCS], float)

    t0 = time.time()
    print(f"Omega={om} sigma/delta*={sf}  {args.n_theta} theta x {len(DELTAS)} delta "
          f"x {args.trials} trials")
    print(f"  reference: plain simulation 0.6325 +- 0.016 (400 trials, paired) | "
          f"memoryless 0.7217 | integral 0.7919 | theta-ensemble 0.7641\n")

    rng = np.random.default_rng(args.seed + 7)
    variants = {("delta kept", "theta stoch"): (True, True),
                ("delta kept", "theta determ"): (True, False),
                ("delta collapsed", "theta stoch"): (False, True),
                ("delta collapsed", "theta determ"): (False, False)}
    depths = {k: [] for k in variants}
    cells = []
    for fc, phi in zip(FUEL_CONCS, phis):
        ts = stage_time(GAMMA0, fc, 2.0)
        K = [[kernel_cell(float(t), float(d), om, int(phi), sigma_counts, ts,
                          args.trials, args.seed + 101 * i + 7 * j + int(fc))
              for j, d in enumerate(DELTAS)] for i, t in enumerate(thetas)]
        # pass 1: the full model, which also measures the mean delta it occupies
        # at each theta -- what the collapsed variants are then forced onto
        d_full, dbar_idx = propagate(K, thetas, int(phi), theta0, True, True,
                                     args.n_traj, rng, record_occupancy=True)
        for k, (kd, st) in variants.items():
            if kd and st:
                depths[k].append(d_full)
                continue
            depths[k].append(propagate(K, thetas, int(phi), theta0, kd, st,
                                       args.n_traj, rng, dbar_idx=dbar_idx))
        cells.append({"phi": int(phi), "kernel": K})
        print(f"  Phi={int(phi):>6} " + "  ".join(
            f"{k[0].split()[1][:4]}/{k[1].split()[1][:5]}={depths[k][-1]:6.2f}"
            for k in variants))

    print(f"\n{'variant':>34} {'exponent':>16} {'depths':>34}")
    res = {}
    for k, v in depths.items():
        p, se = fit(phis, np.array(v))
        res[" x ".join(k)] = {"exponent": p, "se": se, "depths": v}
        print(f"{' x '.join(k):>34} {p:>9.4f}+-{se:<5.3f} "
              + " ".join(f"{x:6.2f}" for x in v))
    print(f"{'plain simulation (measured)':>34} {0.6325:>9.4f}+-{0.016:<5.3f}")

    a = res["delta kept x theta determ"]["exponent"]
    b = res["delta collapsed x theta stoch"]["exponent"]
    print(f"\n  P1: delta-kept/theta-determ = {a:.4f}  vs  "
          f"delta-collapsed/theta-stoch = {b:.4f}")
    print(f"      -> keeping the COORDINATE beats keeping the NOISE"
          if a < b else
          f"      -> P1 REFUTED: keeping the noise beats keeping the coordinate")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"variants": res, "cells": cells},
                                   indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
