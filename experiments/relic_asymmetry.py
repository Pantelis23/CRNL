"""Is the surviving asymmetry dynamical or accidental? -- the antimatter reading.

AM's disagreement reaction is literally an annihilation:

    X + Y  ->  2B        matter + antimatter -> two photons
    2B     ->  X + Y     pair production
    B + X  ->  2X        (no cosmological counterpart -- see below)

and the three ingredients Sakharov requires for a matter excess are all knobs
this rig already has:

  * violation of the conserved number: X - Y is UNCHANGED by annihilation and
    pair production; only the recruitment reactions move it.
  * C/CP violation: the tilt `beta` of `am_asymmetric`, which makes X and Y
    inequivalent while costing no thermodynamic force.
  * departure from equilibrium: gamma < 1 (driven), and the expansion deadline
    1/H of FINDINGS 5.1.

WHERE THE ANALOGY IS NOT ONE. Recruitment `B + X -> 2X` is autocatalysis and has
no counterpart in the standard picture, where a small asymmetry survives
LINEARLY -- annihilation removes matched pairs and whatever excess existed is
left over. Here the excess is AMPLIFIED by an instability. So this is not a model
of baryogenesis; it is the question of what changes when the asymmetry is fed
through a restoring landscape instead of a passive one. That difference is the
whole content, and it should not be dressed up as a cosmology result.

THE QUESTION. Start exactly symmetric (n_X = n_Y). Two things can decide which
species survives:

    the tilt        g/lambda,  g = d(x-y)/dt at the symmetric point = (2 beta/9)(1-gamma)
    shot noise      sigma = sqrt( (D_0/2) / (lambda Omega) ) = sqrt( (1+gamma) / (3(1-2 gamma) Omega) )

using FINDINGS 15's corrected D_0(gamma) = (1+gamma)/9. The relic is DYNAMICAL
when g/lambda > sigma and ACCIDENTAL below it, and the crossover sits at

    beta * sqrt(Omega)  =  (sqrt3/2) * (1-2 gamma)/(1-gamma)      [0.820 at gamma=0.05]

PREDICTION, written before running: P(X survives) collapses onto a single curve
in the variable u = (g/lambda)/sigma, and that curve is the normal CDF Phi(u).
Both parts are parameter-free. The collapse is the structural claim and should
survive even if Phi is only approximate -- the effective-seed picture treats the
accumulated noise as a single Gaussian kick, which is a real approximation.

Computed EXACTLY: the splitting probability from the symmetric start via
`cme.first_passage`, so there is no sampling error to hide a failed collapse.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
from scipy.stats import norm

from crnl.cme import first_passage
from crnl.networks.am_asymmetric import am_asymmetric, beta_critical
from crnl.networks.am_reversible import delta_star, lambda_antisym
from crnl.vectorized import compile_network, gillespie_fast

THETA = 0.80          # "decided" at |n_X - n_Y| >= THETA * delta* * Omega


def seed_ratio(gamma: float, beta: float, omega: float) -> float:
    """u = (g/lambda) / sigma -- the tilt measured in units of the shot noise."""
    lam = lambda_antisym(gamma)
    g = (2.0 * beta / 9.0) * (1.0 - gamma)
    sigma = np.sqrt((1.0 + gamma) / (3.0 * (1.0 - 2.0 * gamma) * omega))
    return float((g / lam) / sigma)


def beta_for_u(gamma: float, u: float, omega: float) -> float:
    """Invert `seed_ratio`, so a sweep can be laid out in u rather than beta."""
    return u / seed_ratio(gamma, 1.0, omega)


def crossover_beta_root_omega(gamma: float) -> float:
    """The value of beta*sqrt(Omega) at which tilt and noise are equal."""
    return float(np.sqrt(3.0) / 2.0 * (1 - 2 * gamma) / (1 - gamma))


def p_survive(gamma: float, beta: float, omega: int) -> float | None:
    """P(X survives | exactly symmetric start), exact."""
    ds = delta_star(gamma)
    nb = int(round(omega * gamma / (1.0 + gamma)))
    if (omega - nb) % 2:                 # n_X = n_Y needs an even remainder
        nb += 1
    half = (omega - nb) // 2
    n0 = np.array([half, half, nb], dtype=np.int64)
    assert n0.sum() == omega and n0[0] == n0[1]
    thr = max(2, int(round(THETA * ds * omega)))

    def absorbing(s, thr=thr):
        return abs(int(s[0]) - int(s[1])) >= thr

    fp = first_passage(am_asymmetric(gamma, beta), int(omega), float(omega),
                       n0, absorbing, None)
    return float(fp["split"]) if fp["valid"] else None


def freeze_out(gamma: float, beta: float, omega: int, hubble: float,
               trials: int, seed: int) -> dict:
    """Outcome shares under an expansion deadline.

    FINDINGS 5.1 proved the expanding SSA is EXACTLY ordinary SSA stopped at
    internal time tau = 1/H (bit-for-bit, 0/300 mismatches), so this runs
    ordinary SSA to t_max = 1/H rather than the expanding integrator. Using the
    expanding path here would be slower and would test the time change again
    rather than the physics.

    PREDICTION, written before running: the deadline decides WHETHER a relic
    forms, not WHICH one. So P(X | decided) stays Phi(u) at every H, while
    P(decided) falls with H and rises with |u|. If instead P(X | decided) moves
    with H, the tilt and the deadline are not separable and the split above is
    wrong.
    """
    ds = delta_star(gamma)
    nb = int(round(omega * gamma / (1.0 + gamma)))
    if (omega - nb) % 2:
        nb += 1
    half = (omega - nb) // 2
    n0 = np.array([half, half, nb], dtype=np.int64)
    thr = max(2, int(round(THETA * ds * omega)))
    comp = compile_network(am_asymmetric(gamma, beta), float(omega))
    rng = np.random.default_rng(seed)
    nx = ny = und = 0
    for _ in range(trials):
        r = gillespie_fast(comp, n0, rng, t_max=1.0 / hubble)
        d = int(r.n_final[0]) - int(r.n_final[1])
        if d >= thr:
            nx += 1
        elif d <= -thr:
            ny += 1
        else:
            und += 1
    dec = nx + ny
    return {"gamma": gamma, "beta": beta, "omega": omega, "hubble": hubble,
            "trials": trials, "p_decided": dec / trials,
            "p_x_given_decided": (nx / dec) if dec else float("nan"),
            "p_undecided": und / trials}


def part_freeze(gamma, omegas, us, hubbles, trials, seed) -> list[dict]:
    rows = []
    for om in omegas:
        print(f"\nOmega={om}   P(decided) / P(X | decided)   [Phi(u) in brackets]")
        print(f"{'u':>6} {'Phi(u)':>8} " + " ".join(f"{'H=' + str(h):>16}" for h in hubbles))
        for u in us:
            b = beta_for_u(gamma, u, om)
            cells = []
            for h in hubbles:
                r = freeze_out(gamma, b, int(om), float(h), trials, seed)
                rows.append(r)
                cells.append(f"{r['p_decided']:>7.3f}/{r['p_x_given_decided']:>8.4f}")
            print(f"{u:>6.2f} {norm.cdf(u):>8.5f} " + " ".join(cells))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.05)
    ap.add_argument("--omegas", type=int, nargs="+", default=[60, 120, 240, 480])
    ap.add_argument("--us", type=float, nargs="+",
                    default=[0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0])
    ap.add_argument("--part", choices=["seed", "freeze", "all"], default="seed")
    ap.add_argument("--hubbles", type=float, nargs="+",
                    default=[0.02, 0.05, 0.10, 0.20])
    ap.add_argument("--trials", type=int, default=800)
    ap.add_argument("--seed", type=int, default=20260727)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/relic_asymmetry.json"))
    args = ap.parse_args()
    if args.part == "freeze":
        rows = part_freeze(args.gamma, args.omegas, args.us, args.hubbles,
                           args.trials, args.seed)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(rows, indent=2))
        print(f"\nwrote {args.out}")
        return

    g = args.gamma
    print(f"gamma={g}  beta_c={beta_critical(g):.4f}  "
          f"crossover at beta*sqrt(Omega) = {crossover_beta_root_omega(g):.4f}")
    print("\nP(X survives) from an exactly symmetric start, exact:")
    print(f"{'u':>6} {'Phi(u)':>8} " +
          " ".join(f"{'Om=' + str(o):>18}" for o in args.omegas))
    print(f"{'':>15} " + " ".join(f"{'beta':>9}{'P(X)':>9}" for _ in args.omegas))

    rows = []
    t0 = time.time()
    for u in args.us:
        cells = []
        for om in args.omegas:
            b = beta_for_u(g, u, om)
            if b >= beta_critical(g):
                cells.append(f"{'past fold':>18}")
                continue
            p = p_survive(g, b, om)
            if p is None:
                cells.append(f"{'invalid':>18}")
                continue
            cells.append(f"{b:>9.5f}{p:>9.5f}")
            rows.append({"gamma": g, "u": u, "omega": om, "beta": b,
                         "p_survive": p, "phi": float(norm.cdf(u)),
                         "beta_root_omega": b * np.sqrt(om)})
        print(f"{u:>6.2f} {norm.cdf(u):>8.5f} " + " ".join(cells))

    # the collapse: spread of P(X) across Omega at fixed u
    print(f"\n{'u':>6} {'P(X) spread across Omega':>26} {'mean':>8} "
          f"{'Phi(u)':>8} {'mean-Phi':>9}")
    for u in args.us:
        ps = [r["p_survive"] for r in rows if r["u"] == u]
        if len(ps) < 2:
            continue
        print(f"{u:>6.2f} {max(ps) - min(ps):>26.5f} {np.mean(ps):>8.5f} "
              f"{norm.cdf(u):>8.5f} {np.mean(ps) - norm.cdf(u):>+9.5f}")

    print(f"\n({time.time() - t0:.0f}s)  A flat column of spreads is the collapse; "
          "a small\nmean-Phi column is the stronger, parameter-free claim.")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
