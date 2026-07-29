"""Base audit: re-derive every load-bearing closed form from the CURRENT code.

The project has grown five new modules and changed one behavioural function
(`information.wall_coefficient`, FINDINGS 15), so the foundations underneath the
early sections are no longer obviously the ones those sections were written
against. The test suite proves the code is self-consistent; this proves it still
agrees with what is PUBLISHED, which is a different question and the one that
rots silently.

Each row states the section, the number as printed in FINDINGS, and the number the
code produces today. Exit status is nonzero if any row disagrees, so this can be
run as a gate rather than read as a report.

Deliberately NOT included: anything whose published value came from sampling.
Those cannot be re-derived exactly and their agreement is the test suite's job.
This file checks only closed forms and exact linear algebra.
"""
from __future__ import annotations

import sys

import numpy as np

CHECKS = []


def check(section, what, published, got, tol=1e-9, note=""):
    ok = abs(float(got) - float(published)) <= tol * max(1.0, abs(published))
    CHECKS.append((section, what, published, float(got), ok, note))


def main() -> int:
    from crnl.cme import first_passage, splitting_probability
    from crnl.networks.am_asymmetric import am_asymmetric
    from crnl.networks.am_fueled import death_waste_fraction
    from crnl.networks.am_reversible import (
        GAMMA_C, am_reversible, cycle_affinity, delta_star, lambda_antisym,
        reverse_pairing,
    )
    from crnl.networks.n_winner_reversible import (
        affinity_critical, breaking_diffusion, diffusion_closed, lambda_breaking,
        lambda_closed,
    )
    from crnl.information import wall_coefficient, wall_coefficient_gain_only
    from experiments.relic_asymmetry import crossover_beta_root_omega

    # -- 9.1: the landscape's price ---------------------------------------
    check("9.1", "gamma_c", 0.5, GAMMA_C)
    net = am_reversible(0.2)
    check("9.1", "A(gamma=0.2) = -3 ln gamma", -3 * np.log(0.2),
          cycle_affinity(net, reverse_pairing(net)))
    check("9.1", "A_c(2) = 3 ln 2", 3 * np.log(2),
          cycle_affinity(am_reversible(0.5), reverse_pairing(am_reversible(0.5))))

    # -- 2 / 12 / 15: the wall coefficient --------------------------------
    check("2", "lambda_antisym(0) = 1/3", 1 / 3, lambda_antisym(0.0))
    check("2", "kappa(0) = 3/2", 1.5, wall_coefficient(0.0))
    check("15", "D_0(2, g) = (1+g)/9  [g=0.35]", (1 + 0.35) / 9,
          breaking_diffusion(2, 0.35))
    check("15", "kappa = lambda/(2 D_0)  [g=0.35]",
          lambda_breaking(2, 0.35) / (2 * breaking_diffusion(2, 0.35)),
          wall_coefficient(0.35))
    check("15", "kappa(0.35) = (3/2)(1-2g)/(1+g)", 1.5 * 0.3 / 1.35,
          wall_coefficient(0.35))
    check("15", "superseded kappa still callable = (3/2)(1-2g)", 1.5 * 0.3,
          wall_coefficient_gain_only(0.35))
    check("15", "kappa/kappa_12 = 1/(1+g)  [g=0.45]", 1 / 1.45,
          wall_coefficient(0.45) / wall_coefficient_gain_only(0.45))

    # -- 2: delta* closed form against the numeric fixed point ------------
    from crnl.networks.am_asymmetric import interior_fixed_points
    pts = interior_fixed_points(0.2, 0.0)
    check("2", "delta*(0.2) vs numeric attractor", delta_star(0.2),
          pts[-1]["x"] - pts[-1]["y"], tol=1e-6)

    # -- 13 / 14: the n-winner laws ---------------------------------------
    for n in (2, 4, 16):
        check("14", f"D_0(n={n}) = (2n-3)/(2n-1)^2",
              (2 * n - 3) / (2 * n - 1) ** 2, diffusion_closed(n))
        check("14", f"lambda(n={n}) = 1/(2n-1)", 1 / (2 * n - 1), lambda_closed(n))
    check("14", "lambda/D_0 -> 1 (n=64)", (2 * 64 - 1) / (2 * 64 - 3),
          lambda_closed(64) / diffusion_closed(64), tol=1e-9)
    check("13", "A_c(2) = 3 ln 2", 3 * np.log(2), affinity_critical(2), tol=1e-6)
    check("13", "A_c(n)/(9 ln n) -> 1 (n=48)", 1.0,
          affinity_critical(48) / (9 * np.log(48)), tol=0.12,
          note="asymptote approached from below; 12% band is the measured gap")

    # -- 16: the tilt costs no thermodynamic force ------------------------
    for beta in (0.0, 0.3, 0.7):
        na = am_asymmetric(0.2, beta)
        check("16", f"A(beta={beta}) = -3 ln gamma", -3 * np.log(0.2),
              cycle_affinity(na, reverse_pairing(na)))

    # -- 18 / 20: the two parameter-free thresholds -----------------------
    check("18", "beta sqrt(Omega) crossover at g=0.05",
          np.sqrt(3) / 2 * 0.9 / 0.95, crossover_beta_root_omega(0.05))
    check("20", "waste fraction at landscape death = 1/3", 1 / 3,
          death_waste_fraction(1.0))

    # -- 21: the general splitting probability matches first_passage ------
    n0, thr = np.array([26, 14, 10]), 20
    absb = (lambda s: abs(int(s[0]) - int(s[1])) >= thr)
    a = first_passage(am_reversible(0.3), 50, 50.0, n0, absb,
                      reverse_pairing(am_reversible(0.3)))
    b = splitting_probability(am_reversible(0.3), 50, 50.0, n0, absb,
                              lambda s: int(s[0]) > int(s[1]))
    check("21", "splitting_probability == first_passage on AM",
          a["split"], b["split"], tol=1e-12)

    w = max(len(c[1]) for c in CHECKS)
    print(f"{'sec':>5}  {'quantity':<{w}}  {'published':>14} {'now':>14}  ok")
    bad = 0
    for sec, what, pub, got, ok, note in CHECKS:
        bad += not ok
        print(f"{sec:>5}  {what:<{w}}  {pub:>14.9g} {got:>14.9g}  "
              f"{'ok' if ok else 'MISMATCH'}" + (f"   [{note}]" if note else ""))
    print(f"\n{len(CHECKS)} checks, {bad} mismatched")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
