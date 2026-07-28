"""The coarse-graining hierarchy: what a simulation is allowed to throw away.

CRNL's whole method is a two-point version of this -- ODE against exact SSA, and
the gap is the subject. This fills in the levels between, so the question "what
can a simulation discard and still get restoration right?" becomes measurable
rather than rhetorical.

    level            noise model                        cost per unit time
    ODE              none                               O(1)
    CLE              Gaussian, variance = propensity    O(1/dt)
    tau-leap         Poisson over a fixed window        O(1/tau)
    SSA              exact jumps, one at a time         O(Omega)
    CME              none -- solves the distribution    O(Omega^2) memory, exact

WHY THIS IS NOT A NUMERICS EXERCISE. Kurtz's theorem says the density process
converges to the mass-action ODE on finite time intervals, and FINDINGS 5.1 leans
on exactly that. It is true, and it does NOT license discarding the molecules for
this observable: restoration is a statement about tails and long times, where the
convergence is not uniform. The error probability vanishes in the limit while
being nonzero at every finite Omega. So the limit theorem cannot tell you what
your simulation may throw away -- only a measurement can.

WHERE THE CLE IS EXPECTED TO FAIL, AND WHY IT IS NOT OBVIOUS. The CLE is the
quadratic truncation of the jump process's large-deviation Hamiltonian,
`H(x,p) = sum_j a_j (exp(p.S_j) - 1)`, keeping `p.S_j + (p.S_j)^2/2`. It therefore
agrees with the exact rate function only to second order in the conjugate
momentum -- i.e. near the saddle, where FINDINGS 15 already checked it: the exact
quasipotential's ridge curvature matches `lambda/(2 D_0)`, the CLE's answer, to
0.1%. The prediction is that this agreement is SPECIFIC to that limit. For a 1-D
birth-death chain the exact barrier is `integral ln(a+/a-)` against the diffusion
approximation's `integral 2(a+-a-)/(a++a-)`, and `ln r > 2(r-1)/(r+1)` for r > 1,
so the exact barrier is LARGER: **the CLE should underestimate the barrier and
overestimate the failure probability, increasingly so as the barrier grows.**

BOTH INTEGRATORS PRESERVE CONSERVATION EXACTLY. A step is always a whole integer
(tau-leap) or real (CLE) combination of stoichiometry columns, each of which
conserves, so the conserved totals cannot drift. Negative counts are handled by
halving the step and retrying rather than by clipping a species at zero -- clipping
would break conservation and, worse, would quietly hand the minority species a
floor it did not earn, which is the failure mode behind three withdrawn results
here. Both integrators report how often they had to retry; a run that retries
constantly is reporting that its step is too big, not a result.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .vectorized import Compiled, propensities_fast

__all__ = ["ApproxResult", "cle_run", "tau_leap_run", "MAX_RETRY"]

#: Step halvings allowed before a step is abandoned. Six is a 64x reduction.
MAX_RETRY = 6


@dataclass
class ApproxResult:
    t_final: float
    n_final: np.ndarray
    steps: int
    retries: int
    hit_budget: bool


def _run(compiled: Compiled, n0, rng, *, step: float, stop, t_max: float,
         max_steps: int, poisson: bool) -> ApproxResult:
    """Shared driver for the two approximate integrators.

    `poisson=True` draws integer firing counts (tau-leaping); False draws the
    Gaussian increment of the chemical Langevin equation and keeps real-valued
    counts. Everything else -- the conservation-safe retry, the stop predicate,
    the budget -- is identical, so a difference between the two levels is the
    noise model and not the harness.
    """
    S = compiled.S.astype(float)
    n = np.array(n0, dtype=float)
    t = 0.0
    steps = retries = 0
    while steps < max_steps and t < t_max:
        h = step
        for attempt in range(MAX_RETRY + 1):
            a = propensities_fast(compiled, n)
            a = np.maximum(a, 0.0)
            if poisson:
                fire = rng.poisson(a * h)
            else:
                mean = a * h
                fire = mean + np.sqrt(mean) * rng.standard_normal(len(a))
            cand = n + S @ fire
            if (cand >= 0.0).all():
                n = cand
                t += h
                break
            h *= 0.5
            retries += 1
        else:
            # could not take a legal step even at 1/64 of the target
            return ApproxResult(t, n, steps, retries, True)
        steps += 1
        if stop is not None and stop(n):
            return ApproxResult(t, n, steps, retries, False)
    return ApproxResult(t, n, steps, retries, steps >= max_steps)


def cle_run(compiled: Compiled, n0, rng, *, dt: float, stop=None,
            t_max: float = np.inf, max_steps: int = 2_000_000) -> ApproxResult:
    """Chemical Langevin equation, Euler-Maruyama, conservation exact.

    dn = sum_j S_j [ a_j dt + sqrt(a_j dt) * N(0,1) ]. Counts stay real-valued --
    rounding them to integers would smuggle in a discreteness the CLE does not
    have, and the point of this level is to measure what its absence costs.
    """
    return _run(compiled, n0, rng, step=dt, stop=stop, t_max=t_max,
                max_steps=max_steps, poisson=False)


def tau_leap_run(compiled: Compiled, n0, rng, *, tau: float, stop=None,
                 t_max: float = np.inf, max_steps: int = 2_000_000) -> ApproxResult:
    """Explicit tau-leaping: Poisson(a_j * tau) firings per window.

    Keeps integer counts and the exact jump structure, and approximates only the
    assumption that propensities are constant across the window -- so it sits
    strictly between the CLE and the SSA in what it discards.
    """
    return _run(compiled, n0, rng, step=tau, stop=stop, t_max=t_max,
                max_steps=max_steps, poisson=True)
