"""Stochastic path: exact Gillespie SSA (design.md §3.1).

This is the lesson, so it is written by hand rather than borrowed. The loop:

    compute propensities a_j
    total a0 = sum a_j
    draw tau  = -ln(u1) / a0
    select reaction j with probability a_j / a0
    apply stoichiometry column j
    advance t += tau
    repeat until a0 == 0 (absorption) or a step budget is exhausted

About thirty lines and worth writing by hand. The units convention that turns
deterministic k into stochastic c, and the homodimer combinatorics, both live
in reactions.py -- this loop just consumes them, so it stays network-agnostic.

RNG discipline (design.md §5): every trial seeds deterministically from
(Omega, trial_index) so any single anomalous trajectory can be replayed exactly.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .reactions import ReactionNetwork


@dataclass
class SSAResult:
    t_final: float
    n_final: np.ndarray  # integer counts, shape (n_species,)
    steps: int
    absorbed: bool  # True if a0 hit 0; False if the step budget ran out
    species: list[str]


def gillespie(
    net: ReactionNetwork,
    n0,
    omega: float,
    rng: np.random.Generator,
    max_steps: int = 10_000_000,
    t_max: float = np.inf,
    record: bool = False,
):
    """Run one exact SSA trajectory to absorption or a budget limit.

    Parameters
    ----------
    n0 : integer molecule counts, shape (n_species,)
    omega : population scale; sets the stochastic constants c = f(k, Omega)
    rng : a numpy Generator, seeded by the caller for replayability
    record : if True, also return arrays of (times, counts) at every step for
        plotting a single trajectory against the ODE overlay (design.md §5:
        divide counts by Omega before comparing).
    """
    S = net.stoichiometry_matrix()
    cs = net.stochastic_constants(omega)
    n = np.array(n0, dtype=np.int64)

    t = 0.0
    steps = 0
    absorbed = False

    if record:
        ts = [0.0]
        ns = [n.copy()]

    while steps < max_steps and t < t_max:
        a = net.propensities(n, cs)
        a0 = a.sum()
        if a0 <= 0.0:
            absorbed = True
            break

        # time to next reaction
        tau = -np.log(rng.random()) / a0
        # select reaction: cumulative search on a fresh uniform draw
        j = int(np.searchsorted(np.cumsum(a), rng.random() * a0))
        # numerical guard: searchsorted can return len(a) if the last draw
        # lands exactly on a0 due to float rounding.
        if j >= net.n_reactions:
            j = net.n_reactions - 1

        n = n + S[:, j].astype(np.int64)
        t += tau
        steps += 1

        if record:
            ts.append(t)
            ns.append(n.copy())

    result = SSAResult(
        t_final=t,
        n_final=n,
        steps=steps,
        absorbed=absorbed,
        species=list(net.species),
    )
    if record:
        return result, np.array(ts), np.array(ns).T  # counts shape (n_species, n_t)
    return result


def seed_for(omega: float, trial: int, base: int = 0) -> np.random.Generator:
    """Deterministic per-trial RNG (design.md §5, RNG discipline).

    Seeded from (Omega, trial, base) so an anomalous trajectory at a given Omega
    and trial index can be replayed bit-for-bit.
    """
    ss = np.random.SeedSequence([int(base), int(round(omega)), int(trial)])
    return np.random.default_rng(ss)
