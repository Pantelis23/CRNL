"""Reversible Approximate Majority, and its closed-form landscape.

AM as written (networks/am.py) is irreversible, so its dissipation is formally
infinite -- there is no free-energy number to report. Adding the reverse of every
reaction, all scaled by a single parameter gamma, makes the cost finite and gives
the drive a thermodynamic meaning.

    forward (rate k)        reverse (rate gamma*k)
    f1: X + Y -> 2B         r1: 2B -> X + Y
    f2: B + X -> 2X         r2: 2X -> B + X
    f3: B + Y -> 2Y         r3: 2Y -> B + Y

Reaction order is fixed as [f1, f2, f3, r1, r2, r3] and documented, but the
reverse pairing is still *derived* rather than hardcoded -- a pairing function
that ignored its argument would silently return garbage for any other network.

All three reverses are homodimers (a single reactant species with coefficient 2).
This is the first network in the project to exercise the homodimer branch of the
units convention (c = 2*gamma*k/Omega, a = c*n(n-1)/2), which reactions.py owns
at engine level precisely because AM never triggered it.

Detailed balance requires the Wegscheider condition k_f1 k_f2 k_f3 = k_r1 k_r2
k_r3, i.e. 1 = gamma**3 for uniform rates -- satisfied only at gamma = 1. So
every gamma < 1 is genuinely driven, and gamma is the single knob.
"""

from __future__ import annotations

import numpy as np

from ..reactions import Reaction, ReactionNetwork

#: The bistability threshold. Above this the landscape has a single minimum and
#: no population size can restore, because there is nothing to restore toward.
#: Derived in docs/superpowers/specs/2026-07-25-dissipation-design.md §2.4.
GAMMA_C = 0.5


def am_reversible(gamma: float, k: float = 1.0) -> ReactionNetwork:
    """Reversible AM with every reverse rate scaled by `gamma`."""
    if gamma < 0:
        raise ValueError(f"gamma must be >= 0, got {gamma}")
    kr = gamma * k
    return ReactionNetwork(
        species=["X", "Y", "B"],
        reactions=[
            Reaction({"X": 1, "Y": 1}, {"B": 2}, k, name="f1:X+Y->2B"),
            Reaction({"B": 1, "X": 1}, {"X": 2}, k, name="f2:B+X->2X"),
            Reaction({"B": 1, "Y": 1}, {"Y": 2}, k, name="f3:B+Y->2Y"),
            Reaction({"B": 2}, {"X": 1, "Y": 1}, kr, name="r1:2B->X+Y"),
            Reaction({"X": 2}, {"B": 1, "X": 1}, kr, name="r2:2X->B+X"),
            Reaction({"Y": 2}, {"B": 1, "Y": 1}, kr, name="r3:2Y->B+Y"),
        ],
        name=f"am-reversible-gamma{gamma}",
    )


def reverse_pairing(net: ReactionNetwork) -> np.ndarray:
    """Index of each reaction's reverse, or -1 where there is none.

    Derived by matching reactants against products, NOT hardcoded: a function
    that returned a fixed table while ignoring `net` would silently mispair any
    other network (e.g. a reversible n-winner). Note that checking S columns
    negate each other is necessary but not sufficient to identify a reverse pair.
    """
    pairing = np.full(net.n_reactions, -1, dtype=np.int64)
    for i, ri in enumerate(net.reactions):
        for j, rj in enumerate(net.reactions):
            if i == j:
                continue
            if ri.reactants == rj.products and ri.products == rj.reactants:
                pairing[i] = j
                break
    return pairing
