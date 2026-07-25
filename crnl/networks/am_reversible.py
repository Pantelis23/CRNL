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

import math

import numpy as np
from scipy.linalg import null_space

from ..reactions import Reaction, ReactionNetwork

#: The bistability threshold. Above this the landscape has a single minimum and
#: no population size can restore, because there is nothing to restore toward.
#: At the symmetric fixed point (1/3, 1/3, 1/3), the transverse (decision) mode
#: of the reduced system has eigenvalue (1 - 2*gamma)/3, which vanishes at
#: gamma = 1/2.
GAMMA_C = 0.5


def am_reversible(gamma: float, k: float = 1.0) -> ReactionNetwork:
    """Reversible AM with every reverse rate scaled by `gamma`."""
    if not math.isfinite(gamma) or gamma < 0:
        raise ValueError(f"gamma must be finite and >= 0, got {gamma}")
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

    If a reaction has more than one candidate reverse, the pairing is
    ambiguous and this raises `ValueError` rather than resolving it arbitrarily
    (e.g. by taking the first match) -- a silent, arbitrary resolution would let
    a downstream entropy-production sum mix fluxes from different pairs and
    report a wrong number with no error raised.
    """
    pairing = np.full(net.n_reactions, -1, dtype=np.int64)
    for i, ri in enumerate(net.reactions):
        matches = [j for j, rj in enumerate(net.reactions)
                   if j != i
                   and ri.reactants == rj.products
                   and ri.products == rj.reactants]
        if len(matches) > 1:
            raise ValueError(
                f"reaction {i} ({ri.name!r}) has {len(matches)} candidate "
                f"reverses {matches}; the pairing is ambiguous"
            )
        if matches:
            pairing[i] = matches[0]
    return pairing


def cycle_affinity(net: ReactionNetwork, pairing: np.ndarray) -> float:
    """Thermodynamic force magnitude around the network's cycle, in units of k_B T.

    A reversible pair (a, b) contributes a stoichiometric vector S[:, a] (taken
    in the `a` orientation, `a` being whichever of the pair has the smaller
    reaction index -- an arbitrary but fixed convention, not a claim about which
    member is chemically "forward"). The real precondition for a single well
    defined affinity is not that some list of reactions happens to sum to zero
    (that can hold by accident for disjoint cycles, or fail by accident because
    of how reactions were ordered) -- it is that the *cycle space* spanned by
    the reversible pairs is exactly one-dimensional. We compute that directly as
    the null space of the matrix of per-pair stoichiometries: dimension 0 means
    the pairs don't close into a cycle at all, dimension >1 means there is more
    than one independent cycle (e.g. several disjoint reversible pairs), and
    either way a single number cannot describe the drive.

    Given a 1-D null space, its generator c (normalised so the smallest nonzero
    magnitude entry is 1, i.e. integer traversal counts) gives the affinity as
    A = sum_i c_i * (ln k_{a_i} - ln k_{b_i}). The sign of c is only defined up
    to an overall flip (the null space doesn't know which way around the cycle
    is "positive"), so the result is oriented so A >= 0: this function returns
    the drive MAGNITUDE, not a signed circulation. For reversible AM with
    uniform rates this is exactly -3*ln(gamma) (positive for gamma < 1) and 0.0
    at gamma = 1.

    Computed from the DETERMINISTIC rate constants. Using stochastic_constants
    instead would be wrong by -3*ln(2) here, because the homodimer reverses carry
    a factor 2 in c that has nothing to do with the thermodynamics.
    """
    pairs = [(j, int(pairing[j])) for j in range(net.n_reactions)
             if pairing[j] > j]
    if not pairs:
        raise ValueError("network has no reversible pairs; affinity is undefined")

    # zero rate anywhere in the cycle -> unbounded drive, deliberately, not via
    # divide-by-zero in np.log below.
    if any(net.reactions[a].k == 0.0 or net.reactions[b].k == 0.0 for a, b in pairs):
        return float("inf")

    S = net.stoichiometry_matrix()
    Sp = S[:, [a for a, _ in pairs]]
    ns = null_space(Sp)
    if ns.shape[1] != 1:
        raise ValueError(
            f"the reversible pairs span a {ns.shape[1]}-dimensional cycle "
            "space, not 1; affinity is not a single number"
        )
    c = ns[:, 0]
    mags = np.abs(c)
    nonzero = mags > 1e-9 * mags.max()
    c = c / mags[nonzero].min()

    A = float(sum(
        ci * (math.log(net.reactions[a].k) - math.log(net.reactions[b].k))
        for ci, (a, b) in zip(c, pairs)
    ))
    return abs(A)
