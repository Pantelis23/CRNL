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


def _check_gamma(gamma: float) -> None:
    """Reject non-finite or negative gamma, the one guard every function here needs."""
    if not math.isfinite(gamma) or gamma < 0:
        raise ValueError(f"gamma must be finite and >= 0, got {gamma}")


def am_reversible(gamma: float, k: float = 1.0) -> ReactionNetwork:
    """Reversible AM with every reverse rate scaled by `gamma`."""
    _check_gamma(gamma)
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


def lambda_antisym(gamma: float) -> float:
    """Restoring gain: the unstable eigenvalue of the decision mode (1,-1).

    At the symmetric point (1/3,1/3,1/3) the Jacobian's (x,y) block has
    off-diagonal -(2+2g)/3 and diagonal -(1+4g)/3, so the (1,-1) mode has
    eigenvalue (1-2g)/3. It is +1/3 at gamma=0 (design.md 2.3's AM saddle) and
    vanishes at gamma = 1/2 -- which is what makes GAMMA_C exact.
    """
    _check_gamma(gamma)
    return (1.0 - 2.0 * gamma) / 3.0


def lambda_sym(gamma: float) -> float:
    """Eigenvalue of the symmetric (1,1) mode at (1/3,1/3,1/3); always stable."""
    _check_gamma(gamma)
    return -(1.0 + 2.0 * gamma)


def delta_star(gamma: float) -> float:
    """Separation x* - y* of the two attractors; 0 at and above GAMMA_C.

    Off-symmetry requires b* = gamma/(1+gamma), which gives
        x* = [1 +- sqrt(1 - 4 g^3/(1-g))] / (2(1+g))
    so the separation is sqrt(disc)/(1+g). The discriminant numerator factors as
    (1-2*gamma)*(2*gamma**2+gamma+1): the quadratic factor has discriminant -7,
    so it is strictly positive for every real gamma, and (1-gamma) > 0 on the
    valid domain gamma < GAMMA_C (gamma >= GAMMA_C is handled by the early
    return below, and _check_gamma has already excluded gamma < 0). So `disc`
    is strictly positive whenever this line runs -- there is no live "disc <= 0"
    case to guard against. The unique real root of the numerator is gamma = 1/2:
    a pitchfork, with delta* ~ (4 sqrt(2)/3) sqrt(g_c - g).
    """
    _check_gamma(gamma)
    if gamma >= GAMMA_C:
        return 0.0
    disc = 1.0 - 4.0 * gamma ** 3 / (1.0 - gamma)
    return float(np.sqrt(disc) / (1.0 + gamma))


def fixed_points(gamma: float) -> list[dict]:
    """All fixed points in the simplex, with a stability field.

    `kind` is a GEOMETRIC label ({"symmetric", "attractor"}, plus "blank" for
    the gamma=0 boundary case below) -- it does not by itself say which points
    are stable. `stable` carries that: at gamma=0.2 the symmetric point has
    reduced eigenvalues (+0.2, -1.4) -- a SADDLE -- while at gamma=0.7 it is
    (-0.133, -2.4), the sole STABLE state. Both get kind="symmetric"; only
    `stable` tells them apart, and that distinction is the entire content of
    the bifurcation. Attractors are always stable=True (they only exist below
    GAMMA_C, where they are the two minima).

    There are THREE fixed points for 0 < gamma < GAMMA_C -- not four. The
    all-blank repeller (0,0,1) of irreversible AM leaves the simplex the
    instant gamma > 0 (there db/dt = -2*gamma), and the second symmetric root
    is negative. So the reversible model has no all-blank outcome for gamma >
    0. At gamma = 0 exactly it IS a fixed point (rhs = 0 there too), so it is
    returned as a fourth point, kind="blank", stable=False (it is a repeller).

    Above GAMMA_C only the symmetric point remains, and it is stable: a single
    minimum, no threshold, nothing to restore toward.
    """
    _check_gamma(gamma)
    out = [{"x": 1 / 3, "y": 1 / 3, "b": 1 / 3, "kind": "symmetric",
            "stable": lambda_antisym(gamma) < 0.0}]
    d = delta_star(gamma)
    if d > 0.0:
        b = gamma / (1.0 + gamma)
        total = 1.0 - b                       # x + y
        hi, lo = (total + d) / 2.0, (total - d) / 2.0
        out.append({"x": hi, "y": lo, "b": b, "kind": "attractor", "stable": True})
        out.append({"x": lo, "y": hi, "b": b, "kind": "attractor", "stable": True})
    if gamma == 0.0:
        out.append({"x": 0.0, "y": 0.0, "b": 1.0, "kind": "blank", "stable": False})
    return out


def theta_decide(gamma: float, frac: float = 0.7) -> float:
    """Decision threshold on delta = (n_X - n_Y)/Omega, scaled to the landscape.

    MUST scale with delta_star(gamma). A fixed threshold (e.g. 0.5) is
    unreachable above gamma ~ 0.417, where delta*(0.49) = 0.187: "deciding" then
    means fluctuating *past* the attractor, which inflates the measured
    dissipation by an order of magnitude and is a protocol artifact, not physics.
    """
    d = delta_star(gamma)
    if d <= 0.0:
        raise ValueError(
            f"no bistable landscape at gamma={gamma} (>= GAMMA_C={GAMMA_C}); "
            "there is no decision to threshold"
        )
    return frac * d


def initial_counts(omega: int, gamma: float, bias_frac: float = 0.2,
                   count_diff: int | None = None) -> np.ndarray:
    """Integer start [n_X, n_Y, n_B] with B empty.

    Follows the project convention (design.md 4, restoration_wall.py): B(0) = 0
    and the committed molecules carry the bias.

    Two ways to set the bias:

    * `count_diff` (PREFERRED when sweeping gamma) -- an explicit integer count
      difference n_X - n_Y, realised exactly.
    * `bias_frac` -- a fraction of delta_star(gamma), which keeps the protocol
      difficulty comparable across gamma in the same spirit as radix_wall.py's
      fixed pairwise margin. But note it is NOT representable on the integer
      lattice: at Omega=60, gamma=0.45 the target 0.0801 realises as 0.0667, 17%
      low, and the error jitters non-monotonically with gamma. Since one molecule
      of bias is worth ~20 k_B of dissipation, a fraction-driven gamma sweep can
      manufacture a fold-back in the headline curve that looks like physics.
      Use it for single-gamma work; use `count_diff` for sweeps.
    """
    if count_diff is None:
        delta0 = bias_frac * delta_star(gamma)
        count_diff = int(round(delta0 * omega))
    if not 0 <= count_diff <= omega:
        raise ValueError(
            f"count_diff={count_diff} is not reachable at omega={omega}")
    if (omega + count_diff) % 2 != 0:
        count_diff += 1          # keep n_X, n_Y integral
    n_x = (omega + count_diff) // 2
    n_y = omega - n_x
    return np.array([n_x, n_y, 0], dtype=np.int64)
