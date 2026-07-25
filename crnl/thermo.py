"""Entropy production for reversible chemical reaction networks.

For a jump n -> n' via reaction rho, the medium (heat) entropy production is

    dS = ln[ a_rho(n) / a_{-rho}(n') ]

with the reverse propensity evaluated at the POST-jump state. That convention is
not cosmetic: evaluating the reverse pre-jump is state-dependent, hits log(0),
and fabricates dissipation at equilibrium.

For reversible AM there is also an exact closed form, because rank(S) = 2 leaves
a single cycle of uniform affinity:

    dS(n -> n' via rho) = ln W(n') - ln W(n) + s_rho * ln(1/gamma)
    W(n) = N! / prod(n_i!)          s_rho = +1 forward, -1 reverse

so along any trajectory

    dS_total = ln[W(n_stop)/W(n_0)] + (A/3) * (M_forward - M_reverse)

Two independently interpretable terms (Hill/Schnakenberg): a boundary term that
depends only on the endpoints, and a cycle term counting net traversals. This is
what `decompose` reports, and it means a simulation needs only an integer counter
-- no logarithms in any hot loop.

`entropy_step` remains the general primitive (the identity above is specific to
this network's symmetry, and fails for asymmetric rate constants). It requires a
COMPLETE reverse pairing and raises rather than dividing by zero.
"""

from __future__ import annotations

import math

import numpy as np

from .reactions import ReactionNetwork


def ln_multinomial(n) -> float:
    """ln W(n) = ln( N! / prod(n_i!) ), the boundary term of the decomposition."""
    n = np.asarray(n, dtype=np.int64)
    total = int(n.sum())
    return float(math.lgamma(total + 1)
                 - sum(math.lgamma(int(c) + 1) for c in n))


def entropy_step(
    net: ReactionNetwork,
    pairing: np.ndarray,
    j: int,
    n_before,
    n_after,
    cs: np.ndarray,
) -> float:
    """Medium entropy production of one jump, ln[a_j(n) / a_reverse(n')].

    The reverse propensity is evaluated AFTER the jump. Raises if reaction j has
    no reverse, or if either propensity is non-positive (which would mean the jump
    was not actually firable, or the pairing is wrong).
    """
    rev = int(pairing[j])
    if rev < 0:
        raise ValueError(
            f"reaction {j} ({net.reactions[j].name!r}) has no reverse; "
            "entropy production is undefined for an incomplete pairing"
        )
    a_fwd = float(net.propensities(np.asarray(n_before), cs)[j])
    a_rev = float(net.propensities(np.asarray(n_after), cs)[rev])
    if a_fwd <= 0.0 or a_rev <= 0.0:
        raise ValueError(
            f"non-positive propensity in entropy_step (forward {a_fwd}, "
            f"reverse {a_rev}); the jump is not firable or the pairing is wrong"
        )
    return math.log(a_fwd / a_rev)


def decompose(n0, n_stop, net_reaction_firings: float, affinity: float,
              boundary: float | None = None) -> dict:
    """Split total entropy production into boundary and cycle contributions.

        dS_total = boundary + (affinity/3) * net_reaction_firings

    `net_reaction_firings` is the NET per-reaction count (forward - reverse)
    summed over all three pairs, exactly as an integer counter in a simulation
    would accumulate it and exactly what `cme.first_passage` returns. The /3
    lives HERE and only here: it converts per-reaction firings into cycle
    traversals, since one traversal is one firing of each of f1, f2, f3.
    Dividing by 3 on both sides makes every result 3x too small, and no test of
    either function alone can detect it -- see
    tests/test_cme.py::test_decomposition_composes_with_first_passage.

    `boundary` should be the EXPECTED ln W difference from the exact solve
    (cme.first_passage returns it). Passing a single representative `n_stop`
    instead is not a harmless approximation: real absorbing states carry a
    substantial blank population, and a hand-built n_B = 0 stopping state gets
    the sign of the boundary term wrong. `n_stop` is therefore only used when
    `boundary` is None, for the single-trajectory case where it is exact.
    """
    if boundary is None:
        if n_stop is None:
            raise ValueError("give either an explicit boundary or a stopping state")
        boundary = ln_multinomial(n_stop) - ln_multinomial(n0)
    cycle = (affinity / 3.0) * float(net_reaction_firings)
    return {"boundary": float(boundary), "cycle": cycle,
            "total": float(boundary) + cycle}
