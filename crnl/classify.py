"""Classifier: attractor-count-agnostic (design.md §3.4).

"End-state" is two different objects in the two engines, and they coincide only
in fortunate networks like AM. Conflating them is the fastest way to make the
word "general" a lie, so this module carries both criteria explicitly.

  * Stochastic absorption -- a configuration where every propensity is zero
    (a0 == 0). The chain physically halts. Cheap and exact. AM has exactly three
    such configurations: (Omega,0,0), (0,Omega,0), (0,0,Omega). Three is a
    theorem about AM, not a constant hardwired here.

  * Deterministic stable end-state -- a stable root of S.v(x) = 0 (Jacobian
    eigenvalues all with negative real part). A settling criterion, not a
    halting one: the ODE never stops, it asymptotes.

  * Dwelling test -- for a general CRN with a stable *interior* fixed point the
    stochastic system never absorbs (a0 never reaches 0); it fluctuates around
    the point forever. There the stochastic classifier must instead check that
    the trajectory has dwelt within a small ball of a known fixed point. AM
    never needs this path (all three attractors are simultaneously zero-
    propensity corners and stable roots), but the engine carries it or it is
    AM-only wearing a general costume.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .deterministic import integrate, jacobian
from .reactions import ReactionNetwork


# --------------------------------------------------------------------------- #
# Stochastic side                                                             #
# --------------------------------------------------------------------------- #

def is_absorbing(net: ReactionNetwork, n, omega: float) -> bool:
    """True if every propensity is zero at counts n (a0 == 0)."""
    cs = net.stochastic_constants(omega)
    return float(net.propensities(np.asarray(n), cs).sum()) <= 0.0


def classify_winner(result, blank: str = "B") -> str:
    """General n-outcome absorption classifier (spec §3.3).

    Returns a committed species name (single winner), "blank" (all committed
    zero -- the (0,...,0,Omega) corner), "coexist" (absorbed with >=2 surviving
    committed species -- impossible in AM, possible in n-winner), or "undecided"
    (budget exhausted, not absorbed -> the dwelling regime).
    """
    idx = {s: i for i, s in enumerate(result.species)}
    committed = [s for s in result.species if s != blank]
    if not result.absorbed:
        return "undecided"
    survivors = [s for s in committed if result.n_final[idx[s]] > 0]
    if len(survivors) == 0:
        return "blank"
    if len(survivors) == 1:
        return survivors[0]
    return "coexist"


def classify_am_outcome(result, species_order=("X", "Y", "B")) -> str:
    """Bin an SSA absorption of AM into {'X', 'Y', 'B', 'undecided'}.

    The third bin matters (design.md §3.4): all-blank (0,0,Omega) is
    stochastically absorbing and genuinely reachable at finite Omega, when the
    last committed pair annihilates via r1 before autocatalysis amplifies a
    lead. The deterministic view calls (0,0,1) a repeller and would never rest
    there; the stochastic view lands there with nonzero probability because it
    can hit the corner exactly at integer count. It gets its own bin.
    """
    w = classify_winner(result, blank="B")
    if w == "blank":
        return "B"
    if w == "coexist":
        return "undecided"  # AM has no coexisting absorbing state; treat as non-decision
    return w


# --------------------------------------------------------------------------- #
# Deterministic side                                                          #
# --------------------------------------------------------------------------- #

@dataclass
class FixedPoint:
    x: np.ndarray
    eigenvalues: np.ndarray
    kind: str  # 'stable', 'unstable', 'saddle', 'marginal'

    @property
    def is_stable(self) -> bool:
        return self.kind == "stable"


def _classify_eigs(eigs: np.ndarray, tol: float = 1e-6) -> str:
    re = eigs.real
    if np.all(re < -tol):
        return "stable"
    if np.all(re > tol):
        return "unstable"
    if np.any(re > tol) and np.any(re < -tol):
        return "saddle"
    return "marginal"


def stoichiometric_basis(net: ReactionNetwork, tol: float = 1e-9) -> np.ndarray:
    """Orthonormal basis (columns) for the stoichiometric subspace Im(S).

    The dynamics never leaves x0 + Im(S); every conservation law lives in the
    orthogonal complement (the left null space of S). Judging stability on the
    full state space would attach a spurious zero eigenvalue per conservation
    law -- e.g. AM's 2->2 reactions make (1,1,1) a left null vector of S, so the
    raw 3x3 Jacobian at every fixed point carries an extra 0 that would demote
    the rails from 'stable' to 'marginal'. Classifying *within* Im(S) recovers
    the reduced 2D picture of design.md §2.3.
    """
    S = net.stoichiometry_matrix()
    U, sv, _ = np.linalg.svd(S, full_matrices=True)
    rank = int(np.sum(sv > tol * max(1.0, sv[0] if sv.size else 1.0)))
    return U[:, :rank]


def classify_fixed_point(
    net: ReactionNetwork, x, tol: float = 1e-6, restrict: bool = True
) -> FixedPoint:
    """Classify a candidate fixed point by the Jacobian eigenvalues.

    With ``restrict=True`` (default) the Jacobian is projected onto the
    stoichiometric subspace before taking eigenvalues, so conservation laws do
    not contaminate the classification. This reproduces the design.md §2.3
    table: rails -> two -1 eigenvalues (stable), interior -> (+1/3, -1) (saddle),
    all-blank -> (+1, +1) (unstable), with no spurious zeros.
    """
    x = np.asarray(x, dtype=float)
    J = jacobian(net, x)
    if restrict:
        Q = stoichiometric_basis(net)
        J = Q.T @ J @ Q
    eigs = np.linalg.eigvals(J)
    return FixedPoint(x, eigs, _classify_eigs(eigs, tol))


def find_stable_endpoint(
    net: ReactionNetwork,
    x0,
    t_span: tuple[float, float] = (0.0, 400.0),
) -> np.ndarray:
    """Where the deterministic flow settles from x0 (its asymptotic attractor)."""
    traj = integrate(net, x0, t_span=t_span)
    return traj.final()


# --------------------------------------------------------------------------- #
# Dwelling test -- the interior-fixed-point path AM never exercises           #
# --------------------------------------------------------------------------- #

def dwells_near(
    trajectory_counts: np.ndarray,
    times: np.ndarray,
    omega: float,
    target: np.ndarray,
    radius: float = 0.05,
    relax_time: float = 5.0,
    dwell_multiple: float = 3.0,
) -> bool:
    """True if a recorded SSA trajectory has stayed within `radius` (in
    concentration units) of `target` for longer than `dwell_multiple` local
    relaxation times.

    This is the general-CRN stochastic settling criterion for interior fixed
    points that never absorb. `trajectory_counts` is (n_species, n_t) integer
    counts; divide by Omega to compare against `target` in concentration units.
    """
    conc = trajectory_counts / omega
    target = np.asarray(target, dtype=float)[:, None]
    inside = np.linalg.norm(conc - target, axis=0) <= radius
    if not inside.any():
        return False
    # find the longest contiguous run that is inside the ball
    longest = 0.0
    run_start = None
    for i, flag in enumerate(inside):
        if flag and run_start is None:
            run_start = times[i]
        elif not flag and run_start is not None:
            longest = max(longest, times[i - 1] - run_start)
            run_start = None
    if run_start is not None:
        longest = max(longest, times[-1] - run_start)
    return longest >= dwell_multiple * relax_time
