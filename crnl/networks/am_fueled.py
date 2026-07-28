"""AM whose drive is a finite fuel, so restoration can run out.

THE GAP THIS CLOSES (THEORIES T10b). Everywhere else in this project `gamma` is a
free parameter held fixed forever: an infinite reservoir, set once and maintained
at no cost. FINDINGS 9 therefore measures what restoration *dissipates* while
nothing ever runs down, and 12.1's depth ceiling is purely noise-limited with no
thermodynamic competitor. A drive that cannot be exhausted is the last piece of
free lunch in the model.

HOW FUEL ENTERS, AND WHY IT IS NOT A HARNESS TRICK. A drive gamma < 1 physically
IS a coupling to fuel hydrolysis, so the honest construction makes the fuel a
reactant:

    f1: X + Y + F -> 2B + W        r1: 2B + W -> X + Y + F     (rate x gamma_inf)
    f2: B + X + F -> 2X + W        r2: 2X + W -> B + X + F     (rate x gamma_inf)
    f3: B + Y + F -> 2Y + W        r3: 2Y + W -> B + Y + F     (rate x gamma_inf)

Every forward firing burns one F and makes one W; every reverse firing undoes it.
`gamma_inf` is the intrinsic ratio of the *uncoupled* reaction (1 = isoenergetic),
and the drive the chemistry actually feels is

    gamma_eff = gamma_inf * w / f

which starts near 0 with a full tank and rises as the tank empties. That is the
mirror of FINDINGS 19's cooling, which drove gamma DOWN.

The alternative -- keeping the 3-species network and letting the integrator adjust
gamma from a running count of firings -- would be the harness doing chemistry, the
failure mode that has already cost this project three withdrawn results. Fuel here
is a species with stoichiometry, and `n_F` is a genuinely independent coordinate:
a complete cycle f1 -> f2 -> f3 returns (X, Y, B) to where it started while
consuming three fuel, so the fixed-gamma model is a PROJECTION that discards a
coordinate which must exist.

TWO CONSERVATION LAWS: X + Y + B is conserved as before, and F + W = Phi.

WHAT IT COSTS STRUCTURALLY. These reactions are 3 -> 3, so the network leaves the
uniform-order-2 class. It is still *uniform* order 3, so `expanding.common_order`
accepts it and FINDINGS 19's reduction survives with lambda = 2H rather than H --
but the ordinary AM scaling does not carry over, and trimolecular steps are a real
idealisation (a physical model would resolve the fuel binding as a separate fast
step). Stated rather than hidden.

VERIFIED EXACTLY: at fixed f and w the (X, Y, B) drift of this network equals
`am_reversible(gamma_inf * w / f)` with time rescaled by f, to 1e-16 over 200
random interior states (`test_reduces_to_am_reversible_at_fixed_fuel`). That is
the anchor, in the style of FINDINGS 19's 0/300 check against `gillespie_expanding`.

PREDICTIONS, written before any run:

  P1  Restoration dies when the waste fraction reaches `1/(1 + 2 gamma_inf)` --
      exactly 1/3 at gamma_inf = 1 -- because that is where gamma_eff crosses
      gamma_c = 1/2 and the landscape stops existing. Parameter-free: no Omega,
      no Phi, no protocol.
  P2  Fuel consumed = net forward reaction firings, which is exactly the quantity
      FINDINGS 9's decomposition already measures. So 9's "cost of remembering"
      stops being a rate and becomes a LIFETIME.
  P3  At fixed fuel CONCENTRATION the fuel ceiling is Omega-INDEPENDENT (budget
      and burn rate are both extensive) while the noise ceiling grows like
      exp(Omega). So beyond a crossover Omega restoration is fuel-limited and
      **more molecules buy nothing** -- the exact mirror of FINDINGS 1's wall.
  P4  The bit is lost BEFORE the formal death point, because gamma_eff rises
      continuously and the barrier degrades all the way up to it. The loss
      fraction w_loss/w_death should approach 1 from below as Omega grows.
"""

from __future__ import annotations

import numpy as np

from ..reactions import Reaction, ReactionNetwork

__all__ = ["am_fueled", "gamma_effective", "death_waste_fraction",
           "initial_counts", "GAMMA_C"]

#: Same bifurcation as am_reversible -- gamma_eff is the same object.
GAMMA_C = 0.5


def am_fueled(gamma_inf: float = 1.0, k: float = 1.0) -> ReactionNetwork:
    """AM with every step coupled to F -> W. Species order (X, Y, B, F, W)."""
    if not np.isfinite(gamma_inf) or gamma_inf <= 0:
        raise ValueError(f"gamma_inf must be finite and > 0, got {gamma_inf}")
    kr = k * gamma_inf
    return ReactionNetwork(
        species=["X", "Y", "B", "F", "W"],
        reactions=[
            Reaction({"X": 1, "Y": 1, "F": 1}, {"B": 2, "W": 1}, k, name="f1"),
            Reaction({"B": 1, "X": 1, "F": 1}, {"X": 2, "W": 1}, k, name="f2"),
            Reaction({"B": 1, "Y": 1, "F": 1}, {"Y": 2, "W": 1}, k, name="f3"),
            Reaction({"B": 2, "W": 1}, {"X": 1, "Y": 1, "F": 1}, kr, name="r1"),
            Reaction({"X": 2, "W": 1}, {"B": 1, "X": 1, "F": 1}, kr, name="r2"),
            Reaction({"Y": 2, "W": 1}, {"B": 1, "Y": 1, "F": 1}, kr, name="r3"),
        ],
        name=f"am-fueled-ginf{gamma_inf}",
    )


def gamma_effective(n_f, n_w, gamma_inf: float = 1.0):
    """The drive the chemistry feels: gamma_inf * w / f.

    Returns inf when the tank is empty -- the reverses then dominate outright and
    there is no landscape, which is a different statement from gamma = 1 and
    should not be silently clipped to it.
    """
    n_f = np.asarray(n_f, dtype=float)
    n_w = np.asarray(n_w, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        g = gamma_inf * n_w / n_f
    return np.where(n_f > 0, g, np.inf)


def death_waste_fraction(gamma_inf: float = 1.0) -> float:
    """Waste fraction at which gamma_eff reaches GAMMA_C and restoration ends.

    gamma_inf * w / (Phi - w) = 1/2  =>  w/Phi = 1/(1 + 2 gamma_inf).
    1/3 at gamma_inf = 1. Depends on nothing else -- not Omega, not Phi.
    """
    return float(1.0 / (1.0 + 2.0 * gamma_inf))


def initial_counts(omega: int, phi: int, waste0: int = 0,
                   committed_bias: float = 0.0, gamma_inf: float = 1.0):
    """Start state (X, Y, B, F, W) at the attractor of the initial gamma_eff.

    `waste0` seeds the tank with some waste so gamma_eff starts positive; with
    waste0 = 0 the drive is infinite (gamma_eff = 0) and the landscape is as deep
    as it can be, which is the natural "full tank" condition.
    """
    from .am_reversible import delta_star
    g0 = float(gamma_effective(phi - waste0, waste0, gamma_inf))
    g0 = min(g0, 0.999) if np.isfinite(g0) else 0.0
    nb = int(round(omega * g0 / (1.0 + g0)))
    rest = omega - nb
    d = delta_star(g0) if g0 < GAMMA_C else 0.0
    sep = int(round((d + committed_bias) * omega))
    if (rest - sep) % 2:
        sep -= 1
    nx = (rest + sep) // 2
    ny = rest - nx
    return np.array([nx, ny, nb, phi - waste0, waste0], dtype=np.int64)
