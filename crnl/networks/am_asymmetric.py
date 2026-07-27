"""Tilted reversible AM: what does an asymmetric landscape cost, and who pays?

Every network in this project so far has been symmetric under relabelling the
symbols, so the two attractors are mirror images and one wall coefficient
describes both. Real chemistry is not like that. This adds a second knob, a tilt
`beta` on the two autocatalytic branches:

    f1: X + Y -> 2B      k                r1: 2B -> X + Y      gamma*k
    f2: B + X -> 2X      k*(1+beta)       r2: 2X -> B + X      gamma*k*(1+beta)
    f3: B + Y -> 2Y      k*(1-beta)       r3: 2Y -> B + Y      gamma*k*(1-beta)

X is simply the better catalyst. Nothing else changes.

THE POINT OF SCALING EACH REVERSE BY ITS OWN FORWARD. Every reversible pair
keeps the ratio k_f/k_r = 1/gamma, so the cycle affinity stays A = -3 ln gamma
for every beta (`test_affinity_is_beta_independent` checks this against the
general `cycle_affinity`, not against this docstring). The Wegscheider product
is k^3 (1-beta^2) forward against gamma^3 k^3 (1-beta^2) reverse, so the
condition is still gamma^3 = 1 and equilibrium is still gamma = 1 alone.

**The tilt is therefore free in the thermodynamic FORCE. It is not free in the
dissipation RATE**, which is force times flux and the flux does move. Those are
different quantities and this project has already published one result that
confused a fixed knob with a fixed effect (FINDINGS 9.2's threshold, 10.3's
control rails). `ep_rate` measures the rate; nothing here should be read as
claiming asymmetry is thermodynamically free.

WHAT THE TILT DOES TO THE LANDSCAPE. It breaks the pitchfork into an imperfect
bifurcation. At beta = 0 there are three interior fixed points below gamma_c:
two attractors and a saddle at (1/3,1/3,1/3). Raising beta deepens the X basin,
shrinks the Y basin, and slides the saddle toward Y until saddle and Y-attractor
annihilate in a saddle-node at `beta_c(gamma)`. Above beta_c the network is
monostable: **it answers X no matter what it is shown.** That is a hard ceiling
on usable tilt, and it is the first thing any asymmetric-landscape claim has to
be checked against -- a "restorer" operating past beta_c restores nothing, it
just reports a constant, which would look like perfect fidelity to any metric
that only ever feeds it X.

WHY THE SYMMETRIC POINT MOVES OFF (1/3,1/3,1/3) IMMEDIATELY. At x = y = b = 1/3,

    d(x-y)/dt = 2*beta*x*(b - gamma*x) = (2*beta/9)*(1 - gamma),

which is nonzero for every beta > 0, gamma < 1. So there is no residual
symmetric state to expand around and the closed forms in `am_reversible`
(delta_star, lambda_antisym) do not carry over -- they are recovered here only
as the beta -> 0 limit, which `test_reduces_to_symmetric` pins numerically.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq, fsolve

from ..reactions import Reaction, ReactionNetwork
from .am_reversible import GAMMA_C

__all__ = [
    "am_asymmetric",
    "drift",
    "interior_fixed_points",
    "beta_critical",
    "basin_boundary",
]


def _check(gamma: float, beta: float) -> None:
    """Reject parameters outside the domain where the network is well posed.

    beta = 1 exactly is excluded, not clipped: it sets k_f3 = 0, which removes
    Y's autocatalysis entirely and with it the reversible pair (f3, r3), so the
    cycle space collapses and the affinity is no longer -3 ln gamma. That is a
    different network, not an extreme case of this one.
    """
    if not math.isfinite(gamma) or gamma < 0:
        raise ValueError(f"gamma must be finite and >= 0, got {gamma}")
    if not math.isfinite(beta) or not (-1.0 < beta < 1.0):
        raise ValueError(
            f"beta must be finite and strictly inside (-1, 1), got {beta}; "
            "at |beta| = 1 one autocatalytic pair vanishes and the cycle "
            "affinity is no longer -3 ln gamma"
        )


def am_asymmetric(gamma: float, beta: float, k: float = 1.0) -> ReactionNetwork:
    """Reversible AM with the autocatalytic branches tilted by `beta`.

    beta > 0 favours X, beta < 0 favours Y, beta = 0 reproduces
    `am_reversible(gamma, k)` reaction for reaction.
    """
    _check(gamma, beta)
    kx, ky = k * (1.0 + beta), k * (1.0 - beta)
    return ReactionNetwork(
        species=["X", "Y", "B"],
        reactions=[
            Reaction({"X": 1, "Y": 1}, {"B": 2}, k, name="f1:X+Y->2B"),
            Reaction({"B": 1, "X": 1}, {"X": 2}, kx, name="f2:B+X->2X"),
            Reaction({"B": 1, "Y": 1}, {"Y": 2}, ky, name="f3:B+Y->2Y"),
            Reaction({"B": 2}, {"X": 1, "Y": 1}, gamma * k, name="r1:2B->X+Y"),
            Reaction({"X": 2}, {"B": 1, "X": 1}, gamma * kx, name="r2:2X->B+X"),
            Reaction({"Y": 2}, {"B": 1, "Y": 1}, gamma * ky, name="r3:2Y->B+Y"),
        ],
        name=f"am-asymmetric-gamma{gamma}-beta{beta}",
    )


def drift(xy, gamma: float, beta: float) -> np.ndarray:
    """Mass-action velocity (dx/dt, dy/dt) on the simplex, with b = 1 - x - y.

    Written out rather than taken from `deterministic.rhs` so that the fixed
    point solves below do not depend on the reaction ordering above; the two are
    cross-checked in `test_drift_matches_engine`.
    """
    x, y = float(xy[0]), float(xy[1])
    b = 1.0 - x - y
    kx, ky = 1.0 + beta, 1.0 - beta
    dx = -x * y + kx * b * x + gamma * b * b - gamma * kx * x * x
    dy = -x * y + ky * b * y + gamma * b * b - gamma * ky * y * y
    return np.array([dx, dy])


def _jac(xy, gamma: float, beta: float) -> np.ndarray:
    """Analytic Jacobian of `drift` in the reduced (x, y) coordinates."""
    x, y = float(xy[0]), float(xy[1])
    b = 1.0 - x - y
    kx, ky = 1.0 + beta, 1.0 - beta
    # d/dx and d/dy of b are both -1
    dxdx = -y + kx * (b - x) - 2.0 * gamma * b - 2.0 * gamma * kx * x
    dxdy = -x - kx * x - 2.0 * gamma * b
    dydx = -y - ky * y - 2.0 * gamma * b
    dydy = -x + ky * (b - y) - 2.0 * gamma * b - 2.0 * gamma * ky * y
    return np.array([[dxdx, dxdy], [dydx, dydy]])


def interior_fixed_points(gamma: float, beta: float,
                          tol: float = 1e-10) -> list[dict]:
    """Every fixed point strictly inside the simplex, with stability.

    Found by multi-start Newton from a grid and then deduplicated, rather than
    by a closed form: with the tilt the fixed-point system is a pair of
    inhomogeneous quadratics whose resultant is a quartic, and picking the right
    roots of that numerically is no more trustworthy than this. The grid is
    dense enough that the count is stable under refinement, which
    `test_fixed_point_count_is_grid_independent` checks by rerunning at double
    resolution.

    `kind` is assigned by stability and position, and deliberately does NOT
    reuse `am_reversible`'s "symmetric" label: with beta > 0 no fixed point sits
    at (1/3,1/3,1/3), so calling one of them symmetric would be false. A saddle
    is whatever has one positive and one negative eigenvalue.
    """
    _check(gamma, beta)
    seeds = []
    for gx in np.linspace(0.02, 0.96, 22):
        for gy in np.linspace(0.02, 0.96, 22):
            if gx + gy < 0.995:
                seeds.append((gx, gy))
    found: list[np.ndarray] = []
    for s in seeds:
        sol, info, ier, _ = fsolve(drift, s, args=(gamma, beta),
                                   fprime=_jac, full_output=True)
        if ier != 1:
            continue
        x, y = float(sol[0]), float(sol[1])
        b = 1.0 - x - y
        if min(x, y, b) <= 1e-7:            # interior only
            continue
        if np.max(np.abs(drift(sol, gamma, beta))) > 1e-9:
            continue
        if any(np.hypot(x - f[0], y - f[1]) < 1e-6 for f in found):
            continue
        found.append(np.array([x, y]))

    out = []
    for f in sorted(found, key=lambda v: v[0]):
        ev = np.linalg.eigvals(_jac(f, gamma, beta))
        npos = int(np.sum(ev.real > tol))
        kind = {0: "attractor", 1: "saddle", 2: "repeller"}.get(npos, "?")
        out.append({"x": float(f[0]), "y": float(f[1]),
                    "b": float(1.0 - f[0] - f[1]),
                    "eigenvalues": np.sort_complex(ev),
                    "kind": kind, "stable": npos == 0})
    return out


def beta_critical(gamma: float, hi: float = 0.999) -> float:
    """Tilt at which the Y attractor and the saddle annihilate.

    Bisection on the interior fixed-point COUNT (3 below, 1 above). Bisecting on
    a count rather than on a smooth residual is the honest move here: near the
    saddle-node the two merging points are separated by O(sqrt(beta_c - beta)),
    so any distance-based criterion resolves the fold to only about the square
    of the bisection tolerance, while the count is exact until the two roots
    actually coincide.

    Returns `nan` when the landscape is already monostable at beta = 0, i.e.
    gamma >= GAMMA_C -- there is no fold to find, and returning 0.0 would be
    read as "any tilt destroys it" rather than "there was nothing there".
    """
    _check(gamma, 0.0)
    if gamma >= GAMMA_C:
        return float("nan")
    if len(interior_fixed_points(gamma, 0.0)) != 3:
        raise RuntimeError(
            f"expected 3 interior fixed points at gamma={gamma}, beta=0"
        )
    if len(interior_fixed_points(gamma, hi)) == 3:
        return float("nan")            # bistable all the way to |beta| -> 1

    def f(b: float) -> float:
        return 1.0 if len(interior_fixed_points(gamma, b)) == 3 else -1.0

    return float(brentq(f, 0.0, hi, xtol=1e-6))


def basin_boundary(gamma: float, beta: float) -> float:
    """The saddle's x - y, the natural decision coordinate's zero crossing.

    At beta = 0 this is 0 by symmetry and "X wins iff x > y". With a tilt the
    dividing line moves, and this returns where to. For beta > 0 it is NEGATIVE:
    at gamma = 0.2, beta = 0.2 the saddle sits at x - y = -0.177, meaning X still
    wins from a starting Y-majority of up to 0.177. That is a systematic bias,
    a different failure mode from the random error the wall protects against, and
    the two are easy to confuse because both show up as "wrong answer".

    Note the bias is carried almost entirely by the SADDLE, not the attractors:
    at gamma = 0.05 the attractors sit at (0.952, 0.000) for every beta tested
    while the boundary moves to -0.141. Reading a tilt off the attractor
    positions would report no tilt at all.

    Returns `nan` when there is no saddle (monostable, beta >= beta_c).
    """
    pts = interior_fixed_points(gamma, beta)
    sad = [p for p in pts if p["kind"] == "saddle"]
    if len(sad) != 1:
        return float("nan")
    return float(sad[0]["x"] - sad[0]["y"])
