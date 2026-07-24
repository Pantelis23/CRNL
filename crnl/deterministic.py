"""Deterministic path: mass-action ODEs (design.md §3.1).

Integrate dx/dt = S . v(x) with scipy's adaptive solver. Integrator internals
are not the lesson, so we borrow them; LSODA is the default (adaptive, switches
to implicit) so a future network with rate constants spanning orders of
magnitude does not explode and masquerade as diverging physics (design.md §5,
stiffness).

Boundary drift (design.md §5): the polynomial RHS can push a concentration
slightly negative near x = 0, exactly where the interesting saddle dynamics
live. We do NOT clamp -- clamping injects mass and breaks conservation at the
most delicate point. Instead we integrate with tight tolerances and *monitor*
the conserved total, treating drift as something to watch, not assume away.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .reactions import ReactionNetwork


@dataclass
class Trajectory:
    t: np.ndarray  # (n_t,)
    x: np.ndarray  # (n_species, n_t) concentrations in [0, 1]-ish
    species: list[str]
    conserved_drift: float  # max |sum(x) - sum(x0)| over the trajectory

    def final(self) -> np.ndarray:
        return self.x[:, -1]

    def as_dict(self) -> dict[str, np.ndarray]:
        return {s: self.x[i] for i, s in enumerate(self.species)}


def integrate(
    net: ReactionNetwork,
    x0,
    t_span: tuple[float, float] = (0.0, 200.0),
    n_eval: int = 400,
    method: str = "LSODA",
    rtol: float = 1e-9,
    atol: float = 1e-12,
) -> Trajectory:
    """Integrate the mass-action ODE from x0 over t_span.

    Returns a Trajectory carrying the conservation drift so callers can see the
    numerical health of the run rather than trust it blindly.
    """
    x0 = np.asarray(x0, dtype=float)
    S = net.stoichiometry_matrix()

    def f(_t, x):
        return S @ net.fluxes(x)

    t_eval = np.linspace(t_span[0], t_span[1], n_eval)
    sol = solve_ivp(
        f, t_span, x0, method=method, t_eval=t_eval, rtol=rtol, atol=atol,
        dense_output=False,
    )
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")

    total0 = x0.sum()
    drift = float(np.max(np.abs(sol.y.sum(axis=0) - total0)))
    return Trajectory(sol.t, sol.y, list(net.species), drift)


def jacobian(net: ReactionNetwork, x) -> np.ndarray:
    """Exact Jacobian of the RHS at x, J = S . G where G_rj = dv_r/dx_j.

    Analytic, not finite-difference. A central difference would be corrupted
    exactly at the boundary (a concentration equal to 0), because probing
    x_j - eps < 0 there gets floored to 0 and halves the one-sided derivative --
    which is precisely where AM's most interesting fixed points sit (the rails
    and the all-blank repeller). The monomial structure gives the derivative in
    closed form and stays network-agnostic:

        v_r = k_r * prod_i x_i ** R_ir
        dv_r/dx_j = k_r * R_jr * x_j ** (R_jr - 1) * prod_{i != j} x_i ** R_ir
                    (0 when R_jr == 0)
    """
    x = np.asarray(x, dtype=float)
    S = net.stoichiometry_matrix()
    R = net.reactant_matrix()  # (n_species, n_reactions)
    ks = np.array([r.k for r in net.reactions])
    n_s, n_r = R.shape

    G = np.zeros((n_r, n_s))  # dv_r / dx_j
    for r in range(n_r):
        for j in range(n_s):
            power = R[j, r]
            if power == 0:
                continue
            term = ks[r] * power
            for i in range(n_s):
                exp = R[i, r] - (1 if i == j else 0)
                if exp == 0:
                    continue
                term *= x[i] ** exp
            G[r, j] = term
    return S @ G
