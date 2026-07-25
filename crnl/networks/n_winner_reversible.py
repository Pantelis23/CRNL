"""Reversible n-winner AM: what drive does an n-symbol landscape cost?

`am_reversible` answered this for n=2: the landscape dies at gamma_c = 1/2, so a
two-symbol landscape costs a minimum affinity A = 3 ln 2 (FINDINGS 9.1). This
generalises it, and the answer is not the obvious one.

STRUCTURE, and three ways n >= 3 differs from n = 2:

  * The network has C(n,2) disagreement reactions X_i + X_j <-> 2B and n
    recruitment reactions B + X_i <-> 2X_i, all with reverse rate gamma*k.
  * The CYCLE SPACE IS NO LONGER ONE-DIMENSIONAL. Counting each reversible pair
    as one edge -- the counting under which 9.1 called AM's cycle space 1-D --
    there are C(n,2) + n edges against rank(S) = n, so the cycle dimension is
    exactly C(n,2): 1 at n=2, then 3, 6, 15 at n = 3, 4, 6. The drive is still
    set by one parameter because gamma is uniform, but the network stops being
    the single-cycle object that made 9.1's argument so short.
    (Counting forward and reverse as separate edges instead gives 4, 9, 16 --
    a different and less useful number. The two must not be mixed.)
  * The ELEMENTARY CYCLE IS STILL THREE REACTIONS, for every n: fire
    X_i + X_j -> 2B, then B + X_i -> 2X_i, then B + X_j -> 2X_j, and every count
    returns to its starting value. So the affinity per cycle stays A = -3 ln gamma
    and the question is entirely about where gamma_c sits.
  * The SYMMETRIC POINT MOVES. At n=2 it is pinned at (1/3, 1/3, 1/3) for every
    gamma. For n >= 3 it is NOT at x = b = 1/(n+1), and it depends on gamma
    (measured at n=3: x = 0.2050 / 0.2269 / 0.2431 at gamma = 0.02 / 0.2 / 0.6).

The measured law is in FINDINGS 13: gamma_c(n) -> n^-3, hence A_c(n) -> 9 ln n.
n = 2 is NOT on that asymptote -- A_c(2) = 3 ln 2 = 2.079 where 9 ln 2 = 6.238 --
so the famous case is the special one, and the asymptotic law is approached from
below over decades in n.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy.optimize import brentq

from ..deterministic import jacobian
from ..reactions import Reaction, ReactionNetwork


def n_winner_reversible(n: int, gamma: float, k: float = 1.0) -> ReactionNetwork:
    """n committed species plus a blank, every reaction reversible.

    Reaction order is documented and fixed: all C(n,2) disagreement pairs in
    `combinations` order, then the n recruitments, then their reverses in the
    same order -- so `reverse_pairing`-style logic (forward j pairs with
    j + n_forward) holds by construction.
    """
    if n < 2:
        raise ValueError(f"n-winner AM needs n >= 2 committed species, got {n}")
    if not np.isfinite(gamma) or gamma < 0.0:
        raise ValueError(f"gamma must be finite and non-negative, got {gamma}")

    committed = [f"X{i + 1}" for i in range(n)]
    species = committed + ["B"]
    forward, reverse = [], []
    for i, j in combinations(range(n), 2):
        xi, xj = committed[i], committed[j]
        forward.append(Reaction({xi: 1, xj: 1}, {"B": 2}, k, name=f"dis:{xi}+{xj}->2B"))
        reverse.append(Reaction({"B": 2}, {xi: 1, xj: 1}, gamma * k,
                                name=f"rev-dis:2B->{xi}+{xj}"))
    for i in range(n):
        xi = committed[i]
        forward.append(Reaction({"B": 1, xi: 1}, {xi: 2}, k, name=f"rec:B+{xi}->2{xi}"))
        reverse.append(Reaction({xi: 2}, {"B": 1, xi: 1}, gamma * k,
                                name=f"rev-rec:2{xi}->B+{xi}"))
    return ReactionNetwork(species=species, reactions=forward + reverse,
                           name=f"n-winner-am-reversible-n{n}-g{gamma}")


def symmetric_state(n: int, gamma: float) -> tuple[float, float]:
    """The permutation-symmetric fixed point (x, b), all x_i = x, b = 1 - n x.

    Setting dx/dt = 0 with b = 1 - n x gives a quadratic:

        gamma (n-1) b^2 + b x - [(n-1) + gamma] x^2 = 0

    At n = 2 this is satisfied by x = 1/3 for EVERY gamma, which is why 9.1's
    symmetric point is gamma-independent. That degeneracy is specific to n = 2.
    """
    if n < 2:
        raise ValueError(f"need n >= 2, got {n}")
    a = gamma * (n - 1) * n ** 2 - n - ((n - 1) + gamma)
    b_ = -2 * gamma * (n - 1) * n + 1
    c = gamma * (n - 1)
    if abs(a) < 1e-14:
        roots = [-c / b_] if abs(b_) > 1e-14 else []
    else:
        disc = b_ * b_ - 4 * a * c
        if disc < 0:
            raise ValueError(f"no real symmetric state at n={n}, gamma={gamma}")
        roots = [(-b_ + np.sqrt(disc)) / (2 * a), (-b_ - np.sqrt(disc)) / (2 * a)]
    for x in roots:
        blank = 1.0 - n * x
        if x > 1e-13 and blank > 1e-13:
            return float(x), float(blank)
    raise ValueError(f"no interior symmetric state at n={n}, gamma={gamma}")


def breaking_mode(n: int) -> np.ndarray:
    """Unit vector for 'one species pulls ahead, the rest share the loss'.

    This is an exact eigenvector of the Jacobian at the symmetric point by
    permutation symmetry (verified numerically to a relative residual of 1e-16),
    so the Rayleigh quotient below is the eigenvalue and not merely a bound.
    """
    v = np.zeros(n + 1)
    v[0] = 1.0
    v[1:n] = -1.0 / (n - 1)
    return v / np.linalg.norm(v)


def lambda_breaking(n: int, gamma: float) -> float:
    """Growth rate of the symmetry-breaking mode. The n-winner lambda_antisym.

    Reduces to the closed form (1 - 2 gamma)/3 at n = 2 (agreement 1e-12).
    Positive => the symmetric state is unstable and a landscape exists;
    negative => a single minimum and nothing to restore toward.
    """
    x, blank = symmetric_state(n, gamma)
    net = n_winner_reversible(n, gamma)
    state = np.concatenate([np.full(n, x), [blank]])
    v = breaking_mode(n)
    return float(v @ (jacobian(net, state) @ v))


def gamma_critical(n: int, lo: float = 1e-14, hi: float = 0.999) -> float:
    """The drive at which the n-symbol landscape disappears.

    Measured to converge on gamma_c ~ n^-3 (local exponent -3.02 at n=256), with
    n=2 giving exactly 1/2. Bracketed rather than iterated from a guess, because
    gamma_c spans 1/2 down to 6e-8 over n = 2..256.
    """
    f_lo, f_hi = lambda_breaking(n, lo), lambda_breaking(n, hi)
    if f_lo * f_hi > 0:
        raise ValueError(
            f"no sign change for n={n} on [{lo}, {hi}]: "
            f"lambda={f_lo:.4g} .. {f_hi:.4g}")
    return float(brentq(lambda g: lambda_breaking(n, g), lo, hi,
                        xtol=1e-16, rtol=8.9e-16))


def affinity_critical(n: int) -> float:
    """Minimum cycle affinity that buys an n-symbol landscape, A_c = -3 ln gamma_c.

    The 3 is the elementary cycle length, which is 3 for every n (see the module
    docstring), so this is directly comparable to 9.1's A_c(2) = 3 ln 2.
    """
    return float(-3.0 * np.log(gamma_critical(n)))


def breaking_diffusion(n: int, gamma: float = 0.0) -> float:
    """Finite-count diffusion D_0 in the symmetry-breaking direction.

    The van Kampen / chemical-Langevin diffusion is D = D_0 / Omega with

        D_0 = sum_r (v . S_r)^2 v_r(x)

    evaluated at the symmetric state, with v the unit breaking mode. At
    gamma = 0 the sum has the closed form

        D_0(n) = (2n - 3) / (2n - 1)^2

    which is 1/9 at n = 2 -- exactly design.md section 9's D = 1/(9 Omega) for
    irreversible AM. Verified against this numeric sum to 7 decimals for
    n = 2..64.

    Why it matters: paired with lambda_breaking(n, 0) = 1/(2n-1), the
    quasipotential barrier c ~ lambda/(2 D_0) becomes

        lambda / D_0 = (2n - 1) / (2n - 3)   ->   1

    so lambda and D_0 vanish at the SAME rate and their ratio saturates. That is
    the mechanism behind FINDINGS 3's measured saturation of c(n), which had no
    explanation when it was measured. See FINDINGS 14 for how far it goes
    quantitatively (it does not go all the way -- a constant factor ~2.3 remains).
    """
    from .n_winner import n_winner

    if gamma == 0.0:
        net = n_winner(n)
        x = 1.0 / (2 * n - 1)
        blank = (n - 1) / (2 * n - 1)
    else:
        net = n_winner_reversible(n, gamma)
        x, blank = symmetric_state(n, gamma)
    state = np.concatenate([np.full(n, x), [blank]])
    S = net.stoichiometry_matrix()
    flux = net.fluxes(state)
    v = breaking_mode(n)
    return float(sum((v @ S[:, r]) ** 2 * flux[r] for r in range(net.n_reactions)))


def predicted_barrier(n: int, delta: float = 0.10) -> float:
    """Quasipotential prediction for the n-winner restoration barrier c(n).

        c = lambda * delta^2 / (2 D_0) = delta^2 (2n-1) / (2 (2n-3))

    -> delta^2 / 2 as n grows. At n = 2 this is 1.5 delta^2, reproducing
    design.md section 9's c(eps) = (3/2) eps^2 exactly.
    """
    return float(lambda_closed(n) * delta ** 2 / (2.0 * diffusion_closed(n)))


def lambda_closed(n: int) -> float:
    """lambda_breaking(n, gamma=0) = 1/(2n-1), in closed form.

    Exists because the numeric route builds the network, and n-winner has
    C(n,2)+n reactions -- at n=4096 that is 8.4 million, which is not a
    calculation, it is an out-of-memory error. Pinned against the numeric route
    for n = 2..64 by tests/test_n_winner_reversible.py.
    """
    return 1.0 / (2 * n - 1)


def diffusion_closed(n: int) -> float:
    """breaking_diffusion(n, gamma=0) = (2n-3)/(2n-1)^2, in closed form.

    Same reason as lambda_closed, and pinned the same way. Equals 1/9 at n=2.
    """
    return (2 * n - 3) / (2 * n - 1) ** 2
