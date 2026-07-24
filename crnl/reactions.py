"""Species + reaction data model.

This module owns two things that every other module depends on and that no
individual network is allowed to override:

  * the *units convention* that maps deterministic rate constants to stochastic
    ones (design.md §3.2), and
  * the *propensity builder* that counts distinct reactant combinations,
    including the homodimer 1/2 (design.md §3.3).

Getting either wrong rescales the noise without touching the deterministic
trajectory -- the most insidious bug class in the project (the ODE still looks
perfect while the comparison underneath it is meaningless). So both live here,
at engine level, and are enforced for all networks at once.

Convention, fixed once (design.md §2.1): the flux of a reaction is

    v_j = k_j * prod_i [X_i] ** reactant_coeff(i, j)

and the stoichiometric coefficient lives in S (products - reactants). So for
2A -> B, v = k[A]^2, S_A = -2, giving d[A]/dt = -2k[A]^2.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class Reaction:
    """A single mass-action reaction.

    ``reactants`` / ``products`` map species name -> integer stoichiometric
    coefficient. ``k`` is the deterministic rate constant in the convention
    of design.md §2.1 (the coefficient lives in the stoichiometry, never
    folded into k).
    """

    reactants: dict[str, int]
    products: dict[str, int]
    k: float
    name: str = ""

    @property
    def order(self) -> int:
        return sum(self.reactants.values())


@dataclass
class ReactionNetwork:
    """A chemical reaction network as data.

    The network derives *both* dynamics -- deterministic (S . v(x)) and
    stochastic (Gillespie propensities) -- from this same object. AM is the
    first network loaded into this engine; it is not the engine.
    """

    species: list[str]
    reactions: list[Reaction]
    name: str = ""
    _index: dict[str, int] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if len(set(self.species)) != len(self.species):
            raise ValueError(f"duplicate species in {self.species}")
        self._index = {s: i for i, s in enumerate(self.species)}
        for r in self.reactions:
            for s in (*r.reactants, *r.products):
                if s not in self._index:
                    raise ValueError(
                        f"reaction {r.name!r} references unknown species {s!r}"
                    )
            if r.k < 0:
                raise ValueError(f"reaction {r.name!r} has negative k={r.k}")

    # -- shape -----------------------------------------------------------------

    @property
    def n_species(self) -> int:
        return len(self.species)

    @property
    def n_reactions(self) -> int:
        return len(self.reactions)

    def index(self, species: str) -> int:
        return self._index[species]

    # -- deterministic pieces --------------------------------------------------

    def stoichiometry_matrix(self) -> np.ndarray:
        """S with shape (n_species, n_reactions): products - reactants."""
        S = np.zeros((self.n_species, self.n_reactions))
        for j, r in enumerate(self.reactions):
            for s, c in r.reactants.items():
                S[self._index[s], j] -= c
            for s, c in r.products.items():
                S[self._index[s], j] += c
        return S

    def reactant_matrix(self) -> np.ndarray:
        """Reactant stoichiometry R with shape (n_species, n_reactions).

        Used to build fluxes as monomials: v_j = k_j * prod_i x_i ** R_ij.
        """
        R = np.zeros((self.n_species, self.n_reactions))
        for j, r in enumerate(self.reactions):
            for s, c in r.reactants.items():
                R[self._index[s], j] += c
        return R

    def fluxes(self, x: np.ndarray) -> np.ndarray:
        """Deterministic mass-action flux vector v(x) (design.md §2.1).

        v_j = k_j * prod_i x_i ** R_ij. The stoichiometric power comes from the
        reactants; k carries no combinatorial factor.
        """
        x = np.asarray(x, dtype=float)
        R = self.reactant_matrix()
        ks = np.array([r.k for r in self.reactions])
        # prod_i x_i ** R_ij  -- guard against 0**0 = 1 (correct here) and any
        # tiny-negative x from integrator drift by flooring the base at 0.
        base = np.clip(x, 0.0, None)[:, None]
        monomials = np.prod(np.where(R > 0, base ** R, 1.0), axis=0)
        return ks * monomials

    def rhs(self, x: np.ndarray) -> np.ndarray:
        """dx/dt = S . v(x)."""
        return self.stoichiometry_matrix() @ self.fluxes(x)

    # -- stochastic pieces -----------------------------------------------------

    def stochastic_constants(self, omega: float) -> np.ndarray:
        """Map deterministic k -> stochastic c for a given population scale Omega.

        design.md §3.2, decided once here so no network gets a vote:

            order          deterministic   stochastic c
            unimolecular   k               k
            heterobimol.   k               k / Omega
            homodimer A+A  k               2k / Omega

        The homodimer factor 2 is *not* a separate convention from the 1/2 in
        ``propensities`` -- it is the same fact from the other side. Propensity
        for A+A is c*n(n-1)/2; each firing changes n_A by -2; the mean drift is
        -2 * c n^2 / 2 = -c n^2, i.e. d[A]/dt = -c*Omega*[A]^2. Matching the
        deterministic -2k[A]^2 forces c = 2k/Omega. Write k/Omega for a
        homodimer and you have silently halved its noise relative to its own ODE.

        Higher orders are handled by the general rule: for a reaction of total
        order m, c = k * s / Omega**(m-1), where s is the product of the
        factorials of the reactant coefficients (s = prod_i coeff_i!). This
        reproduces the table (unimolecular s=1 m=1; heterobimol. s=1 m=2;
        homodimer s=2! =2 m=2) and extends it consistently to trimolecular etc.
        """
        cs = np.empty(self.n_reactions)
        for j, r in enumerate(self.reactions):
            m = r.order
            s = 1
            for coeff in r.reactants.values():
                s *= _factorial(coeff)
            cs[j] = r.k * s / (omega ** (m - 1))
        return cs

    def propensities(self, n: np.ndarray, cs: np.ndarray) -> np.ndarray:
        """Propensity vector a_j = c_j * (# distinct reactant combinations).

        design.md §3.3. The combinatorial count is a *falling-factorial* product
        over reactant species, divided by the symmetry factor of repeated
        reactants:

            hetero A + B : c * n_A * n_B
            homodimer A+A: c * n_A (n_A - 1) / 2   (NOT c * n_A^2)

        The 1/2 here and the 2 in ``stochastic_constants`` must be owned
        together -- implement exactly one and you are off by a factor of 2 in
        the noise. AM has no homodimer and never triggers this, which is exactly
        why it must live at engine level rather than be discovered later.
        """
        n = np.asarray(n)
        a = np.array(cs, dtype=float)
        for j, r in enumerate(self.reactions):
            for s, coeff in r.reactants.items():
                ni = n[self._index[s]]
                # distinct ways to choose `coeff` molecules of species s:
                #   ni * (ni-1) * ... * (ni-coeff+1) / coeff!
                # The /coeff! symmetry factor is cancelled by the coeff! folded
                # into c (see stochastic_constants). So here we use the falling
                # factorial *without* dividing by coeff!, and the factorials in
                # c supply exactly the 1/coeff! that combinations require.
                comb = 1.0
                for d in range(coeff):
                    comb *= ni - d
                a[j] *= comb / _factorial(coeff)
        # For any valid integer state the falling factorial is exactly 0 when
        # ni < coeff (one factor (ni - d) hits 0), so the propensity is 0, never
        # negative -- this clip is a defensive no-op on integer inputs, kept only
        # to guard against non-integer states passed in by exploratory callers.
        return np.clip(a, 0.0, None)


def _factorial(k: int) -> int:
    out = 1
    for i in range(2, k + 1):
        out *= i
    return out
