"""Vectorized propensity / SSA path (spec 2026-07-24 §3.2).

At n~100 the n-winner network has ~5000 reactions; the readable reference loop
in reactions.py is too slow for thousands of trials. This path computes the
IDENTICAL propensity vector with a few NumPy ops. The reference stays the
correctness oracle (see tests). It does NOT invent a new units convention.

Propensity identity used: a_j = c_j * prod_i C(n_i, R_ij), where C is the
binomial coefficient (n_i choose R_ij) = fallingfactorial(n_i, R_ij) / R_ij!.
This equals the reference's `comb / coeff!` per reactant species exactly.
Vectorization: the falling factorial is built in kmax (<=max reactant coeff)
whole-array passes, then scattered into per-reaction products with
np.multiply.at.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import factorial

import numpy as np

from .reactions import ReactionNetwork


@dataclass
class Compiled:
    S: np.ndarray            # (n_species, n_reactions) int64
    cs: np.ndarray           # (n_reactions,) float
    n_species: int
    n_reactions: int
    # reactant structure, flattened over (species, reaction) nonzero entries:
    react_sp: np.ndarray     # species index of each reactant entry
    react_rx: np.ndarray     # reaction index of each reactant entry
    react_coeff: np.ndarray  # reactant coefficient (int) of each entry
    coeff_fact: np.ndarray   # factorial(coeff) per entry (float, precomputed)
    kmax: int                # max reactant coefficient across the network


def compile_network(net: ReactionNetwork, omega: float) -> Compiled:
    R = net.reactant_matrix().astype(np.int64)     # (n_species, n_reactions)
    S = net.stoichiometry_matrix().astype(np.int64)
    cs = np.asarray(net.stochastic_constants(omega), dtype=float)
    sp, rx = np.nonzero(R)
    coeff = R[sp, rx]
    coeff_fact = np.array([float(factorial(int(c))) for c in coeff])
    kmax = int(coeff.max()) if coeff.size else 0
    return Compiled(
        S=S, cs=cs, n_species=net.n_species, n_reactions=net.n_reactions,
        react_sp=sp, react_rx=rx, react_coeff=coeff, coeff_fact=coeff_fact,
        kmax=kmax,
    )


def propensities_fast(compiled: Compiled, n: np.ndarray) -> np.ndarray:
    n = np.asarray(n, dtype=np.int64)
    counts = n[compiled.react_sp].astype(np.float64)   # count per reactant entry
    coeffs = compiled.react_coeff
    # falling factorial: prod_{d=0}^{coeff-1} (count - d), built in kmax passes
    num = np.ones_like(counts)
    for d in range(compiled.kmax):
        num *= np.where(coeffs > d, counts - d, 1.0)
    combs = num / compiled.coeff_fact                  # C(count, coeff)
    a = compiled.cs.copy()
    np.multiply.at(a, compiled.react_rx, combs)        # product per reaction
    return np.clip(a, 0.0, None)
