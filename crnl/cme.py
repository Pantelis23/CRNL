"""Exact chemical master equation on a conserved simplex.

For a network whose reactions all conserve total count, the reachable state space
at total N is the simplex {n : sum(n) = N}, of size C(N + s - 1, s - 1) -- 7381
states for three species at N = 120. That is small enough to solve exactly with
sparse linear algebra, which makes the CME the PRIMARY instrument for this
project's dissipation work rather than a check on sampling.

Measured on this machine (do not replace these with rounder numbers -- an earlier
draft claimed 0.03 s at N=120 and 0.5 s at N=400, which is 7-20x optimistic):

    N=120:  enumerate 0.014 s | generator 0.18 s | stationary 0.20 s | ep_rate 1.7 s
    N=400:  enumerate 0.068 s | generator 2.05 s | stationary 2.45 s | ep_rate 10.6 s

Still decisive against the alternative: a direct SSA measurement of one rare-flip
lifetime at N=120 takes hundreds of hours.

Provides the generator, the stationary distribution, the Schnakenberg entropy
production rate, and first-passage quantities (crnl/cme.py is network-agnostic;
the reversible-AM specifics live in networks/am_reversible.py).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .reactions import ReactionNetwork


def enumerate_states(n_species: int, total: int) -> tuple[np.ndarray, dict]:
    """All non-negative integer vectors of length n_species summing to `total`."""
    if n_species < 1:
        raise ValueError("need at least one species")

    def rec(remaining: int, slots: int):
        if slots == 1:
            yield (remaining,)
            return
        for first in range(remaining + 1):
            for rest in rec(remaining - first, slots - 1):
                yield (first,) + rest

    states = np.array(list(rec(total, n_species)), dtype=np.int64)
    index = {tuple(s): i for i, s in enumerate(states)}
    return states, index


def _transitions(net: ReactionNetwork, total: int, omega: float):
    """Yield (i, j_state, reaction, rate) for every firable jump."""
    states, index = enumerate_states(net.n_species, total)
    S = net.stoichiometry_matrix().astype(np.int64)
    cs = net.stochastic_constants(omega)
    for i, n in enumerate(states):
        a = net.propensities(n, cs)
        for j in range(net.n_reactions):
            if a[j] <= 0.0:
                continue
            n2 = n + S[:, j]
            if n2.min() < 0:
                continue
            yield i, index[tuple(n2)], j, float(a[j])


def generator(net: ReactionNetwork, total: int, omega: float) -> sp.csr_matrix:
    """Master-equation generator Q: Q[i,k] = rate i->k, rows summing to zero."""
    states, _ = enumerate_states(net.n_species, total)
    m = len(states)
    rows, cols, vals = [], [], []
    diag = np.zeros(m)
    for i, k, _j, rate in _transitions(net, total, omega):
        rows.append(i)
        cols.append(k)
        vals.append(rate)
        diag[i] -= rate
    rows.extend(range(m))
    cols.extend(range(m))
    vals.extend(diag.tolist())
    return sp.csr_matrix((vals, (rows, cols)), shape=(m, m))


def stationary(net: ReactionNetwork, total: int, omega: float) -> np.ndarray:
    """Stationary distribution: the normalised left kernel of Q.

    Solved as a linear system with one balance equation replaced by the
    normalisation sum(p) = 1, which is well conditioned at these sizes (measured
    max|pQ| ~ 1e-16 at N=400) and avoids an eigensolver.

    Requires an irreducible chain. For reversible AM that means total >= 3: with
    N <= 1 every state is absorbing (all reactions are bimolecular) and at N = 2
    the chain splits into three closed classes. Without this guard the solve
    returns all-NaN behind nothing but a MatrixRankWarning.
    """
    if total < 3:
        raise ValueError(
            f"total={total} gives a reducible chain (every state absorbing for "
            "N<=1; three closed classes at N=2). Stationary distribution "
            "undefined; use total >= 3."
        )
    Q = generator(net, total, omega)
    m = Q.shape[0]
    A = Q.T.tolil()
    A[0, :] = 1.0                     # replace one equation with sum(p) = 1
    b = np.zeros(m)
    b[0] = 1.0
    p = np.asarray(spla.spsolve(A.tocsr(), b)).ravel()
    if not np.all(np.isfinite(p)):
        raise RuntimeError(
            f"stationary solve failed at total={total} (non-finite result); "
            "the chain is probably reducible"
        )
    # Round-off leaves many tiny negative entries (up to ~89% of them at large N,
    # all at the 1e-16 level). Clip them, but refuse to hide a real negative.
    worst = float(-p.min()) if p.min() < 0 else 0.0
    if worst > 1e-9:
        raise RuntimeError(
            f"stationary solve produced a negative probability of {worst:.3e} "
            "at total={total}; the solve is not trustworthy"
        )
    p = np.maximum(p, 0.0)
    return p / p.sum()


def ep_rate(net: ReactionNetwork, total: int, omega: float,
            pairing: np.ndarray) -> float:
    """Stationary entropy-production rate, sum_n p(n) sum_j a_j(n) ln(a_j/a_rev).

    Equals the Schnakenberg form A*J for a single-cycle network, is non-negative,
    and vanishes exactly at detailed balance.
    """
    from .thermo import entropy_step

    states, _ = enumerate_states(net.n_species, total)
    cs = net.stochastic_constants(omega)
    S = net.stoichiometry_matrix().astype(np.int64)
    p = stationary(net, total, omega)
    sigma = 0.0
    for i, n in enumerate(states):
        if p[i] <= 0.0:
            continue
        a = net.propensities(n, cs)
        for j in range(net.n_reactions):
            if a[j] <= 0.0:
                continue
            n2 = n + S[:, j]
            if n2.min() < 0:
                continue
            if net.propensities(n2, cs)[int(pairing[j])] <= 0.0:
                continue
            sigma += p[i] * a[j] * entropy_step(net, pairing, j, n, n2, cs)
    return float(sigma)
