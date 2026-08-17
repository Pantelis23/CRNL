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

The advantage over SSA is NOT primarily wall-clock, and an earlier version of this
docstring overstated it ("hundreds of hours at N=120"). Measured: SSA runs at ~84,000
steps/s and ~0.4*N steps per unit time, so one flip at N=120, gamma=0.35 costs 5.5
minutes. The real advantages are that this is exact (no sampling error) and that ONE
solve yields the whole first-passage field. SSA only becomes hopeless at strong drive
(58 h per flip at N=120, gamma=0.30; 3e5 h at gamma=0.25) -- which is exactly where
first_passage's own validity guard rejects the solve, so that corner belongs to
neither instrument yet.

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
            f"at total={total}; the solve is not trustworthy"
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


#: Mean first-passage times above this are not trustworthy: the linear solve
#: loses all precision and can return a NEGATIVE time (observed -6.25e16 at
#: N=30, gamma=0.02 with the Part-B absorbing set). Such points are reported
#: invalid, never silently fitted. 1e10 rather than something larger because
#: measured solves degrade well before the sign flips: rows at tau ~ 1e12 came
#: back with ~100% relative residual while still looking finite and positive.
MFPT_MAX = 1e10

#: Reject a solve whose relative residual ||Mx-b||/||b|| exceeds this. Healthy
#: solves measure 1e-14..1e-10; the untrustworthy ones measure 1e-1..1e0, so
#: there is no ambiguity and no tuning to do.
RESIDUAL_MAX = 1e-8


def first_passage(
    net: ReactionNetwork,
    total: int,
    omega: float,
    start,
    is_absorbing,
    pairing: np.ndarray | None = None,
) -> dict:
    """Mean first-passage time and splitting probability to an absorbing set.

    Solves the standard system on the transient states:

        sum_k Q[i,k] * T[k] = -1        (mean time)
        sum_k Q[i,k] * h[k] = 0         (harmonic; h = 1 on the favoured
                                         absorbing states, 0 on the others)

    Returns `valid=False` with the residual when the solve is untrustworthy
    rather than returning a plausible-looking number. Validity is judged on BOTH
    the relative residual and the magnitude: a computed-then-ignored residual is
    how rows with ~100% solve error end up published as fact. See MFPT_MAX and
    RESIDUAL_MAX.
    """
    states, index = enumerate_states(net.n_species, total)
    m = len(states)
    absorbing = np.array([bool(is_absorbing(s)) for s in states])
    if not absorbing.any():
        raise ValueError("no absorbing states: first passage is undefined")
    if absorbing.all():
        raise ValueError("every state is absorbing: nothing to solve")

    # Slice the CSR directly. Converting to lil first is NOT wrong but is
    # catastrophically slow: 48.6 s versus 0.002 s for the identical result at
    # N=400, a 24,000x penalty on the operation this module exists to make cheap.
    Q = generator(net, total, omega)
    trans = np.where(~absorbing)[0]
    tmap = {int(i): r for r, i in enumerate(trans)}
    Qtt = Q[trans][:, trans].tocsr()

    # mean time
    b = -np.ones(len(trans))
    T = spla.spsolve(Qtt, b)
    residual = float(np.linalg.norm(Qtt @ T - b) / np.linalg.norm(b))

    # splitting probability: favour absorbing states with n_X > n_Y
    favoured = np.array([absorbing[i] and states[i][0] > states[i][1]
                         for i in range(m)])
    Qta = Q[trans][:, np.where(absorbing)[0]].tocsr()
    fav_a = favoured[absorbing].astype(float)
    h = spla.spsolve(Qtt, -(Qta @ fav_a))

    start = np.asarray(start, dtype=np.int64)
    si = index[tuple(start)]
    if absorbing[si]:
        return {"mean_time": 0.0, "split": float(favoured[si]),
                "net_reaction_firings": 0.0, "boundary": 0.0,
                "valid": True, "residual": 0.0}
    r = tmap[si]
    mean_time = float(T[r])
    valid = bool(np.isfinite(mean_time)
                 and 0.0 < mean_time <= MFPT_MAX
                 and residual <= RESIDUAL_MAX)

    # Expected NET PER-REACTION FIRINGS (forward - reverse) before absorption.
    # For any quantity accumulated at state-dependent rate w_i, the expectation M
    # solves Qtt @ M = -w -- the same system shape as the mean time (w = 1).
    # Summing w against T instead would be a different object entirely.
    #
    # DO NOT divide by 3 here. `thermo.decompose` applies the affinity/3 that
    # converts per-reaction firings into cycle traversals, and dividing in both
    # places makes every dissipation number exactly 3x too small. The key name
    # says "reaction_firings" so the units cannot be misread.
    if pairing is None:
        net_firings = float("nan")
    else:
        forward = {j for j in range(net.n_reactions) if pairing[j] > j}
        w = np.zeros(len(trans))
        cs = net.stochastic_constants(omega)
        S = net.stoichiometry_matrix().astype(np.int64)
        for i in trans:
            a = net.propensities(states[i], cs)
            acc = 0.0
            for j in range(net.n_reactions):
                if a[j] <= 0.0:
                    continue
                if (states[i] + S[:, j]).min() < 0:
                    continue
                acc += a[j] * (1.0 if j in forward else -1.0)
            w[tmap[int(i)]] = acc
        M = spla.spsolve(Qtt, -w)
        net_firings = float(M[r]) if valid else float("nan")

    # Expected boundary term E[ln W(n_absorbed)] - ln W(n_0), from the same
    # harmonic solve as `split` but with the absorbing states valued by ln W
    # instead of a 0/1 indicator. Building a representative stopping state by
    # hand instead gets the SIGN wrong: the real absorbing states carry a
    # substantial blank population, which a hand-built n_B = 0 state does not.
    from .thermo import ln_multinomial

    lnw_a = np.array([ln_multinomial(states[i]) for i in np.where(absorbing)[0]])
    g = spla.spsolve(Qtt, -(Qta @ lnw_a))
    boundary = float(g[r] - ln_multinomial(start)) if valid else float("nan")

    return {"mean_time": mean_time, "split": float(h[r]),
            "net_reaction_firings": net_firings, "boundary": boundary,
            "valid": valid, "residual": residual}


def splitting_probability(net: ReactionNetwork, total: int, omega: float,
                          start, is_absorbing, is_favoured) -> dict:
    """P(absorbed in the favoured set | start), for an arbitrary favoured set.

    `first_passage` hardcodes its favoured set as `n[0] > n[1]`, which is exactly
    right for two committed species and silently WRONG for more: in a 3-winner
    race a state where X3 has won can still satisfy n_X1 > n_X2 and would be
    scored as a success. This takes the predicate instead, so the caller says what
    winning means.

    Solves the same harmonic system on the transient states,
    `Q_tt h = -Q_ta f`, and returns `valid=False` with the residual rather than a
    plausible-looking number when the solve is untrustworthy.
    """
    states, index = enumerate_states(net.n_species, total)
    absorbing = np.array([bool(is_absorbing(s)) for s in states])
    if not absorbing.any():
        raise ValueError("no absorbing states: splitting probability undefined")
    if absorbing.all():
        raise ValueError("every state is absorbing: nothing to solve")

    Q = generator(net, total, omega)
    trans = np.where(~absorbing)[0]
    tmap = {int(i): r for r, i in enumerate(trans)}
    Qtt = Q[trans][:, trans].tocsr()
    fav = np.array([bool(is_favoured(s)) for s in states])[absorbing].astype(float)
    Qta = Q[trans][:, np.where(absorbing)[0]].tocsr()
    rhs = -(Qta @ fav)
    h = spla.spsolve(Qtt, rhs)
    residual = float(np.linalg.norm(Qtt @ h - rhs) / max(np.linalg.norm(rhs), 1e-300))

    si = index[tuple(np.asarray(start, dtype=np.int64))]
    if absorbing[si]:
        return {"split": float(fav[absorbing[:si + 1].sum() - 1]),
                "valid": True, "residual": 0.0}
    val = float(h[tmap[si]])
    valid = bool(np.isfinite(val) and -1e-9 <= val <= 1 + 1e-9
                 and residual <= RESIDUAL_MAX)
    return {"split": min(max(val, 0.0), 1.0), "valid": valid,
            "residual": residual}


def first_passage_moments(net: ReactionNetwork, total: int, omega: float,
                          start, is_absorbing) -> dict:
    """Mean, second moment and VARIANCE of the first-passage time, exactly.

    `first_passage` already solves `Qtt T = -1` for the mean. The second moment
    satisfies the same system shape with the mean as its source:

        Qtt @ m2 = -2 T

    which is the k = 2 case of the standard recursion `Qtt @ m_k = -k m_{k-1}`.
    Variance is then `m2 - T^2`.

    This exists because the VARIANCE of the decision time is a pure-noise
    observable -- the deterministic limit gives exactly zero -- which makes it
    usable as a reference where a mean is not: a mean is dominated by the drift
    that every approximation retains, so it cannot discriminate between them
    (FINDINGS 25). `mean_time` here is solved identically to `first_passage`'s and
    the two are pinned together by test.
    """
    states, index = enumerate_states(net.n_species, total)
    absorbing = np.array([bool(is_absorbing(s)) for s in states])
    if not absorbing.any():
        raise ValueError("no absorbing states: first passage is undefined")
    if absorbing.all():
        raise ValueError("every state is absorbing: nothing to solve")

    Q = generator(net, total, omega)
    trans = np.where(~absorbing)[0]
    tmap = {int(i): r for r, i in enumerate(trans)}
    Qtt = Q[trans][:, trans].tocsr()

    b = -np.ones(len(trans))
    T = spla.spsolve(Qtt, b)
    res_mean = float(np.linalg.norm(Qtt @ T - b) / np.linalg.norm(b))

    rhs2 = -2.0 * T
    m2 = spla.spsolve(Qtt, rhs2)
    res_m2 = float(np.linalg.norm(Qtt @ m2 - rhs2) / max(np.linalg.norm(rhs2), 1e-300))

    start = np.asarray(start, dtype=np.int64)
    si = index[tuple(start)]
    if absorbing[si]:
        return {"mean_time": 0.0, "second_moment": 0.0, "var_time": 0.0,
                "std_time": 0.0, "valid": True, "residual": 0.0}
    r = tmap[si]
    mean_time, second = float(T[r]), float(m2[r])
    var = second - mean_time ** 2
    residual = max(res_mean, res_m2)
    # A negative variance is not a small number to be clipped: it means the solve
    # is untrustworthy, and reporting sqrt of it would launder that into a plot.
    valid = bool(np.isfinite(mean_time) and np.isfinite(second)
                 and 0.0 < mean_time <= MFPT_MAX
                 and var > 0.0 and residual <= RESIDUAL_MAX)
    return {"mean_time": mean_time, "second_moment": second,
            "var_time": var, "std_time": float(np.sqrt(var)) if var > 0 else float("nan"),
            "valid": valid, "residual": residual}
