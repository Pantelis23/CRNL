"""Chemical freeze-out in an exponentially expanding volume.

An extension that stretches the landscape while the reaction runs. The volume
expands as Omega(t) = Omega0 * exp(H t) (de Sitter-like, constant "Hubble" rate
H). Bimolecular propensities carry the c = k/Omega scaling, so every reaction of
a *purely bimolecular* network slows as 1/Omega(t): relative to the compile-time
Omega0,

    a0(t) = a0_state * exp(-lambda * t),   lambda = (m - 1) * H,

where m is the (common) reaction order -- for AM / n-winner, m = 2 so lambda = H.
Because the future integrated propensity from any state is then *finite*
(a0_now / lambda), the exact next-event waiting time has a closed form whose
non-solvability IS freeze-out: with probability exp(-a0_now / lambda) no further
reaction ever fires and the state is frozen, permanently.

This is the chemical analogue of cosmological freeze-out (the Gamma-vs-H
competition that set the relic dark-matter abundance and primordial helium):
slow expansion -> consensus/equilibrium is reached; fast expansion -> the
decision freezes half-made, a non-equilibrium relic. It shares the *mathematical
structure* of freeze-out, not the astrophysics.

The algorithm is exact for exponential expansion (no constant-rate-between-events
approximation, which would fail exactly at freeze-out). Restricted to networks
whose reactions all share one order so a0(t) is a single exponential; AM and
n-winner qualify.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .vectorized import Compiled, propensities_fast


@dataclass
class ExpandingResult:
    t_final: float          # physical time at termination
    n_final: np.ndarray     # integer counts (conserved for 2->2 networks)
    steps: int
    status: str             # "absorbed" | "frozen" | "budget"
    species: list

    @property
    def frozen(self) -> bool:
        return self.status == "frozen"


def next_event_time(a0_now: float, lam: float, u: float):
    """Exact waiting time for a Poisson rate decaying as a0_now * exp(-lam * s).

    Solves the inversion  integral_0^tau a0_now e^{-lam s} ds = -ln(u).
    Returns tau, or None if the event never fires (freeze-out). u ~ U(0,1).

    lam <= 0 reduces to the ordinary exponential wait -ln(u)/a0_now, so the whole
    method degrades continuously to standard Gillespie as H -> 0.
    """
    target = -np.log(u)
    if lam <= 0.0:
        return target / a0_now
    budget = a0_now / lam            # total future integrated propensity, finite
    if target >= budget:
        return None                  # freeze: the integral never reaches target
    return -np.log1p(-lam * target / a0_now) / lam


def common_order(compiled: Compiled) -> int:
    """The shared reaction order m, or raise if the network is not uniform-order.

    a0(t) is a single exponential only when every reaction scales with Omega the
    same way, i.e. shares one order. (A mix of unimolecular and bimolecular would
    make a0(t) a sum of exponentials -- out of scope for the closed form.)
    """
    order = np.zeros(compiled.n_reactions, dtype=np.int64)
    np.add.at(order, compiled.react_rx, compiled.react_coeff)
    m = int(order.max()) if order.size else 0
    if not np.all(order == m):
        raise ValueError(
            "expanding SSA needs all reactions the same order (a0(t) must be a "
            f"single exponential); got orders {sorted(set(order.tolist()))}"
        )
    return m


def gillespie_expanding(
    compiled: Compiled,
    n0,
    rng: np.random.Generator,
    hubble: float,
    max_steps: int = 20_000_000,
    species=None,
) -> ExpandingResult:
    """Exact SSA in an exponentially expanding volume (see module docstring).

    `compiled` must be compiled at the initial Omega0. `hubble` is H (rate of
    exponential expansion, in units of the reaction rate k). Terminates when the
    state reaches an absorbing corner ("absorbed"), when expansion freezes the
    reaction ("frozen"), or the step budget is hit ("budget").
    """
    if hubble < 0.0:
        raise ValueError(
            "hubble must be >= 0 (de Sitter-like expansion); a contracting volume "
            "(H<0) makes the propensity grow and needs a different inversion"
        )
    m = common_order(compiled)
    lam = (m - 1) * hubble
    S = compiled.S
    n = np.array(n0, dtype=np.int64)
    t = 0.0
    steps = 0
    status = "budget"
    n_rx = compiled.n_reactions

    while steps < max_steps:
        a = propensities_fast(compiled, n)      # propensities at Omega0
        a0_state = float(a.sum())
        if a0_state <= 0.0:
            status = "absorbed"
            break
        a0_now = a0_state * np.exp(-lam * t)    # actual rate at current time
        tau = next_event_time(a0_now, lam, rng.random())
        if tau is None:
            status = "frozen"
            break
        # reaction choice is unaffected by the overall Omega scaling (ratios fixed)
        j = int(np.searchsorted(np.cumsum(a), rng.random() * a0_state))
        if j >= n_rx:
            j = n_rx - 1
        n = n + S[:, j]
        t += tau
        steps += 1

    labels = list(species) if species is not None else [
        f"s{i}" for i in range(compiled.n_species)]
    return ExpandingResult(t, n, steps, status, labels)


def classify_freeze(result: ExpandingResult, x="X", y="Y", blank="B") -> str:
    """Classify a (frozen or absorbed) AM outcome into the decision phases.

        "X" / "Y"    -- resolved: exactly one committed species survives
        "undecided"  -- froze mid-contest: BOTH committed species still present
        "blank"      -- froze/absorbed with all committed gone

    P("undecided") vs H is the freeze-out order parameter: 0 at H=0 (always
    resolves), rising to 1 as expansion outruns consensus.
    """
    idx = {s: i for i, s in enumerate(result.species)}
    nx = result.n_final[idx[x]]
    ny = result.n_final[idx[y]]
    if nx > 0 and ny > 0:
        return "undecided"
    if nx > 0:
        return "X"
    if ny > 0:
        return "Y"
    return "blank"
