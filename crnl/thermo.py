"""Entropy production for reversible chemical reaction networks.

For a jump n -> n' via reaction rho, the medium (heat) entropy production is

    dS = ln[ a_rho(n) / a_{-rho}(n') ]

with the reverse propensity evaluated at the POST-jump state. That convention is
not cosmetic: evaluating the reverse pre-jump is state-dependent, hits log(0),
and fabricates dissipation at equilibrium.

For reversible AM there is also an exact closed form, because rank(S) = 2 leaves
a single cycle of uniform affinity:

    dS(n -> n' via rho) = ln W(n') - ln W(n) + s_rho * ln(1/gamma)
    W(n) = N! / prod(n_i!)          s_rho = +1 forward, -1 reverse

so along any trajectory

    dS_total = ln[W(n_stop)/W(n_0)] + (A/3) * (M_forward - M_reverse)

Two independently interpretable terms (Hill/Schnakenberg): a boundary term that
depends only on the endpoints, and a cycle term counting net traversals. This is
what `decompose` reports, and it means a simulation needs only an integer counter
-- no logarithms in any hot loop.

`entropy_step` remains the general primitive (the identity above is specific to
this network's symmetry, and fails for asymmetric rate constants). It requires a
COMPLETE reverse pairing and raises rather than dividing by zero.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .reactions import ReactionNetwork


def ln_multinomial(n) -> float:
    """ln W(n) = ln( N! / prod(n_i!) ), the boundary term of the decomposition."""
    n = np.asarray(n, dtype=np.int64)
    total = int(n.sum())
    return float(math.lgamma(total + 1)
                 - sum(math.lgamma(int(c) + 1) for c in n))


def entropy_step(
    net: ReactionNetwork,
    pairing: np.ndarray,
    j: int,
    n_before,
    n_after,
    cs: np.ndarray,
) -> float:
    """Medium entropy production of one jump, ln[a_j(n) / a_reverse(n')].

    The reverse propensity is evaluated AFTER the jump. Raises if reaction j has
    no reverse, or if either propensity is non-positive (which would mean the jump
    was not actually firable, or the pairing is wrong).
    """
    rev = int(pairing[j])
    if rev < 0:
        raise ValueError(
            f"reaction {j} ({net.reactions[j].name!r}) has no reverse; "
            "entropy production is undefined for an incomplete pairing"
        )
    a_fwd = float(net.propensities(np.asarray(n_before), cs)[j])
    a_rev = float(net.propensities(np.asarray(n_after), cs)[rev])
    if a_fwd <= 0.0 or a_rev <= 0.0:
        raise ValueError(
            f"non-positive propensity in entropy_step (forward {a_fwd}, "
            f"reverse {a_rev}); the jump is not firable or the pairing is wrong"
        )
    return math.log(a_fwd / a_rev)


def decompose(n0, n_stop, net_reaction_firings: float, affinity: float,
              boundary: float | None = None) -> dict:
    """Split the trajectory's MEDIUM entropy production into two contributions.

        dS_total = boundary + (affinity/3) * net_reaction_firings

    NAMING, because a second-law statement will eventually be written against
    this function: "total" means total over the TRAJECTORY, not total over the
    universe. Every quantity here is MEDIUM (heat) entropy production -- the
    sum of ln[a_rho(n)/a_{-rho}(n')] over jumps. The system term is not computed
    anywhere in this repo. So `total` is NOT the quantity the second law bounds,
    and it can legitimately be negative at a conditioned stopping time (measured:
    -12.0 k_B at Omega=60, gamma=1, from a (1/3,1/3,1/3) start stopped at
    |delta| >= 0.5). Do not "fix" such a result; it is the heat-only statement
    failing, exactly as it should.

    `net_reaction_firings` is the NET per-reaction count (forward - reverse)
    summed over all three pairs, exactly as an integer counter in a simulation
    would accumulate it and exactly what `cme.first_passage` returns. The /3
    lives HERE and only here: it converts per-reaction firings into cycle
    traversals, since one traversal is one firing of each of f1, f2, f3.
    Dividing by 3 on both sides makes every result 3x too small, and no test of
    either function alone can detect it -- see
    tests/test_cme.py::test_decomposition_composes_with_first_passage.

    `boundary` should be the EXPECTED ln W difference from the exact solve
    (cme.first_passage returns it). Passing a single representative `n_stop`
    instead is not a harmless approximation: real absorbing states carry a
    substantial blank population, and a hand-built n_B = 0 stopping state gets
    the sign of the boundary term wrong. `n_stop` is therefore only used when
    `boundary` is None, for the single-trajectory case where it is exact.
    """
    if boundary is None:
        if n_stop is None:
            raise ValueError("give either an explicit boundary or a stopping state")
        boundary = ln_multinomial(n_stop) - ln_multinomial(n0)
    cycle = (affinity / 3.0) * float(net_reaction_firings)
    return {"boundary": float(boundary), "cycle": cycle,
            "total": float(boundary) + cycle}


class FlipCounter:
    """Schmitt-trigger flip detection on the decision coordinate.

    A "flip" is the signal committing to the OPPOSITE side, not a zero
    crossing: an unhysteretic counter on sign(delta) counts every thermal
    wobble near delta=0 and overestimates the rate by orders of magnitude.

    CONVENTION, and it matters by a factor of 2: this counts ONE-WAY crossings,
    so in a symmetric bistable system `flips / T -> 1/tau`, where tau is the
    one-sided mean first-passage time that cme.first_passage computes. A round
    trip counts as TWO. (1/(2 tau) is the ROUND-TRIP rate; using it makes a
    measured tau exactly half the exact one, which was verified to sit far
    outside sampling noise -- 0.37..0.58 with 52..82 flips per point.)

    The arm MUST scale with the landscape. A fixed arm of 0.3 never arms above
    gamma = 0.473 (delta_star(0.49) = 0.187), so the counter would read zero
    forever in exactly the range where SSA is affordable -- and read it as a
    physical result rather than a protocol failure. Callers pass
    0.7 * delta_star(gamma).

    `side` starts at 0 (unarmed): the first arrival at either arm arms the
    trigger WITHOUT counting a flip, so a run starting mid-landscape does not
    manufacture one.
    """

    def __init__(self, arm: float, side: int = 0) -> None:
        if not arm > 0.0:
            raise ValueError(f"flip arm must be positive, got {arm}")
        self.arm = float(arm)
        self.side = int(side)
        self.flips = 0

    def update(self, delta: float) -> None:
        if delta >= self.arm:
            if self.side == -1:
                self.flips += 1
            self.side = +1
        elif delta <= -self.arm:
            if self.side == +1:
                self.flips += 1
            self.side = -1


@dataclass
class InstrumentedResult:
    """One instrumented trajectory.

    `net_firings` is the NET PER-REACTION count (forward - reverse) summed over
    all pairs -- the same units as cme.first_passage's `net_reaction_firings`,
    and what thermo.decompose expects. NOT divided by 3 here.
    """

    t_final: float
    n_final: np.ndarray
    steps: int
    absorbed: bool
    stopped: bool
    species: list[str]
    net_firings: int
    flips: int


def gillespie_instrumented(
    compiled,
    n0,
    rng: np.random.Generator,
    pairing: np.ndarray,
    *,
    stop=None,
    flip_arm: float | None = None,
    flip_index: tuple[int, int] = (0, 1),
    omega: float | None = None,
    max_steps: int = 10_000_000,
    t_max: float = np.inf,
    halt_before_tmax: bool = False,
    species=None,
) -> InstrumentedResult:
    """Exact SSA with an entropy-production counter, a stop predicate, and flips.

    Deliberately a COPY of vectorized.gillespie_fast, not a refactor of it. The
    fast loop is the verified hot path used by every other experiment here;
    adding branches for three new callers was measured at 3x the stated scope,
    and the tempting shortcut -- reusing the next iteration's propensity vector
    to price the current jump -- silently drops the FINAL jump's contribution, a
    systematic -A/3 bias that breaks the exact identity. The cost of copying is
    drift, and tests/test_thermo_ssa.py::test_instrumented_matches_fast_bit_for_bit
    pins it.

    No logarithms are evaluated in this loop: entropy production comes from an
    integer counter plus the closed form in the module docstring.

    `stop(n) -> bool` halts at the first satisfying state (checked on the
    initial state too). `flip_arm` enables flip counting on
    delta = (n[i] - n[j])/omega for flip_index = (i, j), and requires `omega`.

    `halt_before_tmax` matters whenever the counter is compared against a
    fixed-time-window expectation. gillespie_fast's convention -- inherited here
    so the two loops stay bit-identical -- is to APPLY the jump that crosses
    t_max, so net_firings includes one reaction beyond the window. Measured bias
    at gamma=0.3, Omega=30, t_max=2: +0.62 firings, +6.7%, which is z = +13 at
    20000 trials but sits inside 4 SEM at 800 -- i.e. a comparison against an
    exact <M(t)> passes at low trial counts and FAILS as you add trials, sending
    the reader after the innocent quadrature. Set True to stop at the last jump
    that fits entirely within t_max. Default False preserves the bit-for-bit
    identity pinned by test_instrumented_matches_fast_bit_for_bit.
    """
    from .vectorized import propensities_fast

    pairing = np.asarray(pairing)
    if (pairing < 0).any():
        raise ValueError(
            "incomplete reverse pairing: entropy production is undefined. "
            "An unpaired reaction would be silently counted as a REVERSE "
            "firing (pairing[j] > j is False for -1), corrupting net_firings."
        )
    if flip_arm is not None and omega is None:
        raise ValueError("flip counting needs omega to form delta = (n_i - n_j)/omega")

    S = compiled.S
    n = np.array(n0, dtype=np.int64)
    t = 0.0
    steps = 0
    absorbed = False
    stopped = False
    n_rx = compiled.n_reactions

    # Forward half of each reversible pair, by the same convention as
    # cme.first_passage: j is forward iff its partner has a higher index.
    forward = np.array([pairing[j] > j for j in range(n_rx)], dtype=bool)
    net_firings = 0

    fi, fj = flip_index
    fc = FlipCounter(flip_arm) if flip_arm is not None else None
    if fc is not None:
        fc.update(float(int(n[fi]) - int(n[fj])) / float(omega))

    if stop is not None and stop(n):
        stopped = True

    while not stopped and steps < max_steps and t < t_max:
        a = propensities_fast(compiled, n)
        a0 = a.sum()
        if a0 <= 0.0:
            absorbed = True
            break
        tau = -np.log(rng.random()) / a0
        j = int(np.searchsorted(np.cumsum(a), rng.random() * a0))
        if j >= n_rx:
            j = n_rx - 1
        if halt_before_tmax and t + tau > t_max:
            t = t_max
            break
        n = n + S[:, j]
        t += tau
        steps += 1
        net_firings += 1 if forward[j] else -1
        if fc is not None:
            fc.update(float(int(n[fi]) - int(n[fj])) / float(omega))
        if stop is not None and stop(n):
            stopped = True

    labels = list(species) if species is not None else [
        f"s{i}" for i in range(compiled.n_species)]
    return InstrumentedResult(
        t_final=t, n_final=n, steps=steps, absorbed=absorbed, stopped=stopped,
        species=labels, net_firings=int(net_firings),
        flips=(fc.flips if fc is not None else 0),
    )
