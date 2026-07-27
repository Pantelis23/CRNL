"""Expansion that COOLS, so the landscape itself changes as the volume grows.

WHAT THIS FIXES. Every reaction in AM is 2 -> 2, so under expansion every
propensity scales by the identical factor and the *ratios* never move. gamma,
delta*(gamma), kappa(gamma), beta and gamma_c are all invariant under expansion:
the landscape is frozen and only the clock slows. That is why FINDINGS 5.1's
expanding SSA reduced *exactly* to ordinary SSA stopped at internal time 1/H --
it could not have come out any other way, and `expanding.common_order` enforces
the uniform-order property that makes it true.

Real freeze-out is not a clock running out. It is the EQUILIBRIUM MOVING: as the
volume expands the medium cools, the reverse (pair-production-like) reactions are
suppressed, and the departure from equilibrium is an event at a definite time
rather than a standing condition. The old model has no temperature, so it cannot
represent that, and in particular it has no relic ABUNDANCE -- only a relic sign.

THE MODEL. Adiabatic expansion gives T proportional to Omega(t)^-w (w = 1/3 for a
radiation-like medium, 2/3 for a non-relativistic one). The reverse rates carry an
activation energy, so gamma = exp(-dE/T) and

    gamma(t) = gamma0 ** exp(w * H * t)

which falls monotonically from gamma0 to 0. Forward rates are untouched; only the
balance moves, which is the whole point.

THE USEFUL CHANGE OF VARIABLE. In FINDINGS 5.1's internal time
tau = (1 - e^{-Ht})/H the dilution disappears from the forward reactions, and
with s = H*tau in [0, 1) the drive becomes

    gamma(s) = gamma0 ** ((1 - s)^-w)

**which does not depend on H at all.** So the cooling SCHEDULE is universal in s,
and H enters only as an overall rate: the same sweep of gamma from gamma0 down
through gamma_c to 0 always happens, and H decides how many reactions fit inside
it. That is the competition the fixed-gamma model cannot have -- cooling deepens
the landscape while dilution starves it -- and it is why this is not a time
change. `test_reduces_to_expanding_at_w_zero` pins the w = 0 case back onto
`gillespie_expanding` so the generalisation is anchored on the known result.

WHAT IS DELIBERATELY NOT MODELLED. The forward rates have no temperature
dependence, so this is the minimal change that lets the balance move, not a
thermochemistry. And the medium is still well mixed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .vectorized import Compiled, propensities_fast

__all__ = ["CoolingSchedule", "CoolingResult", "gillespie_cooling",
           "reverse_mask"]


def reverse_mask(net, pairing: np.ndarray) -> np.ndarray:
    """True for the reverse member of each reversible pair.

    The lower index of a pair is taken as forward -- the same arbitrary but fixed
    convention `cycle_affinity` uses. Reactions with no partner are treated as
    forward (not cooled), which is correct for an irreversible network and would
    be wrong for one whose only reverse is unpaired; there is no such network here
    and the assertion below would catch it.
    """
    mask = np.zeros(net.n_reactions, dtype=bool)
    for j in range(net.n_reactions):
        if pairing[j] >= 0 and pairing[j] < j:
            mask[j] = True
    assert mask.sum() * 2 == int((pairing >= 0).sum()), (
        "reverse_mask must select exactly half of the paired reactions")
    return mask


@dataclass(frozen=True)
class CoolingSchedule:
    """gamma(s)/gamma0 on the rescaled internal clock s = H*tau in [0, 1).

    Precomputes the cumulative integral of g(s) once, so drawing an event is a
    monotone interpolation rather than a quadrature.
    """
    gamma0: float
    w: float
    s: np.ndarray
    g: np.ndarray
    cum: np.ndarray                      # integral of g from 0 to s

    @staticmethod
    def build(gamma0: float, w: float, n_grid: int = 200_001,
              s_max: float = 1.0 - 1e-9) -> "CoolingSchedule":
        if not (0.0 < gamma0 < 1.0):
            raise ValueError(f"gamma0 must be in (0,1), got {gamma0}")
        if w < 0.0:
            raise ValueError(f"w must be >= 0, got {w}")
        s = np.linspace(0.0, s_max, n_grid)
        if w == 0.0:
            g = np.ones_like(s)
        else:
            # gamma(s)/gamma0 = gamma0**((1-s)^-w - 1); the exponent diverges as
            # s -> 1, so this underflows to 0 rather than overflowing.
            expo = np.power(1.0 - s, -w) - 1.0
            with np.errstate(over="ignore", under="ignore"):
                g = np.exp(expo * np.log(gamma0))
            g[~np.isfinite(g)] = 0.0
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1])
                                               * np.diff(s))])
        return CoolingSchedule(float(gamma0), float(w), s, g, cum)

    def g_at(self, s: float) -> float:
        return float(np.interp(s, self.s, self.g))

    def cum_at(self, s: float) -> float:
        return float(np.interp(s, self.s, self.cum))

    def gamma_at(self, s: float) -> float:
        return self.gamma0 * self.g_at(s)

    def s_of_gamma(self, gamma: float) -> float:
        """Where the drive passes a given value -- e.g. gamma_c, the bifurcation."""
        if self.w == 0.0:
            return float("nan")
        target = gamma / self.gamma0
        if target >= 1.0:
            return 0.0
        idx = np.searchsorted(-self.g, -target)      # g is decreasing
        if idx >= len(self.s):
            return float("nan")
        return float(self.s[idx])


@dataclass
class CoolingResult:
    s_final: float
    n_final: np.ndarray
    steps: int
    status: str                          # "frozen" | "absorbed" | "budget"
    gamma_final: float
    species: list

    @property
    def frozen(self) -> bool:
        return self.status == "frozen"


def gillespie_cooling(
    compiled: Compiled,
    n0,
    rng: np.random.Generator,
    hubble: float,
    schedule: CoolingSchedule,
    is_reverse: np.ndarray,
    max_steps: int = 20_000_000,
    species=None,
) -> CoolingResult:
    """Exact SSA in an expanding volume whose cooling suppresses the reverses.

    `compiled` must be compiled at Omega0 with the network at gamma0; the reverse
    reactions are then scaled by g(s) = gamma(s)/gamma0 as the run proceeds.

    Between events the state is fixed, so the integrated rate to rescaled time s is

        F * (s - s0)  +  R * (cum(s) - cum(s0))   =   E * H,     E = -ln U

    with F and R the forward and reverse propensity sums at Omega0. This is exact:
    no thinning, no rejection. The horizon is s = 1 (internal time 1/H); if the
    remaining integral cannot reach E the run is FROZEN, which is the same
    freeze-out criterion `expanding.gillespie_expanding` uses.
    """
    if hubble <= 0.0:
        raise ValueError("hubble must be > 0; use gillespie_fast for H = 0")
    is_reverse = np.asarray(is_reverse, dtype=bool)
    n = np.array(n0, dtype=np.int64)
    S = compiled.S
    s = 0.0
    steps = 0
    status = "budget"

    while steps < max_steps:
        a = propensities_fast(compiled, n)
        fwd = float(a[~is_reverse].sum())
        rev = float(a[is_reverse].sum())
        if fwd + rev <= 0.0:
            status = "absorbed"
            break
        target = -np.log(rng.random()) * hubble
        cum0 = schedule.cum_at(s)
        # remaining integral out to the horizon
        remaining = fwd * (1.0 - s) + rev * (schedule.cum[-1] - cum0)
        if remaining < target:
            status = "frozen"
            break
        # invert F*(s-s0) + R*(cum(s)-cum0) = target, monotone increasing in s.
        # Binary search over the grid: building the whole accumulated array here
        # instead is O(n_grid) PER EVENT and made a 2000-step trial cost 4e8
        # operations. Only the probed indices are ever evaluated.
        gs_arr, cum_arr = schedule.s, schedule.cum
        lo_i = int(np.searchsorted(gs_arr, s))
        hi_i = len(gs_arr) - 1
        while lo_i < hi_i:
            mid = (lo_i + hi_i) // 2
            if fwd * (gs_arr[mid] - s) + rev * (cum_arr[mid] - cum0) < target:
                lo_i = mid + 1
            else:
                hi_i = mid
        idx = max(lo_i, 1)
        a_lo = fwd * (gs_arr[idx - 1] - s) + rev * (cum_arr[idx - 1] - cum0)
        a_hi = fwd * (gs_arr[idx] - s) + rev * (cum_arr[idx] - cum0)
        frac = 0.0 if a_hi <= a_lo else (target - a_lo) / (a_hi - a_lo)
        s = float(gs_arr[idx - 1] + frac * (gs_arr[idx] - gs_arr[idx - 1]))
        if s >= schedule.s[-1]:
            status = "frozen"
            break
        # the reaction MIX is time-dependent -- this is the whole point
        gs = schedule.g_at(s)
        w = np.where(is_reverse, a * gs, a)
        tot = float(w.sum())
        if tot <= 0.0:
            status = "absorbed"
            break
        j = int(np.searchsorted(np.cumsum(w), rng.random() * tot))
        n = n + S[:, min(j, compiled.n_reactions - 1)]
        steps += 1

    labels = list(species) if species is not None else [
        f"s{i}" for i in range(compiled.n_species)]
    return CoolingResult(s, n, steps, status, schedule.gamma_at(s), labels)
