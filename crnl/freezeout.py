"""Freeze-out as a time change: the expanding SSA *is* ordinary SSA with a clock.

`crnl/expanding.py` runs the SSA in a volume Omega(t) = Omega0 e^{Ht}, where a
purely bimolecular network's total propensity decays as a0_state * e^{-lam t}
with lam = (m-1)H = H. That construction has an exact reduction which the
original module notices but does not use:

    define the INTERNAL CLOCK   tau(t) = (1 - e^{-lam t}) / lam,   dtau = e^{-lam t} dt

then the next-event condition

    -ln u = int_0^{tau_step} a0_state e^{-lam (t+s)} ds
          = a0_state * [e^{-lam t} - e^{-lam (t+tau_step)}] / lam
          = a0_state * Delta tau

is *exactly* the ordinary Gillespie increment `Delta tau = -ln(u)/a0_state`, and
reaction choice is untouched by the overall scaling. Since tau(inf) = 1/lam,

    **the expanding SSA at rate H is ordinary SSA stopped at internal time 1/H.**

Freeze-out is not a new dynamics: it is the ordinary chain having failed to
absorb within a finite time budget `tau_max = 1/H`. `test_freezeout.py` checks
the equivalence bit-for-bit (identical seeds give identical state sequences).

Two things follow.

1. **One run measures every H.** Record the observable along the internal clock
   once and read it at tau = 1/H for as many H as you like. The old sweep paid
   for each H separately.

2. **There is no critical rate.** D(H, Omega) is the consensus-time distribution
   of ordinary AM read at 1/H, and the consensus time from a symmetric start
   diverges logarithmically in Omega: d = x - y obeys d' = b*d exactly, b -> 1/3,
   and d(0) ~ Omega^{-1/2} is shot noise, so tau* = (3/2) ln Omega + O(1).
   Hence H*(Omega) -> 0 like 1/((3/2) ln Omega). See FINDINGS.md Sec.5.1.

This module provides the fast instrument: an SSA vectorised *across trials* that
accumulates observables on an internal-time grid. It is an independent
implementation of the same Markov chain as `vectorized.gillespie_fast` (own
propensity code, own reaction-selection code), so agreement between them is a
real cross-check rather than a tautology.
"""

from __future__ import annotations

import numpy as np

from .vectorized import Compiled


def bimolecular_pairs(compiled: Compiled) -> tuple[np.ndarray, np.ndarray]:
    """Per-reaction (species_a, species_b) for a network of A + B -> ... only.

    Every reaction must have exactly two reactant entries, each with coefficient
    1 (no homodimers). AM and n-winner AM both satisfy this. Raising here rather
    than silently mis-computing propensities is the point: the fast kernel's
    whole speed advantage comes from assuming this shape.
    """
    n_rx = compiled.n_reactions
    sp_a = np.full(n_rx, -1, dtype=np.int64)
    sp_b = np.full(n_rx, -1, dtype=np.int64)
    for sp, rx, coeff in zip(compiled.react_sp, compiled.react_rx,
                             compiled.react_coeff):
        if coeff != 1:
            raise ValueError(
                f"reaction {rx} has a reactant with coefficient {coeff}; the "
                "fast freeze-out kernel handles only A + B -> ... reactions")
        if sp_a[rx] < 0:
            sp_a[rx] = sp
        elif sp_b[rx] < 0:
            sp_b[rx] = sp
        else:
            raise ValueError(f"reaction {rx} has more than two reactant species")
    if np.any(sp_b < 0):
        bad = int(np.flatnonzero(sp_b < 0)[0])
        raise ValueError(
            f"reaction {bad} is not bimolecular with two distinct reactants; "
            "the fast freeze-out kernel requires A + B -> ...")
    return sp_a, sp_b


def _scatter_fill(acc_list, val_list, k, hi, n_grid):
    """Add val to acc[k_i : hi_i] for every trial i, vectorised.

    `acc_list` are (n_grid,) accumulators, `val_list` the matching per-trial
    values held constant over that grid span (the state does not change between
    events, so the observable is piecewise constant on the internal clock).
    """
    need = hi - k
    m = need > 0
    if not m.any():
        return
    counts = need[m]
    total = int(counts.sum())
    starts = np.repeat(k[m], counts)
    within = np.arange(total) - np.repeat(np.cumsum(counts) - counts, counts)
    off = starts + within
    for acc, val in zip(acc_list, val_list):
        acc += np.bincount(off, weights=np.repeat(val[m], counts),
                           minlength=n_grid)


def internal_clock_sweep(
    compiled: Compiled,
    n0,
    taus,
    trials: int,
    rng: np.random.Generator,
    observables,
    max_iters: int = 200_000_000,
) -> dict:
    """SSA vectorised across trials, sampled on an internal-time grid.

    Runs `trials` independent copies of the chain from `n0` and, for each grid
    point tau_k, accumulates `observables(state)` evaluated at the state holding
    at time tau_k (i.e. after the last event with t <= tau_k). Because the
    expanding SSA is this chain read at tau = 1/H (see module docstring), the
    returned curves give D(H) for every H simultaneously.

    `observables(n)` takes an (n_trials, n_species) int array and returns a list
    of (n_trials,) float arrays. Returns their means on the grid, plus the SEM
    of the first one and the fraction of trials already absorbed.
    """
    taus = np.asarray(taus, dtype=float)
    if taus.ndim != 1 or taus.size == 0:
        raise ValueError("taus must be a non-empty 1-D array")
    if np.any(np.diff(taus) <= 0):
        raise ValueError("taus must be strictly increasing")
    n_grid = taus.size
    sp_a, sp_b = bimolecular_pairs(compiled)
    cs = compiled.cs
    S = compiled.S.T.astype(np.int64)          # (n_reactions, n_species)

    n = np.tile(np.asarray(n0, dtype=np.int64), (trials, 1))
    probe = observables(n)
    n_obs = len(probe)
    acc = [np.zeros(n_grid) for _ in range(n_obs)]
    acc_sq = np.zeros(n_grid)                   # for the SEM of observable 0
    acc_absorbed = np.zeros(n_grid)

    t = np.zeros(trials)
    k = np.zeros(trials, dtype=np.int64)
    iters = 0

    while n.shape[0]:
        if iters >= max_iters:
            raise RuntimeError(
                f"internal_clock_sweep hit max_iters={max_iters} with "
                f"{n.shape[0]} trials still short of tau={taus[-1]}")
        iters += 1
        a = cs[None, :] * n[:, sp_a] * n[:, sp_b]
        a0 = a.sum(axis=1)
        absorbed = a0 <= 0.0

        vals = observables(n)
        safe_a0 = np.where(absorbed, 1.0, a0)
        step = -np.log(rng.random(n.shape[0])) / safe_a0
        t_next = np.where(absorbed, np.inf, t + step)
        hi = np.searchsorted(taus, t_next, side="left")

        _scatter_fill(acc, vals, k, hi, n_grid)
        _scatter_fill([acc_sq], [vals[0] ** 2], k, hi, n_grid)
        _scatter_fill([acc_absorbed], [absorbed.astype(float)], k, hi, n_grid)
        k = hi

        keep = (~absorbed) & (k < n_grid)
        if not keep.all():
            n, t, k = n[keep], t[keep], k[keep]
            a, a0, t_next = a[keep], a0[keep], t_next[keep]
            if n.shape[0] == 0:
                break

        u = rng.random(n.shape[0]) * a0
        j = (np.cumsum(a, axis=1) < u[:, None]).sum(axis=1)
        np.clip(j, 0, compiled.n_reactions - 1, out=j)
        n = n + S[j]
        t = t_next

    mean0 = acc[0] / trials
    var0 = np.maximum(acc_sq / trials - mean0 ** 2, 0.0)
    return {
        "taus": taus,
        "means": [a_ / trials for a_ in acc],
        "sem0": np.sqrt(var0 / trials),
        "p_absorbed": acc_absorbed / trials,
        "trials": trials,
        "iters": iters,
    }


def am_observables(n: np.ndarray) -> list:
    """AM readouts matching experiments/expansion.run_point exactly.

    order  = |X-Y|/(X+Y), and 1.0 when X+Y = 0 (the all-blank corner counts as
             "decided", which is run_point's convention -- kept so the two
             instruments are comparable, not because it is the only choice).
    relic  = min(X,Y)/Omega_total
    """
    x = n[:, 0].astype(float)
    y = n[:, 1].astype(float)
    s = x + y
    order = np.where(s > 0, np.abs(x - y) / np.where(s > 0, s, 1.0), 1.0)
    total = n.sum(axis=1).astype(float)          # conserved, but read per trial
    relic = np.minimum(x, y) / total
    return [order, relic]


def n_winner_observables(n_committed: int, convention: str = "dominance"):
    """Order parameter for n-winner AM, on an internal-time grid.

    "dominance" is `experiments/expansion_radix.py`'s convention, kept identical
    so the curves are directly comparable with FINDINGS Sec.6:

        D = (max_i x_i - 1/n) / (1 - 1/n),   x_i = share of the committed pool,
        and D = 0 for the all-blank state (fully undecided).

    "gap" is `(top - runner_up) / sum(committed)`, which reduces exactly to AM's
    `|X-Y|/(X+Y)` at n_committed = 2. Second observable returned is the surviving
    committed-species count, the relic richness Sec.6 reports.
    """
    if convention not in ("dominance", "gap"):
        raise ValueError(f"unknown convention {convention!r}")
    n_c = n_committed

    def obs(n: np.ndarray) -> list:
        c = n[:, :n_c].astype(float)
        s = c.sum(axis=1)
        alive = s > 0
        safe = np.where(alive, s, 1.0)
        if convention == "dominance":
            top = c.max(axis=1) / safe
            order = np.where(alive, (top - 1.0 / n_c) / (1.0 - 1.0 / n_c), 0.0)
        else:
            srt = np.sort(c, axis=1)
            order = np.where(alive, (srt[:, -1] - srt[:, -2]) / safe, 1.0)
        return [order, (n[:, :n_c] > 0).sum(axis=1).astype(float)]
    return obs


def deterministic_times(omega: float, level: float = 0.5,
                        seed_scale: float = 1.0) -> dict:
    """Third route to the same law: no stochastics at all, only the ODE.

    Start the mass-action ODE at the SSA's own start, `(1/2, 1/2, 0)`, displaced
    along the decision axis by `delta0 = seed_scale * Omega^{-1/2}` -- the
    *quenched* stand-in for shot noise, which fixes the initial asymmetry at its
    typical size and then follows the deterministic flow. Returns the time to
    reach order parameter `level` and the time for the loser to fall to one
    molecule (`y = 1/Omega`, the absorption proxy).

    The two logarithms this exposes are the whole content of the law:
    `d(t_level)/d(ln Omega) -> 3/2` (symmetry breaking at rate lambda = 1/3 from
    an `Omega^{-1/2}` seed) and `d(t_clear)/d(ln Omega) -> 5/2` (that, plus
    clearing the last molecules off the rail at unit rate). Both approach their
    limits **from above**, which is why fits over a narrow Omega range read a
    slope well over 3/2 -- see FINDINGS.md Sec.5.1.

    `seed_scale` is the ONE thing this route cannot predict: the amplitude of the
    effective seed. It multiplies delta0 and therefore only shifts every time by
    the same `3 ln(seed_scale)`, leaving both slopes untouched. Calibrate it
    against one measured point if you want the offsets to line up; the slopes are
    the prediction.
    """
    from scipy.integrate import solve_ivp

    from .networks import approximate_majority

    net = approximate_majority()
    S = net.stoichiometry_matrix()
    if seed_scale <= 0.0:
        raise ValueError("seed_scale must be positive")
    delta0 = seed_scale * float(omega) ** -0.5
    if delta0 >= 1.0:
        raise ValueError(f"seed delta0={delta0:.3g} exceeds the simplex")

    def rhs(_t, x):
        return S @ net.fluxes(x)

    def hit_level(_t, x):
        s = x[0] + x[1]
        return (x[0] - x[1]) / s - level if s > 0 else 1.0 - level

    def hit_one(_t, x):
        return x[1] - 1.0 / float(omega)

    hit_level.terminal = False
    hit_one.terminal = True
    hit_one.direction = -1
    x0 = [0.5 * (1.0 + delta0), 0.5 * (1.0 - delta0), 0.0]
    sol = solve_ivp(rhs, (0.0, 400.0), x0,
                    rtol=1e-11, atol=1e-14, events=[hit_level, hit_one])
    if not sol.success:
        raise RuntimeError(f"ODE failed at Omega={omega}: {sol.message}")
    if not len(sol.t_events[0]) or not len(sol.t_events[1]):
        raise RuntimeError(f"events not reached at Omega={omega}, level={level}")
    return {"t_level": float(sol.t_events[0][0]),
            "t_clear": float(sol.t_events[1][0])}


def am_state_index(total: int):
    """Triangular indexing of the AM simplex {(x, y, b) : x+y+b = total}.

    Returns (xs, ys, index_of) where index_of(x, y) is a vectorised lookup. The
    ordering is x-major, matching `cme.enumerate_states` for the (X, Y, B)
    species order, so the two generators can be compared entry by entry.
    """
    n = total
    xs = np.concatenate([np.full(n - x + 1, x) for x in range(n + 1)])
    ys = np.concatenate([np.arange(n - x + 1) for x in range(n + 1)])
    offset = np.cumsum(np.concatenate([[0], np.arange(n + 1, 0, -1)]))[:n + 1]

    def index_of(x, y):
        return offset[x] + y

    return xs, ys, index_of


def am_generator(total: int):
    """Exact CME generator for irreversible AM at population `total`, vectorised.

    `cme.generator` is network-agnostic and loops in Python (2 s at N=400); this
    builds the same matrix from index arithmetic in milliseconds, which is what
    makes an exact D(tau) affordable at N in the thousands. Verified against
    `cme.generator` entry-for-entry in the tests.
    """
    import scipy.sparse as sp

    n = total
    xs, ys, index_of = am_state_index(n)
    bs = n - xs - ys
    i = np.arange(len(xs))
    c = 1.0 / n

    rows, cols, vals = [], [], []
    # r1: X + Y -> 2B
    m = (xs > 0) & (ys > 0)
    rows.append(i[m]); cols.append(index_of(xs[m] - 1, ys[m] - 1))
    vals.append(c * xs[m] * ys[m])
    # r2: B + X -> 2X
    m = (bs > 0) & (xs > 0)
    rows.append(i[m]); cols.append(index_of(xs[m] + 1, ys[m]))
    vals.append(c * bs[m] * xs[m])
    # r3: B + Y -> 2Y
    m = (bs > 0) & (ys > 0)
    rows.append(i[m]); cols.append(index_of(xs[m], ys[m] + 1))
    vals.append(c * bs[m] * ys[m])

    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    vals = np.concatenate(vals)
    m = len(xs)
    Q = sp.coo_matrix((vals, (rows, cols)), shape=(m, m)).tocsr()
    diag = np.asarray(Q.sum(axis=1)).ravel()
    return (Q - sp.diags(diag)).tocsr(), xs, ys


def am_order_exact(total: int, tau_max: float, n_grid: int = 401):
    """Exact D(tau) for AM from a symmetric start -- no sampling, no SSA.

    Integrates dp/dtau = Q^T p on the conserved simplex and returns the mean
    order parameter |X-Y|/(X+Y) on a uniform tau grid, using the same convention
    as `am_observables` (all-blank counts as 1.0). Because the expanding SSA is
    this chain read at tau = 1/H, this is an exact freeze-out curve.
    """
    from scipy.sparse.linalg import expm_multiply

    Q, xs, ys = am_generator(total)
    _, _, index_of = am_state_index(total)
    p0 = np.zeros(Q.shape[0])
    p0[index_of(total // 2, total - total // 2)] = 1.0
    s = (xs + ys).astype(float)
    order = np.where(s > 0, np.abs(xs - ys) / np.where(s > 0, s, 1.0), 1.0)
    P = expm_multiply(Q.T.tocsc(), p0, start=0.0, stop=tau_max, num=n_grid,
                      endpoint=True)
    taus = np.linspace(0.0, tau_max, n_grid)
    return taus, P @ order


def crossing_tau(taus, curve, level: float) -> float:
    """Internal time at which a rising curve crosses `level` (log-interpolated).

    Interpolates in log tau, matching freezeout_scaling.crossing's log
    interpolation in H, so the two are the same estimator of the same point.
    Returns nan if the curve does not cross.
    """
    taus = np.asarray(taus, dtype=float)
    curve = np.asarray(curve, dtype=float)
    for i in range(len(taus) - 1):
        lo, hi = curve[i], curve[i + 1]
        if lo <= level <= hi and hi > lo:
            f = (level - lo) / (hi - lo)
            return float(np.exp(np.log(taus[i])
                                + f * (np.log(taus[i + 1]) - np.log(taus[i]))))
    return float("nan")
