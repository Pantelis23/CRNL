"""Exact cascade arithmetic: what one restoring stage costs, and what it buys.

FINDINGS 7 measured WHY restoration matters -- a restoring cascade carries a bit to
depth 45 where analog passthrough is a coin flip by 22 -- but priced nothing, because
irreversible AM's dissipation is formally infinite. With the reversible network there
is a finite number per stage, and this module computes it exactly.

METHOD. A stage is: seed a fresh vessel from the previous stage's transmitted output
(B = 0, per design.md 4), let it run for a FIXED time t_stage, and read out the
composition it actually reached. No threshold, no sign(), no renormalisation -- the
stage emits the delta the chemistry produced. Then

    K[i, o]  = P(output count-diff o | seeded input i) after t_stage
    M[i]     = <net forward firings> over the stage from input i
    boundary = <lnW(out)> - lnW(in)
    cost[i]  = thermo.decompose(...)  =  boundary + (A/3) * M

and the cascade is a matrix product over depth. No trials, no sampling error.

THREE DESIGN POINTS, each of which a previous version got wrong.

1. THE CONTROL'S DYNAMIC RANGE IS AN AXIS, NOT A DETAIL.  A passive channel confined
   to +-1 while its noise is sigma_ch = 0.35*delta_star(gamma) has an ABSOLUTE
   dynamic range and LANDSCAPE-scaled noise. That mismatch grows with gamma and
   manufactures results: it made the control's depth-30 fidelity swing 0.538 -> 0.692
   across gamma while a control railed to the chemistry's own +-delta_star is flat to
   0.003 (0.5273/0.5253/0.5286/0.5299). An earlier headline -- "restoration needs a
   minimum Omega" -- was entirely this artifact: at gamma=0.30, Omega=30 the chemistry
   LOSES to the +-1 control (0.5640 vs 0.6002) and WINS against the matched-rail
   control (0.5640 vs 0.5220). So `control` takes an explicit `rail`, and callers are
   expected to report more than one. Never quote a single control.

2. <M> IS COMPUTED EXACTLY, NOT BY QUADRATURE.  Augmenting the generator with one row
   integrates the firing rate along the trajectory analytically:

       d/dt [p; m] = [[Q^T, 0], [w^T, 0]] [p; m],  m(0) = 0  =>  m(t) = int_0^t <w> ds

   Simpson's rule over expm_multiply time points was measured at 5.8% max relative
   error at t_stage=8, 240% at t=32 and 538% at t=150 -- always UNDERSTATING the cost,
   and worsening along exactly the axis this experiment sweeps. It also costs ngrid
   times the memory. There is no ngrid here.

3. TWO ALPHABETS, NOT ONE.  n_X - n_Y parity is NOT conserved (B+X->2X shifts it by 1)
   while a B=0 seed forces n_X - n_Y = Omega (mod 2). So seeded inputs are
   parity-restricted (Omega+1 of them, spanning -Omega..Omega in steps of 2) and
   outputs are not (2*Omega+1). About HALF the output mass lands on the opposite
   parity, so collapsing the two alphabets without renormalising loses half the mass
   per stage and drives fidelity to 0 like 0.5^depth -- not to 0.5, as an earlier note
   claimed. Row sums are asserted, not trusted.

Noise is in LANDSCAPE UNITS: sigma_ch = noise_frac * delta_star(gamma). delta_star(0)
= 1 EXACTLY, so FINDINGS 7 (absolute sigma against rails at +-1) is the gamma->0
member of this family.

The comparison to 7 is DIRECTIONAL, NOT QUANTITATIVE, and the depth-1 point shows why.
Both apply the channel before each readout, but 7 starts from a deliberately weak
s_init = 0.3 against sigma = 0.35 and so begins at Phi(0.3/0.35) = 0.80, whereas a
stage here starts from the previous stage's own rail (delta* = 0.952 at gamma=0.05,
sigma = 0.333) and begins at 0.9975. Both are defensible; they are different
experiments. `start_delta` exists so 7's convention can be reproduced when wanted.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.stats import norm

from .cme import enumerate_states, generator
from .networks.am_reversible import (
    GAMMA_C,
    am_reversible,
    cycle_affinity,
    delta_star,
    reverse_pairing,
)
from .thermo import decompose, ln_multinomial


def _check_bistable(gamma: float) -> None:
    if not 0.0 < gamma < GAMMA_C:
        raise ValueError(
            f"no bistable landscape at gamma={gamma} (need 0 < gamma < {GAMMA_C}); "
            "there is nothing to restore toward, and sigma_ch = noise_frac * "
            "delta_star(gamma) would be zero"
        )


def input_alphabet(omega: int) -> np.ndarray:
    """Count-differences reachable by a B=0 seed: Omega+1 values, step 2.

    NOT Omega/2 + 1 -- an earlier draft miscounted this, which put every memory
    estimate 2x low and made a test assert 16 where the answer is 31 at Omega=30.
    """
    return np.array([d for d in range(-omega, omega + 1) if (omega + d) % 2 == 0])


def stage_kernel(gamma: float, omega: int, t_stage: float,
                 chunk: int = 64) -> dict:
    """Exact per-stage map for reversible AM run for a fixed time t_stage.

    Returns din, dout, K (input x output), M (exact expected net forward
    firings), boundary, affinity. `chunk` bounds peak memory by propagating the
    input basis in column blocks; results are identical (measured max deviation
    2.8e-14).
    """
    _check_bistable(gamma)
    if not t_stage > 0.0:
        raise ValueError(f"t_stage must be positive, got {t_stage}")

    net = am_reversible(gamma)
    pairing = reverse_pairing(net)
    states, index = enumerate_states(3, omega)
    n_states = len(states)
    Q = generator(net, omega, float(omega))
    cs = net.stochastic_constants(float(omega))

    fwd = [j for j in range(net.n_reactions) if pairing[j] > j]
    rev = [j for j in range(net.n_reactions) if pairing[j] < j]
    w = np.array([net.propensities(n, cs)[fwd].sum()
                  - net.propensities(n, cs)[rev].sum() for n in states])
    lnW = np.array([ln_multinomial(n) for n in states])
    diff = np.array([int(n[0]) - int(n[1]) for n in states])

    din = input_alphabet(omega)
    dout = np.arange(-omega, omega + 1)
    oidx = {int(d): i for i, d in enumerate(dout)}

    # augmented generator: the extra row accumulates int_0^t <w> ds exactly
    Qa = sp.bmat([[Q.T.tocsr(), sp.csr_matrix((n_states, 1))],
                  [sp.csr_matrix(w.reshape(1, -1)), sp.csr_matrix((1, 1))]],
                 format="csc")

    K = np.zeros((len(din), len(dout)))
    M = np.empty(len(din))
    boundary = np.empty(len(din))
    lnW_in = np.array([ln_multinomial(
        np.array([(omega + int(d)) // 2, omega - (omega + int(d)) // 2, 0]))
        for d in din])

    for lo in range(0, len(din), chunk):
        block = din[lo:lo + chunk]
        cols = np.zeros((n_states + 1, len(block)))
        for c, d in enumerate(block):
            n_x = (omega + int(d)) // 2
            cols[index[(n_x, omega - n_x, 0)], c] = 1.0
        out = spla.expm_multiply(Qa, cols, start=0.0, stop=t_stage,
                                 num=2, endpoint=True)[-1]
        pend = out[:n_states, :]                       # (states, block)
        M[lo:lo + len(block)] = out[n_states, :]
        boundary[lo:lo + len(block)] = (lnW @ pend) - lnW_in[lo:lo + len(block)]
        for s in range(n_states):
            K[lo:lo + len(block), oidx[int(diff[s])]] += pend[s]

    rows = K.sum(axis=1)
    if not np.allclose(rows, 1.0, rtol=0.0, atol=1e-8):
        raise AssertionError(
            f"kernel rows do not sum to 1 (min {rows.min():.9f}, max "
            f"{rows.max():.9f}) -- output mass was lost, almost certainly the "
            "parity/alphabet error this module documents")
    return {"din": din, "dout": dout, "K": K, "M": M, "boundary": boundary,
            "affinity": cycle_affinity(net, pairing)}


def channel_matrix(src, dst, omega: int, sigma: float) -> np.ndarray:
    """C[s, d] = P(next seeded input dst[d] | current value src[s]).

    Gaussian on delta with sd sigma, requantised onto the parity-restricted
    input lattice by nearest-neighbour bin edges; the tails fall into the end
    bins, which is the right physics (a channel cannot push delta past the
    lattice).
    """
    gs, gd = np.asarray(src) / omega, np.asarray(dst) / omega
    edges = np.concatenate(([-np.inf], (gd[:-1] + gd[1:]) / 2.0, [np.inf]))
    C = np.empty((len(gs), len(gd)))
    for s, g0 in enumerate(gs):
        C[s] = np.diff(norm.cdf(edges, loc=g0, scale=sigma))
    return C


def _fidelity(p, grid) -> tuple[float, float]:
    """(half-credit, strict) probability the bit is still +1.

    Both are reported because the half-credit convention is cosmetically
    load-bearing near gamma_c: at gamma=0.45 the strict count is 0.4938 and
    half-credit rounds it to the tidier-looking 0.5005.
    """
    pos = float(p[grid > 0].sum())
    tie = float(p[grid == 0].sum())
    return pos + 0.5 * tie, pos


def run_cascade(gamma: float, omega: int, t_stage: float, depth: int,
                noise_frac: float = 0.35, chunk: int = 64,
                start_delta: float | None = None) -> dict:
    """P(bit correct) and cumulative dissipation at each depth, exactly.

    Stage order is CHANNEL THEN CHEMISTRY, matching cascade.py (FINDINGS 7),
    which adds noise before each restoring step. Reading out before the first
    channel step instead reports a trivial p_correct[0] = 1.0.

    `start_delta` defaults to the chemistry's own rail delta_star(gamma) -- the
    previous stage's clean output. Pass 0.3 to reproduce FINDINGS 7's weak-signal
    convention; the depth-1 value is Phi(start_delta/sigma_ch) either way (0.9975
    from the rail at gamma=0.05, 0.80 from 0.3 against sigma=0.35).
    """
    _check_bistable(gamma)
    d_star = delta_star(gamma)
    sigma = noise_frac * d_star
    k = stage_kernel(gamma, omega, t_stage, chunk)
    din, dout, K, M = k["din"], k["dout"], k["K"], k["M"]
    A = k["affinity"]

    # the /3 belongs to thermo.decompose and nowhere else, so build the
    # per-input stage cost through it rather than writing (A/3)*M inline
    stage_cost = np.array([
        decompose(None, None, float(M[i]), A, boundary=float(k["boundary"][i]))["total"]
        for i in range(len(din))])

    C_in = channel_matrix(dout, din, omega, sigma)     # output -> next input
    C_start = channel_matrix(din, din, omega, sigma)   # start   -> first input

    s0 = d_star if start_delta is None else float(start_delta)
    p = np.zeros(len(din))
    p[np.argmin(np.abs(din / omega - s0))] = 1.0

    go = dout / omega
    p_correct, p_strict, cum, total = [], [], [], 0.0
    first = True
    for _ in range(depth):
        p = p @ (C_start if first else C_in)           # channel corrupts
        first = False
        total += float(p @ stage_cost)
        q = p @ K                                      # chemistry restores
        half, strict = _fidelity(q, go)
        p_correct.append(half)
        p_strict.append(strict)
        cum.append(total)
        p = q
    return {"gamma": gamma, "omega": omega, "t_stage": t_stage,
            "delta_star": d_star, "sigma_ch": sigma, "affinity": A,
            "noise_frac": noise_frac, "start_delta": s0,
            "p_correct": p_correct, "p_correct_strict": p_strict,
            "cum_ds": cum, "ds_per_stage": cum[-1] / depth,
            "ds_final_stage": cum[-1] - (cum[-2] if depth > 1 else 0.0)}


def run_control(omega: int, depth: int, sigma: float, start_delta: float,
                rail: float = 1.0) -> dict:
    """Passive channel, no restoring stage: costs zero and loses the bit.

    `rail` is the control's dynamic range as a fraction of full scale, and it
    MUST be treated as an axis (see this module's docstring, point 1). rail=1.0
    gives the channel the whole lattice while the chemistry is confined to
    +-delta_star(gamma); rail=delta_star(gamma) matches them. The two disagree by
    up to 0.16 in depth-30 fidelity and they disagree MONOTONICALLY IN GAMMA, so
    picking one silently is how a control-side effect gets reported as chemistry.

    Deliberately NOT a free sign() limiter: a sign() is itself a restoring
    element, and handing the harness one for free assumes away the thesis.
    """
    din = input_alphabet(omega)
    keep = np.abs(din / omega) <= rail + 1e-12
    din = din[keep]
    if len(din) < 2:
        raise ValueError(f"rail={rail} leaves {len(din)} lattice sites at Omega={omega}")
    C = channel_matrix(din, din, omega, sigma)
    p = np.zeros(len(din))
    p[np.argmin(np.abs(din / omega - start_delta))] = 1.0
    g = din / omega
    half_l, strict_l = [], []
    for _ in range(depth):
        p = p @ C
        half, strict = _fidelity(p, g)
        half_l.append(half)
        strict_l.append(strict)
    return {"rail": rail, "sigma_ch": sigma, "p_correct": half_l,
            "p_correct_strict": strict_l, "cum_ds": [0.0] * depth,
            "ds_per_stage": 0.0}


#: Fidelity difference below which a chemistry-vs-control comparison is called a
#: tie. Both arms decay toward 0.5 with depth, so at depth 300 the raw
#: comparison returns verdicts on differences of 6e-07 -- noise dressed as a
#: result. 0.01 is ~3x the worst start-quantisation shift measured (0.003 per
#: lattice site at Omega=30).
TIE_BAND = 0.01


def verdict(chem: float, ctrl: float, band: float = TIE_BAND) -> str:
    if chem > ctrl + band:
        return "wins"
    if chem < ctrl - band:
        return "loses"
    return "tie"
