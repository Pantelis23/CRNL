"""What a cascade costs per bit it actually delivers.

Every comparison in FINDINGS 10 needed a CONTROL, and the control's dynamic range
turned out to be a free parameter that had already manufactured one withdrawn
claim. This module removes the comparator entirely.

Send b = +-1 equiprobably as delta_0 = +-delta_star(gamma). The cascade is
symmetric, so after D stages the two output distributions p_+ and p_- are mirror
images, and

    I(b ; X_D) = H( (p_+ + p_-)/2 ) - H(p_+)          [bits]

is exactly what survived. Divide the cumulative dissipation by it and the result
is k_B T per bit delivered -- no control, no rail convention, no tie band, and
directly comparable to k_B T ln 2.

THE DEPTH IS NOT OPTIONAL, and this is the third time this project has met the
same trap. Asked at DEPTH 1 the question is degenerate: a stage with t_stage ->
0 does nothing, costs nothing, and still scores well because a single channel
application barely damages a bit sitting on a rail. Measured at gamma=0.15,
Omega=30: 0.89 kT/bit at t_stage=0.05 versus 20.2 at t_stage=16 -- the cheapest
"restoration event" is the one that does not restore. The two earlier versions
of this failure were a stop predicate that fired on the initial state (83-96% of
stages ran zero reactions) and a control whose rails were free to differ from the
chemistry's.

The degeneracy disappears once the bit must survive a depth at which a passive
channel would lose it. There, doing nothing scores infinitely badly: at depth 30
the same t_stage=0.05 costs 5493 kT/bit and t_stage=1 costs 155984, while an
interior optimum in t_stage appears. So the well-posed question is

    "what is the cheapest way to deliver one bit to depth D?"

and D is part of the question, not a nuisance parameter. `cost_per_bit` therefore
requires a depth and returns the whole depth profile, so the D-dependence is
always visible rather than collapsed to a single headline number.
"""

from __future__ import annotations

import numpy as np

from .cascade_exact import channel_matrix, stage_kernel
from .networks.am_reversible import delta_star
from .thermo import decompose


def shannon_bits(p) -> float:
    """H(p) in bits, ignoring zero-probability outcomes."""
    p = np.asarray(p, dtype=float)
    p = p[p > 0.0]
    return float(-(p * np.log2(p)).sum())


def mutual_information_bits(p_plus, p_minus) -> float:
    """I(b ; X) for an equiprobable binary input, in bits.

    I = H(mixture) - <H(conditional)>. Written without assuming symmetry so it
    stays correct if an asymmetric network is ever pushed through it (the
    symmetric case is the one measured here, where H(p_+) == H(p_-)).
    """
    p_plus = np.asarray(p_plus, dtype=float)
    p_minus = np.asarray(p_minus, dtype=float)
    mixture = 0.5 * (p_plus + p_minus)
    info = shannon_bits(mixture) - 0.5 * (shannon_bits(p_plus)
                                          + shannon_bits(p_minus))
    # I >= 0 always. When the bit is completely destroyed the two conditionals
    # coincide and the difference of two ~equal entropies lands at float noise:
    # measured -8.88e-16 (4*eps) at gamma=0.45, t_stage=64. Clamping keeps the
    # sign honest, and the guard distinguishes that noise from a real defect --
    # without it two cells of identical physics reported `inf` while their
    # neighbours reported a large finite number.
    if info < -1e-9:
        raise ValueError(
            f"mutual information came out {info}, which is not float noise; "
            "the conditional distributions are probably misaligned")
    return max(info, 0.0)


def cost_per_bit(gamma: float, omega: int, t_stage: float, depth: int,
                 noise_frac: float = 0.35, chunk: int = 64) -> list[dict]:
    """Exact I(b;X_d) and cumulative dissipation at every depth d = 1..depth.

    Returns one row per depth with `I_bits`, `ds` (cumulative k_B T) and
    `kT_per_bit`. Read the profile, not just the last row: cost per bit grows
    superlinearly with depth because information decays exponentially while cost
    accumulates linearly.
    """
    if depth < 1:
        raise ValueError(f"depth must be at least 1, got {depth}")
    d_star = delta_star(gamma)
    sigma = noise_frac * d_star
    k = stage_kernel(gamma, omega, t_stage, chunk)
    din, dout, K, M, A = k["din"], k["dout"], k["K"], k["M"], k["affinity"]

    stage_cost = np.array([
        decompose(None, None, float(M[i]), A,
                  boundary=float(k["boundary"][i]))["total"]
        for i in range(len(din))])
    C_in = channel_matrix(dout, din, omega, sigma)
    C_start = channel_matrix(din, din, omega, sigma)

    def propagate(sign: int):
        p = np.zeros(len(din))
        p[np.argmin(np.abs(din / omega - sign * d_star))] = 1.0
        total = 0.0
        out = []
        first = True
        for _ in range(depth):
            p = p @ (C_start if first else C_in)
            first = False
            total += float(p @ stage_cost)
            q = p @ K
            out.append((q.copy(), total))
            p = q
        return out

    plus, minus = propagate(+1), propagate(-1)
    rows = []
    for d in range(depth):
        p_p, total = plus[d]
        p_m, _ = minus[d]
        info = mutual_information_bits(p_p, p_m)
        rows.append({
            "depth": d + 1, "I_bits": info, "ds": total,
            "kT_per_bit": (total / info) if info > 1e-12 else float("inf"),
        })
    return rows
