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

    # The start quantises onto the input lattice, and at small Omega that shift
    # is large and SYSTEMATIC: at gamma=0.05 (delta*=0.9521) every Omega <= 20
    # snaps to the full rail 1.0 (+5.0%) while Omega=30 lands at 0.9333 (-2.0%).
    # Since a further-out start is an easier start, this flatters small Omega --
    # precisely the regime where an efficiency optimum is being looked for. It is
    # reported per row rather than buried.
    start_idx = int(np.argmin(np.abs(din / omega - d_star)))
    realised_start = float(din[start_idx] / omega)

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
            "realised_start": realised_start, "delta_star": d_star,
            "start_shift_frac": (realised_start - d_star) / d_star,
        })
    return rows


def flip_probability(gamma: float, omega: int, t_stage: float,
                     noise_frac: float = 0.35, depth: int = 40,
                     chunk: int = 64) -> float:
    """Per-stage probability that the cascade loses the bit.

    Extracted from the exponential decay of I(b;X_D) rather than measured
    directly: for a binary symmetric chain with per-stage flip probability p,
    `I ~ (1-2p)^{2D}`, so the slope of ln I against depth gives p without ever
    needing a threshold or a control. Returns NaN if the information falls below
    numerical resolution too quickly to fit.
    """
    prof = cost_per_bit(gamma, omega, t_stage, depth, noise_frac, chunk)
    info = np.array([r["I_bits"] for r in prof])
    d = np.arange(1, len(info) + 1)
    usable = info > 1e-9
    if usable.sum() < 4:
        return float("nan")
    slope = np.polyfit(d[usable], np.log(info[usable]), 1)[0]
    return float(0.5 * (1.0 - np.exp(slope / 2.0)))


def wall_coefficient(gamma: float) -> float:
    """kappa(gamma): the restoration-wall coefficient, P(err) ~ exp(-kappa*e^2*Omega).

    design.md 9 derives kappa = 3/2 for IRREVERSIBLE AM, from the saddle's
    restoring gain lambda = 1/3 against finite-count diffusion D = 1/(9 Omega).
    At finite gamma the gain is lambda_antisym(gamma) = (1-2gamma)/3, so the
    barrier scales with it:

        kappa(gamma) = (3/2) * (1 - 2 gamma) = (9/2) * lambda_antisym(gamma)

    Using the gamma=0 value everywhere fails badly as gamma -> gamma_c: the
    collapse of `predicted_exponent` degrades from R^2 = 0.99 at gamma=0.05 to
    0.60 at gamma=0.45, and pooling all gamma gives 0.69. With this correction
    the pooled collapse is 0.93.
    """
    return 1.5 * (1.0 - 2.0 * gamma)


def predicted_exponent(gamma: float, omega: int,
                       noise_frac: float = 0.35) -> float:
    """Saddle-point prediction for -ln p, with NO fitted parameter.

    A flip needs the channel to displace the state from the rail +delta* to some
    delta, and then finite-count noise to carry it the rest of the way. The two
    costs add in the exponent:

        f(delta) = (delta* - delta)^2 / (2 sigma^2)  +  kappa * Omega * delta^2
                    channel (Gaussian)                  restoration wall

    Minimising over delta gives delta = delta*/(1 + 2 kappa Omega sigma^2) and

        -ln p  ~  kappa Omega delta*^2 / (1 + 2 kappa Omega sigma^2)

    which interpolates BOTH regimes with one expression:

      * 2 kappa Omega sigma^2 << 1 -> kappa Omega delta*^2, the restoration wall
        of FINDINGS 1-2: exponential in Omega, the population-limited side;
      * 2 kappa Omega sigma^2 >> 1 -> delta*^2 / (2 sigma^2), independent of
        Omega: the channel-limited floor, where more molecules buy nothing.

    The crossover sits at Omega_x = 1 / (2 kappa sigma^2).
    """
    d_star = delta_star(gamma)
    sigma = noise_frac * d_star
    kappa = wall_coefficient(gamma)
    return float(kappa * omega * d_star ** 2
                 / (1.0 + 2.0 * kappa * omega * sigma ** 2))


def crossover_omega(gamma: float, noise_frac: float = 0.35) -> float:
    """Population at which the restoration wall gives way to the channel floor."""
    sigma = noise_frac * delta_star(gamma)
    return float(1.0 / (2.0 * wall_coefficient(gamma) * sigma ** 2))
