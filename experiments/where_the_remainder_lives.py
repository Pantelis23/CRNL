"""§117 (T-CASC-aj) -- where the closure's remainder lives, factor by factor.

After §116 the closure of §103/§107 -- with stage 2's rate replaced by §115's exact Q-process gap
and the measured side re-run with matched seeding -- sits FLAT at a 21% offset over Omega = 14-40
(model/measured contam/pure 1.203 ... 1.226) and RISES beyond (1.300 at 55, 1.442 at 70).

The closure is a product of three factors, and each one has a measured counterpart in the same
matched-seed joint solve, so the error can be split with nothing substituted:

    contam/pure  =  [ los0 * p_t ]  /  [ (1 - los0) * los1 ]

    los0  = P(stage 1 low at t)                 model: pi_low(r3) * (1 - exp(-k0 t))
    p_t   = P(stage 2 low | stage 1 low)        model: 0.9376, ONE value from Omega = 30
    los1  = P(stage 2 low | stage 1 high)       model: pi * (1 - exp(-k2 t))

T-CASC-aj's suspect was p_t: one constant calibrated mid-range. A second suspect is visible in
§116's own numbers, and it is an inconsistency of my own making: §116 paired the EXACT conditioned
rate k2 with the PINNED occupancy pi. Both belong to the same Q-process; the consistent two-state
pair is (pi_cond, lambda1 - lambda0). pi_cond/pi_pinned runs 1.03 (Omega = 14) -> 1.38 (70), which
grows with Omega exactly where the remainder rises.

PREDICTIONS. NOT BLIND on the second point, and said so: I read pi_cond/pi_pinned off §116's JSON
before writing this, and divided §116.5's matched ratios by it in my head (giving ~1.17 -> 1.04).

  P1  WIRING. The three model/measured factor ratios multiply back to §116.5's matched-seed column,
      to 1e-6 (the eigensolver floor, T-NUM-a -- not tighter; §114.8(d), §116.4).
  P2  THE RISE. It lives in los1, and it is the pinned occupancy: los1's model/measured ratio with
      pinned pi rises with Omega and with pi_cond it does not rise. If the rise instead lives in
      p_t, T-CASC-aj's original suspect was right and the pinned pi is innocent.
  P3  THE FLAT OFFSET, a genuine guess: it lives in p_t, because 0.9376 was measured under §101's
      DEFAULT seed and everything else here is matched. If it lives in los0 instead, the stage-1
      two-state form is what carries it.
  P4  THE CONSISTENT CLOSURE. With (pi_cond, exact gap) for stage 2 and p_t, los0 left as modelled,
      report the ratio across Omega and its span. No threshold: read the cells (rule 18, and §116
      printed two wrong verdict clauses off a span criterion).

WHAT THIS CANNOT SETTLE. The exact gap and pi_cond come from a two-stage stationary eigenproblem, so
a closure built on them is no longer §103's "no joint solve anywhere" -- it trades a time-evolved
joint solve for a stationary joint eigenproblem. That is reported as a cost, not hidden.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib

import numpy as np

from experiments.depth_compounding import R2, R3
from experiments.escape_accounts_for_it import escape_rate
from experiments.free_upstream_depth import solve
from experiments.the_corrected_closure import P_TRANSMIT_MEASURED, pi_low

T = 2.0


def measured_factors(om, t=T):
    """los0, p_t, los1 from ONE matched-seed free-chain joint solve."""
    p, ref, dims, strides, walled = solve(om, 2, 0, t, matched_seed=True)
    m = dims[0]
    g = p.reshape(m, dims[1])
    s1lo = np.arange(m) < R2 * om
    s2lo = np.arange(dims[1]) < R2 * om
    los0 = float(g[s1lo, :].sum())
    both = float(g[np.ix_(s1lo, s2lo)].sum())
    pure = float(g[np.ix_(~s1lo, s2lo)].sum())
    return {"los0": los0, "p_t": both / los0, "los1": pure / (1 - los0),
            "contam": both, "pure": pure}


def p5_window(rows):
    """KILL TEST for the account of what remains (rule 17), criterion fixed before running.
    Suspect: the remaining 3-10% is the two-state form's finite-window transient on each stage.
    If so, at t = 8 both deficits shrink toward 1. If stage 1's ~5% stays put, it is not a
    window effect -- the next suspect is the definition itself: "below the saddle at t" counts
    saddle-crossers that return, which a two-state rate does not."""
    print("\nP5 -- kill test: the same factors at t = 8")
    out = []
    for r in rows:
        om = r["omega"]
        f = measured_factors(om, 8.0)
        k0 = escape_rate(om, R3)
        a = pi_low(om, R3) * (1 - math.exp(-k0 * 8.0)) / f["los0"]
        c = r["pi_cond"] * (1 - math.exp(-r["k_exact"] * 8.0)) / f["los1"]
        out.append({"omega": om, "r_los0_t8": a, "r_los1_cond_t8": c,
                    "r_los0_t2": r["r_los0"], "r_los1_cond_t2": r["r_los1_cond"]})
        print(f"   Om={om:>3}  los0: t=2 {r['r_los0']:.4f} -> t=8 {a:.4f}"
              f"    los1 (cond): t=2 {r['r_los1_cond']:.4f} -> t=8 {c:.4f}", flush=True)
    return out


def p6_long_window(rows):
    """Written before running. P5 refuted the account of stage 2's residual: its conditioned
    factor goes 0.990 -> 0.849 at Omega = 70 between t = 2 and t = 8, so the t = 2 agreement was a
    time-dependent ratio passing through 1. For a two-state Q-process pi_cond * lambda IS the
    forward flux, so at long t model/measured should return to 1 -- IF the drift is a seed
    transient (stage 2 starts in the pinned-upstream QSD, not the conditioned one).

      P6a  it is a transient: the ratio turns and converges back toward 1 by t = 32-64.
      P6b  if instead it converges to a constant clearly below 1 (e.g. ~0.85), the identity
           pi_cond * lambda = measured flux is broken, and the free chain is not measuring what
           the Q-process describes -- the next suspect being that the FREE chain conditions on
           stage 1 high at the END, not on stage 1 never having failed.
    """
    print("\nP6 -- stage 2's factor over a long window (Omega = 30, 70), one chained evolution")
    import scipy.sparse.linalg as spla
    from experiments.free_upstream_depth import build_gen, seed_matched
    out = []
    byom = {r["omega"]: r for r in rows}
    for om in (30, 70):
        # Omega = 70's generator is ~5x larger and each 2-unit step costs minutes; it is taken
        # to t = 16 only, which the t = 8 column of P5 already shows to be near its plateau.
        marks = (2.0, 4.0, 8.0, 16.0, 32.0, 64.0) if om == 30 else (2.0, 4.0, 8.0, 16.0)
        r = byom[om]
        Q, ref, dims, strides, cap, walled = build_gen(om, 2, 0)
        p = seed_matched(om, ref, dims, strides, walled)
        QT = (Q.T * 2.0).tocsr()
        seq, t = [], 0.0
        s1lo = np.arange(dims[0]) < R2 * om
        s2lo = np.arange(dims[1]) < R2 * om
        while t < marks[-1] - 1e-9:
            p = spla.expm_multiply(QT, p)          # advance by exactly 2
            t += 2.0
            if any(abs(t - mk) < 1e-9 for mk in marks):
                g = p.reshape(dims[0], dims[1])
                los0 = float(g[s1lo, :].sum())
                los1 = float(g[np.ix_(~s1lo, s2lo)].sum()) / (1 - los0)
                model = r["pi_cond"] * (1 - math.exp(-r["k_exact"] * t))
                seq.append({"t": t, "ratio": model / los1, "los1": los1})
                print(f"      Om={om} t={t:g}: model/measured = {model / los1:.4f}", flush=True)
        out.append({"omega": om, "seq": seq})
        print(f"   Om={om:>3}: " + "  ".join(f"t={x['t']:g}: {x['ratio']:.4f}" for x in seq),
              flush=True)
    return out


def p7_absorbing(rows, om=30):
    """Written before running. P6 at Omega = 30 plateaus at model/measured ~0.836 (t = 16-32):
    the FREE chain's stage 2 falls ~16% faster than the Q-process's (pi_cond, lambda) predicts.
    Suspect (P6b): the free chain counts stage-1 NEAR-MISSES -- dips below the saddle that recover
    -- whose effect on stage 2 survives, while the Q-process kills those trajectories.

    Test: the ABSORBING chain is exactly what the Q-process describes (stage 1 killed at its
    saddle). Evolve it from the same seeds and track P(stage 2 low | stage 1 alive).
      P7a  it converges to ratio 1 (up to a seed transient): the Q-process pair is right, and the
           free chain's 16% is a genuine near-miss channel the Q-process excludes by construction.
      P7b  it ALSO plateaus near 0.84: the (pi_cond, lambda) pair is wrong, and the near-miss
           account is refuted.
    """
    import scipy.sparse.linalg as spla
    from experiments.the_survival_term import _joint_alive
    from experiments.two_axes_for_the_position import Elem, PUBLISHED_RAILS
    r = {x["omega"]: x for x in rows}[om]
    e = Elem("published", *PUBLISHED_RAILS, cap_mult=1.25)
    Qa = _joint_alive(e, om)
    m = e.cap(om) + 1
    alive = np.arange(m)[np.arange(m) >= e.r2 * om]
    up, pi1 = e.stage1_stationary(om)
    mus, _ = e.operating_points(om, 2)
    keep, w2 = e.qsd(om, mus[0])
    g = np.zeros((len(alive), m))
    for a, wa in zip(up, pi1):
        g[list(alive).index(a), keep] = wa * w2
    p = g.ravel()
    QT = (Qa.T * 2.0).tocsc()
    s2lo = np.tile(np.arange(m) < e.r2 * om, len(alive))
    print(f"\nP7 -- the ABSORBING chain (what the Q-process describes), Omega = {om}")
    seq, t = [], 0.0
    for _ in range(16):
        p = spla.expm_multiply(QT, p)
        t += 2.0
        if t in (2.0, 4.0, 8.0, 16.0, 32.0):
            cond = float(p[s2lo].sum() / p.sum())
            model = r["pi_cond"] * (1 - math.exp(-r["k_exact"] * t))
            seq.append({"t": t, "ratio": model / cond})
            print(f"      t={t:g}: model/measured = {model / cond:.4f}", flush=True)
    return {"omega": om, "seq": seq}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/where_the_remainder_lives.json"))
    args = ap.parse_args()
    prev = json.loads(pathlib.Path("results/the_residual_had_a_sign.json").read_text())
    rows116 = {r["omega"]: r for r in prev["rows"]}
    seed116 = {r["omega"]: r["ratio_exact_matched"] for r in prev["seed_kill"]["rows"]}

    print(f"{'Om':>4}{'los0 m/meas':>12}{'p_t meas':>10}{'p_t m/meas':>11}"
          f"{'los1 pin':>10}{'los1 cond':>11}{'product':>9}{'§116.5':>8}{'consistent':>12}")
    rows = []
    for om in sorted(rows116):
        r = rows116[om]
        f = measured_factors(om)
        k0 = escape_rate(om, R3)
        los0_m = pi_low(om, R3) * (1 - math.exp(-k0 * T))
        los1_pin = r["pi_pinned"] * (1 - math.exp(-r["k_exact"] * T))
        los1_cond = r["pi_cond"] * (1 - math.exp(-r["k_exact"] * T))
        a = los0_m / f["los0"]
        b = P_TRANSMIT_MEASURED / f["p_t"]
        c_pin = los1_pin / f["los1"]
        c_cond = los1_cond / f["los1"]
        d = (1 - f["los0"]) / (1 - los0_m)               # the (1 - los0) factor, inverted
        prod = a * b * d / c_pin
        consistent = a * b * d / c_cond
        rows.append({"omega": om, **f, "los0_model": los0_m, "pi_cond": r["pi_cond"],
                     "k_exact": r["k_exact"], "r_los0": a, "r_pt": b,
                     "r_los1_pinned": c_pin, "r_los1_cond": c_cond, "r_one_minus_los0": d,
                     "product": prod, "seed116": seed116[om], "consistent": consistent})
        print(f"{om:>4}{a:>12.4f}{f['p_t']:>10.4f}{b:>11.4f}{c_pin:>10.4f}{c_cond:>11.4f}"
              f"{prod:>9.4f}{seed116[om]:>8.4f}{consistent:>12.4f}", flush=True)

    w = max(abs(r["product"] / r["seed116"] - 1) for r in rows)
    print(f"\nP1 -- factors multiply back to §116.5: worst rel {w:.1e}"
          f" -> {'HOLDS' if w < 1e-5 else 'FAILS'}")
    for key, lab in (("r_los0", "los0"), ("r_pt", "p_t"), ("r_los1_pinned", "los1 (pinned pi)"),
                     ("r_los1_cond", "los1 (conditioned pi)")):
        v = [r[key] for r in rows]
        print(f"   {lab:>22} model/measured: " + " ".join(f"{x:.3f}" for x in v)
              + f"   span {max(v)/min(v):.3f}x")
    cs = [r["consistent"] for r in rows]
    print(f"P4 -- consistent closure: " + " ".join(f"{x:.3f}" for x in cs)
          + f"   span {max(cs)/min(cs):.3f}x")
    p5 = p5_window(rows)
    p6 = p6_long_window(rows)
    p7 = p7_absorbing(rows)
    args.out.write_text(json.dumps({"rows": rows, "p1_worst": w, "p5": p5, "p6": p6, "p7": p7},
                                   indent=2, default=float))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
