"""§104 (T-CASC-t) -- predicting p_transmit, the last empirical number in the chain.

§103 closed the composition arc on single-element quantities except for ONE input:
p_transmit, taken from §100's measurement (0.9376 at t0 = 2.0). It is not a constant of the
coupling -- §100.2 measured it running 0.7254 -> 0.9830 as the window grows 16x, and it must
reach 1 as t -> infinity, because it is an ENDPOINT CO-OCCURRENCE and given enough time a
chain whose upstream fell is a chain that has failed. So a predicted version has to be a
function of the window, not a number.

THE STRUCTURE. Stage 1 falls at some time s in [0, t]; stage 2 then has only (t - s) left to
follow it down. With the upstream escape rate k1 tiny (3.8e-3 at Omega = 30, so k1*t << 1
over every window here) the fall time is near-uniform on [0, t], and

    p_transmit(t)  =  (1/t) INT_0^t [1 - exp(-k_low (t - s))] ds
                   =  1 - (1 - exp(-k_low t)) / (k_low t)

with ONE rate: k_low, the rate at which a stage whose input has collapsed crosses its own
saddle. That is a 1-D quantity -- the same pinned generator §102 used, evaluated at the LOW
rail instead of the high one.

WHAT WAS SCOUTED, and is therefore not a prediction (rule 2). I inverted that closed form on
each of §100.2's five measured windows separately BEFORE writing this file. It gives
k_low = 7.071, 8.771, 8.013, 7.463, 7.353 for t = 0.5, 1, 2, 4, 8 -- a spread of 1.24x across
a 16x change in window. So the one-parameter FORM is already known to work; what is not known,
and is the test below, is whether an INDEPENDENTLY COMPUTED rate lands in that range.

PREDICTIONS.

  P1  WIRING. With the upstream pinned at its low rail the downstream must be MONOSTABLE-low:
      no high rail, so F(r1) is undefined and the descent has no barrier to cross. If the
      pinned landscape still has three roots the instrument is not measuring what the name
      says.

  P2  THE TEST. The descent rate computed from the exact 1-D chain -- mean first-passage from
      the high rail down to the saddle, upstream pinned at r1 -- must land inside
      [7.07, 8.77], the interval the five windows independently imply. **The verdict is
      whether it lands in that interval, not whether it matches any single window**, because
      the windows themselves disagree by 1.24x and a tighter gate would be testing noise in
      the near-uniform-fall approximation rather than the rate (rule 20).

  P3  THE CURVE. Feeding that computed rate through the closed form must reproduce all five
      measured p_transmit values. I predict the worst residual is at the SHORTEST window,
      because that is where the near-uniform-fall approximation is weakest -- there stage 2's
      own independent escape and stage 1's chance of returning are least negligible relative
      to the transmitted signal. Residuals are reported per window; no aggregate.

  P4  OMEGA-INSENSITIVITY, which is what makes this a prediction rather than a refit. With the
      upstream pinned low the descent is a deterministic slide down a monostable landscape, so
      k_low is a MACROSCOPIC rate: it must be nearly independent of Omega, unlike every escape
      rate in §102 which varies by three orders. Checked at Omega = 14, 30, 55. If k_low
      instead scales with Omega it is not a descent rate and P2 passing would be a
      coincidence.

WHAT THIS DOES NOT CLAIM. p_transmit as measured is an endpoint co-occurrence, so this
predicts a quantity that is itself a proxy (§100.2 says so). Reproducing it makes §103's model
parameter-free on its own terms; it does not make the proxy a transmission probability.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

import experiments.chemical_cascade as cc
from experiments.depth_compounding import C, R1, R2, R3

# §100.2, Omega = 30 -- the target, not an input.
MEASURED = {0.5: 0.7254, 1.0: 0.8860, 2.0: 0.9376, 4.0: 0.9665, 8.0: 0.9830}
IMPLIED = (7.071, 8.771)          # the interval the five windows imply, scouted


def pinned_roots(x_up):
    """Positive roots of the downstream drift with the upstream pinned at x_up."""
    return cc.downstream_roots(x_up, C, R3, "hill")


def descent_rate(om, x_up=R1, cap_mult=1.25):
    """1/MFPT from the high rail down to the saddle, upstream pinned at x_up.

    Exact for the 1-D birth-death chain: absorbing at the saddle, reflecting at the cap.
    """
    cap = int(np.ceil(cap_mult * R3 * om))
    # "low" is counts < R2*om (last_low), so the first absorbing state is one BELOW
    # ceil(R2*om). The first version absorbed one site high; the fix lowers k_low.
    a = int(np.ceil(R2 * om)) - 1           # absorbing boundary: "stage is low"
    b = int(round(R3 * om))                 # start at the high rail
    idx = np.arange(a, cap + 1)
    n = len(idx)
    A = np.zeros((n, n))
    rhs = -np.ones(n)
    for i, s in enumerate(idx):
        if s == a:
            A[i, i] = 1.0
            rhs[i] = 0.0
            continue
        lam, mu = cc.rates_stage(float(s), x_up * om, om, C, R3, False, "hill")
        if s == cap:
            lam = 0.0
        A[i, i] = -(lam + mu)
        if s + 1 <= cap:
            A[i, i + 1] = lam
        if s - 1 >= a:
            A[i, i - 1] = mu
    T = np.linalg.solve(A, rhs)
    mfpt = float(T[list(idx).index(b)])
    return 1.0 / mfpt, mfpt


def p_transmit(k, t):
    """Near-uniform fall time on [0, t], then a descent at rate k."""
    return float(1.0 - (1.0 - np.exp(-k * t)) / (k * t))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/predicting_transmission.json"))
    args = ap.parse_args()

    print("P1 -- with the upstream pinned at its low rail, is the downstream monostable?")
    for x, name in ((R3, "high rail"), (R2, "saddle"), (R1, "low rail")):
        r = pinned_roots(x)
        print(f"  x_up = {x:.4f} ({name:>9}): {len(r)} positive root(s)"
              f"   {np.array2string(np.asarray(r), precision=4)}")
    mono = len(pinned_roots(R1)) < 3

    print(f"\nP2 -- the descent rate, computed; the windows imply [{IMPLIED[0]:.3f},"
          f" {IMPLIED[1]:.3f}]")
    k30, mfpt30 = descent_rate(30)
    print(f"  Omega = 30:  MFPT = {mfpt30:.5f}   k_low = {k30:.4f}"
          f"   {'INSIDE' if IMPLIED[0] <= k30 <= IMPLIED[1] else 'OUTSIDE'} the interval")

    print("\nP3 -- the curve, with nothing fitted")
    print(f"{'t0':>6}{'predicted':>12}{'measured':>11}{'residual':>11}")
    rows = []
    for t, m in sorted(MEASURED.items()):
        p = p_transmit(k30, t)
        rows.append({"t": t, "pred": p, "meas": m, "rel": (p - m) / m})
        print(f"{t:>6.1f}{p:>12.4f}{m:>11.4f}{(p-m)/m:>10.2%}")
    worst = max(rows, key=lambda r: abs(r["rel"]))
    print(f"  worst residual at t0 = {worst['t']}  ({worst['rel']:+.2%})"
          f"   -- P3 predicted the shortest window")

    print("\nP4 -- is k_low a macroscopic rate? (an escape rate would move 3 orders)")
    ks = {}
    for om in (14, 30, 55):
        k, _ = descent_rate(om)
        ks[om] = k
        print(f"  Omega = {om:3d}:  k_low = {k:.4f}")
    spread = max(ks.values()) / min(ks.values())
    print(f"  spread over Omega 14-55: {spread:.4f}x")

    args.out.write_text(json.dumps(
        {"monostable_at_low_rail": mono, "k_low": ks, "curve": rows,
         "implied_interval": IMPLIED}, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
