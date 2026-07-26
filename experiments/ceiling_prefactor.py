"""Is §12.1's depth-ceiling "factor ≈3" a genuine, σ-independent prefactor?

This is the second-axis test that ejected §14 from the (now dissolved) Q7 cluster,
applied to the one member that had never been run: §12.1 predicts a cascade depth
ceiling `D_max ≈ exp(δ*²/2σ²)/4` from a saddle point and measures about 3× that.
If the shortfall is a Laplace prefactor it must not depend on the channel noise σ.

§12.1's own data cannot answer this. Read at a *fixed* Ω the factor drifts
3.00 → 2.41 at Ω=64 but is flat 3.00 → 3.33 at Ω=128, because the death depth is
still converging in Ω and converges more slowly the smaller σ is (the ceiling is
higher, so the crossing is further away). So each σ is extrapolated to Ω→∞ the way
§12.1 did at σ/δ* = 0.45: the increments shrink geometrically, so the tail is
summed with the last observed ratio.

Prediction, written before running: the converged factor DRIFTS with σ. A saddle
point discards the Gaussian fluctuations, and the curvatures that a Kramers
prefactor is built from are σ-dependent here. A constant would be the surprise.

Result (§12.1, FINDINGS): it drifts 3.07 → 4.05, so it is not a prefactor — and
freeing the coefficient in the exponent fits at R² = 0.9998, which locates the whole
residual in the exponent (k = 1.0695 against the predicted 1) rather than in an
amplitude. Forcing k = 1 is what manufactures the apparent drifting "prefactor".

Runtime: ~10 minutes; the σ = 0.28 row dominates it.
"""
import time

import numpy as np

from crnl.information import depth_at_information
from crnl.networks.am_reversible import delta_star

GAMMA = 0.05
DS = delta_star(GAMMA)
OMEGAS = [16, 32, 64, 128, 192]
SIGMAS = [0.45, 0.38, 0.32, 0.28]   # in units of delta*


def converged_depth(depths):
    """Sum the geometric tail of a decelerating sequence, as §12.1 did.

    Returns (D_inf, ratio). Falls back to the last value when the increments
    are not decelerating, in which case no extrapolation is defensible.
    """
    incs = [depths[i] - depths[i - 1] for i in range(1, len(depths))]
    if len(incs) < 2 or incs[-1] <= 0:
        return depths[-1], None
    ratios = [incs[i] / incs[i - 1] for i in range(1, len(incs)) if incs[i - 1] > 0]
    r = ratios[-1] if ratios else 0.0
    if not (0.0 < r < 1.0):
        return depths[-1], r
    return depths[-1] + incs[-1] * r / (1.0 - r), r


def main():
    print(f"gamma={GAMMA}, delta*={DS:.4f}")
    header = (f"{'sig/d*':>7} " + " ".join(f"{'Om=' + str(o):>9}" for o in OMEGAS)
              + f" {'D_inf':>8} {'r':>6} {'pred':>9} {'factor':>7} {'s':>6}")
    print(header)
    rows = []
    for nf in SIGMAS:
        t0 = time.time()
        # generous cap: the ceiling itself grows like exp(1/2sig^2)
        max_depth = int(min(6000, 40 * np.exp(1.0 / (2 * nf * nf)) + 60))
        depths = [depth_at_information(GAMMA, om, 16.0, 0.5, max_depth, nf)
                  for om in OMEGAS]
        el = time.time() - t0
        if any(d is None for d in depths):
            shown = " ".join(f"{d:>9.2f}" if d is not None else f"{'>max':>9}"
                             for d in depths)
            print(f"{nf:>7.2f} {shown}   (not all resolved) {el:>6.0f}")
            continue
        dinf, r = converged_depth(depths)
        pred = np.exp(DS ** 2 / (2 * (nf * DS) ** 2)) / 4.0
        rows.append((nf, dinf, pred))
        print(f"{nf:>7.2f} " + " ".join(f"{d:>9.2f}" for d in depths)
              + f" {dinf:>8.2f} {r if r else 0:>6.2f} {pred:>9.2f}"
              f" {dinf / pred:>7.2f} {el:>6.0f}")

    if len(rows) < 3:
        return
    sig = np.array([r[0] for r in rows])
    dinf = np.array([r[1] for r in rows])
    facs = dinf / np.array([r[2] for r in rows])
    print("\nconverged factor across sigma: " + ", ".join(f"{f:.2f}" for f in facs))
    print(f"spread {max(facs) / min(facs):.2f}x  "
          "-- a genuine prefactor is flat (cf. §14's ejected 'constant', 1.7x)")

    # Locate the residual: free the coefficient in the exponent.
    x = 1.0 / (2 * sig ** 2)                 # = delta*^2 / 2 sigma^2
    k, c = np.polyfit(x, np.log(dinf), 1)
    resid = np.log(dinf) - (k * x + c)
    r2 = 1 - resid.var() / np.log(dinf).var()
    print(f"\nD_max = {np.exp(c):.3f} * exp({k:.4f} * d*^2/2sig^2)")
    print(f"  saddle point predicts k = 1 and prefactor 0.25")
    print(f"  R^2 = {r2:.6f}   residuals(ln) = "
          + ", ".join(f"{v:+.4f}" for v in resid))
    print(f"  forcing k=1 reproduces the 'prefactor' as "
          + ", ".join(f"{v:.2f}" for v in 4 * np.exp(c) * np.exp((k - 1) * x)))


if __name__ == "__main__":
    main()
