"""T-CASC-b: is §72's transfer about the ELEMENT, or was it never about the element?

§72 measured the saturated depth ceiling on Schloegl and found the ratio to
exp(Delta^2/2 sigma^2)/4 landing on AM's -- 2.71/3.04/3.75 against 3.00/3.38/3.33 -- and read
it as **the first quantity in this project to transfer across substrates.**

**There is a nastier reading, and §72 did not test it.** The formula contains no chemistry at
all: it is rail separation against channel noise. §72's own reduction makes the point sharply
-- the chemistry enters ONLY through `p_cross(n)`, the probability of committing to the wrong
rail from state n, and that function is monotone, near 1 at the low rail and near 0 at the
high one. **If any monotone commitment function with the same rail geometry gives the same
ratio, then the ceiling was never a property of the element and §72's headline is deflated**
from "transfers across substrates" to "was never about the substrate".

This is a direct attack on a result committed one section ago, which is why it is worth
running: §72 is exactly the kind of agreeable finding this project has had to withdraw before
(§67's affinity-floor pattern died in §68, one commit later).

THE TEST. Hold the rails, the channel and the depth criterion fixed, and swap ONLY the
commitment function:

  (a) **Schloegl exact** -- §72's baseline, the birth-death splitting probability.
  (b) **Langevin double well** -- U(x) = a(x - r1)^2 (x - r3)^2, splitting probability by the
      exact 1-D formula P(low first | x) = int_x^r3 e^{2U/D} / int_r1^r3 e^{2U/D}. Different
      physics, no chemistry, no counting noise, same geometry. Its "Omega" is 1/D.
  (c) **A bare sigmoid** -- p_cross(x) = 1/(1 + exp(k (x - r2))), carrying NO dynamics
      whatsoever, only the same monotone shape and the same midpoint. k is set so its slope at
      the saddle matches (a)'s, which is the only thing that could plausibly matter.
  (d) **A step function** -- the degenerate limit: commit to whichever rail is nearer. No
      element at all.

PREDICTIONS, written before running, and the honest one is uncomfortable.

  P1  GATE. (a) reproduces §72's numbers exactly (same code path, same Omega, same rails).
      If not, this file is not testing what §72 measured.
  P2  **THE TEST. PREDICTED: (c) and (d) give ratios in the same 2.7-3.8 band as (a).** The
      formula has no chemistry in it, so the ceiling is most likely set by the channel and the
      threshold, with the element contributing only the location of the midpoint. **That would
      DEFLATE §72** to "the ceiling is a property of the readout geometry", which is a weaker
      and more honest claim than "the first quantity that transfers".
  P3  **WHAT WOULD SAVE §72.** If (c) and (d) land clearly outside (a)'s band -- say the step
      function gives a ratio near 1, or above 6 -- then the commitment function's SHAPE
      matters, the element is doing real work, and §72's reading stands. The discriminator is
      stated before the numbers so it cannot be chosen afterwards.
  P4  **RULE 9.** Sweep sigma/Delta across the same three values §72 used, not one. A
      difference that appears at one channel width and not the others is not a difference.
  P5  If (b) -- a genuine but non-chemical dynamics -- agrees with (a) while (c)/(d) do not,
      that is the most informative outcome available: it would say the ceiling needs a real
      escape problem but not a chemical one, which is a substrate-independence claim with
      actual content.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from experiments.cascade_saturated import d_max_saturated, p_cross

AM_RATIO = {0.45: 9.0 / 3.0, 0.35: 50.0 / 14.8, 0.28: 489.0 / 147.0}


def grid(r1, r3, cap_mult, omega):
    cap = int(np.ceil(cap_mult * r3 * omega))
    return np.arange(cap + 1) / omega, cap


def pc_schlogl(x, r1, r2, r3, omega, cap_mult=1.6):
    pc, _, _ = p_cross(omega, r1, r2, r3, cap_mult)
    return pc


def pc_langevin(x, r1, r2, r3, omega, cap_mult=1.6, a=1.0):
    """Exact splitting probability for dx = -U'(x)dt + sqrt(2D)dW, U = a(x-r1)^2 (x-r3)^2.

    D is set so the barrier in units of D matches the chemical one's scale: D = 1/omega,
    which is the standard weak-noise correspondence (noise ~ 1/Omega).
    """
    # For dX = -U'(X)dt + sqrt(2D)dW the scale density is exp(U/D), NOT exp(2U/D):
    # the first version had that factor of 2 and, combined with a plain cumsum over an
    # integrand spanning e^9450, produced 0/0 and a silently clipped nonsense column.
    # Done here in logs with a suffix log-sum-exp, as §72's p_cross fix does.
    D = 1.0 / omega
    U = a * (x - r1) ** 2 * (x - r3) ** 2
    lo = int(np.argmin(np.abs(x - r1)))
    hi = int(np.argmin(np.abs(x - r3)))
    w = U / D
    seg = w[lo:hi + 1]
    m = seg.max()
    e = np.exp(seg - m)
    suff = np.cumsum(e[::-1])[::-1]
    out = np.ones_like(x)
    out[lo:hi + 1] = suff / suff[0]
    out[hi:] = 0.0
    out[:lo + 1] = 1.0
    if not np.all(np.isfinite(out)):
        return None
    return np.clip(out, 0.0, 1.0)


def pc_sigmoid(x, r1, r2, r3, omega, slope):
    return 1.0 / (1.0 + np.exp(slope * (x - r2)))


def pc_step(x, r1, r2, r3, omega):
    return (x < r2).astype(float)


def eps_from(pc, x, r1, r3, sigma):
    g_hi = np.exp(-0.5 * ((x - r3) / sigma) ** 2)
    g_lo = np.exp(-0.5 * ((x - r1) / sigma) ** 2)
    g_hi /= g_hi.sum()
    g_lo /= g_lo.sum()
    return float(g_hi @ pc), float(g_lo @ (1.0 - pc))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r1", type=float, default=0.1)
    ap.add_argument("--r2", type=float, default=1.0)
    ap.add_argument("--r3", type=float, default=1.9)
    ap.add_argument("--omega", type=int, default=14400)
    ap.add_argument("--fracs", type=float, nargs="+", default=[0.45, 0.35, 0.28])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/ceiling_is_it_the_element.json"))
    args = ap.parse_args()
    D = (args.r3 - args.r1) / 2.0
    x, cap = grid(args.r1, args.r3, 1.6, args.omega)

    pcs = {"a schlogl exact": pc_schlogl(x, args.r1, args.r2, args.r3, args.omega)}
    # match the sigmoid's slope at the saddle to Schloegl's, the only shape parameter that
    # could plausibly matter; everything else about it is arbitrary by construction.
    sad = int(np.argmin(np.abs(x - args.r2)))
    dsl = float(-(pcs["a schlogl exact"][sad + 1] - pcs["a schlogl exact"][sad - 1])
                / (x[sad + 1] - x[sad - 1]))
    lang = pc_langevin(x, args.r1, args.r2, args.r3, args.omega)
    if lang is not None:
        pcs["b langevin well"] = lang
    pcs["c sigmoid"] = pc_sigmoid(x, args.r1, args.r2, args.r3, args.omega, 4.0 * dsl)
    pcs["d step"] = pc_step(x, args.r1, args.r2, args.r3, args.omega)

    print(f"Omega = {args.omega}, rails {args.r1}/{args.r3}, Delta = {D}")
    print(f"Schloegl slope at the saddle: {dsl:.3f} per unit x  (sigmoid matched to it)")
    print(f"\n{'commitment function':>22}" + "".join(f"{f'f={f}':>14}" for f in args.fracs))
    rows = []
    for name, pc in pcs.items():
        cells = []
        for f in args.fracs:
            e_hi, e_lo = eps_from(pc, x, args.r1, args.r3, f * D)
            dm = d_max_saturated(e_hi, e_lo)
            pred = float(np.exp(1.0 / (2 * f ** 2)) / 4.0)
            ratio = dm / pred if dm else None
            cells.append(ratio)
            rows.append({"pc": name, "frac": f, "e_hi": e_hi, "e_lo": e_lo,
                         "dmax": dm, "ratio": ratio})
        print(f"{name:>22}" + "".join(f"{r:>14.2f}" if r else f"{'--':>14}" for r in cells))
    print(f"{'AM (published)':>22}"
          + "".join(f"{AM_RATIO[f]:>14.2f}" for f in args.fracs))

    print(f"\n=== P2/P3: does the ELEMENT matter, or only the readout geometry?")
    base = [r["ratio"] for r in rows if r["pc"] == "a schlogl exact" and r["ratio"]]
    band = (min(base), max(base))
    print(f"  Schloegl's band across channels: {band[0]:.2f}..{band[1]:.2f}")
    for name in ("b langevin well", "c sigmoid", "d step"):
        v = [r["ratio"] for r in rows if r["pc"] == name and r["ratio"]]
        if not v:
            continue
        inside = all(band[0] * 0.75 <= r <= band[1] * 1.33 for r in v)
        print(f"  {name:>18}: {min(v):.2f}..{max(v):.2f}"
              f"  -> {'INSIDE Schloegl band' if inside else 'OUTSIDE'}")
    triv = [r["ratio"] for r in rows if r["pc"] == "d step" and r["ratio"]]
    if triv:
        inside = all(band[0] * 0.75 <= r <= band[1] * 1.33 for r in triv)
        print(f"\n  -> {'DEFLATES §72: a step function with NO dynamics reproduces the ceiling, so it is a property of the readout geometry, not of the element' if inside else 'STANDS: the commitment function shape matters, so the element is doing real work and §72s reading survives'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"rows": rows, "slope": dsl, "Delta": D},
                                   indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
