"""T-CASC-a: does the DEPTH CEILING transfer, when the thermodynamic costs did not?

Everything measured in §62-§70 is a single element, and the founding claim is about
COMPOSITION: a restoring gate lets you build deep, because error does not accumulate. §12.1
priced that for AM and found a hard ceiling set by the INTER-STAGE channel, not by molecule
count -- at large Omega the per-stage error saturates at an Omega-independent floor and

    D_max ~ exp(delta*^2 / 2 sigma^2) / 4

with delta* the element's rail half-separation and sigma the channel noise. Measured as the
depth at which the mutual information falls through 0.5, AM ran about 3x that prediction
(9 against 3.0, 44-50 against 14.8, 355-489 against 147).

**THE PRIOR HERE IS THE OPPOSITE OF §67's AND §68's, AND THAT IS THE POINT.** Those two asked
whether THERMODYNAMIC quantities transfer between substrates, and neither did: the affinity
floor is element-specific (§68: A_c(p) = 2 ln[(p+1)/(p-1)] runs 2.20 down to 0.50), and the
cost per e-fold has no counterpart on a chemostatted element at all (§67). But **the depth
ceiling is not a thermodynamic quantity.** It is information-theoretic -- von Neumann (1956)
and Evans-Schulman's noisy-circuit depth bounds -- and its formula mentions only the geometry
of the two rails against the channel noise. Nothing in it refers to closed bookkeeping, to
affinity, or to how the element is driven.

**So this one SHOULD transfer, and a failure would be far more surprising than §67's or
§68's.** Writing the prediction that way is the whole test: three prior transfer failures make
"does not transfer" the cheap answer, and it must not be the automatic one.

SUBSTRATE. §67/§69's Schloegl element: `pX <-> (p+1)X`, `0 <-> X`, chemostatted, one dynamic
species, no exchange symmetry, and 1-D so the stage kernel is an exact matrix exponential
rather than a 2-D CME solve. Rails at r1 and r3, so the analogue of delta* is the rail
half-separation Delta = (r3 - r1)/2, and the channel adds Gaussian noise of width
sigma = f * Delta in the SAME coordinate -- matching §12's convention exactly, which is what
makes the comparison a comparison (rule 11).

Stage order is CHANNEL THEN CHEMISTRY, matching §7 and `cascade_exact.run_cascade`; reading out
before the first channel step would report a trivial I = 1 at depth 0.

PREDICTIONS, written before running.

  P1  GATE. The stage kernel is a proper stochastic map (rows sum to 1 to 1e-12) and, run for
      a long time from a rail, leaves the rail occupied -- i.e. the element restores at all.
      If the kernel does not preserve the rails there is no restoration to compose and nothing
      below counts.
  P2  GATE. With the channel switched OFF (sigma = 0), the mutual information must NOT decay
      to 0.5 within the depths swept: a ceiling that appears without channel noise would be
      the kernel leaking, not composition. **This is the rule-10 guard** -- the harness must
      not manufacture the very decay being measured.
      **SECOND VERSION: the first tested `d_max(I) is None` and printed HOLDS on a run whose
      I was 0.000000.** d_max returns None both when I never falls through 0.5 and when it
      started below 0.5 and so never CROSSED it -- the two opposite cases the gate exists to
      separate. It is now `min(I) > 0.5` over the noiseless run, which only the intended case
      can satisfy. Caught by reading the output against P1, which had failed on the same run.
  P3  **THE TEST, absolute (rule 16).** Measure D_max as the depth at which I falls through
      0.5, at several sigma/Delta, and compare against **the same formula AM was measured
      against**, exp(Delta^2/2 sigma^2)/4. Report the RATIO measured/predicted, because AM's
      own ratio is about 3 and not 1 -- **the transferable claim is the ratio, not the raw
      depth.** If Schloegl's ratio matches AM's, the ceiling transfers and is the first
      quantity in this project that does.
  P4  **RULE 9, an axis I did not choose.** Sweep Omega and the autocatalytic order p. §12's
      claim is that D_max is Omega-INDEPENDENT at large Omega; if Schloegl's ceiling keeps
      rising with Omega, the ceiling is not the same object and P3's ratio is meaningless.
  P6  **P3 IS GATED ON P4, added after the first full run.** That run printed
      "THE DEPTH CEILING TRANSFERS" for f = 0.45 while P4 had already reported the ceiling
      still moving with Omega -- i.e. a verdict announced from a comparison its own
      precondition had just invalidated. A verdict must not print when the gate it depends on
      has failed; this is the same defect class as a gate on the wrong quantity (§64) and a
      gate applied in one place but not another (§66).
  P5  **VERDICT RULE unit-tested on engineered data before this runs** (§66's convention,
      which has caught a defect in four consecutive sections). Both "transfers" and "does not"
      must be reachable, and the comparison is against AM's PUBLISHED ratios, not a fresh AM
      run whose noise differs.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
from scipy.linalg import expm

from experiments.affinity_floor_family import degenerate_consts

# §12.1's published AM numbers: predicted D_max and the measured depth at Omega = 128
AM_PREDICTED = {0.45: 3.0, 0.35: 14.8, 0.28: 147.0}
AM_MEASURED = {0.45: 9.0, 0.35: 50.0, 0.28: 489.0}


def schlogl_consts(r1, r2, r3, k1r=1.0):
    e1, e2, e3 = r1 + r2 + r3, r1 * r2 + r1 * r3 + r2 * r3, r1 * r2 * r3
    return k1r * e1, k1r, k1r * e3, k1r * e2


def rates(omega, c, cap):
    k1a, k1r, k2b, k2r = c
    n = np.arange(0, cap + 1, dtype=float)
    lam = k1a * n * (n - 1.0) / omega + k2b * omega
    mu = (k1r * n * (n - 1.0) * (n - 2.0) / omega ** 2 + k2r * n)
    return np.maximum(lam, 0.0), np.maximum(mu, 0.0)


def stage_kernel(omega, r1, r2, r3, t_stage, cap_mult=2.0):
    """Exact per-stage map exp(Q t) for the 1-D chain. Rows are output distributions."""
    c = schlogl_consts(r1, r2, r3)
    cap = int(np.ceil(cap_mult * r3 * omega))
    lam, mu = rates(omega, c, cap)
    Q = np.zeros((cap + 1, cap + 1))
    idx = np.arange(cap + 1)
    Q[idx[:-1], idx[:-1] + 1] = lam[:-1]
    Q[idx[1:], idx[1:] - 1] = mu[1:]
    Q[idx, idx] = -Q.sum(axis=1)
    return expm(Q * t_stage), cap


def channel(cap, omega, sigma, r3):
    """Gaussian noise of width sigma (in concentration units) on the readout, then re-seed."""
    x = np.arange(cap + 1) / omega
    d = x[:, None] - x[None, :]
    C = np.exp(-0.5 * (d / sigma) ** 2)
    return C / C.sum(axis=1, keepdims=True)


def mutual_information(p_hi, p_lo, r2, omega):
    """I(input bit ; sign of the output about the saddle), inputs equiprobable."""
    n_sad = r2 * omega
    ns = np.arange(len(p_hi))
    up = ns > n_sad
    a = float(p_hi[up].sum())          # P(read high | started high)
    b = float(p_lo[up].sum())          # P(read high | started low)
    out = np.array([0.5 * (a + b), 1.0 - 0.5 * (a + b)])

    def h(v):
        v = v[v > 0]
        return float(-(v * np.log2(v)).sum())

    return h(out) - 0.5 * (h(np.array([a, 1 - a])) + h(np.array([b, 1 - b])))


def run(omega, r1, r2, r3, t_stage, depth, sigma, cap_mult=2.0):
    K, cap = stage_kernel(omega, r1, r2, r3, t_stage, cap_mult)
    C = (channel(cap, omega, sigma, r3) if sigma > 0
         else np.eye(cap + 1))
    step = C @ K                                   # CHANNEL THEN CHEMISTRY (§7's order)
    p_hi = np.zeros(cap + 1); p_hi[int(round(r3 * omega))] = 1.0
    p_lo = np.zeros(cap + 1); p_lo[int(round(r1 * omega))] = 1.0
    Is = []
    for _ in range(depth):
        p_hi = p_hi @ step
        p_lo = p_lo @ step
        Is.append(mutual_information(p_hi, p_lo, r2, omega))
    return np.array(Is), K


def d_max(Is, level=0.5):
    """First depth where I falls through `level`, linearly interpolated. None if never."""
    for i in range(len(Is) - 1):
        if Is[i] >= level > Is[i + 1]:
            t = (Is[i] - level) / (Is[i] - Is[i + 1])
            return float(i + 1 + t)
    return None


def verdict(ratio_s, ratio_am, tol=0.5):
    """Reachable both ways. Compares RATIOS, since AM's own ratio is ~3, not 1."""
    if ratio_s is None:
        return "none", "no ceiling reached within the depths swept"
    rel = abs(ratio_s - ratio_am) / abs(ratio_am)
    if rel <= tol:
        return "transfers", (f"Schloegl's measured/predicted = {ratio_s:.2f} against AM's "
                             f"{ratio_am:.2f} -- {100*rel:.0f}%, inside {100*tol:.0f}%. "
                             f"THE DEPTH CEILING TRANSFERS, unlike every thermodynamic "
                             f"quantity tried.")
    return "does-not", (f"Schloegl's measured/predicted = {ratio_s:.2f} against AM's "
                        f"{ratio_am:.2f} -- {100*rel:.0f}%, outside {100*tol:.0f}%. "
                        f"Does not transfer.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--r1", type=float, default=0.1)
    ap.add_argument("--r2", type=float, default=1.0)
    ap.add_argument("--r3", type=float, default=1.9)
    ap.add_argument("--omegas", type=int, nargs="+", default=[400, 600])
    ap.add_argument("--fracs", type=float, nargs="+", default=[0.45, 0.35, 0.28])
    ap.add_argument("--t", type=float, default=2.0)
    ap.add_argument("--depth", type=int, default=700)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/cascade_schlogl.json"))
    args = ap.parse_args()
    D = (args.r3 - args.r1) / 2.0
    print(f"Schloegl rails at {args.r1} and {args.r3}; Delta = {D}; saddle {args.r2}")

    print("\n=== P1 GATE: is the stage kernel stochastic, and does it hold a rail?")
    K, cap = stage_kernel(args.omegas[0], args.r1, args.r2, args.r3, args.t)
    rowsum = float(np.abs(K.sum(axis=1) - 1.0).max())
    p = np.zeros(cap + 1); p[int(round(args.r3 * args.omegas[0]))] = 1.0
    for _ in range(20):
        p = p @ K
    held = float(p[np.arange(cap + 1) > args.r2 * args.omegas[0]].sum())
    print(f"  |row sums - 1| max {rowsum:.2e}; mass still above the saddle after 20 "
          f"noiseless stages: {held:.6f}")
    print(f"  -> P1 {'HOLDS' if rowsum < 1e-12 and held > 0.99 else 'FAILS'}")

    print("\n=== P2 GATE (rule 10): with the channel OFF, does I survive?")
    Is0, _ = run(args.omegas[0], args.r1, args.r2, args.r3, args.t, args.depth, 0.0)
    print(f"  I after {args.depth} noiseless stages: {Is0[-1]:.6f}"
          f"  (min over the run {Is0.min():.6f})")
    p2 = bool(Is0.min() > 0.5)          # NOT `d_max is None`: see P2's second version
    print(f"  -> P2 {'HOLDS: no ceiling without channel noise' if p2 else f'FAILS: I falls to {Is0.min():.4f} with the channel OFF -- the kernel leaks and is manufacturing the decay'}")

    print(f"\n=== P3/P4: the ceiling against exp(Delta^2/2 sigma^2)/4")
    print(f"{'f=s/D':>7}{'Om':>6}{'sigma':>9}{'predicted':>11}{'measured':>10}{'ratio':>8}")
    rows = []
    for f in args.fracs:
        sigma = f * D
        pred = float(np.exp(D ** 2 / (2 * sigma ** 2)) / 4.0)
        for om in args.omegas:
            Is, _ = run(om, args.r1, args.r2, args.r3, args.t, args.depth, sigma)
            dm = d_max(Is)
            rows.append({"frac": f, "omega": om, "sigma": sigma, "pred": pred,
                         "dmax": dm, "ratio": (dm / pred) if dm else None,
                         "I_end": float(Is[-1])})
            print(f"{f:>7.2f}{om:>6}{sigma:>9.4f}{pred:>11.2f}"
                  + (f"{dm:>10.2f}{dm/pred:>8.2f}" if dm else f"{'>depth':>10}{'--':>8}"))

    print(f"\n=== P4 (rule 9): is the ceiling Omega-independent, as §12 claims for AM?")
    for f in args.fracs:
        sel = [r for r in rows if r["frac"] == f and r["dmax"]]
        if len(sel) >= 2:
            v = [r["dmax"] for r in sel]
            print(f"  f={f}: D_max = " + ", ".join(f"{x:.1f}" for x in v)
                  + f"  spread {100*(max(v)-min(v))/np.mean(v):.1f}%")
    print(f"  -> {'Omega-independent, so it is the same object §12 measured' if all(np.ptp([r['dmax'] for r in rows if r['frac']==f and r['dmax']] or [0]) / max(np.mean([r['dmax'] for r in rows if r['frac']==f and r['dmax']] or [1]), 1e-9) < 0.25 for f in args.fracs) else 'MOVES with Omega: not the same object, and P3s ratio is not comparable'}")

    sat = {}
    for f in args.fracs:
        v = [r["dmax"] for r in rows if r["frac"] == f and r["dmax"]]
        sat[f] = (len(v) >= 2 and np.ptp(v) / max(np.mean(v), 1e-9) < 0.10)

    print(f"\n=== P3 VERDICT: ratios, against AM's published ratios")
    print(f"  (P6: printed ONLY where P4's Omega-saturation gate passed)")
    for f in args.fracs:
        if not sat.get(f):
            v = [r["dmax"] for r in rows if r["frac"] == f and r["dmax"]]
            print(f"  f={f}: WITHHELD -- D_max still moving with Omega"
                  + (f" ({', '.join(f'{x:.1f}' for x in v)})" if v else "")
                  + ", so the ratio is not comparable to AM's saturated value")
            continue
        sel = [r for r in rows if r["frac"] == f and r["ratio"]]
        if not sel or f not in AM_MEASURED:
            continue
        rs = float(np.median([r["ratio"] for r in sel]))
        ram = AM_MEASURED[f] / AM_PREDICTED[f]
        code, msg = verdict(rs, ram)
        print(f"  f={f}: {msg}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"rows": rows, "Delta": D,
                                    "I_noiseless_end": float(Is0[-1])},
                                   indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
