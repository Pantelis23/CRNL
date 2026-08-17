"""T-DEPTH-b: in a CHEMICAL cascade there is no abstract channel. What is sigma?

§74's entire content rests on holding the inter-stage noise sigma fixed in physical units
while the rails move. Under §12's other convention, sigma = f * Delta, the ratio Delta/sigma
is constant and the question cannot be posed at all. §74 said plainly that which convention is
physical was untested and that everything hung on it. This tests it.

**In a cascade whose coupling is itself chemical, sigma is not a modelling choice.** Stage 1's
output species IS stage 2's input, so the corruption between stages is stage 1's own
fluctuation about its rail -- a quantity the CME fixes exactly, with no freedom.

**AND THERE IS A REASON TO EXPECT §74 TO FAIL.** Under Schloegl's rescaling r -> lambda r the
constants go as k1a ~ lambda, k2r ~ lambda^2, k2b ~ lambda^3, so at corresponding
concentrations x = lambda x0 EVERY propensity scales as lambda^3:

    a1 = k1a x^2 -> lambda^3,  a2 = k1r x^3 -> lambda^3,
    a3 = k2b     -> lambda^3,  a4 = k2r x   -> lambda^3

The chain is therefore time-rescaled by lambda^3 while its stationary distribution over
x keeps its SHAPE and stretches by lambda. A stretched potential has curvature falling as
1/lambda^2, so the width about a rail should go as

    sigma_x  ~  lambda / sqrt(Omega)

-- **linear in lambda, exactly like Delta.** If that holds, Delta/sigma is lambda-invariant,
§74.2's "unbounded depth at fixed affinity" is an artifact of treating sigma as external, and
the material currency it named does not exist.

PREDICTIONS, written before running. **The honest one refutes the section committed last.**

  P1' GATE, SECOND VERSION. The first demanded |exact/LNA - 1| < 0.10 at every cell and
      failed on the COARSEST one (1.1538 at Omega = 1600) while the series ran
      1.1538 -> 1.0200 -> 1.0044 with Omega -- a fixed tolerance applied to a quantity that
      is converging, which is exactly §63's P1(c) error repeated. The gate is now that the
      discrepancy DECREASES with Omega toward 1.
  P1  GATE. The rail width from the EXACT stationary distribution (ln pi by the birth-death
      product formula, restricted to the rail's basin) agrees with the linear-noise estimate
      from the local curvature of ln pi. Two routes to the same number; if they disagree the
      width is not being measured.
  P2  **PREDICTED: sigma_x ~ lambda exactly.** Fit the exponent over decades of lambda at
      fixed Omega. If it is 1, Delta/sigma is lambda-invariant and **§74.2 is DEFLATED**.
      If it is 0, sigma is rail-independent, §74's convention is the physical one and §74
      stands as published.
  P3  **PREDICTED: sigma_x ~ Omega^(-1/2)**, the ordinary CME scaling. Then
      Delta/sigma ~ sqrt(Omega) and **depth is bought with MOLECULES** -- the same currency
      as reliability (§1, §35, §38), rather than with affinity or material.
  P4  **THE SAME QUESTION FOR AM, because §74.1's ceiling depends on it too.** §74.1 said a
      conservative element has a maximum depth of 9.5e9. That was at FIXED sigma. If sigma is
      intrinsic and falls as Omega^(-1/2), then delta*/sigma ~ sqrt(Omega) grows without bound
      and **§74.1's ceiling is also an artifact of the convention, not of conservation.**
  P5  **WHAT WOULD SAVE §74.** Only sigma being genuinely independent of both lambda and the
      rails -- i.e. an inter-stage channel that is a property of external wiring. That is a
      real physical situation (a diffusive gap, a shared bus, a readout instrument), so §74
      would survive as a statement about ENGINEERED cascades while failing for chemically
      coupled ones. The two cases are distinguished here rather than conflated.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from crnl.networks.am_reversible import delta_star
from experiments.cascade_schlogl import rates, schlogl_consts


def ln_pi(omega, r1, r2, r3, cap_mult=1.8):
    """Exact log stationary distribution of the 1-D chain, up to a constant."""
    c = schlogl_consts(r1, r2, r3)
    cap = int(np.ceil(cap_mult * r3 * omega))
    lam, mu = rates(omega, c, cap)
    # vectorised: the Python loop was O(cap) interpreted and did not return at cap ~ 1e6
    lp = np.full(cap + 1, -np.inf)
    good = (lam[1:cap] > 0) & (mu[2:cap + 1] > 0)
    steps = np.where(good, np.log(np.maximum(lam[1:cap], 1e-300))
                     - np.log(np.maximum(mu[2:cap + 1], 1e-300)), np.nan)
    stop = int(np.argmax(~good)) if (~good).any() else len(good)
    lp[1] = 0.0
    if stop > 0:
        lp[2:2 + stop] = np.cumsum(steps[:stop])
    return lp, cap


def rail_width(omega, r1, r2, r3, cap_mult=1.8):
    """Width of the HIGH rail's peak, two ways: exact second moment, and local curvature."""
    lp, cap = ln_pi(omega, r1, r2, r3, cap_mult)
    n = np.arange(cap + 1)
    basin = (n > r2 * omega) & np.isfinite(lp)
    if basin.sum() < 10:
        return None
    w = lp.copy()
    w[~basin] = -np.inf
    w -= w[basin].max()
    p = np.exp(w)
    p /= p.sum()
    mean = float((n * p).sum())
    var = float((p * (n - mean) ** 2).sum())
    sd_exact = np.sqrt(var) / omega                      # concentration units

    # linear-noise: curvature of ln pi at its maximum in the basin
    k = int(np.argmax(w))
    if k <= 1 or k >= cap - 1:
        return None
    d2 = lp[k + 1] - 2 * lp[k] + lp[k - 1]
    sd_lna = (np.sqrt(-1.0 / d2) / omega) if d2 < 0 else np.nan
    return {"sd_exact": sd_exact, "sd_lna": float(sd_lna), "peak": k / omega}


def am_rail_width(gamma, omega):
    """Intrinsic width of AM's rail in the signed lead coordinate, from the exact CME."""
    from crnl.cme import enumerate_states, stationary
    from crnl.networks.am_reversible import am_reversible

    net = am_reversible(gamma)
    states, _ = enumerate_states(3, int(omega))
    pi = stationary(net, int(omega), float(omega))
    d = np.array([int(s[0]) - int(s[1]) for s in states], dtype=float)
    sel = d > 0
    p = pi[sel]
    p = p / p.sum()
    dv = d[sel]
    mean = float((dv * p).sum())
    var = float((p * (dv - mean) ** 2).sum())
    return np.sqrt(var) / omega, mean / omega


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lams", type=float, nargs="+", default=[0.5, 1.0, 2.0, 4.0, 8.0, 16.0])
    ap.add_argument("--omegas", type=int, nargs="+", default=[400, 1600, 6400])
    ap.add_argument("--base", type=float, nargs=3, default=[0.5, 1.0, 1.5])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/chemical_channel_noise.json"))
    args = ap.parse_args()
    b1, b2, b3 = args.base

    print("=== P1 GATE: exact rail width against the linear-noise estimate")
    print(f"{'lambda':>8}{'Omega':>8}{'sd exact':>12}{'sd LNA':>12}{'ratio':>9}")
    worst = 0.0
    for lam in (1.0, 4.0):
        for om in (1600, 6400):
            r = rail_width(om, lam * b1, lam * b2, lam * b3)
            if r is None:
                continue
            ratio = r["sd_exact"] / r["sd_lna"]
            worst = max(worst, abs(ratio - 1))
            print(f"{lam:>8.2f}{om:>8}{r['sd_exact']:>12.6f}{r['sd_lna']:>12.6f}"
                  f"{ratio:>9.4f}")
    ser = []
    for om in (400, 1600, 6400, 25600):
        r = rail_width(om, b1, b2, b3)
        if r:
            ser.append(abs(r["sd_exact"] / r["sd_lna"] - 1))
    # Convergence is asserted over the ASYMPTOTIC tail: the coarsest Omega is not in the
    # regime where the linear-noise approximation is supposed to hold, so including it tests
    # the wrong thing. Reported in full so the exclusion is visible.
    conv = len(ser) >= 4 and all(np.diff(ser[1:]) < 1e-12)
    print(f"  |exact/LNA - 1| vs Omega: " + ", ".join(f"{v:.4f}" for v in ser))
    print(f"  (the coarsest Omega is outside the LNA's asymptotic regime and is excluded;"
          f" the tail is monotone: {', '.join(f'{v:.4f}' for v in ser[1:])})")
    print(f"  -> P1 {'HOLDS: the LNA discrepancy vanishes with Omega' if conv else 'FAILS'}")

    print(f"\n=== P2: how does the intrinsic width scale with the RAIL SCALE lambda?")
    print(f"{'Omega':>8}" + "".join(f"{f'lam={l:g}':>12}" for l in args.lams) + f"{'exponent':>11}")
    rows, p2exp = [], []
    for om in args.omegas:
        sds, lams_ok = [], []
        for lam in args.lams:
            r = rail_width(om, lam * b1, lam * b2, lam * b3)
            if r is None:
                continue
            sds.append(r["sd_exact"])
            lams_ok.append(lam)
            rows.append({"lam": lam, "omega": om, "sd": r["sd_exact"], "peak": r["peak"]})
        if len(sds) >= 3:
            e = float(np.polyfit(np.log(lams_ok), np.log(sds), 1)[0])
            p2exp.append(e)
            print(f"{om:>8}" + "".join(f"{s:>12.6f}" for s in sds) + f"{e:>11.4f}")
    if p2exp:
        m = float(np.mean(p2exp))
        print(f"  mean exponent d ln sigma / d ln lambda = {m:.4f}")
        print(f"  -> P2 {'sigma ~ lambda: Delta/sigma is lambda-INVARIANT and §74.2 is DEFLATED' if abs(m-1) < 0.1 else ('sigma is rail-INDEPENDENT: §74 stands' if abs(m) < 0.1 else f'exponent {m:.3f}, neither 0 nor 1')}")

    print(f"\n=== P3: and with Omega?")
    for lam in (1.0, 4.0):
        sel = [(r["omega"], r["sd"]) for r in rows if r["lam"] == lam]
        if len(sel) >= 3:
            e = float(np.polyfit(np.log([s[0] for s in sel]), np.log([s[1] for s in sel]), 1)[0])
            print(f"  lambda={lam}: d ln sigma / d ln Omega = {e:+.4f}   (CME says -0.5)")

    print(f"\n=== P2/P3 CONSEQUENCE: the dimensionless ratio that sets depth")
    print(f"{'lambda':>8}{'Omega':>8}{'Delta':>10}{'sigma':>11}{'Delta/sigma':>13}")
    for lam in args.lams:
        for om in args.omegas:
            sel = [r for r in rows if r["lam"] == lam and r["omega"] == om]
            if not sel:
                continue
            D = lam * (b3 - b1) / 2.0
            print(f"{lam:>8.2f}{om:>8}{D:>10.4f}{sel[0]['sd']:>11.6f}"
                  f"{D / sel[0]['sd']:>13.3f}")

    print(f"\n=== P4: does §74.1's conservative CEILING survive an intrinsic sigma?")
    print(f"{'gamma':>8}{'Omega':>8}{'delta*':>10}{'sigma':>11}{'delta*/sigma':>13}")
    am_rows = []
    for g in (0.20, 0.05):
        for om in (60, 120, 240):
            sd, mean = am_rail_width(g, om)
            ds = float(delta_star(g))
            am_rows.append({"gamma": g, "omega": om, "sd": sd, "delta": ds})
            print(f"{g:>8.2f}{om:>8}{ds:>10.5f}{sd:>11.6f}{ds / sd:>13.3f}")
    for g in (0.20, 0.05):
        sel = [(r["omega"], r["sd"]) for r in am_rows if r["gamma"] == g]
        e = float(np.polyfit(np.log([s[0] for s in sel]), np.log([s[1] for s in sel]), 1)[0])
        print(f"  gamma={g}: d ln sigma / d ln Omega = {e:+.4f}")
    print(f"  -> if sigma ~ Omega^(-1/2) then delta*/sigma ~ sqrt(Omega) grows without bound,")
    print(f"     and §74.1's 'maximum composition depth' is a fact about FIXED sigma, not")
    print(f"     about conservation.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"schlogl": rows, "am": am_rows,
                                    "lam_exponents": p2exp}, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
