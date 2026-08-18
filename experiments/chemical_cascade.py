"""The founding claim's untested half: a cascade whose stages are coupled by CHEMISTRY

Every multi-stage cascade this project has built couples stages through `channel()` -- run for
time t, read the output through a Gaussian channel, RE-SEED the next stage (§12, §71, §72). §75
argued the physical cascade has no readout at all: stage i's output SPECIES is stage i+1's input.
§76-§79 then computed epsilon for ONE element and used D_max = c*/epsilon analytically.
**The multi-stage chemically-coupled chain was never built.**

TWO COUPLINGS, because "drive the next stage with the previous one" does not name a mechanism:

    SOURCE     stage i+1's influx is driven:      lam = k1a n(n-1)/Om + (k2b/r3) n_up
    CATALYTIC  stage i+1's autocatalysis is:      lam = (k1a/r3)(n_up/Om) n(n-1)/Om + k2b Om

Both are exactly neutral when the upstream sits at r3 -- they reproduce the isolated element.

**WHAT THE FIRST VERSION OF THIS FILE GOT WRONG, recorded because it is the point.** It gated only
NEUTRALITY: with the upstream pinned at r3, is the downstream the same element? Both schemes pass
that. It then measured error accumulation under SOURCE coupling and found it strongly sublinear --
adding a stage cost only 0.095 of an isolated element's error -- which reads like restoration
filtering upstream errors, the founding claim confirmed. **It is an artifact.** Decomposing the
joint law exactly gives P(stage 2 low | stage 1 low) = 0.0017, 0.0086, 0.026, 0.064, 0.136 at
t = 1..16: **the flip does not propagate at all.** Starving the influx by r1/r3 = 21x leaves the
downstream still bistable, because the autocatalysis carries the landscape and the source term is
a small correction. Under SOURCE coupling these are three nearly independent elements, and "error
does not accumulate" is true because nothing is connected.

**A NEUTRALITY GATE IS HALF A GATE.** It tests the null condition -- that the coupling does
nothing when the upstream is correct. A cascade also has to TRANSMIT. That is P1(b) below and it
is the gate the first version lacked.

PREDICTIONS.

  P1  GATE, and it takes both halves.
      (a) NEUTRALITY: with the upstream pinned at r3, the downstream must be the isolated element
          to machine precision, in propensities and in roots.
      (b) TRANSMISSION: with the upstream pinned at r1, the downstream must LOSE its high rail --
          the drift must become monostable low. A scheme that stays bistable does not carry a
          signal. **Predicted: SOURCE fails this and CATALYTIC passes** -- found, not guessed,
          and recorded above.
  P2  **DOES ERROR ACCUMULATE LINEARLY IN DEPTH, under a coupling that transmits?** P_flip of the
      last stage against D, exactly, from the joint CME. **Predicted: close to the union bound
      1-(1-eps)^D**, because an upstream flip lasts ~exp(A*Omega) and propagates in O(1/|f'|) =
      0.15 here. Sublinear would mean the chain genuinely filters, which is the founding claim.
  P3  **IS THE CHAIN'S PER-STAGE EPSILON THE ISOLATED ELEMENT'S?** Under SOURCE coupling the
      answer was already measured: stage 2's own escape (conditioned on the upstream staying
      high) ran 1.06-1.09 times the isolated element's. **Predicted: larger again under CATALYTIC
      coupling**, and by more, since the upstream now multiplies the autocatalytic term and its
      fluctuation enters the barrier directly rather than as an additive source.
  P4  **THE DECOMPOSITION, which is what caught the artifact.** Report P(down low | up low) --
      the transmission probability -- alongside every accumulation number. An accumulation result
      without it cannot distinguish "restoration filters errors" from "the stages are not
      connected".
  P5  **THE STAGE TIME, measured not imported.** §80 used t = 2.0 from §71/§72's readout protocol;
      a continuously coupled chain has no such clock. Report the propagation time from |f'(r3)|
      and the time for P(down low | up low) to saturate.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from experiments.cascade_schlogl import schlogl_consts

RAILS = (0.15, 1.0, 3.1827)
SCHEMES = ("source", "catalytic", "hill")
HILL_N, HILL_K = 4.0, 1.0        # cooperative coupling, half-max AT THE SADDLE


def rates_stage(n_self, n_up, omega, c, r3, first, scheme):
    k1a, k1r, k2b, k2r = c
    mu = k1r * n_self * (n_self - 1.0) * (n_self - 2.0) / omega ** 2 + k2r * n_self
    auto = k1a * n_self * (n_self - 1.0) / omega
    if first:
        lam = auto + k2b * omega
    elif scheme == "source":
        lam = auto + (k2b / r3) * n_up
    elif scheme == "catalytic":
        lam = (n_up / omega / r3) * auto + k2b * omega
    else:
        lam = hill(n_up / omega, r3) * auto + k2b * omega
    return max(lam, 0.0), max(mu, 0.0)


def hill(x_up, r3):
    """Saturating transfer, normalised so h(r3) = 1 -- the ingredient a transistor has."""
    f = lambda z: z ** HILL_N / (HILL_K ** HILL_N + z ** HILL_N)
    return f(max(x_up, 0.0)) / f(r3)


def downstream_roots(x_up, c, r3, scheme):
    """Positive roots of the downstream drift with the upstream PINNED at x_up."""
    k1a, k1r, k2b, k2r = c
    if scheme == "source":
        coeff = [-k1r, k1a, -k2r, (k2b / r3) * x_up]
    elif scheme == "catalytic":
        coeff = [-k1r, (k1a / r3) * x_up, -k2r, k2b]
    else:
        coeff = [-k1r, k1a * hill(x_up, r3), -k2r, k2b]
    r = np.roots(coeff)
    return np.sort([z.real for z in r if abs(z.imag) < 1e-9 and z.real > 1e-12])


def build(D, omega, scheme, rails=RAILS, cap_mult=1.25):
    r1, r2, r3 = rails
    c = schlogl_consts(r1, r2, r3)
    cap = int(np.ceil(cap_mult * r3 * omega))
    m = cap + 1
    N = m ** D
    strides = [m ** (D - 1 - i) for i in range(D)]
    rows, cols, vals = [], [], []
    for idx in range(N):
        ns = [(idx // strides[i]) % m for i in range(D)]
        tot = 0.0
        for i in range(D):
            lam, mu = rates_stage(ns[i], ns[i - 1] if i else 0, omega, c, r3, i == 0, scheme)
            if ns[i] < cap and lam > 0:
                rows.append(idx); cols.append(idx + strides[i]); vals.append(lam); tot += lam
            if ns[i] > 0 and mu > 0:
                rows.append(idx); cols.append(idx - strides[i]); vals.append(mu); tot += mu
        rows.append(idx); cols.append(idx); vals.append(-tot)
    return sp.csr_matrix((vals, (rows, cols)), shape=(N, N)), cap, m, strides, c


def seed_high(D, omega, m, strides, r3):
    p = np.zeros(m ** D)
    p[sum(int(round(r3 * omega)) * s for s in strides)] = 1.0
    return p


def masks(D, m, strides, omega, r2):
    idx = np.arange(m ** D)
    return [((idx // strides[k]) % m) < r2 * omega for k in range(D)]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=30)
    ap.add_argument("--times", type=float, nargs="+", default=[1, 2, 4, 8, 16])
    ap.add_argument("--depths", type=int, nargs="+", default=[1, 2])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/chemical_cascade.json"))
    args = ap.parse_args()
    r1, r2, r3 = RAILS
    c = schlogl_consts(r1, r2, r3)
    k1a, k1r, k2b, k2r = c
    om = args.omega
    out = {"rails": list(RAILS), "omega": om}

    print("=== P1(a) NEUTRALITY and P1(b) TRANSMISSION -- a neutrality gate is half a gate")
    print(f"{'scheme':>11}{'upstream':>11}{'downstream roots':>34}{'verdict':>26}")
    gate = {}
    for scheme in SCHEMES:
        worst = 0.0
        for nn in (5, 30, 60, 90):
            a, b = rates_stage(float(nn), r3 * om, om, c, r3, False, scheme)
            a0, b0 = rates_stage(float(nn), 0.0, om, c, r3, True, scheme)
            worst = max(worst, abs(a - a0), abs(b - b0))
        hi = downstream_roots(r3, c, r3, scheme)
        lo = downstream_roots(r1, c, r3, scheme)
        transmits = len(lo) == 1 and lo[0] < r2
        gate[scheme] = {"neutral_err": worst, "roots_hi": hi.tolist(),
                        "roots_lo": lo.tolist(), "transmits": bool(transmits)}
        print(f"{scheme:>11}{'at r3':>11}{str(np.round(hi, 4)):>34}"
              f"{('neutral, |err| %.0e' % worst):>26}")
        print(f"{'':>11}{'at r1':>11}{str(np.round(lo, 4)):>34}"
              f"{('TRANSMITS' if transmits else 'still bistable: NO SIGNAL'):>26}")
    print(f"\n  NOISE MARGIN: the upstream value at which the downstream loses its high rail")
    print(f"{'scheme':>11}{'x_up critical':>15}{'(r3-xc)/r3':>13}{'in rail widths':>16}")
    from experiments.derive_eta import schlogl_V
    sd_up = float(np.sqrt(schlogl_V(r1, r2, r3) / om))
    for scheme in SCHEMES:
        xs = np.linspace(r1, r3, 4001)
        xc = None
        for x in xs[::-1]:
            if len(downstream_roots(x, c, r3, scheme)) < 3:
                xc = float(x)
                break
        if xc is None:
            print(f"{scheme:>11}{'never loses it':>15}{'--':>13}{'--':>16}")
            gate[scheme]["x_crit"] = None
            continue
        gate[scheme]["x_crit"] = xc
        print(f"{scheme:>11}{xc:>15.4f}{(r3-xc)/r3:>13.4f}{(r3-xc)/sd_up:>16.2f}")
    print(f"  (upstream rail width sigma = sqrt(V/Omega) = {sd_up:.4f} at Omega = {om})")
    print(f"  **A cascade needs the downstream to survive ordinary upstream fluctuation.**")
    print(f"  A margin of a few sigma is a noise margin; a fraction of a sigma is not.")
    out["gate"] = gate
    ok_a = all(g["neutral_err"] < 1e-9 for g in gate.values())
    print(f"  -> P1(a) {'HOLDS for both schemes' if ok_a else 'FAILS'}")
    print(f"  -> P1(b) source {'transmits' if gate['source']['transmits'] else 'DOES NOT TRANSMIT'}"
          f";  catalytic {'transmits' if gate['catalytic']['transmits'] else 'DOES NOT TRANSMIT'}")
    live = [s for s in SCHEMES if gate[s]["transmits"]]
    print(f"  -> only {live} can be called a cascade. The other is measured anyway, as the"
          f" control that shows what a non-transmitting chain looks like.")

    print(f"\n=== P2/P3/P4: exact joint CME at Omega = {om}")
    fp = (2 * k1a * r3 - 3 * k1r * r3 ** 2 - k2r)
    print(f"  P5: |f'(r3)| = {abs(fp):.3f}, so one stage relaxes in {1/abs(fp):.3f} time units;"
          f" §80 used t = 2.0, i.e. {2*abs(fp):.0f} relaxation times -- ample to propagate")
    res = {}
    for scheme in SCHEMES:
        print(f"\n  --- {scheme} coupling")
        print(f"{'t':>7}{'P(s1 lo)':>12}{'P(sD lo)':>12}{'union 1-(1-e)^D':>17}"
              f"{'P(sD lo|s1 lo)':>16}{'own/eps1':>10}")
        rows = []
        for D in args.depths:
            if D < 2:
                continue
            Q, cap, m, strides, _ = build(D, om, scheme)
            p = seed_high(D, om, m, strides, r3)
            mk = masks(D, m, strides, om, r2)
            tprev = 0.0
            for t in args.times:
                p = spla.expm_multiply(Q.T * (t - tprev), p)
                tprev = t
                P1 = float(p[mk[0]].sum())
                PD = float(p[mk[-1]].sum())
                both = float(p[mk[-1] & mk[0]].sum())
                own = float(p[mk[-1] & ~mk[0]].sum())
                rows.append({"D": D, "t": t, "P1": P1, "PD": PD, "cond": both / P1,
                             "own_over_eps1": own / P1,
                             "union": 1 - (1 - P1) ** D})
                print(f"{t:>7.1f}{P1:>12.4e}{PD:>12.4e}{1-(1-P1)**D:>17.4e}"
                      f"{both/P1:>16.4f}{own/P1:>10.4f}")
        res[scheme] = rows
    out["runs"] = res

    print("\n=== P4: the decomposition is the verdict")
    for scheme in SCHEMES:
        rr = [r for r in res[scheme] if r["D"] == max(args.depths)]
        if not rr:
            continue
        cond = [r["cond"] for r in rr]
        ratio = [r["PD"] / r["union"] for r in rr]
        print(f"  {scheme:>10}: P(down low | up low) = "
              + ", ".join(f"{v:.4f}" for v in cond))
        print(f"  {'':>10}  measured/union       = " + ", ".join(f"{v:.3f}" for v in ratio))
        if max(cond) < 0.25:
            print(f"  {'':>10}  -> the signal does not propagate: any 'no accumulation' here is")
            print(f"  {'':>10}     the stages being disconnected, NOT restoration filtering.")
        else:
            print(f"  {'':>10}  -> the signal propagates, so the accumulation number is about"
                  f" composition")

    print("\n=== P6 (rule 17): is the NOISE MARGIN the controlling variable, or the scheme?")
    print("    Sweep the Hill exponent AND its half-max, which move the margin without")
    print("    touching the element, the rails, or Omega. If margin is the variable, both")
    print("    knobs must collapse onto one curve.")
    print(f"{'n':>5}{'K':>6}{'x_crit':>9}{'margin/sigma':>14}{'own/eps1':>11}{'meas/union':>12}")
    global HILL_N, HILL_K
    n0, k0 = HILL_N, HILL_K
    p6 = []
    # The two knobs must span OVERLAPPING margins, or "do they lie on one curve" is never
    # tested -- the first sweep put both K variants outside the n range entirely, which is
    # rule 19's extrapolation trap wearing a different hat.
    combos = [(1.0, 1.0), (2.0, 1.0), (3.0, 1.0), (4.0, 1.0), (6.0, 1.0), (8.0, 1.0),
              (12.0, 1.0), (4.0, 0.6), (4.0, 0.8), (4.0, 1.2), (4.0, 1.4), (4.0, 1.6),
              (6.0, 1.3), (8.0, 1.5)]
    for n, K in combos:
        HILL_N, HILL_K = n, K
        xs = np.linspace(r1, r3, 4001)
        xc = next((float(x) for x in xs[::-1]
                   if len(downstream_roots(x, c, r3, "hill")) < 3), None)
        if xc is None:
            print(f"{n:>5}{K:>6}   never loses the high rail")
            continue
        Q, cap, m, strides, _ = build(2, om, "hill")
        q = spla.expm_multiply(Q.T * 2.0, seed_high(2, om, m, strides, r3))
        mk = masks(2, m, strides, om, r2)
        P1 = float(q[mk[0]].sum()); PD = float(q[mk[-1]].sum())
        own = float(q[mk[-1] & ~mk[0]].sum())
        p6.append({"n": n, "K": K, "x_crit": xc, "margin": (r3 - xc) / sd_up,
                   "own": own / P1, "vs_union": PD / (1 - (1 - P1) ** 2)})
        print(f"{n:>5}{K:>6}{xc:>9.4f}{(r3-xc)/sd_up:>14.2f}{own/P1:>11.3f}"
              f"{PD/(1-(1-P1)**2):>12.3f}")
    HILL_N, HILL_K = n0, k0
    out["p6"] = p6
    if len(p6) >= 8:
        mg = np.array([z["margin"] for z in p6])
        ow = np.log(np.array([z["own"] for z in p6]))
        a, b = np.polyfit(mg, ow, 1)
        resid = ow - (a * mg + b)
        isK = np.array([z["K"] != 1.0 for z in p6])
        rms_n = float(np.sqrt((resid[~isK] ** 2).mean()))
        rms_K = float(np.sqrt((resid[isK] ** 2).mean()))
        overlap = (mg[isK].min() < mg[~isK].max()) and (mg[isK].max() > mg[~isK].min())
        srt = sorted(p6, key=lambda z: z["margin"])
        print(f"  log(own/eps1) vs margin: slope {a:.3f}/sigma, RMS residual"
              f" {float(np.sqrt((resid**2).mean())):.3f} in log")
        print(f"    n-knob RMS {rms_n:.3f}   K-knob RMS {rms_K:.3f}"
              f"   ranges overlap: {overlap}")
        print(f"    margin {srt[0]['margin']:.2f} sigma -> {srt[0]['own']:.1f}x;"
              f"  margin {srt[-1]['margin']:.2f} sigma -> {srt[-1]['own']:.2f}x")
        ok6 = overlap and rms_K < 3 * max(rms_n, 0.02) and abs(a) > 0.5
        print(f"  -> P6 {'HOLDS: both knobs collapse onto ONE curve in the margin -- moving the margin by changing the exponent or the half-max gives the same per-stage error, so the MARGIN is the controlling variable and not the scheme parameters' if ok6 else 'FAILS: the two knobs do not collapse, so the margin is a correlate and not the control'}")
        print("     (strict monotonicity is NOT the test here and the first version used it,")
        print("      failing on a tied pair 0.8% apart in margin -- and with the K points")
        print("      outside the n range, where the collapse was never actually tested.)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
