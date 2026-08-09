"""T-COST-i: is §47's 8% shortfall physics, or endpoints and finite Omega?

§47 predicted the T_det/MFPT gap with no fitted constant and got pred/gap = 0.921 (rho)
and 0.931 (gamma). The shortfall was PREDICTED IN ADVANCE as next-order quasi-steady
state -- and §47 then refused to attribute it there, because it does not behave like one:
pred/gap is 0.936 where the gap is 0.150 and 0.780 where it is 0.051, **larger shortfall
at smaller gap**, the wrong direction for an O(gap) correction.

**TWO INSTRUMENT EFFECTS ARE CANDIDATES AND BOTH WERE MEASURED BEFORE THESE PREDICTIONS
WERE WRITTEN.**

**(a) The endpoints do not match -- rule 11, literally.** `T_det = int ddelta/mu` is taken
over the NOMINAL range [eps*delta*, theta*delta*], while the CME first-passage runs from
the LATTICE start d0/Omega to the LATTICE threshold thr/Omega. Measured at Omega = 300:

    cell                eps_nom  eps_real   theta_nom  theta_real
    gamma=0.20           0.3500    0.3470      0.8000      0.8002
    gamma=0.07           0.3500    0.3498      0.8000      0.7995
    gamma=0.35           0.3500    0.3514      0.8000      0.8025
    gamma=0.20, rho=32   0.3500    0.3497      0.8000      0.8010

**The sign tracks the shortfall.** gamma = 0.35 has the largest overshoot (start LATER
than nominal, so the true traversal is SHORTER than T_det computes, inflating the gap) and
also the worst pred/gap, 0.780. gamma = 0.07 is nearly matched and sits at 0.936.

**(b) The additive signature.** The absolute shortfall `gap - pred` clusters near +0.003
across cells whose gaps differ 7x -- rho = 8, 16, 32 give 0.0031, 0.0029, 0.0029 at gaps of
0.039, 0.022, 0.013. **A smooth physical correction cannot produce nearly identical
absolute values at wildly different gaps; a shared lattice rounding can.**

PREDICTIONS, written before running:

  P1  GATE. Recomputing `T_det` between the REALIZED lattice endpoints must change it by
      the amount the endpoint offsets imply and no more -- report both integrals. If
      matched and nominal T_det agree to better than 0.1% everywhere, effect (a) is too
      small to matter and P2 is dead on arrival.
  P2  THE TEST, and it is a PER-CELL prediction rather than an average. Matching endpoints
      must help **gamma = 0.35 most** (realized eps 0.3514, the largest overshoot, worst
      pred/gap 0.780) and **gamma = 0.07 least** (realized eps 0.3498, nearly matched,
      pred/gap 0.936). An improvement that is uniform across cells would mean something
      else changed and would NOT support the endpoint account.
  P3  With endpoints matched, sweep Omega = 300/500/700/1000. **pred/gap -> 1.** Any
      residual Omega-dependence is the genuine absorption effect -- §39.1's candidate (iv),
      which §39.2 left live for the time and which nothing has yet measured.
  P4  THE SIGN IS FORCED. Absorption at a threshold selects the leading edge of the packet,
      so the MFPT sits BELOW the mean arrival, so the measured gap exceeds the pure lag and
      **pred/gap must approach 1 FROM BELOW.** Convergence from above would contradict the
      absorption picture and would mean the lag model overshoots.
  P5  REFUTING OUTCOME. If matched endpoints do not reduce the per-cell scatter, the
      lattice is not the cause. If pred/gap then fails to converge to 1 with Omega, the lag
      model is systematically off by a measured factor and T-COST-i stays open with a
      number attached rather than a story.
  P6  If P2 and P3 both hold, **the lag model is exact** and §39.2's mechanism is confirmed
      in absolute terms: `T_det/MFPT - 1 = <eps>_time` with no constant, the finite-Omega
      residual being absorption bias. §47's "8% shortfall" would then be withdrawn as
      instrument, exactly as §47 withdrew §46's.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.cme import first_passage
from crnl.networks.am_reversible import reverse_pairing
from experiments.arrhenius_optimum import am_rho, delta_star_rho
from experiments.cost_absolute import sigma_and_mu
from experiments.lag_absolute import predict
from experiments.slaving_axis import sep_of, slaved


def traversal(net, d_lo, d_hi, pairing, n=1201):
    xs = np.linspace(d_lo, d_hi, n)
    inv = []
    for x in xs:
        sm = sigma_and_mu(net, float(x), pairing)
        if sm is None or sm[1] <= 0:
            return None
        inv.append(1.0 / sm[1])
    return float(np.trapezoid(inv, xs))


def cell(gamma, rho, omega, eps_frac, theta, h):
    ds = delta_star_rho(gamma, rho)
    if ds <= 0:
        return None
    net = am_rho(gamma, rho)
    pairing = reverse_pairing(net)

    st = slaved(net, eps_frac * ds)
    if st is None:
        return None
    nb = int(round(st[2] * omega))
    rest = omega - nb
    d0 = max(1, int(round(eps_frac * ds * omega)))
    if (rest - d0) % 2:
        d0 -= 1
    thr = max(2, int(round(theta * ds * omega)))
    if thr <= d0 or rest - d0 < 0:
        return None
    n0 = np.array([(rest + d0) // 2, (rest - d0) // 2, nb], dtype=np.int64)

    fp = first_passage(net, int(omega), float(omega), n0,
                       lambda s, t=thr: abs(int(s[0]) - int(s[1])) >= t)
    mfpt = float(fp["mean_time"])
    if not np.isfinite(mfpt) or mfpt <= 0:
        return None

    lo_nom, hi_nom = eps_frac * ds, theta * ds
    lo_real, hi_real = d0 / omega, thr / omega
    t_nom = traversal(net, lo_nom, hi_nom, pairing)
    t_real = traversal(net, lo_real, hi_real, pairing)
    if t_nom is None or t_real is None:
        return None

    p_nom = predict(net, lo_nom, hi_nom, pairing, h)
    p_real = predict(net, lo_real, hi_real, pairing, h)
    if p_nom is None or p_real is None:
        return None

    g_nom = t_nom / mfpt - 1.0
    g_real = t_real / mfpt - 1.0
    return {"gamma": gamma, "rho": rho, "omega": omega, "sep": float(sep_of(net)[0]),
            "eps_real": d0 / omega / ds, "theta_real": thr / omega / ds,
            "t_nom": t_nom, "t_real": t_real, "mfpt": mfpt,
            "gap_nom": g_nom, "gap_real": g_real,
            "pred_nom": p_nom, "pred_real": p_real,
            "ratio_nom": p_nom / g_nom if g_nom else float("nan"),
            "ratio_real": p_real / g_real if g_real else float("nan")}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cells", type=str, nargs="+",
                    default=["0.07:1", "0.20:1", "0.35:1", "0.20:0.5", "0.20:32"])
    ap.add_argument("--omegas", type=int, nargs="+", default=[300, 500, 700, 1000])
    ap.add_argument("--eps", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--h", type=float, default=1e-4)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/lag_endpoints.json"))
    args = ap.parse_args()

    t0 = time.time()
    cells = [(float(c.split(":")[0]), float(c.split(":")[1])) for c in args.cells]

    print("=== P1/P2: does matching the endpoints to the lattice help, and per cell?")
    print(f"{'gamma':>6}{'rho':>6}{'Om':>6}{'eps_real':>10}{'T_nom':>9}{'T_real':>9}"
          f"{'dT%':>7}{'pred/gap nom':>14}{'pred/gap real':>15}")
    rows = []
    for g, r in cells:
        for om in args.omegas:
            c = cell(g, r, om, args.eps, args.theta, args.h)
            if c is None:
                print(f"{g:>6.2f}{r:>6.1f}{om:>6}   dropped")
                continue
            rows.append(c)
            dt = 100 * (c["t_real"] - c["t_nom"]) / c["t_nom"]
            print(f"{g:>6.2f}{r:>6.1f}{om:>6}{c['eps_real']:>10.4f}"
                  f"{c['t_nom']:>9.4f}{c['t_real']:>9.4f}{dt:>7.2f}"
                  f"{c['ratio_nom']:>14.3f}{c['ratio_real']:>15.3f}")
        print()

    at300 = [c for c in rows if c["omega"] == args.omegas[0]]
    print("=== P2: the per-cell prediction -- gamma=0.35 most, gamma=0.07 least")
    for c in at300:
        if c["rho"] == 1.0:
            imp = abs(c["ratio_nom"] - 1) - abs(c["ratio_real"] - 1)
            print(f"  gamma={c['gamma']:.2f}: |1-ratio| {abs(c['ratio_nom']-1):.4f}"
                  f" -> {abs(c['ratio_real']-1):.4f}   improvement {imp:+.4f}")
    ga = {c["gamma"]: c for c in at300 if c["rho"] == 1.0}
    if 0.35 in ga and 0.07 in ga:
        i35 = abs(ga[0.35]["ratio_nom"] - 1) - abs(ga[0.35]["ratio_real"] - 1)
        i07 = abs(ga[0.07]["ratio_nom"] - 1) - abs(ga[0.07]["ratio_real"] - 1)
        print(f"  -> gamma=0.35 improves {i35:+.4f}, gamma=0.07 improves {i07:+.4f}"
              f"   {'AS PREDICTED' if i35 > i07 else 'NOT as predicted'}")

    print(f"\n=== P3/P4: with endpoints matched, does pred/gap -> 1 with Omega?")
    print(f"{'gamma':>6}{'rho':>6}" + "".join(f"{f'Om={o}':>12}" for o in args.omegas))
    for g, r in cells:
        line = f"{g:>6.2f}{r:>6.1f}"
        for om in args.omegas:
            m = [c for c in rows if c["gamma"] == g and c["rho"] == r
                 and c["omega"] == om]
            line += f"{m[0]['ratio_real']:>12.4f}" if m else f"{'--':>12}"
        print(line)

    for om in args.omegas:
        v = np.array([c["ratio_real"] for c in rows if c["omega"] == om])
        if v.size:
            print(f"  Omega={om:>5}: mean {v.mean():.4f}  spread"
                  f" {100*(v.max()-v.min())/v.mean():>5.1f}%"
                  f"  ({'below 1' if v.mean() < 1 else 'ABOVE 1'})")

    print(f"\n=== P4: is the approach FROM BELOW, as absorption requires?")
    means = [np.mean([c["ratio_real"] for c in rows if c["omega"] == om])
             for om in args.omegas]
    print("  " + "  ".join(f"Om={o}:{m:.4f}" for o, m in zip(args.omegas, means)))
    below = all(m < 1.0 for m in means)
    rising = means[-1] > means[0]
    print(f"  all below 1: {below};  rising with Omega: {rising}"
          f"  -> {'consistent with absorption bias' if below and rising else 'NOT the absorption signature'}")

    print(f"\n=== P5/P6: verdict")
    last = np.array([c["ratio_real"] for c in rows if c["omega"] == args.omegas[-1]])
    nom_last = np.array([c["ratio_nom"] for c in rows if c["omega"] == args.omegas[-1]])
    if last.size:
        print(f"  at Omega={args.omegas[-1]}: pred/gap nominal endpoints"
              f" {nom_last.mean():.4f}, matched endpoints {last.mean():.4f}")
        if abs(last.mean() - 1) < 0.05:
            print("  -> P6: the lag model is EXACT within 5%. §47's 8% shortfall was")
            print("     instrument -- endpoints and finite Omega -- and is withdrawn.")
        else:
            print(f"  -> P5: a {100*(1-last.mean()):.1f}% shortfall survives both fixes.")
            print("     T-COST-i stays open with a measured value, not a story.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
