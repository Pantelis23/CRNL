"""T-COST-h: the slow-manifold lag, predicted ABSOLUTELY and with no fitted parameter

§39.2 identified the `T_det/MFPT` gap as the slow-manifold lag and measured it as C/sep,
finding C axis-dependent. §46 then found the 1/sep SCALING axis-dependent too: it holds on
the T axis (12.3%) and the gamma axis (22.2%) and fails on the rho axis (72.1%), where
gap*sep is U-shaped rather than constant. Three sections have now FITTED this gap. Rule 16
says a model that is only ever fitted is never tested.

**SINGULAR PERTURBATION GIVES IT IN CLOSED FORM, AND THE FITTED CONSTANT DISAPPEARS.**
Work in (delta, s) with s = x + y, b = 1 - s, and write

    mu(delta, s) = d(delta)/dt ,      nu(delta, s) = ds/dt

The slaved manifold is nu(delta, s*(delta)) = 0. As delta advances the manifold moves and
the pool lags. Near the manifold nu ~ (dnu/ds)(s - s*), and in quasi-steady state the lag
Delta = s - s* satisfies (ds*/ddelta)*mu = (dnu/ds)*Delta, so

    Delta = (ds*/ddelta) * mu / (dnu/ds)

The lag changes the signal drift by (dmu/ds)*Delta, and **mu cancels**:

    eps(delta) = (dmu/ds)*Delta / mu = (dmu/ds) * (ds*/ddelta) / (dnu/ds)          (*)

a dimensionless property of the vector field alone. Then with T = int ddelta/mu,

    **T_det/MFPT - 1  ~  <eps>_time  =  int eps (ddelta/mu) / int (ddelta/mu)**

**No constant. No fit.** Both sides are computable, so the model can be tested in absolute
terms against a quantity §39.2 obtained exactly -- which is precisely what rule 16 demands
and what §22.2's withdrawn convolution failed to do.

Note (*) is NOT |lambda_slow|/|lambda_fast|. The reduced slow eigenvalue is
lambda_red = dmu/ddelta + (dmu/ds)(ds*/ddelta), so the lag term is lambda_red - dmu/ddelta,
and it coincides with the eigenvalue ratio only when the manifold's motion happens to be
set by the slow eigenvalue. **That is the whole reason it can succeed where 1/sep failed.**

PREDICTIONS, written before running:

  P1  GATE, and rule 13 applies: (*) is built from three finite differences, so the step h
      is the approximation's own numerical parameter and is a SECOND AXIS. Convergence in h
      is checked WITHIN a cell before any comparison BETWEEN cells. If eps moves by more
      than 1% between h and h/2, nothing below is admissible.
  P2  THE TEST, absolute. `<eps>_time / (T_det/MFPT - 1)` ~ 1 on all three axes. A ratio
      near 1 means the lag picture is right and the fitted C of §39.2 was the average of
      (*) over whatever path that axis happened to trace.
  P3  THE DISCRIMINATING REQUIREMENT. It must work where 1/sep WORKED (gamma) **and** where
      1/sep FAILED (rho, 72.1% spread). Reproducing gamma alone is no advance -- 1/sep
      already does that. **rho is the test.**
  P4  SHAPE, not just magnitude. §46 found gap*sep U-shaped in rho, bottoming near 0.545 at
      rho = 2-3 and rising to 0.946 and 1.082 at the ends. The prediction must reproduce
      that shape per cell, not merely the average. Reported cell by cell.
  P5  CONSISTENCY with what already works. On the gamma axis, where 1/sep holds with
      C = 0.5963, the prediction must satisfy <eps>_time * sep ~ 0.6 -- i.e. it must
      EXPLAIN the fitted constant rather than merely coexist with it.
  P6  THE REFUTING OUTCOME, named in advance. If <eps>_time misses on rho as badly as 1/sep
      did, the lag PICTURE -- not just its 1/sep proxy -- is wrong on that axis, and
      §39.2's identification of the mechanism needs qualifying on a third axis. Given that
      the gap is 8-13% and quasi-steady state is a leading-order approximation in exactly
      that small parameter, **a systematic 10-20% shortfall is expected and would NOT be a
      refutation**; a factor of 2, or the wrong shape in rho, would be.

=============================================================================
SECOND PASS -- THE PREDICTION MAY HAVE CAUGHT AN ERROR IN §46'S MEASUREMENT.

The first pass gives pred/gap = 0.921 (rho) and 0.931 (gamma), a systematic ~8% shortfall
of exactly the predicted sign and size, and <eps>*sep = 0.627 on the gamma axis against the
0.5963 that §46 FITTED. But P4 failed in a revealing way:

    predicted eps*sep   0.664 0.646 0.637 0.631 0.633 0.641 0.648 0.657 0.662 0.672 0.677
    measured  gap*sep   0.845 0.788 0.657 0.737 0.745 0.571 0.685 0.539 0.720 0.773 0.871

**The prediction is nearly FLAT -- 7% across the whole rho axis -- while the measurement
scatters 62%.** The model says the 1/sep law SHOULD hold on rho, with C ~ 0.65. §46 said it
fails there, spread 72.1%.

**One of the two is wrong, and rule 13 says look at the instrument first.** §46 computed
that spread at Omega = 200 and never checked Omega-convergence of the GAP -- while its own
table shows Omega = 200 and 300 differing by 16% in a single cell (0.0805 vs 0.0938 at
rho = 1), which alone moves gap*sep by 17%. §39.2 established Omega-convergence of this
quantity only over Omega = 400-800. **I checked convergence in h for the prediction and
never checked it in Omega for the measurement -- rule 13 applied to one side only.**

  P7  The measured gap is NOT converged at Omega = 200-300. At Omega = 300/500/700 the
      per-cell gap settles and gap*sep tightens substantially toward the predicted ~0.65.
      Cells at rho = 3 and rho = 6, which produced §46's two low outliers (0.571, 0.539),
      should move UP the most.
  P8  If the scatter dissolves, **§46's headline is withdrawn**: the 1/sep scaling does not
      fail on rho, it was measured at insufficient Omega. The absolute prediction would
      then have predicted a measurement error, which is a stronger outcome than agreeing
      with a correct one. Rule 14 -- a withdrawal is a claim, and gets the same scrutiny.
  P9  If the scatter SURVIVES at Omega = 700, §46 stands, (*) is missing something real on
      the rho axis, and the flatness of the prediction is the thing to explain.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
from scipy.optimize import brentq

from crnl.networks.am_reversible import reverse_pairing
from experiments.arrhenius_optimum import am_rho, cost_cell, delta_star_rho
from experiments.cost_absolute import sigma_and_mu
from experiments.slaving_axis import sep_of


def field(net, delta, s):
    """(mu, nu) = (d delta/dt, d s/dt) at the given (delta, s)."""
    x, y, b = 0.5 * (s + delta), 0.5 * (s - delta), 1.0 - s
    if min(x, y, b) < 0:
        return None
    v = net.stoichiometry_matrix() @ net.fluxes(np.array([x, y, b]))
    return float(v[0] - v[1]), float(v[0] + v[1])


def s_star(net, delta):
    """Pool value on the slaved manifold: nu(delta, s*) = 0."""
    def f(s):
        r = field(net, delta, s)
        return np.nan if r is None else r[1]
    lo, hi = abs(delta) + 1e-9, 1.0 - 1e-9
    a, b = f(lo), f(hi)
    if not (np.isfinite(a) and np.isfinite(b)) or a * b > 0:
        return None
    return float(brentq(f, lo, hi, xtol=1e-14))


def eps_at(net, delta, h):
    """The dimensionless lag correction (*), by central differences with step h."""
    s0 = s_star(net, delta)
    if s0 is None:
        return None
    sp, sm = s_star(net, delta + h), s_star(net, delta - h)
    if sp is None or sm is None:
        return None
    ds_dd = (sp - sm) / (2 * h)

    fp, fm = field(net, delta, s0 + h), field(net, delta, s0 - h)
    if fp is None or fm is None:
        return None
    dmu_ds = (fp[0] - fm[0]) / (2 * h)
    dnu_ds = (fp[1] - fm[1]) / (2 * h)
    if abs(dnu_ds) < 1e-12:
        return None
    return float(dmu_ds * ds_dd / dnu_ds)


def predict(net, d_lo, d_hi, pairing, h, n=241):
    """<eps>_time, weighted by dt = ddelta/mu along the traversal."""
    xs = np.linspace(d_lo, d_hi, n)
    e, w = [], []
    for x in xs:
        ev = eps_at(net, float(x), h)
        sm = sigma_and_mu(net, float(x), pairing)
        if ev is None or sm is None or sm[1] <= 0:
            return None
        e.append(ev)
        w.append(1.0 / sm[1])
    e, w = np.array(e), np.array(w)
    return float(np.trapezoid(e * w, xs) / np.trapezoid(w, xs))


def cell(gamma, rho, omega, eps_frac, theta, h):
    ds = delta_star_rho(gamma, rho)
    if ds <= 0:
        return None
    net = am_rho(gamma, rho)
    pairing = reverse_pairing(net)
    d_lo, d_hi = eps_frac * ds, theta * ds
    pred = predict(net, d_lo, d_hi, pairing, h)
    if pred is None:
        return None
    c = cost_cell(gamma, rho, omega, eps_frac, theta)
    if c is None or not np.isfinite(c.get("mean_time", np.nan)):
        return None
    xs = np.linspace(d_lo, d_hi, 601)
    inv_mu = []
    for x in xs:
        sm = sigma_and_mu(net, float(x), pairing)
        if sm is None or sm[1] <= 0:
            return None
        inv_mu.append(1.0 / sm[1])
    t_det = float(np.trapezoid(inv_mu, xs))
    gap = t_det / c["mean_time"] - 1.0
    return {"gamma": gamma, "rho": rho, "omega": omega, "sep": float(sep_of(net)[0]),
            "gap": float(gap), "pred": pred,
            "pred_over_gap": float(pred / gap) if gap != 0 else float("nan"),
            "t_det": t_det, "mfpt": float(c["mean_time"])}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rhos", type=float, nargs="+",
                    default=[0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 16.0, 32.0])
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.07, 0.12, 0.20, 0.28, 0.35])
    ap.add_argument("--gamma", type=float, default=0.20)
    ap.add_argument("--omega", type=int, default=300)
    ap.add_argument("--eps", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--h", type=float, default=1e-4)
    ap.add_argument("--conv-omegas", type=int, nargs="+", default=[300, 500, 700])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/lag_absolute.json"))
    args = ap.parse_args()

    t0 = time.time()
    eps_f, th = args.eps, args.theta

    print("=== P1 GATE (rule 13): is eps converged in its own step h?")
    net = am_rho(args.gamma, 1.0)
    pairing = reverse_pairing(net)
    ds = delta_star_rho(args.gamma, 1.0)
    print(f"{'h':>10}{'<eps>_time':>14}{'rel change':>13}")
    prev = None
    for h in (4e-4, 2e-4, 1e-4, 5e-5, 2.5e-5):
        p = predict(net, eps_f * ds, th * ds, pairing, h)
        rc = abs(p - prev) / abs(p) if prev is not None else float("nan")
        print(f"{h:>10.1e}{p:>14.6f}{rc:>13.2e}")
        prev = p
    conv = rc < 0.01
    print(f"  -> P1 {'HOLDS' if conv else 'FAILS'} (last relative change {rc:.2e})")

    print(f"\n=== P2/P3/P4: absolute prediction vs measurement -- rho axis"
          f" (where 1/sep FAILED, 72.1%)")
    print(f"{'rho':>7}{'sep':>9}{'gap':>10}{'<eps>':>10}{'pred/gap':>11}{'gap*sep':>10}")
    rrows = []
    for r in args.rhos:
        c = cell(args.gamma, r, args.omega, eps_f, th, args.h)
        if c is None:
            print(f"{r:>7.2f}   dropped")
            continue
        rrows.append(c)
        print(f"{r:>7.2f}{c['sep']:>9.3f}{c['gap']:>10.4f}{c['pred']:>10.4f}"
              f"{c['pred_over_gap']:>11.3f}{c['gap']*c['sep']:>10.4f}")

    print(f"\n=== gamma axis (where 1/sep HELD, 22.2%)")
    print(f"{'gamma':>7}{'sep':>9}{'gap':>10}{'<eps>':>10}{'pred/gap':>11}{'eps*sep':>10}")
    grows = []
    for g in args.gammas:
        c = cell(g, 1.0, args.omega, eps_f, th, args.h)
        if c is None:
            print(f"{g:>7.2f}   dropped")
            continue
        grows.append(c)
        print(f"{g:>7.2f}{c['sep']:>9.3f}{c['gap']:>10.4f}{c['pred']:>10.4f}"
              f"{c['pred_over_gap']:>11.3f}{c['pred']*c['sep']:>10.4f}")

    print(f"\n=== P2: is the absolute prediction right?")
    for name, rows in (("rho", rrows), ("gamma", grows)):
        v = np.array([c["pred_over_gap"] for c in rows if np.isfinite(c["pred_over_gap"])])
        if v.size:
            print(f"  {name:>6} axis: pred/gap mean {v.mean():.3f},"
                  f" range {v.min():.3f}..{v.max():.3f},"
                  f" spread {100*(v.max()-v.min())/abs(v.mean()):.1f}%")

    print(f"\n=== P3: the discriminating comparison -- constancy of pred/gap vs 1/sep's")
    print(f"  1/sep gave: T axis 12.3%, gamma axis 22.2%, rho axis 72.1% (§46)")
    for name, rows in (("rho", rrows), ("gamma", grows)):
        v = np.array([c["pred_over_gap"] for c in rows if np.isfinite(c["pred_over_gap"])])
        if v.size > 2:
            sp = 100 * (v.max() - v.min()) / abs(v.mean())
            base = 72.1 if name == "rho" else 22.2
            print(f"  {name:>6} axis: pred/gap spread {sp:.1f}% against 1/sep's {base}%"
                  f"  -> {'BETTER' if sp < base else 'no better'}")

    print(f"\n=== P4: does it reproduce the U-shape in rho, cell by cell?")
    if rrows:
        gs = np.array([c["gap"] * c["sep"] for c in rrows])
        ps = np.array([c["pred"] * c["sep"] for c in rrows])
        i, j = int(np.argmin(gs)), int(np.argmin(ps))
        print(f"  measured gap*sep minimum at rho = {rrows[i]['rho']}"
              f"   predicted <eps>*sep minimum at rho = {rrows[j]['rho']}"
              f"  -> {'SAME cell' if i == j else 'different cells'}")
        print(f"  measured  gap*sep: " + " ".join(f"{v:.3f}" for v in gs))
        print(f"  predicted eps*sep: " + " ".join(f"{v:.3f}" for v in ps))

    print(f"\n=== P5: does it EXPLAIN §46's fitted gamma-axis constant 0.5963?")
    if grows:
        v = np.array([c["pred"] * c["sep"] for c in grows])
        print(f"  <eps>*sep on the gamma axis: {v.min():.3f}..{v.max():.3f},"
              f" mean {v.mean():.3f}   against the fitted 0.5963")

    print(f"\n=== P6: verdict")
    allv = np.array([c["pred_over_gap"] for c in rrows + grows
                     if np.isfinite(c["pred_over_gap"])])
    if allv.size:
        worst = max(abs(allv.max() - 1), abs(1 - allv.min()))
        if worst < 0.25:
            print(f"  pred/gap within {100*worst:.0f}% of 1 everywhere -> the lag picture")
            print(f"  is CONFIRMED in absolute terms, and §39.2's fitted C is explained.")
        elif allv.mean() > 0.5 and allv.mean() < 2.0:
            print(f"  pred/gap averages {allv.mean():.2f} but spans {allv.min():.2f}"
                  f"..{allv.max():.2f} -- right scale, wrong detail.")
        else:
            print(f"  pred/gap = {allv.mean():.2f} on average -> the lag picture is")
            print(f"  REFUTED in absolute terms. §39.2's mechanism needs qualifying.")

    # ===== SECOND PASS: rule 13 on the MEASUREMENT, not just the prediction ==========
    print("\n" + "=" * 78)
    print("P7/P8/P9: is the measured gap converged in Omega? (§46 never checked)")
    conv_rhos = [0.5, 1.5, 3.0, 6.0, 32.0]
    omegas = args.conv_omegas
    print(f"{'rho':>7}{'sep':>9}" + "".join(f"{f'gap*sep(Om={o})':>17}" for o in omegas)
          + f"{'predicted':>11}")
    crows = []
    for r in conv_rhos:
        ds = delta_star_rho(args.gamma, r)
        net = am_rho(args.gamma, r)
        pairing = reverse_pairing(net)
        sep = float(sep_of(net)[0])
        pred = predict(net, eps_f * ds, th * ds, pairing, args.h)
        line = f"{r:>7.2f}{sep:>9.3f}"
        for om in omegas:
            c = cell(args.gamma, r, om, eps_f, th, args.h)
            if c is None:
                line += f"{'--':>17}"
                continue
            crows.append(c)
            line += f"{c['gap']*c['sep']:>17.4f}"
        print(line + f"{pred*sep:>11.4f}")

    print(f"\n=== P7: does the scatter tighten with Omega?")
    for om in omegas:
        v = np.array([c["gap"] * c["sep"] for c in crows if c["omega"] == om])
        if v.size > 2:
            print(f"  Omega={om:>4}: gap*sep spans {v.min():.4f}..{v.max():.4f}"
                  f"   spread {100*(v.max()-v.min())/v.mean():>5.1f}%"
                  f"   mean {v.mean():.4f}")
    pv = np.array([c["pred"] * c["sep"] for c in crows if c["omega"] == omegas[-1]])
    if pv.size:
        print(f"  predicted: spans {pv.min():.4f}..{pv.max():.4f}"
              f"   spread {100*(pv.max()-pv.min())/pv.mean():>5.1f}%   mean {pv.mean():.4f}")

    print(f"\n=== P8/P9: verdict on §46's rho-axis claim")
    hi = np.array([c["gap"] * c["sep"] for c in crows if c["omega"] == omegas[-1]])
    lo = np.array([c["gap"] * c["sep"] for c in crows if c["omega"] == omegas[0]])
    if hi.size > 2 and lo.size > 2:
        s_hi = 100 * (hi.max() - hi.min()) / hi.mean()
        s_lo = 100 * (lo.max() - lo.min()) / lo.mean()
        print(f"  spread at Omega={omegas[0]}: {s_lo:.1f}%   at Omega={omegas[-1]}:"
              f" {s_hi:.1f}%")
        if s_hi < 0.6 * s_lo and s_hi < 35:
            print("  -> P8: the scatter DISSOLVES with Omega. §46's 'the 1/sep scaling")
            print("     fails on rho' is WITHDRAWN -- it was measured at insufficient")
            print("     Omega. The absolute prediction predicted a measurement error.")
        else:
            print("  -> P9: the scatter SURVIVES. §46 stands and (*) is missing something")
            print("     real on the rho axis; the flat prediction is what needs explaining.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"rho_rows": rrows, "gamma_rows": grows,
                                    "conv_rows": crows,
                                    "h": args.h, "omega": args.omega},
                                   indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
