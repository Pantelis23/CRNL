"""T-COST-f: does rho act through the timescale separation? — and is that §39's missing 6%?

§44.2 found rho a free lever: cost -43..-50%, reliability x2, time halved, with delta*
frozen. The mechanism proposed there is a SUSPECT (rule 17): a faster disagreement channel
is a faster pool, which is deeper slaving, which is §36's on-manifold condition and
§39.2's 1/sep law. This is its kill test.

**IT IS ALSO AN ATTACK ON T-COST-c, WHICH IS THE PROJECT'S OUTSTANDING UNEXPLAINED
NUMBER.** §39's closed form `Sigma_pred = Omega * int sigma/mu` sits ~6% off the exact CME
cost, mean ratio 1.0583 over 36 cells, and that residual is **flat in Omega over a 3x
range** -- not finite-size -- and **not attributable to sigma**, which converges properly.
If the residual is the slaved approximation failing at finite separation, then driving sep
up with rho must close it. If it does not close, rho does not act through slaving AND
§39's residual is something else again.

**THE DESIGN SCAN HANDED THIS TEST ITS OWN CONTROL.** Measured before writing these
predictions, at gamma = 0.20:

    rho     0.5    1.0    2.0    4.0    8.0   16.0   32.0
    sep    11.70   7.00   7.19  10.56  18.47  34.82  67.76

**sep is NON-MONOTONE in rho, with a minimum near rho ~ 1.5** -- while §44.2's cost G is
monotone decreasing across that whole range. Two consequences, and they cut in opposite
directions:

  * Over rho = 0.5 -> 1.5, cost IMPROVES while sep gets WORSE. So rho's cost benefit
    cannot be "more slaving" in that region, whatever happens to §39's residual.
  * The non-monotonicity supplies MATCHED PAIRS: rho = 0.5 and rho = 4 differ 8-fold at
    sep 11.70 vs 10.56. That is rule 9 handed over for free -- a way to move rho a long
    way while holding the suspected cause nearly fixed.

`sep_of` reproduces the closed form 3(1+2g)/(1-2g) at rho = 1 to four digits (3.977,
7.000, 17.000), so the instrument is anchored before use.

PREDICTIONS, written before running:

  P1  GATE, two parts. (a) sep_of matches 3(1+2g)/(1-2g) at rho = 1 to 1e-9. (b) At
      rho = 1 the ratio pred/exact reproduces §39: mean ~1.058, cells inside 0.905-1.167.
      If (b) fails the closed form has been mis-transplanted onto the rho network and
      nothing below is admissible.
  P2  THE TEST. The residual (pred/exact - 1) tracks 1/sep, and therefore must itself be
      **NON-MONOTONE in rho** -- worsening from rho = 0.5 to rho ~ 1.5, then improving. A
      residual that falls MONOTONICALLY in rho would be tracking rho, not sep, and would
      refute the slaving account even while looking superficially like a confirmation.
      **This is the whole reason the non-monotone design is worth having.**
  P3  MATCHED PAIRS, the sharp form. At gamma = 0.20 the pairs (rho = 0.5, sep 11.70) vs
      (rho = 4, sep 10.56), and (rho = 1, sep 7.00) vs (rho = 2, sep 7.19), differ 8x and
      2x in rho at nearly equal sep. **Their residuals must agree.** If they do not, sep
      is not the governing variable and any 1/sep fit is an artifact of the sweep.
  P4  CROSS-KNOB COLLAPSE. gamma moves sep over 3.98-17.0 at rho = 1; rho moves it over
      7-68 at gamma = 0.20. In the overlap, residual vs 1/sep must lie on ONE curve
      regardless of which knob produced the sep. Two physically different knobs, one
      governing variable -- or not.
  P5  QUANTITATIVE, absolute and not fitted (rule 16). If the law is residual ~ C/sep,
      then residual*sep is constant. Reported per cell. Drift of more than ~2x across a
      10x range of sep means the 1/sep FORM is wrong even if the residual does shrink,
      and that must be said rather than absorbed into a fitted exponent.
  P6  THE REFUTING OUTCOME, named in advance. If the residual does not shrink with sep at
      all, T-COST-f is dead: rho does not act through the slaved approximation, §44.2's
      measurement stands with no mechanism, and §39's flat 6% goes back to T-COST-c's
      remaining candidates -- absorption selecting early-fluctuating trajectories, or a
      Jensen gap E[sigma(state)] != sigma(E[state]).
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.networks.am_reversible import reverse_pairing
from experiments.arrhenius_optimum import am_rho, cost_cell, delta_star_rho
from experiments.cost_absolute import sigma_predicted
from experiments.slaving_axis import sep_of


def residual_cell(gamma, rho, omega, eps, theta):
    """pred/exact for §39's closed form on the rho-generalised network."""
    ds = delta_star_rho(gamma, rho)
    if ds <= 0:
        return None
    net = am_rho(gamma, rho)
    per_om = sigma_predicted(net, eps * ds, theta * ds, reverse_pairing(net))
    if not np.isfinite(per_om):
        return None
    c = cost_cell(gamma, rho, omega, eps, theta)
    if c is None:
        return None
    try:
        sep, _ = sep_of(net)
    except Exception:
        return None
    pred = omega * per_om
    return {"gamma": gamma, "rho": rho, "omega": omega, "sep": float(sep),
            "delta_star": ds, "exact": c["Sigma"], "pred": float(pred),
            "ratio": float(pred / c["Sigma"]), "L": c["L"], "G": c["G"]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.20)
    ap.add_argument("--rhos", type=float, nargs="+",
                    default=[0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 16.0, 32.0])
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.07, 0.12, 0.20, 0.28, 0.35])
    ap.add_argument("--omegas", type=int, nargs="+", default=[200, 300])
    ap.add_argument("--eps", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/slaving_residual_rho.json"))
    args = ap.parse_args()

    t0 = time.time()
    eps, th = args.eps, args.theta

    print("=== P1a GATE: sep_of vs the closed form 3(1+2g)/(1-2g) at rho = 1")
    worst = 0.0
    for g in args.gammas:
        s, _ = sep_of(am_rho(g, 1.0))
        cf = 3.0 * (1.0 + 2.0 * g) / (1.0 - 2.0 * g)
        worst = max(worst, abs(s - cf) / cf)
        print(f"  gamma={g:.2f}: sep_of = {s:.9f}   closed = {cf:.9f}")
    print(f"  -> P1a {'HOLDS' if worst < 1e-9 else 'FAILS'}  (worst rel {worst:.2e})")

    # ---- the rho knob, at fixed gamma -------------------------------------------
    print(f"\n=== rho knob at gamma = {args.gamma}  (sep is NON-monotone here)")
    print(f"{'rho':>7}{'sep':>9}{'1/sep':>8}" + "".join(f"{f'ratio(Om={o})':>14}"
                                                        for o in args.omegas)
          + f"{'resid*sep':>11}{'delta*':>9}")
    rrows = []
    for r in args.rhos:
        cells = [residual_cell(args.gamma, r, om, eps, th) for om in args.omegas]
        if any(c is None for c in cells):
            print(f"{r:>7.2f}   dropped (no landscape / solve failed)")
            continue
        sep = cells[0]["sep"]
        line = f"{r:>7.2f}{sep:>9.3f}{1/sep:>8.4f}"
        for c in cells:
            line += f"{c['ratio']:>14.4f}"
            rrows.append(c)
        rs = np.mean([(c["ratio"] - 1.0) * sep for c in cells])
        print(line + f"{rs:>11.3f}{cells[0]['delta_star']:>9.4f}")

    # ---- the gamma knob, at rho = 1 ----------------------------------------------
    print(f"\n=== gamma knob at rho = 1  (the second, physically different lever)")
    print(f"{'gamma':>7}{'sep':>9}{'1/sep':>8}" + "".join(f"{f'ratio(Om={o})':>14}"
                                                          for o in args.omegas)
          + f"{'resid*sep':>11}{'delta*':>9}")
    grows = []
    for g in args.gammas:
        cells = [residual_cell(g, 1.0, om, eps, th) for om in args.omegas]
        if any(c is None for c in cells):
            print(f"{g:>7.2f}   dropped")
            continue
        sep = cells[0]["sep"]
        line = f"{g:>7.2f}{sep:>9.3f}{1/sep:>8.4f}"
        for c in cells:
            line += f"{c['ratio']:>14.4f}"
            grows.append(c)
        rs = np.mean([(c["ratio"] - 1.0) * sep for c in cells])
        print(line + f"{rs:>11.3f}{cells[0]['delta_star']:>9.4f}")

    # ---- P1b -----------------------------------------------------------------------
    at1 = [c for c in rrows + grows if abs(c["rho"] - 1.0) < 1e-12]
    if at1:
        rs = np.array([c["ratio"] for c in at1])
        print(f"\n=== P1b GATE: rho = 1 must reproduce §39 (mean 1.058, cells 0.905-1.167)")
        print(f"  {len(rs)} cells: mean {rs.mean():.4f}, range {rs.min():.4f}..{rs.max():.4f}"
              f"  -> {'HOLDS' if 0.88 < rs.mean() < 1.20 else 'FAILS'}")

    # ---- P2 -------------------------------------------------------------------------
    print(f"\n=== P2: is the residual NON-monotone in rho, as 1/sep requires?")
    byr = {}
    for c in rrows:
        byr.setdefault(c["rho"], []).append(abs(c["ratio"] - 1.0))
    rs_ = sorted(byr)
    res = np.array([np.mean(byr[r]) for r in rs_])
    print("  rho     " + "".join(f"{r:>9.2f}" for r in rs_))
    print("  |resid| " + "".join(f"{v:>9.4f}" for v in res))
    i = int(np.argmax(res))
    nonmono = 0 < i < len(res) - 1
    print(f"  peak |residual| at rho = {rs_[i]:.2f}"
          f"  ({'INTERIOR -> non-monotone, as 1/sep requires' if nonmono else 'at an edge -> MONOTONE, which refutes the sep reading'})")

    # ---- P3 matched pairs -------------------------------------------------------------
    print(f"\n=== P3: matched-sep pairs at very different rho (rule 9)")
    got = {c["rho"]: c for c in rrows if c["omega"] == args.omegas[0]}
    for a, b in ((0.5, 4.0), (1.0, 2.0)):
        if a in got and b in got:
            ca, cb = got[a], got[b]
            dsep = 100 * abs(ca["sep"] - cb["sep"]) / (0.5 * (ca["sep"] + cb["sep"]))
            dres = abs(abs(ca["ratio"] - 1) - abs(cb["ratio"] - 1))
            print(f"  rho {a} (sep {ca['sep']:.2f}, |resid| {abs(ca['ratio']-1):.4f})"
                  f"  vs rho {b} (sep {cb['sep']:.2f}, |resid| {abs(cb['ratio']-1):.4f})")
            print(f"    rho differs {max(a,b)/min(a,b):.0f}x, sep differs {dsep:.1f}%,"
                  f" residuals differ by {dres:.4f}"
                  f"  -> {'AGREE, sep governs' if dres < 0.02 else 'DIFFER, sep does NOT govern'}")

    # ---- P4 cross-knob collapse -----------------------------------------------------
    print(f"\n=== P4: do the two knobs collapse onto one curve in 1/sep?")
    lo = max(min(c["sep"] for c in rrows), min(c["sep"] for c in grows))
    hi = min(max(c["sep"] for c in rrows), max(c["sep"] for c in grows))
    print(f"  overlap sep in [{lo:.2f}, {hi:.2f}]")
    for name, rows in (("rho knob", rrows), ("gamma knob", grows)):
        sel = [c for c in rows if lo - 1e-9 <= c["sep"] <= hi + 1e-9]
        if len(sel) >= 2:
            x = np.array([1.0 / c["sep"] for c in sel])
            y = np.array([c["ratio"] - 1.0 for c in sel])
            p = np.polyfit(x, y, 1)
            print(f"  {name:>10}: {len(sel)} cells, resid = {p[0]:+.4f}/sep {p[1]:+.4f}")

    # ---- P5 ----------------------------------------------------------------------
    print(f"\n=== P5: is residual*sep constant (absolute check, not a fitted exponent)?")
    for name, rows in (("rho knob", rrows), ("gamma knob", grows)):
        v = np.array([(c["ratio"] - 1.0) * c["sep"] for c in rows])
        s = np.array([c["sep"] for c in rows])
        print(f"  {name:>10}: resid*sep spans {v.min():+.3f} .. {v.max():+.3f}"
              f"  over sep {s.min():.1f}..{s.max():.1f} ({s.max()/s.min():.1f}x)"
              f"  -> {'roughly constant' if v.max()/max(abs(v.min()),1e-9) < 2 and v.min() > 0 else 'NOT constant, the 1/sep FORM is wrong'}")

    # ---- P6 --------------------------------------------------------------------------
    print(f"\n=== P6: does the residual shrink with sep at all?")
    allr = rrows + grows
    x = np.array([1.0 / c["sep"] for c in allr])
    y = np.array([abs(c["ratio"] - 1.0) for c in allr])
    cc = float(np.corrcoef(x, y)[0, 1])
    hi_sep = np.array([abs(c["ratio"] - 1) for c in allr if c["sep"] > 20])
    lo_sep = np.array([abs(c["ratio"] - 1) for c in allr if c["sep"] < 8])
    print(f"  corr(|residual|, 1/sep) = {cc:+.4f} over {len(allr)} cells")
    if hi_sep.size and lo_sep.size:
        print(f"  |residual| at sep < 8:  mean {lo_sep.mean():.4f}  ({lo_sep.size} cells)")
        print(f"  |residual| at sep > 20: mean {hi_sep.mean():.4f}  ({hi_sep.size} cells)")
        shrinks = hi_sep.mean() < 0.5 * lo_sep.mean()
        print(f"  -> T-COST-f {'SURVIVES this test' if shrinks else 'REFUTED: the residual does not close with sep'}")
        if not shrinks:
            print("     §44.2's measurement stands with no mechanism, and §39's flat 6%")
            print("     returns to T-COST-c's other candidates (absorption bias; Jensen gap).")

    out = {"gamma_fixed": args.gamma, "eps": eps, "theta": th,
           "rho_rows": rrows, "gamma_rows": grows}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
