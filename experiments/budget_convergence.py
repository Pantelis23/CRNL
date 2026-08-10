"""Does §50's ~21% survive more cells and larger Omega? — verifying before theorising

§50 concluded that a complete 1-D account supplies ~21% of the absorption coefficient,
leaving ~79% two-dimensional. **That rests on two cells selected POST-HOC** -- the criterion
was "last two Omega agree within 1%", chosen after seeing which cells qualified -- and one
of the two, gamma = 0.07, has a contaminated top point: its cme_coeff reads 2.61, 2.49, 2.60
then jumps to 3.65 at Omega = 1000, which §49 attributed to a numerical floor.

So the 21% is really one clean cell plus one suspect one, and §51 has since eliminated two
candidate mechanisms for the remaining 79% without touching the measurement itself. **Rule 8
says verify before building further, and two mechanism proposals have already failed on this
quantity.** This adds Omega = 1400 and three new cells rather than a third mechanism.

PREDICTIONS, written before running:

  P1  At Omega = 1400 more cells reach the convergence criterion, and it is applied
      UNCHANGED from §50 (last two Omega within 1%) rather than retuned.
  P2  THE TEST. bd/cme stays near 0.21 in the newly converged cells. **If the newly
      converged cells give something materially different, §50's 21% was a two-cell
      accident** and the "79% is two-dimensional" claim loses its basis.
  P3  gamma = 0.07's Omega = 1000 outlier resolves one way or the other: either it rejoins
      the 2.5-2.6 trend, confirming §49's numerical-floor diagnosis, or the trend genuinely
      rises and §49's diagnosis was wrong.
  P4  The three new cells are chosen to span the axes already swept -- a new gamma, a new
      rho, and a cell with both moved -- so that any drift in bd/cme is attributable.
  P5  If bd/cme turns out to VARY systematically across cells rather than clustering, then
      "21%" is not a constant of the chemistry and §50's headline should be restated as a
      per-cell quantity. That is a real outcome and is reported as one.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from experiments.birthdeath_absorption import cell as bd_cell
from experiments.lag_endpoints import cell as cme_cell


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cells", type=str, nargs="+",
                    default=["0.07:1", "0.20:1", "0.35:1", "0.20:0.5", "0.20:32",
                             "0.12:1", "0.20:4", "0.28:8"])
    ap.add_argument("--omegas", type=int, nargs="+", default=[300, 500, 700, 1000, 1400])
    ap.add_argument("--eps", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--h", type=float, default=1e-4)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/budget_convergence.json"))
    args = ap.parse_args()

    t0 = time.time()
    cells = [(float(c.split(":")[0]), float(c.split(":")[1])) for c in args.cells]

    print("cme_coeff = (gap - <eps>)*Omega, the 2-D absorption coefficient")
    print(f"{'gamma':>6}{'rho':>6}" + "".join(f"{f'Om={o}':>10}" for o in args.omegas)
          + f"{'last-two':>10}{'conv?':>7}")
    rows, cme = [], {}
    for g, r in cells:
        line, vals = f"{g:>6.2f}{r:>6.1f}", []
        for om in args.omegas:
            try:
                c = cme_cell(g, r, om, args.eps, args.theta, args.h)
            except Exception:
                c = None
            if c is None:
                line += f"{'--':>10}"
                continue
            v = (c["gap_real"] - c["pred_real"]) * om
            vals.append(v)
            rows.append({"gamma": g, "rho": r, "omega": om, "cme_coeff": v,
                         "ratio_real": c["ratio_real"]})
            line += f"{v:>10.3f}"
        if len(vals) >= 2:
            lt = abs(vals[-1] - vals[-2]) / abs(vals[-1])
            cme[(g, r)] = (vals[-1], lt < 0.01, vals)
            line += f"{100*lt:>9.1f}%{'YES' if lt < 0.01 else 'no':>7}"
        print(line)

    print(f"\n=== the 1-D chain, and the budget")
    print(f"{'gamma':>6}{'rho':>6}{'bd':>9}{'cme':>9}{'bd/cme':>9}{'cme conv?':>11}")
    out = []
    for g, r in cells:
        if (g, r) not in cme:
            continue
        b = bd_cell(g, r, 2000, args.eps, args.theta)
        if b is None:
            continue
        cv, conv, _ = cme[(g, r)]
        out.append({"gamma": g, "rho": r, "bd": b["bd_coeff"], "cme": cv,
                    "bd_over_cme": b["bd_coeff"] / cv, "cme_converged": bool(conv)})
        print(f"{g:>6.2f}{r:>6.1f}{b['bd_coeff']:>9.3f}{cv:>9.3f}"
              f"{b['bd_coeff']/cv:>9.3f}{'YES' if conv else 'no':>11}")

    conv = [c for c in out if c["cme_converged"]]
    print(f"\n=== P1/P2: converged cells (criterion unchanged from §50)")
    print(f"  {len(conv)} of {len(out)} cells converged, against 2 of 5 in §50")
    if conv:
        v = np.array([c["bd_over_cme"] for c in conv])
        for c in conv:
            print(f"  gamma={c['gamma']:.2f} rho={c['rho']:<5} bd/cme = {c['bd_over_cme']:.3f}")
        print(f"  bd/cme {v.min():.3f}..{v.max():.3f}, mean {v.mean():.3f},"
              f" spread {100*(v.max()-v.min())/v.mean():.1f}%")
        near = abs(v.mean() - 0.215) < 0.04
        print(f"  -> P2 {'HOLDS, §50 21% survives' if near else 'FAILS -- §50 was a two-cell accident'}")
        if v.size > 2 and 100 * (v.max() - v.min()) / v.mean() > 30:
            print("  -> P5: bd/cme VARIES across cells; 21% is not a constant of the")
            print("     chemistry and §50's headline must be restated per cell.")

    print(f"\n=== P3: the gamma = 0.07 outlier")
    if (0.07, 1.0) in cme:
        _, _, vals = cme[(0.07, 1.0)]
        print(f"  series: " + " ".join(f"{v:.2f}" for v in vals))
        if len(vals) >= 5:
            print(f"  Omega=1000 was 3.65 against a 2.5-2.6 trend;"
                  f" Omega=1400 gives {vals[-1]:.2f}"
                  f"  -> {'rejoins the trend, §49 was right' if vals[-1] < 3.2 else 'the trend genuinely rises, §49 was wrong'}")

    print(f"\n=== unconverged cells, reported not hidden")
    for c in out:
        if not c["cme_converged"]:
            print(f"  gamma={c['gamma']:.2f} rho={c['rho']:<5} bd/cme = {c['bd_over_cme']:.3f}"
                  f"  (cme still moving)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"sweep": rows, "budget": out}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
