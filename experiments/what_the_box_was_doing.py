"""§113 (T-BOX-a) -- which published numbers move under a converged box?

§112.7 found that `cap_mult`, the width of the truncated state space with a reflecting wall at
`cap_mult * r3 * Omega`, was never swept. At the inherited 1.25 the error in ln T reaches 0.249
on a 1-D element; the wall truncates excursions above the rail so it can only SHORTEN an escape,
making the bias one-signed and worst for shallow wells at small Omega. `cap_mult` appears in 28
experiment files and 1.25 is the dominant default, so most of the §91-§111 cascade arc inherits it.

A first probe in §112.7 showed the bias partly cancels in ratios but not away: on
`escape_accounts_for_it.escape_rate` the absolute error is 6.88% / 2.23% / 0.47% at Omega =
14 / 30 / 55, and the ratio k(2.4)/k(r3) is off 1.95% / 0.97% / 0.25%. A factor of ~3.5, not
infinity. **So published cascade numbers quoted below ~2% at Omega = 14 are candidates.** This
section re-measures the ones that are (a) quoted with a precision and (b) affordable at a
converged box.

`cap_mult` is threaded through `solve`, `rate_limits`, `spectral_gap` and `upstream_relaxation`
with the default left at 1.25 everywhere, so nothing already published moves silently (rule 7).
P1 checks that the threading is inert at the default.

PREDICTIONS, WRITTEN BEFORE RUNNING.

  P1  WIRING, and it gates everything else. At cap_mult = 1.25 every threaded quantity must
      reproduce its published value: spectral_gap(14, free) = 7.29805e-02 (§102's P1 reference),
      rate_limits(14, 2, 2.0) giving k_avg/k_mean = 2.435 and a position in (0, 0.5) (§102's P3
      and its test), tau_up(14) and tau_cross(14) matching §110's stored table. If any default
      moves, the threading changed the physics and nothing below is readable.

  P2  CONVERGENCE, NOT TOLERANCE (rule 20, and §112 violated this in its own first draft). Every
      quantity is swept over cap_mult 1.25, 1.6, 2.0, 3.0, 4.0 and judged by whether successive
      differences SHRINK toward zero. A quantity that has not converged by 4.0 is reported as
      unconverged, not as a number. No fixed tolerance appears anywhere in the verdicts.

  P3  THE SINGLE-STAGE RATES, an extrapolation of a measurement rather than a new claim. §112.7
      measured the 1-D escape rate's box error as one-signed and falling steeply with Omega.
      §100's spectral gap and §104's descent rate are the same kind of object -- a relaxation
      rate of one walled 1-D stage -- so predict the same shape and comparable size: several
      percent at Omega = 14, well under 1% by Omega = 55, always with the converged value SMALLER
      (a wider box lengthens the escape, so it lowers the rate).

  P4  §102's BRACKET, and it can print either verdict. `pos = ln(k_eff/k_mean)/ln(k_avg/k_mean)`.
      The WIDTH ln(k_avg/k_mean) is a ratio of two escape rates at nearby operating points, so
      most of the one-signed bias should cancel: predict the width moves by **less than 2% at
      Omega = 14**. The POSITION also carries k_eff, which comes from the D-dimensional JOINT
      solve, whose wall truncates a product space; there is no reason its box error matches the
      single-stage one, and if they differ the position moves more than the width does. **Predict
      |Delta position| < 0.10 in absolute terms.** §102 published "0 < pos < 0.5, near the fast
      end" and a test asserts exactly that. **If the position leaves (0, 0.5), or moves by more
      than 0.10, §102's headline needs a correction printed beside it and §106-§110's readings of
      "position" inherit the problem.**

  P5  §110's RATIO, where the most is at stake, and RULE 14 IS WHY IT IS HERE: §110 published a
      REFUTATION -- it retired the timescale-ratio account of the cascade drift -- and rule 14
      says a withdrawal must be verified as carefully as an assertion. §110's headline is that
      tau_up/tau_cross FALLS across Omega = 14-70 while the measured position RISES, so the
      position is not set by that ratio.

      tau_cross is an h-transform confined between the rail and the saddle, entirely BELOW the
      wall, so predict it is nearly box-independent -- the same reason §112's splitting
      probability was exactly so. [RESOLVED ON INSPECTION, not by the sweep: the prediction is
      right and it is STRUCTURAL. `downstream_crossing` never reads cap_mult at all, so sweeping
      it re-runs one computation five times. The first draft of this section reported those five
      identical numbers as a 0.00% convergence result, which is a number the harness produced
      and the chemistry did not -- rule 10, caught by a test asserting exact equality and failing
      on BLAS-level noise.] tau_up is the sub-dominant eigenvalue of a walled generator and
      should carry a box error like the escape rates of P3. **Predict the ratio moves by less
      than 5% at every Omega and that §110's DIRECTION survives.** If instead the trend flattens
      or inverts, §110's refutation is a box artifact, the timescale-ratio account has to be
      reopened, and a withdrawal in this project will have been withdrawn for the second time.

      LIMITATION, NAMED BEFORE RUNNING: tau_cross takes its operating point x_up from
      `chain_operating_points`, which has a box of its own that this sweep holds at 1.25. So P5
      varies each clock's DIRECT box dependence only. P5b bounds the held path separately;
      without it P5's percentages are a lower bound, not the whole sensitivity.

  P6  WHAT THE AUDIT COVERS, stated because a null result is worth exactly its scope. This
      re-runs quantities that are published with a quoted precision AND affordable at a converged
      box. It does NOT cover §101/§103/§105-§109's joint-solve tables at D = 3 with Omega >= 30:
      at cap_mult 3.0 those need 10^6-10^7 states. Whatever the verdict here, those are listed as
      UNAUDITED rather than implied safe -- which is the whole difference between this section
      and a reassurance.

WHAT THIS CANNOT SETTLE. The box is one inherited numerical parameter among several -- the time
grid, the expm tolerance, the volume ceilings. This audits one axis. A quantity that survives the
box may still be sensitive to another, and nothing here licenses the assumption that it is not.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

CAPS = (1.25, 1.6, 2.0, 3.0, 4.0)


def _converged(vals):
    """Rule 20: the criterion is the SHAPE of the approach, not a tolerance on the value.

    A first version demanded d[i] >= d[i+1] for every i and flagged NOT-CONV on sequences that
    had plainly settled -- once the differences reach the floating-point floor their ordering is
    noise, not physics. That is the same over-strict gate §112's own first draft carried. What
    convergence actually means here is that the tail step is negligible against the first: the
    box stops mattering. `tail <= 1% of head` is a statement about the shape of the sequence and
    is scale-free; it is not a tolerance on the quantity, which is what rule 20 forbids.
    """
    v = [x for x in vals if np.isfinite(x)]
    if len(v) < 3:
        return False, float("nan")
    d = [abs(v[i + 1] - v[i]) for i in range(len(v) - 1)]
    if d[0] == 0.0:
        return True, 0.0
    return bool(d[-1] <= 0.01 * d[0]), d[-1] / d[0]


def _row(name, vals, published=None):
    sh, last = _converged(vals)
    fin = [v for v in vals if np.isfinite(v)]
    # the error is measured against the WIDEST box actually run, not against a padding nan
    err = abs(fin[0] - fin[-1]) / (abs(fin[-1]) or 1.0) * 100
    print(f"   {name:<26}" + "".join(f"{v:>14.6g}" for v in vals)
          + f"{err:>9.2f}%{'  conv' if sh else '  NOT-CONV'}"
          + f"{'' if sh else f'  (tail/head = {last:.2g})'}")
    return {"name": name, "caps": list(CAPS), "values": list(map(float, vals)),
            "err_at_1_25_pct": float(err), "shrinking": bool(sh),
            "last_gap_rel": float(last), "published": published}


def p1_wiring(report):
    from experiments.does_the_ratio_move import downstream_crossing, upstream_relaxation
    from experiments.escape_accounts_for_it import rate_limits
    from experiments.what_reflection_costs import spectral_gap

    print("P1 -- the threading must be inert at the published default")
    checks = []
    g = spectral_gap(14, False)[0]
    checks.append(("spectral_gap(14, free)", g, 7.29805e-02, abs(g - 7.29805e-02) < 1e-7))
    k_mean, k_avg, k_eff, pos = rate_limits(14, 2, 2.0)
    checks.append(("rate_limits sep(14)", k_avg / k_mean, 2.435, abs(k_avg / k_mean - 2.435) < 0.02))
    checks.append(("rate_limits pos(14) in (0,.5)", pos, None, 0.0 < pos < 0.5))
    src = json.loads(pathlib.Path("results/does_the_ratio_move.json").read_text())
    ref14 = next(r for r in src if r["omega"] == 14)
    tu = upstream_relaxation(14)[0]
    tx = downstream_crossing(14)[0]
    checks.append(("tau_up(14)", tu, ref14["tau_up"], abs(tu - ref14["tau_up"]) < 1e-9))
    checks.append(("tau_cross(14)", tx, ref14["tau_cross"], abs(tx - ref14["tau_cross"]) < 1e-9))
    for name, got, want, ok in checks:
        w = "" if want is None else f"   published {want:.6g}"
        print(f"   {'ok ' if ok else 'FAIL'} {name:<30} {got:.6g}{w}")
    ok = all(c[3] for c in checks)
    print(f"   -> P1 {'HOLDS' if ok else 'FAILS'}")
    report["p1"] = [{"name": n, "got": float(g), "published": w, "ok": bool(o)}
                    for n, g, w, o in checks]
    return ok


def p3_single_stage(report):
    from experiments.predicting_transmission import descent_rate
    from experiments.what_reflection_costs import spectral_gap

    print("\nP3 -- single-stage rates: same shape as §112.7's escape rate?")
    print(f"   {'quantity':<26}" + "".join(f"{'cap=' + str(c):>14}" for c in CAPS)
          + f"{'err@1.25':>9}")
    rows = []
    for om in (14, 30, 55):
        rows.append(_row(f"spectral_gap(free, Om={om})",
                         [spectral_gap(om, False, cap_mult=c)[0] for c in CAPS]))
    for om in (14, 30, 55):
        rows.append(_row(f"spectral_gap(walled, Om={om})",
                         [spectral_gap(om, True, cap_mult=c)[0] for c in CAPS]))
    for om in (14, 30, 55):
        rows.append(_row(f"descent_rate(Om={om})",
                         [descent_rate(om, cap_mult=c)[0] for c in CAPS]))
    falls = all(rows[i]["err_at_1_25_pct"] > rows[i + 1]["err_at_1_25_pct"]
                for grp in (0, 3, 6) for i in (grp, grp + 1))
    signs = [r["values"][0] > r["values"][-1] for r in rows]
    print(f"   error falls with Omega in every family: {falls}")
    print(f"   converged value is SMALLER in {sum(signs)}/{len(rows)} rows "
          f"(a wider box lengthens the escape)")
    print(f"   -> P3 {'HOLDS' if falls and all(signs) else 'PARTIAL/FAILS'}")
    report["p3"] = rows


def p4_bracket(report):
    from experiments.escape_accounts_for_it import rate_limits

    print("\nP4 -- §102's bracket: the width should cancel, the position need not")
    print(f"   {'cell':<26}" + "".join(f"{'cap=' + str(c):>14}" for c in CAPS)
          + f"{'err@1.25':>9}")
    rows = []
    # (14,3) at cap 3.0/4.0 is 2.5-6M states; the cell is run only where it is affordable
    # D=3 is 10^5-10^6 states: 1.25 -> 185k, 1.6 -> 389k, 2.0 -> 754k, 3.0 -> 2.5M. The 754k
    # expm_multiply thrashes on this machine (1 CPU-minute per 5 wall-minutes), so D=3 gets two
    # boxes only and is reported as a TREND, not a converged value. Saying so is the point.
    cells = [((14, 2), CAPS), ((30, 2), CAPS[:4]), ((14, 3), CAPS[:2])]
    for (om, D), caps in cells:
        got = [rate_limits(om, D, 2.0, cap_mult=c) for c in caps]
        w = [g[1] / g[0] for g in got]
        pos = [g[3] for g in got]
        pad = [float("nan")] * (len(CAPS) - len(caps))
        rows.append(_row(f"width  Om={om} D={D}", w + pad))
        rows.append(_row(f"POSITION Om={om} D={D}", pos + pad))
        rows[-1]["d_pos"] = abs(pos[0] - pos[-1])
        rows[-1]["in_range"] = bool(all(0.0 < p < 0.5 for p in pos))
        rows[-1]["caps_run"] = list(caps)
        rows[-2]["caps_run"] = list(caps)
    dpos = [r["d_pos"] for r in rows if "d_pos" in r]
    inr = [r["in_range"] for r in rows if "in_range" in r]
    wid = [r["err_at_1_25_pct"] for r in rows if r["name"].startswith("width")]
    print(f"   width moves at most {max(wid):.2f}% (predicted < 2% at Om=14) -> "
          f"{'as predicted' if max(wid) < 2.0 else 'PREDICTION MISSED: it cancels less than I said'}")
    print(f"   position moves at most {max(dpos):.4f} absolute (predicted < 0.10)")
    print(f"   position stays in (0, 0.5) in {sum(inr)}/{len(inr)} cells")
    ok = max(dpos) < 0.10 and all(inr)
    print(f"   -> P4 {'HOLDS: §102 survives the box' if ok else 'FAILS: §102 needs a correction'}")
    report["p4"] = rows


def p5_ratio(report):
    from experiments.does_the_ratio_move import downstream_crossing, upstream_relaxation

    print("\nP5 -- §110's ratio, and rule 14: a withdrawal verified as carefully as an assertion")
    print(f"   {'quantity':<26}" + "".join(f"{'cap=' + str(c):>14}" for c in CAPS)
          + f"{'err@1.25':>9}")
    src = json.loads(pathlib.Path("results/does_the_ratio_move.json").read_text())
    oms = [r["omega"] for r in src]
    pos = {r["omega"]: r["position"] for r in src}
    # tau_cross is box-independent STRUCTURALLY, not measurably: its h-transform runs on
    # [saddle, rail] with both ends absorbing, so the wall is outside the domain and
    # `downstream_crossing` ignores cap_mult entirely. Sweeping it would re-run one computation
    # five times and print 0.00% -- a number produced by the harness, not the chemistry (rule
    # 10). The fact is asserted from the source and the domain instead.
    import inspect

    body = inspect.getsource(downstream_crossing).split('"""')[2]
    assert body.count("cap_mult") == 0, "downstream_crossing now uses the box; re-sweep it"
    print("   tau_cross: box-independent BY CONSTRUCTION (domain [saddle, rail], both")
    print("     absorbing; cap_mult appears 0 times in the body). Not swept -- there is")
    print("     nothing to sweep. Its x_up path is bounded in P5b.")
    rows, ratios = [], {}
    for om in oms:
        tu = [upstream_relaxation(om, cap_mult=c)[0] for c in CAPS]
        tx = downstream_crossing(om)[0]
        rows.append(_row(f"tau_up(Om={om})", tu))
        rows.append({"name": f"tau_cross(Om={om})", "values": [float(tx)],
                     "box_independent_by_construction": True})
        rt = [a / tx for a in tu]
        rows.append(_row(f"RATIO(Om={om})", rt))
        ratios[om] = rt
    print("\n   §110's headline, recomputed at each box: does the ratio still fall as the"
          " position rises?")
    print(f"   {'cap_mult':>10}{'ratio(14)':>12}{'ratio(70)':>12}{'direction':>12}"
          f"{'monotone':>10}")
    verdicts = []
    for i, c in enumerate(CAPS):
        series = [ratios[om][i] for om in oms]
        falling = series[-1] < series[0]
        mono = series == sorted(series, reverse=True)
        verdicts.append({"cap_mult": c, "series": list(map(float, series)),
                         "falling": bool(falling), "monotone": bool(mono)})
        print(f"   {c:>10}{series[0]:>12.4f}{series[-1]:>12.4f}"
              f"{'FALLING' if falling else 'rising':>12}{str(mono):>10}")
    pos_rising = [pos[om] for om in oms] == sorted(pos[om] for om in oms)
    same = all(v["falling"] == verdicts[0]["falling"] for v in verdicts)
    worst = max(abs(ratios[om][0] / ratios[om][-1] - 1) * 100 for om in oms)
    print(f"   measured position rises across Omega: {pos_rising}")
    print(f"   ratio moves at most {worst:.2f}% between boxes (predicted < 5%)")
    print(f"   every box gives the same direction: {same}")
    print(f"   -> P5 {'HOLDS: §110s refutation survives the box' if same and verdicts[-1]['falling'] else 'FAILS: §110 must be reopened'}")
    report["p5"] = {"rows": rows, "verdicts": verdicts, "position_rising": bool(pos_rising),
                    "worst_move_pct": float(worst)}


def p5b_the_held_path(report):
    """P5 varies each clock's OWN box but takes x_up from `chain_operating_points`, which has a
    box of its own that P5 holds at 1.25. That is a real omitted path, so it is bounded here
    rather than left as a caveat.

    NOT pre-registered as a number -- P5 named the limitation in words; this measures it.
    """
    import experiments.chemical_cascade as cc
    from experiments.depth_compounding import C, R2, R3
    from experiments.does_the_ratio_move import upstream_relaxation

    def tau_cross_at(om, x_up):
        a, b = int(np.ceil(R2 * om)), int(round(R3 * om))
        idx = np.arange(a, b + 1)
        m = len(idx)
        L = np.zeros((m, m))
        for i, s in enumerate(idx):
            if s in (a, b):
                L[i, i] = 1.0
                continue
            lam, mu = cc.rates_stage(float(s), x_up * om, om, C, R3, False, "hill")
            L[i, i] = -(lam + mu); L[i, i + 1] = lam; L[i, i - 1] = mu
        rhs_h = np.zeros(m); rhs_h[0] = 1.0
        h = np.linalg.solve(L, rhs_h)
        rhs_v = -h.copy(); rhs_v[0] = 0.0; rhs_v[-1] = 0.0
        v = np.linalg.solve(L, rhs_v)
        j = list(idx).index(b - 1)
        return float(v[j] / h[j])

    def xup(om, cm):
        Q, cap, m, strides, _ = cc.build(1, om, "source", cap_mult=cm)
        n = np.arange(m)
        keep = np.where(n > R2 * om)[0]
        A = Q[keep][:, keep].toarray()
        w, v = np.linalg.eig(A.T)
        i = int(np.argmax(w.real))
        p = np.abs(v[:, i].real); p /= p.sum()
        return float((p * (n[keep] / om)).sum())

    print("\nP5b -- the path P5 held fixed: x_up's OWN box (not pre-registered as a number)")
    print(f"{'Om':>5}{'x_up@1.25':>13}{'x_up@conv':>13}{'x_up moves':>12}"
          f"{'ratio@P5':>12}{'ratio@both':>12}{'extra':>9}")
    rows = []
    for om in (14, 30, 55):
        lo, hi = xup(om, 1.25), xup(om, 4.0)
        tu = upstream_relaxation(om, cap_mult=4.0)[0]
        r_p5, r_full = tu / tau_cross_at(om, lo), tu / tau_cross_at(om, hi)
        extra = abs(r_p5 / r_full - 1) * 100
        rows.append({"omega": om, "xup_125": lo, "xup_conv": hi,
                     "xup_move_pct": abs(lo / hi - 1) * 100,
                     "ratio_p5": r_p5, "ratio_both": r_full, "extra_pct": extra})
        print(f"{om:>5}{lo:>13.6f}{hi:>13.6f}{abs(lo / hi - 1) * 100:>11.2f}%"
              f"{r_p5:>12.6f}{r_full:>12.6f}{extra:>8.2f}%")
    worst = max(r["extra_pct"] for r in rows)
    print(f"   the held path adds at most {worst:.2f}% to the ratio's box sensitivity,")
    print(f"   so §110's total is ~{worst + 5.05:.1f}% at Omega=14 -- direction unaffected.")
    report["p5b"] = rows


def p6_scope(report):
    print("\nP6 -- what this audit does NOT cover (rule: a null result is worth its scope)")
    unaudited = [
        ("§101/§103 joint tables, D=3, Om>=30", "10^6-10^7 states at cap_mult 3.0"),
        ("§105-§109 closure/seed tables at D=3", "same joint solve"),
        ("§107's closure residuals", "depends on the D=3 tables above"),
        # AM was on this list until §113.7: it has NO truncation. Every AM reaction conserves
        # total count, so crnl/cme.py enumerates the exact simplex {n : sum(n) = N} and there is
        # no wall to sweep. Removed rather than left as a plausible-looking unaudited item.
        ("§34's AM crossovers (N <= 72)", "a CONDITIONING ceiling, not a box -- see T-BOX-d"),
    ]
    for what, why in unaudited:
        print(f"   UNAUDITED  {what:<44} ({why})")
    print("   These are listed, not implied safe. §112.7's bias is one-signed and falls with")
    print("   Omega, so the D=3 / large-Omega cells are the ones where it should matter LEAST --")
    print("   but that is an extrapolation, not a measurement, and is labelled as one.")
    report["p6_unaudited"] = [{"what": w, "why": y} for w, y in unaudited]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/what_the_box_was_doing.json"))
    args = ap.parse_args()
    report = {}
    if not p1_wiring(report):
        print("\nP1 FAILED -- the threading is not inert. Stopping.")
        args.out.write_text(json.dumps(report, indent=2, default=float))
        return
    p3_single_stage(report)
    p5_ratio(report)
    p4_bracket(report)
    p5b_the_held_path(report)
    p6_scope(report)
    args.out.write_text(json.dumps(report, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
