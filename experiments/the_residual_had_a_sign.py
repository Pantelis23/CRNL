"""§116 (T-CASC-x, T-NUM-c) -- §107's residual drift, with stage 2's rate replaced by the exact one.

§107 applied every known correction to §103's closure and was left with a residual that drifts by
1.83x across Omega = 14-70 (model/measured 0.790 -> 1.449 with the geometric-mean rate). It then
ruled out the one remaining approximation -- the fast/frozen averaging of stage 2's rate -- on a
SIGN argument:

    "§102.1's position moves *toward* the fast end as Omega grows (0.2200 -> 0.1072 over
     Omega = 14 -> 55), so the geometric mean should become *more* accurate at large Omega. The
     measured ratio does the opposite ... **The position argument has the wrong sign for the
     residual.**"

That argument sent §108, §109 and §110 after a mechanism outside the averaging family, and
T-CASC-x ended with "no candidate left". §115 found the position was mis-estimated: under the
exact Q-process gap it moves the OTHER way, toward FROZEN (0.449 -> 0.996 at cap 2.0), and sits
above the geometric (fast) limit at every Omega with the separation growing. **So the sign flips**:
a closure built on the fast-limit rate SHOULD get worse with Omega, and it does. Whether the
magnitude matches is the test, and it is absolute (rule 16) -- nothing below is fitted.

The closure's stage-2 factor is `los[1] = pi_low * (1 - exp(-k2 t))` with k2 = exp(<ln k>). This
replaces k2 by the exact gap lambda1 - lambda0 of the stage-1-alive joint generator (§115.3) and
changes nothing else.

PREDICTIONS. NOT BLIND, AND SAID SO: before writing this file I estimated from §115's stored
positions (exact vs fast-limit) that the swap takes the ratio's span from 1.83x to ~1.26x, with
every ratio moving DOWN and the large-Omega end moving most. The numbers below are computed; the
estimate was not, and agreement would be agreement with arithmetic I had already half done.

  P1  WIRING. With k2 = geometric_rate and pinned pi_low this reproduces §107's stored "+ geometric"
      column (0.790 ... 1.449) to 1e-9.
  P2  THE SWAP. Every ratio moves DOWN (k_exact > k_geom at every Omega), the large-Omega end most.
  P3  THE VERDICT, and it can print either. §107's drift is 1.83x. If after the swap the span is
      below 1.35x, stage 2's rate was the dominant missing factor and §107's sign argument --
      and the "no candidate left" that followed -- was the misestimated position. If the span is
      still above 1.5x, the rate is not what drives the residual and T-CASC-x stays as open as
      §110 left it. Between the two, partial, and reported as such.
  P4  RULE 15. §115 found the pinned pi_low is also biased against the conditioned occupancy.
      Report the closure with BOTH replaced (exact rate, conditioned pi_low) as a second variant.
      No prediction for it: it changes the two-state bookkeeping as well as the rate, and the
      closure's pi_low multiplies a factor that is not in the exact-gap definition.

WHAT THIS CANNOT SETTLE. The measured side is §106's free-chain solve with §101's default seed,
which §109 showed carries seeding transients of its own. This asks whether the closure's error
tracks stage 2's rate; it does not certify the measured ratio. And a level offset that survives
the swap is not explained by it.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib

import numpy as np

from experiments.depth_compounding import R3
from experiments.chain_without_a_joint_solve import chain_operating_points
from experiments.escape_accounts_for_it import escape_rate
from experiments.the_corrected_closure import P_TRANSMIT_MEASURED, geometric_rate, pi_low
from experiments.the_survival_term import conditioned_gap, conditioned_qsd_low
from experiments.two_axes_for_the_position import PUBLISHED_RAILS, Elem

T = 2.0


def closure_with(om, k2, pi2):
    """§107's closure (indexing + two-state, measured p_transmit), stage 2's (rate, pi) given."""
    k0 = escape_rate(om, R3)                 # stage 1 is driven by the rail (§106.3)
    los0 = pi_low(om, R3) * (1 - math.exp(-k0 * T))
    los1 = pi2 * (1 - math.exp(-k2 * T))
    return (los0 * P_TRANSMIT_MEASURED) / ((1 - los0) * los1)


def audit(rep):
    """T-NUM-c: which §102-§110 CONCLUSIONS survive the exact rate. Bookkeeping over stored
    numbers, not a new measurement: every value here was computed in §113 or §115 and seen.
    The criteria are the published claims' own words, fixed before tallying."""
    box = json.loads(pathlib.Path("results/what_the_box_was_doing.json").read_text())
    ex = json.loads(pathlib.Path("results/the_survival_term.json").read_text())
    ratio = {v["cap_mult"]: v["series"] for v in box["p5"]["verdicts"]}
    pos2 = [r["pos_exact_cap2"] for r in ex["p10"]["v1"]]
    pos125 = [r["pos_exact"] for r in ex["p9"]["published"]]
    out = {}
    print("\nT-NUM-c -- published conclusions against the exact rate")
    s110 = all(ratio[c] == sorted(ratio[c], reverse=True) for c in (1.25, 2.0)) \
        and pos2 == sorted(pos2) and pos125 == sorted(pos125)
    out["s110_refutation"] = bool(s110)
    print(f"   §110 'the position rises while tau_up/tau_cross falls': "
          f"{'SURVIVES' if s110 else 'FAILS'} (ratio {ratio[2.0][0]:.3f} -> {ratio[2.0][-1]:.3f};"
          f" exact position {pos2[0]:.3f} -> {pos2[-1]:.3f}, both monotone, both boxes)")
    inside = all(0 < p < 1 for p in pos2 + pos125)
    out["s109_exits_bracket"] = not inside
    print(f"   §109/§110 'and exits the bracket': {'WITHDRAWN' if inside else 'stands'}"
          f" (max exact position {max(pos2):.3f})")
    out["s102_fast_end_Om14"] = pos2[0] < 0.5
    out["s102_fast_end_Om30"] = pos2[2] < 0.5
    print(f"   §102.1 'near the fast end' (position < 0.5): Om=14 "
          f"{'holds' if pos2[0] < 0.5 else 'fails'} ({pos2[0]:.3f}),"
          f" Om=30 {'holds' if pos2[2] < 0.5 else 'FAILS'} ({pos2[2]:.3f}); D=3 unaudited")
    out["s107_toward_fast_end"] = pos2 == sorted(pos2, reverse=True)
    print(f"   §107 'the position moves toward the fast end as Omega grows': "
          f"{'stands' if out['s107_toward_fast_end'] else 'WITHDRAWN -- it moves toward FROZEN'}")
    rep["audit"] = out


def seed_kill_test(rep):
    """T-CASC-ai's kill test, criterion fixed before running. §106's measured contam/pure used
    §101's default seed (stage 1 reflected law, stage 2 a delta at its rail), which §109 showed
    carries stage-2 transients that shrink with Omega -- a large-Omega rise in the direction of
    the U. Re-measure with matched QSD seeding and recompute the exact-rate column.
    Criterion: if max/min of the exact-rate ratio across Omega falls below 1.10, the U was the
    measurement's seed. Otherwise the remainder persists and the seed is exonerated."""
    from experiments.free_upstream_depth import channel_split, solve
    print("\nT-CASC-ai kill test -- §106's measured ratio re-measured with matched seeding")
    out = []
    for r in rep["rows"]:
        om = r["omega"]
        p, ref, dims, strides, walled = solve(om, 2, 0, T, matched_seed=True)
        tot, pure, contam = channel_split(p, om, dims, strides, ref, walled)
        meas_m = contam / pure
        mus, _ = chain_operating_points(om, 2)
        model = closure_with(om, r["k_exact"], r["pi_pinned"])
        out.append({"omega": om, "meas_default": r["rebuilt"], "meas_matched": meas_m,
                    "ratio_exact_matched": model / meas_m})
        print(f"   Om={om:>3}  measured contam/pure (matched) {meas_m:.4e}"
              f"   exact-rate ratio: default seed {r['ratio_exact']:.4f}"
              f"  matched seed {model / meas_m:.4f}", flush=True)
    ys = [o["ratio_exact_matched"] for o in out]
    span = max(ys) / min(ys)
    verdict = span < 1.10
    # The first version printed "the remainder PERSISTS; the seed is exonerated" off the span
    # alone. The first clause was right and the second false: the seed removes the whole low-Omega
    # arm and the level (§116.5). A span merges the two arms, so they are now read separately.
    lo = [o["ratio_exact_matched"] for o in out if o["omega"] <= 40]
    hi = [o["ratio_exact_matched"] for o in out if o["omega"] >= 40]
    print(f"   span {span:.3f}x -> {'the U WAS the seed' if verdict else 'the remainder PERSISTS'}")
    print(f"   by arm: Omega <= 40 spans {max(lo)/min(lo):.3f}x (the seed's arm), "
          f"Omega >= 40 rises {hi[-1]/hi[0]:.3f}x (what survives)")
    rep["seed_kill"] = {"rows": out, "span": span, "u_was_seed": bool(verdict)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/the_residual_had_a_sign.json"))
    args = ap.parse_args()

    meas = {r["omega"]: r["meas_ratio"] for r in
            json.loads(pathlib.Path("results/where_the_expansion_frays.json").read_text())["cells"]}
    pub = {r["omega"]: r["ratios"][3] for r in
           json.loads(pathlib.Path("results/the_corrected_closure.json").read_text())}
    e = Elem("published", *PUBLISHED_RAILS, cap_mult=1.25)     # the box §106/§107 used

    print(f"{'Om':>5}{'§107 (geom)':>13}{'rebuilt':>10}{'k_geom':>12}{'k_exact':>12}"
          f"{'ratio':>8}{'exact rate':>12}{'+ cond pi':>11}")
    rows = []
    for om in sorted(pub):
        mus, _ = chain_operating_points(om, 2)
        pin = pi_low(om, mus[0])
        kg = geometric_rate(om, mus[0])
        l0, l1, _ = conditioned_gap(e, om)
        kx = l1 - l0
        pc = conditioned_qsd_low(e, om)
        r_geom = closure_with(om, kg, pin) / meas[om]
        r_exact = closure_with(om, kx, pin) / meas[om]
        r_both = closure_with(om, kx, pc) / meas[om]
        rows.append({"omega": om, "published": pub[om], "rebuilt": r_geom, "k_geom": kg,
                     "k_exact": kx, "pi_pinned": pin, "pi_cond": pc,
                     "ratio_exact": r_exact, "ratio_both": r_both})
        print(f"{om:>5}{pub[om]:>13.4f}{r_geom:>10.4f}{kg:>12.4e}{kx:>12.4e}{kx/kg:>8.3f}"
              f"{r_exact:>12.4f}{r_both:>11.4f}", flush=True)

    w = max(abs(r["rebuilt"] / r["published"] - 1) for r in rows)
    print(f"\nP1 -- rebuilt vs §107's stored column: worst rel {w:.1e}"
          f"  -> {'HOLDS' if w < 1e-9 else 'FAILS'}")
    down = all(r["ratio_exact"] < r["rebuilt"] for r in rows)
    print(f"P2 -- every ratio moves down: {down}")
    span = lambda key: max(r[key] for r in rows) / min(r[key] for r in rows)
    s0, s1, s2 = span("rebuilt"), span("ratio_exact"), span("ratio_both")
    print(f"P3 -- span across Omega: geometric {s0:.3f}x  ->  exact rate {s1:.3f}x"
          f"   (both replaced: {s2:.3f}x)")
    verdict = ("the rate WAS the dominant missing factor" if s1 < 1.35 else
               "the rate is NOT what drives the residual" if s1 > 1.5 else
               "PARTIAL -- the rate accounts for some of the drift, not all")
    print(f"   -> {verdict}")
    mono = [r["ratio_exact"] for r in rows]
    print(f"   exact-rate ratios by Omega: " + ", ".join(f"{x:.3f}" for x in mono)
          + f"   (monotone: {mono == sorted(mono)})")

    rep = {"rows": rows}
    audit(rep)
    seed_kill_test(rep)
    args.out.write_text(json.dumps({**rep, "span_geom": s0, "span_exact": s1,
                                    "span_both": s2, "p1_worst": w, "p2_all_down": down},
                                   indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
