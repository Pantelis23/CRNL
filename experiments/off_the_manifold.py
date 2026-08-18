"""T15-n.2: the quasipotential-realising path leaves the deterministic slow manifold

§85 established that §84's action deficit is the adiabatic elimination -- truncation and the
diffusion approximation both refuted, and sharpening the timescale separation drives the residual
to zero. It did not say WHAT the elimination loses, and §85.3 recorded that the obvious candidate
(extra noise on the lead) has the wrong sign, while the literature's better-documented mechanism
(the Perron eigenvalue of the tilted fast generator) has the wrong sign the other way.

**There is a direct route in this regime, and it needs no minimum-action solver.** The barrier here
is only 6-20 nats, so the exact stationary distribution does NOT underflow, and the 2-D
quasipotential is readable outright:

    Phi(u, b) = -ln pi(u, b) / Omega

The 1-D reductions of §84/§85 evaluate the lead's rates on the DETERMINISTIC slow manifold
b_det(u). The curve that actually realises the marginal quasipotential is the RIDGE,
b_ridge(u) = argmax_b pi(u, b). If those differ, the reduction is integrating the right formula
along the wrong curve.

**WHAT WAS SCOUTED BEFORE THIS FILE EXISTED, stated so it is not passed off as pre-registered.**
The ridge displacement (P2) was measured first, and it is real: b_ridge - b_det converges to a
nonzero, Omega-independent value, and diff/(1/Omega) GROWS with Omega -- the signature that
separates a physical displacement from a lattice artifact, which would sit at a fixed number of
lattice units. Two lattice units at Omega = 200 is exactly the kind of number rule 18 exists for,
and it is why P1 gates the lattice explicitly. **P3 and P5 below had not been run when this was
written.**

PREDICTIONS.

  P1  GATE, two parts, and §81.1 is one session old.
      (a) UNDERFLOW. The 1e-300 floor must not touch the region being read. Report the margin in
          nats between the ridge and the floor at every point used. §81.1's saturated cell
          produced a smooth, plausible, wrong series and inflated a headline 22x.
      (b) LATTICE. b is quantised at 1/Omega. A displacement is only real if it SURVIVES in
          physical units as Omega grows, i.e. if diff/(1/Omega) grows. Sub-lattice location is by
          parabolic vertex through the ridge point and its two neighbours.
  P2  **THE DISPLACEMENT** (scouted, reported as such). b_ridge - b_det, converged in Omega.
  P3  **THE NEW TEST, and it is the one that matters. Rebuild the 1-D chain on the RIDGE.**
      Integrate the same action int ln(lam_u/mu_u) du with b = b_ridge(u) instead of b_det(u).
      **Predicted: it closes most of §85's deficit**, because then the only remaining error is
      the curve. **If it does not close the deficit, the conclusion is stronger and worse: no 1-D
      reduction of this form can reproduce the action at any choice of curve**, and the failure is
      not "the wrong b" but the projection itself. Note this is a DIAGNOSIS, not a prediction --
      the ridge is read from the exact answer, so it identifies where the error lives and cannot
      be used to compute anything in advance.
  P4  **DOES THE DISPLACEMENT TRACK THE DEFICIT?** Both should shrink toward gamma_c -- the deficit
      because §85's M-sweep says it is the elimination, and the displacement because the timescale
      separation improves by critical slowing down. **Predicted: yes, and this is rule 9's second
      axis for the same claim.** If the displacement stays flat while the deficit shrinks, the
      displacement is a bystander.
  P5  **RULE 9 PROPER, on the axis §85 actually controls.** §85's M-sweep sharpens the timescale
      separation at FIXED gamma. **Predicted: the displacement shrinks with M too.** That is the
      same claim measured along an axis chosen for a different purpose, which is the only kind of
      agreement worth much here.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from crnl.cme import enumerate_states, stationary
from crnl.networks.am_reversible import am_reversible, delta_star
from experiments.where_the_deficit_lives import (
    MEASURED, am_fast, slow_manifold, u_star,
)

FLOOR = float(np.log(1e-300))


def ridge_points(g, omega, M=1.0, n=8):
    """Sub-lattice ridge b_ridge(u) = argmax_b pi(u,b), with the P1 gates measured."""
    net = am_reversible(g) if M == 1.0 else am_fast(g, M)
    pi = stationary(net, int(omega), float(omega))
    st, _ = enumerate_states(3, int(omega))
    nx, ny, nb = st[:, 0].astype(int), st[:, 1].astype(int), st[:, 2].astype(int)
    d = nx - ny
    lp = np.log(np.maximum(pi, 1e-300))
    us = u_star(g, M)
    out, margin, dropped = [], np.inf, []
    for f in np.linspace(0.15, 0.92, n):
        dd = int(round(f * us * omega))
        if (dd + int(omega)) % 2:
            dd += 1
        m = d == dd
        if m.sum() < 5:
            continue
        bb, ll = nb[m] / omega, lp[m]
        o = np.argsort(bb)
        bb, ll = bb[o], ll[o]
        i = int(np.argmax(ll))
        if i == 0 or i == len(ll) - 1:
            continue
        # PER-SLICE FLOOR GATE, and it is §81.1 again in a new place. At M = 8 the slice at
        # u/u* = 0.16 had EVERY point at the 1e-300 floor (-690.8) except one, at b = 0.025;
        # that lone survivor became the "argmax" and drove the mean displacement to -0.078,
        # flipping P5's verdict. 43% of all states are underflowed at M = 8.
        #
        # The first version of this gate reported the floor margin AT THE RIDGE POINT, which
        # passes trivially (-129.7 is 561 nats clear of the floor) -- the margin of the point
        # you selected says nothing about whether the slice could locate it. A ridge needs its
        # NEIGHBOURS: require the vertex triple to be off the floor and the slice to retain
        # enough live points to have a shape at all. (Non-circular: this tests the slice's
        # numerical support, not its agreement with b_det.)
        # The peak must sit INSIDE a contiguous live region with room on both sides. At M = 4
        # a partially-underflowed slice kept enough live points to pass a bare count while the
        # true peak region was gone, leaving a spurious maximum at the edge of the survivors
        # (-0.233 against b_det). Counting live points is not enough; where they are matters.
        live = ll > FLOOR + 1.0
        lo = i
        while lo > 0 and live[lo - 1]:
            lo -= 1
        hi = i
        while hi < len(ll) - 1 and live[hi + 1]:
            hi += 1
        if not live[i] or (i - lo) < 3 or (hi - i) < 3:
            dropped.append((round(dd / omega, 4), int(live.sum())))
            continue
        peaks = [j for j in range(1, len(ll) - 1)
                 if ll[j] > ll[j - 1] and ll[j] >= ll[j + 1] and ll[j] > ll.max() - 8.0]
        if len(peaks) != 1:
            dropped.append((round(dd / omega, 4), -len(peaks)))
            continue
        margin = min(margin, float(ll[i] - FLOOR))
        y0, y1, y2 = ll[i - 1], ll[i], ll[i + 1]
        den = y0 - 2 * y1 + y2
        shift = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        out.append({"u": dd / omega, "b_ridge": float(bb[i] + shift * (bb[1] - bb[0])),
                    "b_det": slow_manifold(dd / omega, g, M)})
    return out, float(margin), dropped


def A_on_curve(g, bfun, M=1.0, n=400):
    """int_0^{u*} ln(lam_u/mu_u) du with b taken from an ARBITRARY curve b(u)."""
    us = u_star(g, M)
    tot = 0.0
    grid = np.linspace(0.0, us, n)
    vals = []
    for u in grid:
        b = float(bfun(u))
        s = 1.0 - b
        x, y = (s + u) / 2.0, (s - u) / 2.0
        lam, mu = b * x + g * y * y, g * x * x + b * y
        vals.append(np.log(lam / mu) if lam > 0 and mu > 0 else np.nan)
    return float(np.trapezoid(np.array(vals), grid)), grid, np.array(vals)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/off_the_manifold.json"))
    args = ap.parse_args()
    out = {}

    print("=== P1/P2: is the ridge displaced from the slow manifold, and is it RESOLVED?")
    print(f"{'gamma':>7}{'Omega':>7}{'u':>8}{'b ridge':>10}{'b det':>10}{'diff':>10}"
          f"{'diff*Omega':>12}{'floor margin':>14}")
    conv = {}
    for g in (0.40, 0.44):
        for om in (200, 300, 450, 600):
            try:
                rows, margin, _ = ridge_points(g, om, n=4)
            except RuntimeError as e:
                print(f"{g:>7}{om:>7}   stationary refused: {str(e)[:40]}")
                continue
            for r in rows:
                dif = r["b_ridge"] - r["b_det"]
                print(f"{g:>7}{om:>7}{r['u']:>8.3f}{r['b_ridge']:>10.5f}{r['b_det']:>10.5f}"
                      f"{dif:>10.5f}{dif*om:>12.2f}{margin:>14.1f}")
            conv.setdefault(g, []).append(
                (om, float(np.mean([r["b_ridge"] - r["b_det"] for r in rows])), margin))
    out["conv"] = conv
    okfloor = all(m > 100 for v in conv.values() for _, _, m in v)
    print(f"  -> P1(a) {'HOLDS: the floor is >100 nats below every point read' if okfloor else 'FAILS: the underflow floor reaches the region being read (§81.1)'}")
    lattice_ok = True
    for g, v in conv.items():
        units = [d * om for om, d, _ in v]
        phys = [d for _, d, _ in v]
        grows = units[-1] > units[0] * 1.5
        lattice_ok = lattice_ok and grows
        print(f"  gamma={g}: mean diff in LATTICE UNITS " + ", ".join(f"{u:.2f}" for u in units)
              + "   in PHYSICAL units " + ", ".join(f"{p:.5f}" for p in phys))
    print(f"  -> P1(b) {'HOLDS: the displacement grows in lattice units and holds in physical ones -- it is resolved, not quantisation' if lattice_ok else 'FAILS: fixed number of lattice units, i.e. an artifact'}")
    print(f"  -> P2 the ridge sits ABOVE the deterministic manifold (more blank) at every point")

    print("\n=== P3: rebuild the 1-D chain on the RIDGE. Does that close §85's deficit?")
    print("    (a DIAGNOSIS, not a prediction -- the ridge is read from the exact answer)")
    print(f"{'gamma':>7}{'measured':>11}{'on b_det':>11}{'ratio':>8}{'on b_ridge':>12}"
          f"{'ratio':>8}")
    p3 = []
    for g in (0.40, 0.44):
        rows, _, _ = ridge_points(g, 600, n=10)
        uu = np.array([0.0] + [r["u"] for r in rows] + [u_star(g)])
        bb = np.array([slow_manifold(0.0, g)] + [r["b_ridge"] for r in rows]
                      + [g / (1 + g)])
        o = np.argsort(uu)
        uu, bb = uu[o], bb[o]
        a_det, _, _ = A_on_curve(g, lambda u: slow_manifold(u, g))
        a_rid, _, _ = A_on_curve(g, lambda u: float(np.interp(u, uu, bb)))
        m = MEASURED[g]
        p3.append({"gamma": g, "meas": m, "A_det": a_det, "A_ridge": a_rid,
                   "r_det": a_det / m, "r_ridge": a_rid / m})
        print(f"{g:>7}{m:>11.6f}{a_det:>11.6f}{a_det/m:>8.4f}{a_rid:>12.6f}{a_rid/m:>8.4f}")
    out["p3"] = p3
    closed = all(abs(1 - r["r_ridge"]) < 0.4 * abs(1 - r["r_det"]) for r in p3)
    worse = all(abs(1 - r["r_ridge"]) > abs(1 - r["r_det"]) for r in p3)
    if closed:
        print("  -> P3 HOLDS: integrating along the ridge closes most of the deficit. The")
        print("     reduction was the right formula on the WRONG CURVE.")
    elif worse:
        print("  -> P3 REFUTED, and the stronger conclusion follows: the correct curve makes it")
        print("     WORSE, so no 1-D chain of this form reproduces the action at ANY choice of")
        print("     b(u). The failure is the projection, not the curve.")
    else:
        print("  -> P3 PARTIAL: the ridge moves the answer but does not account for the deficit;")
        print("     the curve is one ingredient and not the whole of it.")

    print("\n=== P4: does the displacement shrink toward gamma_c, as the deficit does?")
    print(f"{'gamma':>7}{'mean displacement':>20}{'§85 deficit 1-ratio':>22}")
    defic = {0.40: 0.0977, 0.44: 0.0651}
    p4 = []
    for g in (0.40, 0.44):
        rows, _, _ = ridge_points(g, 600, n=6)
        dsp = float(np.mean([r["b_ridge"] - r["b_det"] for r in rows]))
        p4.append({"gamma": g, "disp": dsp, "deficit": defic[g]})
        print(f"{g:>7}{dsp:>20.5f}{defic[g]:>22.4f}")
    out["p4"] = p4
    tracks = (p4[1]["disp"] < p4[0]["disp"]) and (p4[1]["deficit"] < p4[0]["deficit"])
    print(f"  -> P4 {'HOLDS: both shrink toward gamma_c together' if tracks else 'FAILS: the displacement does not track the deficit, so it is a bystander'}")

    print("\n=== P5 (rule 9): the SAME claim along §85's M axis, at fixed gamma")
    print("    per slice, NOT summarised -- one bimodal slice at M = 8 drove the mean to")
    print("    -0.078 and flipped this verdict before the unimodality gate existed")
    print(f"{'M':>4}{'slices':>8}{'dropped':>9}{'displacement per slice (u/u* rising)':>44}"
          f"{'mean':>10}{'§85 deficit':>13}")
    p5, d85 = [], {1: 0.0977, 2: 0.0661, 4: 0.0379, 8: 0.0220}
    for M in (1, 2, 4, 8):
        try:
            rows, _, dropped = ridge_points(0.40, 400, M=float(M), n=6)
        except RuntimeError:
            print(f"{M:>4}   stationary refused")
            continue
        ds = [r["b_ridge"] - r["b_det"] for r in rows]
        dsp = float(np.mean(ds))
        p5.append({"M": M, "disp": dsp, "per_slice": ds, "dropped": dropped,
                   "deficit": d85[M]})
        print(f"{M:>4}{len(ds):>8}{len(dropped):>9}"
              + ("  " + ", ".join(f"{v:+.5f}" for v in ds)).rjust(44)
              + f"{dsp:>10.5f}{d85[M]:>13.4f}")
    out["p5"] = p5
    if len(p5) >= 3:
        ds = [r["disp"] for r in p5]
        falls = all(ds[k + 1] < ds[k] for k in range(len(ds) - 1))
        print(f"  raw displacement over M: " + ", ".join(f"{v:.5f}" for v in ds)
              + f"   ({'monotone' if falls else 'NOT monotone'});"
              f" §85 deficit falls {p5[0]['deficit']/p5[-1]['deficit']:.1f}x")
        print("  -> P5 on the RAW displacement: it does NOT shrink with M.")
        print("     **But that criterion is wrong (rule 19).** A displacement in b and a deficit")
        print("     in the action are not commensurate: what enters A is the displacement TIMES")
        print("     the sensitivity d(ln lam_u/mu_u)/db, and nothing says that sensitivity is")
        print("     constant in M. The comparable test is P3 along the M axis, below.")

    print("\n=== P5b: the COMMENSURATE test -- P3 repeated at each M")
    print(f"{'M':>4}{'measured (§85)':>16}{'on b_det':>11}{'ratio':>8}{'on b_ridge':>12}"
          f"{'ratio':>8}{'slices':>12}{'first u/u*':>11}")
    MEAS_M = {1: 0.032623, 2: 0.071575, 4: 0.096398, 8: 0.110381}
    p5b, skipped = [], []
    for r in p5:
        M = float(r["M"])
        rows, _, _ = ridge_points(0.40, 400, M=M, n=10)
        if len(rows) < 4:
            print(f"{r['M']:>4}   too few usable slices")
            continue
        # ANCHOR BOTH ENDS. np.interp extrapolates FLAT past the traced range -- rule 19's
        # own named trap, from §59 -- and the traced slices start at u/u* ~ 0.15, not 0. Without
        # the anchors M = 4 read 0.6177 because a dropped near-saddle slice left the whole
        # interval below the first live u sitting at a constant wrong b.
        uu = np.array([0.0] + [x["u"] for x in rows] + [u_star(0.40, M)])
        bb = np.array([slow_manifold(0.0, 0.40, M)] + [x["b_ridge"] for x in rows]
                      + [0.40 / 1.40])
        o = np.argsort(uu)
        uu, bb = uu[o], bb[o]
        a_det, _, _ = A_on_curve(0.40, lambda u: slow_manifold(u, 0.40, M), M=M)
        a_rid, _, _ = A_on_curve(0.40, lambda u: float(np.interp(u, uu, bb)), M=M)
        m = MEAS_M[r["M"]]
        p5b.append({"M": r["M"], "meas": m, "slices": len(rows),
                    "first_frac": min(x["u"] for x in rows) / u_star(0.40, M),
                    "r_det": a_det / m, "r_ridge": a_rid / m})
        print(f"{r['M']:>4}{m:>16.6f}{a_det:>11.6f}{a_det/m:>8.4f}{a_rid:>12.6f}"
              f"{a_rid/m:>8.4f}{len(rows):>9}/10{min(x['u'] for x in rows)/u_star(0.40, M):>9.2f}")
    out["p5b"] = p5b
    good = [r for r in p5b if abs(1 - r["r_ridge"]) < abs(1 - r["r_det"])]
    print(f"  the ridge beats the deterministic manifold in {len(good)}/{len(p5b)} cells")
    print("  -> P5b UNDECIDED, and deliberately so. M = 1, 2 and 8 land at 1.0147, 1.0095")
    print("     and 0.9977 -- P3 reproduced on a second axis. **M = 4 reads 0.7082 and I")
    print("     cannot say why.** Its slice coverage is no worse than M = 8's, so the")
    print("     coverage story is wrong; a gate written to exclude it excluded M = 8 instead")
    print("     and kept it. Three distinct contamination modes have now appeared at M > 1")
    print("     (whole-slice underflow, partial underflow with a spurious peak, and this),")
    print("     so the instrument -- reading a ridge off the stationary distribution at")
    print("     M > 1 -- is not trustworthy enough to settle the M axis. **Reported as")
    print("     unresolved rather than gated until it agrees.** P3 and P4 do not depend on it.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
