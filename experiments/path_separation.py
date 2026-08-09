"""T-COST-f2: is the governing variable the PATH separation, not the point separation?

§45 found the §39 residual tracking `sep` in shape -- it peaks at rho = 1.5, exactly where
sep bottoms out, while cost is monotone -- but failing every sharp test: the 8x matched-sep
pair disagreed 3x, the rho and gamma knobs gave different curves in 1/sep (slopes 0.441 vs
1.147, opposite intercepts), and resid*sep drifted 4.5x.

`sep_of` takes the eigenvalue ratio at the SYMMETRIC POINT (x = y, b = b*). The traversal
runs over delta in [eps*delta*, theta*delta*], away from it. **A separation measured at one
point need not represent the separation along the path**, and the two knobs deform the
manifold differently -- which would explain both failures without abandoning slaving.

**THE SLOW EIGENVALUE CHANGES SIGN MID-PATH, so the averaging convention is forced, not
chosen.** Measured before writing these predictions, at gamma = 0.20:

    delta/delta*     0.35    0.45    0.55    0.65    0.75    0.80
    sep (rho=0.5)   18.8    31.3   188.2    37.0    15.3    11.5
    sep (rho=4.0)   17.3    29.5   220.4    33.0    14.1    10.7

sep diverges near delta/delta* ~ 0.57 where the drift mu peaks and its derivative crosses
zero. **The arithmetic mean of sep is therefore meaningless and only a HARMONIC mean is
defined** -- which is also what the physics wants, since corrections of size 1/sep add
along the path.

**A STRUCTURAL FACT THAT MAY MATTER MORE THAN THE TEST.** Over rho = 0.5 -> 4 the fast
eigenvalue goes -1.132 -> -3.831 (3.38x) and the slow one 0.0602 -> 0.2210 (3.67x): **both
scale together**, which is precisely why the point ratio matched. rho = 0.5 -> 4 is close
to a uniform time rescale, and §44's P1a proved Sigma is EXACTLY invariant under a uniform
rescale -- so under a perfectly uniform one the residual could not move at all. The 3x
residual difference must come from the ~8% departure from uniformity, or from something
that is not a rate ratio. Both absolute eigenvalues are reported so a later section can use
them; no new variable is fitted here.

THREE CONVENTIONS, all reported (rule 15 -- report every candidate, not the flattering one):

    1/sep_path = < 1/sep >_delta                        uniform in delta
    1/sep_time = int (1/sep) dt / int dt,  dt = ddelta/mu     time-weighted
    1/sep_cost = int (1/sep)(sigma/mu) ddelta / int (sigma/mu) ddelta   cost-weighted

**sep_cost is nominated as primary IN ADVANCE**, because the quantity being explained is a
COST residual and the error should be weighted where the cost accumulates.

PREDICTIONS, written before running:

  P1  GATE. The path separations must differ materially from the point value (7.000 at
      gamma = 0.20, rho = 1). **If any convention reproduces the point value to within 1%,
      the path/point distinction is empty and this experiment has no content.**
  P2  THE TEST, and it is specific. The 8x matched pair (rho = 0.5, rho = 4) has point-seps
      11.70 and 10.56 -- 10% apart -- but residuals 3x apart. **The correct path convention
      must UN-MATCH that pair**, separating their path-seps by something like the ~3x their
      residuals differ by. A convention that leaves them matched has not explained anything.
  P3  COLLAPSE. In 1/sep_path the rho and gamma knobs must fall on ONE curve over the
      overlap. Quantitatively: their fitted slopes, 2.6x apart in 1/sep_point, must come
      within ~30%. That is the pass mark, fixed now.
  P4  ABSOLUTE FORM (rule 16). residual * sep_path constant. Point sep drifted 4.5x; the
      path version must do materially better -- under 2x -- or the 1/sep form is still
      wrong however the average is taken.
  P5  THE REFUTING OUTCOME, named in advance. **If all three conventions leave the knobs
      split, slaving is refuted as rho's mechanism.** T-COST-f closes negative, §44.2's
      lever keeps its measurement permanently without an account, and §39's residual goes
      back to T-COST-c's remaining candidates. Given that both eigenvalues were just shown
      to scale together, **this is the outcome I expect** -- the point ratio and the path
      ratio are unlikely to differ enough to manufacture a 3x.
  P6  If the three conventions disagree about the verdict, the quantity is UNRESOLVED and
      is reported as such rather than by picking the one that works.

=============================================================================
SECOND PASS -- THE FIRST PASS ASKED THE WRONG QUESTION OF THE WRONG QUANTITY.

The first pass ran and refuted P2/P3/P4 under all four conventions; those results stand
below and in FINDINGS. **But re-reading §39.1 and §39.2 before writing them up shows the
framing was wrong in two ways, and the corrected test is different.**

  1. **T-COST-c is CLOSED, not open.** §39.1 showed the entire cost residual is a TIME
     residual -- `Sigma_pred/Sigma_exact` tracks `T_det/MFPT` at correlation +0.9513 -- and
     explicitly WITHDREW candidates (ii) the off-manifold path and (iii) the Jensen gap as
     explanations of the cost. So there are no "remaining candidates" to hand anything back
     to, and the first pass's P5 text says there are.
  2. **§39.2 already established the 1/sep law AND already recorded that its coefficient
     does not transfer between axes**, in as many words and under rule 9: the T-axis value
     (T_det/MFPT - 1)*sep = 0.6465 +- 12.3% predicts 16.2% at gamma = 0.07 against 16%
     measured, but 5.4% at gamma = 0.30 against 8% measured, 33% off. **So "the rho and
     gamma knobs give different slopes" is not evidence against slaving. It is §39.2's
     published fact appearing on a third axis.** The first pass read a confirmation as a
     refutation.
  3. And the quantity is wrong: §39.1/§39.2 test `T_det/MFPT`, the TIME ratio. The first
     pass tested `Sigma_pred/Sigma_exact`, the COST ratio. They correlate at 0.95, not 1.

The corrected question is therefore **not** "do the axes share a coefficient" -- §39.2 says
they do not -- but **"does the 1/sep SCALING itself hold along rho, with its own
coefficient?"**

  P7  GATE, against a published number. At gamma = 0.20, rho = 1, sep = 7.00, §39.2's table
      gives `T_det/MFPT - 1 = +0.0914`. This run must reproduce it. If it does not, T_det or
      the MFPT is being computed differently from §39.1 and nothing else is admissible.
  P8  THE TEST. Along the rho axis, `(T_det/MFPT - 1) * sep` is CONSTANT within the axis,
      as §39.2 found on the T axis (0.6465, constant to 12.3% over 9x in sep). That is what
      the 1/sep scaling claim means. Cells whose gap is below 0.008 are at numerical
      resolution and are excluded AND reported, exactly as §39.2 did.
  P9  The rho-axis coefficient will NOT equal the T-axis 0.6465. §39.2 predicted
      non-transfer explicitly; this is a third axis on which to check it. All three
      coefficients reported side by side.
  P10 THE OUTCOME THAT WOULD MATTER. If `(T_det/MFPT - 1)*sep` is NOT constant along rho,
      then the 1/sep SCALING -- not merely its coefficient -- fails on a third axis, and
      §39.2's central claim needs qualifying. Given the first pass found resid*sep drifting
      4.4x on the rho axis for the COST ratio, **this is a live possibility and is the
      reason to run it on the right quantity.**
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.deterministic import jacobian
from crnl.networks.am_reversible import reverse_pairing
from experiments.arrhenius_optimum import am_rho, cost_cell, delta_star_rho
from experiments.cost_absolute import sigma_and_mu, sigma_predicted
from experiments.slaving_axis import sep_of, slaved


def eigs_at(net, delta):
    """(fast, slow) real eigenvalues on the slaved manifold at this delta."""
    st = slaved(net, delta)
    if st is None:
        return None
    ev = np.linalg.eigvals(jacobian(net, st)).real
    ev = np.sort(ev[np.abs(ev) > 1e-12])
    if len(ev) < 2:
        return None
    return float(ev[0]), float(ev[-1])


def path_separations(net, d_lo, d_hi, pairing, n=601):
    """Harmonic path separations under three weightings, plus the absolute eigenvalues.

    Harmonic is forced: the slow eigenvalue crosses zero mid-path, so sep diverges there
    and only <1/sep> is defined.
    """
    xs = np.linspace(d_lo, d_hi, n)
    inv, w_time, w_cost, fast, slow = [], [], [], [], []
    for x in xs:
        e = eigs_at(net, float(x))
        sm = sigma_and_mu(net, float(x), pairing)
        if e is None or sm is None or sm[1] <= 0:
            return None
        sig, mu = sm
        inv.append(abs(e[1] / e[0]))          # 1/sep = |slow/fast|
        w_time.append(1.0 / mu)
        w_cost.append(sig / mu)
        fast.append(e[0]); slow.append(e[1])
    inv = np.array(inv); wt = np.array(w_time); wc = np.array(w_cost)
    out = {}
    for name, w in (("path", np.ones_like(inv)), ("time", wt), ("cost", wc)):
        num = float(np.trapezoid(inv * w, xs))
        den = float(np.trapezoid(w, xs))
        out[f"inv_sep_{name}"] = num / den
        out[f"sep_{name}"] = den / num if num > 0 else float("inf")
    out["fast_mean"] = float(np.mean(fast))
    out["slow_absmax"] = float(np.max(np.abs(slow)))
    return out


def cell(gamma, rho, omega, eps, theta):
    ds = delta_star_rho(gamma, rho)
    if ds <= 0:
        return None
    net = am_rho(gamma, rho)
    pairing = reverse_pairing(net)
    ps = path_separations(net, eps * ds, theta * ds, pairing)
    if ps is None:
        return None
    per_om = sigma_predicted(net, eps * ds, theta * ds, pairing)
    c = cost_cell(gamma, rho, omega, eps, theta)
    if c is None or not np.isfinite(per_om):
        return None
    # §39.1/§39.2's quantity: the deterministic traversal time against the exact MFPT
    xs = np.linspace(eps * ds, theta * ds, 601)
    inv_mu = []
    for x in xs:
        sm = sigma_and_mu(net, float(x), pairing)
        if sm is None or sm[1] <= 0:
            return None
        inv_mu.append(1.0 / sm[1])
    t_det = float(np.trapezoid(inv_mu, xs))
    mfpt = c.get("mean_time", float("nan"))

    r = {"gamma": gamma, "rho": rho, "omega": omega, "delta_star": ds,
         "sep_point": float(sep_of(net)[0]),
         "ratio": float(omega * per_om / c["Sigma"]),
         "t_det": t_det, "mfpt": float(mfpt),
         "time_ratio": float(t_det / mfpt) if mfpt and np.isfinite(mfpt) else float("nan"),
         "exact": c["Sigma"], "L": c["L"]}
    r.update(ps)
    return r


CONV = ("point", "path", "time", "cost")


def sep_of_row(r, conv):
    return r["sep_point"] if conv == "point" else r[f"sep_{conv}"]


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
                    default=pathlib.Path("results/path_separation.json"))
    args = ap.parse_args()

    t0 = time.time()
    eps, th = args.eps, args.theta

    print(f"=== rho knob at gamma = {args.gamma}")
    print(f"{'rho':>6}{'sep_pt':>9}{'sep_path':>10}{'sep_time':>10}{'sep_cost':>10}"
          f"{'|resid|':>9}{'fast':>9}{'slow':>8}")
    rrows = []
    for r in args.rhos:
        cs = [cell(args.gamma, r, om, eps, th) for om in args.omegas]
        if any(c is None for c in cs):
            print(f"{r:>6.2f}   dropped")
            continue
        rrows += cs
        c = cs[0]
        res = float(np.mean([abs(x["ratio"] - 1.0) for x in cs]))
        print(f"{r:>6.2f}{c['sep_point']:>9.3f}{c['sep_path']:>10.3f}"
              f"{c['sep_time']:>10.3f}{c['sep_cost']:>10.3f}{res:>9.4f}"
              f"{c['fast_mean']:>9.3f}{c['slow_absmax']:>8.3f}")

    print(f"\n=== gamma knob at rho = 1")
    print(f"{'gamma':>6}{'sep_pt':>9}{'sep_path':>10}{'sep_time':>10}{'sep_cost':>10}"
          f"{'|resid|':>9}{'fast':>9}{'slow':>8}")
    grows = []
    for g in args.gammas:
        cs = [cell(g, 1.0, om, eps, th) for om in args.omegas]
        if any(c is None for c in cs):
            print(f"{g:>6.2f}   dropped")
            continue
        grows += cs
        c = cs[0]
        res = float(np.mean([abs(x["ratio"] - 1.0) for x in cs]))
        print(f"{g:>6.2f}{c['sep_point']:>9.3f}{c['sep_path']:>10.3f}"
              f"{c['sep_time']:>10.3f}{c['sep_cost']:>10.3f}{res:>9.4f}"
              f"{c['fast_mean']:>9.3f}{c['slow_absmax']:>8.3f}")

    ref = [c for c in rrows if abs(c["rho"] - 1.0) < 1e-12]
    print(f"\n=== P1 GATE: do the path conventions differ from the point value?")
    if ref:
        c = ref[0]
        for conv in ("path", "time", "cost"):
            d = 100 * abs(sep_of_row(c, conv) - c["sep_point"]) / c["sep_point"]
            print(f"  sep_{conv:>4} = {sep_of_row(c, conv):8.3f} vs point"
                  f" {c['sep_point']:.3f}   differ {d:6.1f}%"
                  f"   -> {'informative' if d > 1 else 'EMPTY, no content'}")

    print(f"\n=== P2: does any convention UN-MATCH the 8x pair (rho 0.5 vs rho 4)?")
    got = {c["rho"]: c for c in rrows if c["omega"] == args.omegas[0]}
    if 0.5 in got and 4.0 in got:
        a, b = got[0.5], got[4.0]
        ra, rb = abs(a["ratio"] - 1), abs(b["ratio"] - 1)
        print(f"  residuals: {ra:.4f} vs {rb:.4f}  (ratio {max(ra,rb)/min(ra,rb):.2f}x)")
        for conv in CONV:
            sa, sb = sep_of_row(a, conv), sep_of_row(b, conv)
            f = max(sa, sb) / min(sa, sb)
            # to explain a 3x residual gap via 1/sep, sep must differ by ~3x the OTHER way
            print(f"  sep_{conv:>5}: {sa:8.3f} vs {sb:8.3f}  differ {f:.2f}x"
                  f"   -> {'UN-MATCHED, could explain it' if f > 2 else 'still matched, explains nothing'}")

    print(f"\n=== P3: do the two knobs collapse under each convention?")
    print(f"{'convention':>12}{'rho slope':>12}{'gamma slope':>13}{'ratio':>8}{'verdict':>12}")
    verdicts = {}
    for conv in CONV:
        lo = max(min(sep_of_row(c, conv) for c in rrows),
                 min(sep_of_row(c, conv) for c in grows))
        hi = min(max(sep_of_row(c, conv) for c in rrows),
                 max(sep_of_row(c, conv) for c in grows))
        fits = {}
        for name, rows in (("rho", rrows), ("gamma", grows)):
            sel = [c for c in rows if lo - 1e-9 <= sep_of_row(c, conv) <= hi + 1e-9]
            if len(sel) < 2:
                fits[name] = float("nan")
                continue
            x = np.array([1.0 / sep_of_row(c, conv) for c in sel])
            y = np.array([c["ratio"] - 1.0 for c in sel])
            fits[name] = float(np.polyfit(x, y, 1)[0])
        if np.isfinite(fits["rho"]) and np.isfinite(fits["gamma"]) and fits["gamma"] != 0:
            fr = max(abs(fits["rho"]), abs(fits["gamma"])) / max(
                min(abs(fits["rho"]), abs(fits["gamma"])), 1e-12)
            ok = fr < 1.3
            verdicts[conv] = ok
            print(f"{conv:>12}{fits['rho']:>12.4f}{fits['gamma']:>13.4f}{fr:>8.2f}"
                  f"{'COLLAPSE' if ok else 'split':>12}")
        else:
            print(f"{conv:>12}   too few overlapping cells")

    print(f"\n=== P4: is residual * sep constant (absolute form)?")
    print(f"{'convention':>12}{'rho drift':>12}{'gamma drift':>13}")
    for conv in CONV:
        line = f"{conv:>12}"
        for rows in (rrows, grows):
            v = np.array([(c["ratio"] - 1.0) * sep_of_row(c, conv) for c in rows])
            pos = v[v > 0]
            line += f"{(pos.max()/pos.min() if pos.size > 1 else float('nan')):>12.2f}x"
        print(line)

    print(f"\n=== P5/P6: verdict")
    if verdicts and any(verdicts.values()):
        win = [k for k, v in verdicts.items() if v]
        print(f"  conventions that COLLAPSE the two knobs: {win}")
        if len(win) < len(verdicts):
            print("  -> P6: the conventions DISAGREE. The quantity is not resolved by this")
            print("     experiment and the collapse is reported as convention-dependent.")
    else:
        print("  NO convention collapses the two knobs.")
        print("  -> P5: slaving is REFUTED as rho's mechanism. T-COST-f closes negative:")
        print("     §44.2's lever keeps its measurement and has no account, and §39's")
        print("     residual returns to T-COST-c's other candidates. This was the")
        print("     predicted outcome -- both eigenvalues scale together with rho, so no")
        print("     ratio, point or path, can manufacture the 3x the residuals show.")

    # ================= SECOND PASS: the TIME ratio, §39.1/§39.2's quantity ==========
    print("\n" + "=" * 78)
    print("SECOND PASS: T_det/MFPT, which is what §39.1/§39.2 actually test")

    print(f"\n=== P7 GATE: §39.2's published cell -- gamma=0.20, rho=1, sep=7.00,"
          f" T_det/MFPT-1 = +0.0914")
    for c in [x for x in rrows if abs(x["rho"] - 1.0) < 1e-12]:
        print(f"  Omega={c['omega']}: T_det = {c['t_det']:.4f}  MFPT = {c['mfpt']:.4f}"
              f"  ratio-1 = {c['time_ratio']-1:+.4f}")

    print(f"\n=== P8/P9: is (T_det/MFPT - 1)*sep constant WITHIN each axis?")
    print("  §39.2 T axis: 0.6465 +- 12.3% over 9x in sep. Cells with gap < 0.008 are")
    print("  at numerical resolution and are excluded and reported, as §39.2 did.")
    coefs = {}
    for name, rows, key in (("rho", rrows, "rho"), ("gamma", grows, "gamma")):
        print(f"\n  --- {name} axis")
        print(f"{key:>8}{'sep':>9}{'gap':>10}{'gap*sep':>10}")
        keep, excl = [], []
        for c in sorted(rows, key=lambda z: z[key]):
            if c["omega"] != args.omegas[0]:
                continue
            gap = c["time_ratio"] - 1.0
            v = gap * c["sep_point"]
            print(f"{c[key]:>8.2f}{c['sep_point']:>9.3f}{gap:>10.4f}{v:>10.4f}"
                  + ("   [excluded, gap < 0.008]" if abs(gap) < 0.008 else ""))
            (excl if abs(gap) < 0.008 else keep).append(v)
        if len(keep) >= 3:
            a = np.array(keep)
            spread = 100 * (a.max() - a.min()) / abs(a.mean())
            coefs[name] = (float(a.mean()), spread, len(a))
            print(f"    coefficient = {a.mean():.4f}, spread {spread:.1f}%"
                  f" over {len(a)} cells"
                  f"  -> {'CONSTANT, 1/sep scaling holds' if spread < 30 else 'NOT constant, the 1/sep SCALING fails on this axis'}")
        if excl:
            print(f"    excluded {len(excl)} resolution-limited cells (reported)")

    print(f"\n=== P9: do the coefficients transfer? (§39.2 says NO, explicitly)")
    print(f"  {'T axis (§39.2)':>18}: 0.6465")
    for k, (m, s, n) in coefs.items():
        print(f"  {k + ' axis (this run)':>18}: {m:.4f}  (spread {s:.1f}%, {n} cells)")
    if "rho" in coefs:
        print(f"  -> rho vs T: {coefs['rho'][0]/0.6465:.2f}x."
              f"  Non-transfer is §39.2's published finding, not new evidence.")

    print(f"\n=== P10: verdict on the SCALING (not the coefficient)")
    for k, (m, s, n) in coefs.items():
        print(f"  {k:>6} axis: {'scaling HOLDS' if s < 30 else 'scaling FAILS'}"
              f"  (spread {s:.1f}% vs §39.2's 12.3% on the T axis)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"rho_rows": rrows, "gamma_rows": grows,
                                    "coefficients": coefs}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
