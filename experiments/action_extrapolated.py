"""T15-m: is §63.2's nu = 1.95 real, or is it an unconverged estimator? — a self-audit

§63.2 reported nu = 1.9496 +- 0.0026 with the pitchfork's 2 excluded at 9.7x the scatter,
using A_eff = the local slope of -ln|lambda_A| in Omega, gated on the DRIFT between the last
two slopes being under 2%. **The gate was the wrong quantity.** If A_eff(Wbar) = A + c/Wbar,
the drift between consecutive estimates is c(1/Wbar1 - 1/Wbar2) while the REMAINING error is
c/Wbar2, and for these Omega ladders the second is **4x the first** -- measured, not assumed:

    gamma    0.38   0.41   0.43   0.45   0.46
    drift    0.03%  0.12%  0.20%  1.49%  4.33%
    remaining 0.14%  0.48%  0.79%  5.95% 17.31%

So §63.2 admitted a point carrying ~6% error at gamma = 0.45 and ~0.1% at gamma = 0.38.
**That bias is systematic in gamma and signed** -- A_eff approaches A from above, and the
excess grows toward gamma_c -- which flattens ln A against ln(gamma_c - gamma) and therefore
LOWERS the fitted nu. §63.2's 1.95 is exactly what that bias would produce.

And the kill test T15-m shipped with is not reachable. It proposed pushing to
gamma = 0.47-0.49 directly; at gamma = 0.47 the slope still drifts 4.2% at Omega = 1000
(17.8s, 5e5 states), so 1% needs Omega ~ 4000 and ~8e6 states. **That estimate was made from
wall-clock at Omega = 640 without asking what Omega the convergence actually required**, and
it is withdrawn here rather than left standing.

PREDICTIONS, written before running.

  P1  GATE. A_eff must be LINEAR in 1/Wbar, which is the correction form the whole repair
      assumes. Report the fit residual per gamma; any gamma whose A_eff is not linear in
      1/Wbar to 1% of its range is excluded and counted, because then the extrapolation is
      not justified and a second-order form is being papered over.
  P2  **THE TEST, and its DIRECTION is predicted, which is the strong form.** After
      extrapolating A_eff -> A, **nu must move UP** from 1.9496, because the bias it removes
      is positive and grows toward gamma_c. A nu that moves DOWN or stays put refutes the
      diagnosis above, and then §63.2's 1.95 needs a different explanation and stands
      meanwhile.
  P3  **PREDICTED VALUE: nu = 2.02 +- 0.03.** Estimated before running from the measured
      remaining errors: +0.06 in ln A at gamma = 0.45 against +0.001 at gamma = 0.38, over
      a span of 0.88 in ln(gamma_c - gamma), biases the slope by -0.067. So 1.9496 + 0.067.
      **This is a number to be hit or missed, not a direction to be confirmed.**
  P4  **VERDICT (rule 19), and it must be able to print three things.** (a) nu consistent
      with 2 within the scatter across nested windows -- the pitchfork, and §63.2 withdrawn;
      (b) nu flat at some other value with 2 outside the scatter -- §63.2's conclusion
      survives with a corrected number; (c) nu still drifting with the window -- unresolved,
      and neither is claimed. The criterion compares trend against scatter, as §63's
      corrected rule does, and (c) is a real outcome rather than a failure to be rounded away.
  P5  **RULE 9, an axis I did not choose.** Re-extrapolate from a DISJOINT Omega ladder. If
      A depends on which Omega were used to reach it, the extrapolation is not converged
      either and P2's answer is not entitled to three digits.
  SECOND PASS, before re-running, on two verdict rules of my own that misfired (rule 3 keeps
  them; rule 19 is why they are here at all):

    (i) P1's residual was normalised by the RANGE of the four A_eff being fitted. That range
        shrinks toward zero as convergence improves, so the criterion punishes exactly the
        gamma it should pass: gamma = 0.43 shifts by 0.62% and scored 0.41, gamma = 0.46
        shifts by 25% and scored 0.019. It rejected all 8. **A residual must be normalised
        by the quantity being estimated, not by the spread of the estimates**, so it is now
        max|pred - A_eff| / |A|.
    (ii) The verdict helper could not print outcome (c) at all. It tested
        `trend > scatter` where trend = |v_first - 2| - |v_last - 2|; for any sequence
        MONOTONE toward 2 those are identically equal, so a drifting sequence always fell
        through to "(a) nu = 2 within the scatter" -- and it duly did, on 1.8522, 1.9857,
        1.9909. A rule with an unreachable branch is not a three-way test. Replaced by:
        estimate = the narrowest window's nu, uncertainty u = |v_last - v_prev|; (a) if
        |v_last - 2| <= u; else (c) if the sequence is monotone toward 2; else (b). Checked
        against §63.2's own numbers, where it reproduces that section's verdict (b).

  P6  **RULE 14: a correction is a claim, so it gets an INDEPENDENT instrument.** The action
      is also the quasipotential barrier, and -ln(pi)/Omega -> V for the CME whether or not
      detailed balance holds. So compute A from the STATIONARY DISTRIBUTION -- the log-ratio
      of pi at the rail to pi on the symmetric line, divided by Omega -- which involves no
      eigenvalue, no antisymmetric block and no local slope. If the two disagree on nu, the
      correction is not established and says so.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.cme import enumerate_states, stationary
from crnl.networks.am_reversible import am_reversible
from experiments.threshold_sharpness import lambda_A

GAMMA_C = 0.5


def action_series(gamma, omegas, floor=1e-12):
    """A_eff between consecutive Omega, and the Wbar each belongs to."""
    oms, lns = [], []
    for om in omegas:
        lam = lambda_A(am_reversible(float(gamma)), om)
        if lam is None or lam >= 0 or abs(lam) < floor:
            continue
        oms.append(om)
        lns.append(-np.log(abs(lam)))
    if len(oms) < 3:
        return None
    a = [(lns[i + 1] - lns[i]) / (oms[i + 1] - oms[i]) for i in range(len(oms) - 1)]
    wb = [(oms[i] + oms[i + 1]) / 2.0 for i in range(len(oms) - 1)]
    return np.array(wb), np.array(a), oms


def extrapolate(wb, a, n_last=4):
    """A_eff = A + c/Wbar, fitted on the last n_last points. Returns (A, resid_frac)."""
    w, y = wb[-n_last:], a[-n_last:]
    if len(w) < 3:
        return None, None
    c, A = np.polyfit(1.0 / w, y, 1)
    pred = A + c / w
    # relative to the estimate, NOT to the spread of the points being fitted
    return float(A), float(np.abs(pred - y).max() / max(abs(A), 1e-30))


def action_from_stationary(gamma, omega):
    """P6: the quasipotential barrier from pi alone -- no eigenvalue, no block."""
    net = am_reversible(float(gamma))
    states, index = enumerate_states(3, int(omega))
    pi = stationary(net, int(omega), float(omega))
    d = np.array([int(s[0]) - int(s[1]) for s in states])
    on_sym = d == 0
    if not on_sym.any():
        return None
    peak = float(np.log(pi.max()))
    saddle = float(np.log(pi[on_sym].max()))
    return (peak - saddle) / omega


def fit_nu(pts, lo):
    sel = [(g, A) for g, A in pts if g >= lo]
    if len(sel) < 3:
        return None
    x = np.log(GAMMA_C - np.array([g for g, _ in sel]))
    y = np.log(np.array([A for _, A in sel]))
    nu, c = np.polyfit(x, y, 1)
    return float(nu), float(np.exp(c)), len(sel)


def verdict(nus, label):
    """Three-way, and all three branches must be reachable (see the SECOND PASS note)."""
    v = np.array([n for _, n in nus])
    est = float(v[-1])                       # the narrowest window
    u = float(abs(v[-1] - v[-2]))            # how much it is still moving there
    gap = abs(est - 2.0)
    mono = bool(np.all(np.diff(np.abs(v - 2.0)) <= 1e-12))
    print(f"  {label}: " + ", ".join(f"{n:.4f}" for _, n in nus))
    print(f"    narrowest window {est:.4f}, still moving by {u:.4f},"
          f" |nu - 2| = {gap:.4f} ({gap/max(u,1e-9):.1f}x that movement),"
          f" {'monotone toward 2' if mono else 'not monotone'}")
    if gap <= u:
        out = "a"
        print(f"    -> (a) nu = 2 within its own residual movement: THE PITCHFORK,"
              f" and §63.2's 1.9496 is withdrawn")
    elif mono:
        out = "c"
        print(f"    -> (c) nu is DRIFTING toward 2, at {est:.4f} and still moving:"
              f" consistent with 2 but NOT arrived; nothing is rounded")
    else:
        out = "b"
        print(f"    -> (b) nu FLAT at {est:.4f}, 2 excluded at {gap/max(u,1e-9):.1f}x"
              f" the residual movement: §63.2's conclusion survives")
    return {"values": [float(n) for _, n in nus], "est": est, "u": u, "gap": gap,
            "monotone": mono, "outcome": out}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.24, 0.28, 0.32, 0.35, 0.38, 0.41, 0.43, 0.45, 0.46])
    ap.add_argument("--ladderA", type=int, nargs="+",
                    default=[40, 60, 90, 120, 160, 200, 260, 320, 400, 500])
    ap.add_argument("--ladderB", type=int, nargs="+",
                    default=[50, 75, 110, 150, 210, 290, 380, 480])
    ap.add_argument("--statomegas", type=int, nargs="+", default=[60, 100, 150, 220])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/action_extrapolated.json"))
    args = ap.parse_args()
    t0 = time.time()

    print("=== P1 GATE: is A_eff linear in 1/Wbar, as the extrapolation assumes?")
    print(f"{'gamma':>7}{'A_eff(last)':>13}{'A extrap':>12}{'shift%':>9}"
          f"{'resid/A':>13}{'used':>6}")
    ptsA, rows = [], []
    for g in args.gammas:
        r = action_series(g, args.ladderA)
        if r is None:
            print(f"{g:>7.3f}   too few resolvable Omega")
            continue
        wb, a, oms = r
        A, resid = extrapolate(wb, a)
        if A is None:
            continue
        shift = (a[-1] - A) / max(abs(A), 1e-30)
        used = resid is not None and resid < 0.02 and A > 0
        if used:
            ptsA.append((g, A))
        rows.append({"gamma": g, "A_last": float(a[-1]), "A": A, "resid": resid,
                     "shift": float(shift), "used": bool(used),
                     "omegas": oms, "slopes": a.tolist()})
        print(f"{g:>7.3f}{a[-1]:>13.7f}{A:>12.7f}{100*shift:>9.2f}"
              f"{resid:>13.4f}{'yes' if used else 'NO':>6}")
    print(f"  -> P1: {len(ptsA)} of {len(rows)} gamma usable")

    print(f"\n=== P2/P3: nu after extrapolation (§63.2 had 1.9496; predicted 2.02 +- 0.03)")
    print(f"{'window':>22}{'n':>4}{'nu':>10}{'amplitude':>12}")
    nusA = []
    for lo in (0.24, 0.30, 0.35, 0.38, 0.41):
        f = fit_nu(ptsA, lo)
        if f is None:
            continue
        nu, amp, n = f
        nusA.append((lo, nu))
        hi = max(g for g, _ in ptsA)
        print(f"{f'[{lo:.2f}, {hi:.2f}]':>22}{n:>4}{nu:>10.4f}{amp:>12.4f}")
    resA = verdict(nusA, "nu across nested windows") if len(nusA) >= 3 else None
    if nusA:
        mv = np.mean([n for _, n in nusA])
        print(f"  P2 direction: nu moved {mv - 1.9496:+.4f} from §63.2's 1.9496"
              f" -> {'UP, as predicted' if mv > 1.9496 else 'DOWN, REFUTING the diagnosis'}")
        print(f"  P3 absolute: predicted 2.02 +- 0.03, got {mv:.4f}"
              f"  -> {'HIT' if abs(mv - 2.02) <= 0.03 else 'MISSED'}")

    print(f"\n=== P5 (rule 9): the same thing from a DISJOINT Omega ladder")
    ptsB = []
    for g in args.gammas:
        r = action_series(g, args.ladderB)
        if r is None:
            continue
        A, resid = extrapolate(r[0], r[1])
        if A is not None and resid is not None and resid < 0.02 and A > 0:
            ptsB.append((g, A))
    common = [g for g, _ in ptsA if any(abs(g - h) < 1e-12 for h, _ in ptsB)]
    print(f"  ladder A {args.ladderA[0]}..{args.ladderA[-1]},"
          f" ladder B {args.ladderB[0]}..{args.ladderB[-1]};"
          f" {len(common)} gamma in both")
    if common:
        da = dict(ptsA); db = dict(ptsB)
        worst = max(abs(da[g] - db[g]) / abs(da[g]) for g in common)
        print(f"  worst relative disagreement in A: {100*worst:.2f}%")
        nusB = [(lo, fit_nu(ptsB, lo)[0]) for lo in (0.24, 0.30, 0.35, 0.38, 0.41)
                if fit_nu(ptsB, lo)]
        if len(nusB) >= 3:
            verdict(nusB, "nu from ladder B")

    print(f"\n=== P6 (rule 14): the action from the STATIONARY DISTRIBUTION alone")
    print(f"{'gamma':>7}" + "".join(f"{f'W={o}':>12}" for o in args.statomegas))
    ptsS = []
    for g in args.gammas:
        vals = []
        for om in args.statomegas:
            try:
                vals.append(action_from_stationary(g, om))
            except Exception:
                vals.append(None)
        good = [(om, v) for om, v in zip(args.statomegas, vals) if v and v > 0]
        print(f"{g:>7.3f}" + "".join(f"{v:>12.6f}" if v else f"{'--':>12}" for v in vals))
        if len(good) >= 3:
            o = np.array([x for x, _ in good], float)
            y = np.array([v for _, v in good], float)
            c, A = np.polyfit(1.0 / o, y, 1)
            if A > 0:
                ptsS.append((g, float(A)))
    nusS = [(lo, fit_nu(ptsS, lo)[0]) for lo in (0.24, 0.30, 0.35)
            if fit_nu(ptsS, lo)]
    if len(nusS) >= 3:
        verdict(nusS, "nu from pi alone")
        vS = np.mean([n for _, n in nusS])
        vA = np.mean([n for _, n in nusA]) if nusA else float("nan")
        print(f"  eigenvalue route {vA:.4f} vs stationary route {vS:.4f}"
              f"  -> {'AGREE, the correction is established' if abs(vA-vS) < 0.08 else 'DISAGREE: the correction is NOT established (rule 14)'}")
    else:
        print("  too few gamma with a converged stationary action")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(
        {"rows": rows, "nu_A": nusA, "nu_B": locals().get("nusB", []),
         "nu_stat": nusS, "verdictA": resA,
         "pts": {"A": ptsA, "B": ptsB, "stat": ptsS}}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
