"""T16-b: a closed form for §33's crossover, tested in absolute terms.

§33 measured that re-merging beats the single-tank hold below a crossover in Omega
that moves with tau, and located it at Omega ~ 19 for k = 3 at tau = t_relax. That was
a measured boundary with no theory. The exponent count that §33 already validated
(P5's integer structure) determines it, and the derivation is short.

With m = ceil((k+1)/2) the merged round fails when m of k flip, so for tau << T

    L_hold(Omega)    = T(k*Omega)
    L_remerge(Omega) = T(Omega)^m / ( C(k,m) * tau^(m-1) )

and writing ln T(N) = c*N + a, the two are equal at

    **Omega_x = [ (m-1)*(a - ln tau) - ln C(k,m) ] / [ c * (k - m) ]**

**FOR EVERY ODD k, m-1 = k-m = (k-1)/2 EXACTLY**, so the ratio (m-1)/(k-m) is 1 and

    Omega_x = (a - ln tau)/c  -  ln C(k,m) / ( c*(k-1)/2 )

**The leading term contains no k at all.** That is a structural prediction and a
surprising one, because the win MARGIN at fixed Omega differs by more than 2x between
k = 3 and k = 5 (§33: 0.482 against 0.269 at Omega = 14) -- the margins differ, the
crossings nearly coincide.

WHY THIS IS AN ABSOLUTE TEST AND NOT A FIT. `c` and `a` are obtained from the HOLD
protocol alone -- a straight-line fit to ln T(N) over the exactly-solvable sizes. **No
crossover measurement enters the prediction at any point**, and the combinatorial term
is fixed by k. §30.1 is the standing reminder of why this matters: a formula that
reproduces the sweep it was built from proves nothing, and that one died to the first
sweep aimed at it. In-sample at gamma = 0.30 the form already reproduces §33 to 0.8%
(k=3) and 0.1% (k=5) against ratios extrapolated by an independent route; this
experiment asks whether it holds where it was never fitted.

THE CROSSOVER IS MEASURED CONTINUOUSLY, by interpolating ln(L_hold/L_remerge) to its
zero -- **not** as "the largest integer Omega at which re-merge still wins". That
quantised version is failure pattern 2 in THEORIES §4, which has already cost this
project a result twice (0.474796 and 0.495227 rounding identically), and §33's own
table reports exactly that quantised form. Quote the continuous quantity; round only
at display.

PREDICTIONS, written before running:

  P1  The closed form holds with the odd-k identity, so the predicted crossover is
      nearly k-INDEPENDENT: spread across k = 3, 5, 7 under about 5% at fixed gamma
      and tau, despite the win margins differing by 2x or more.
  P2  THE ABSOLUTE TEST. Predicted Omega_x agrees with the measured continuous
      crossover within about 10% at every reachable (gamma, k). **A known systematic
      is expected and its DIRECTION is predicted:** ln T(N) is not exactly linear --
      Kramers carries a power-law prefactor -- so a straight-line fit absorbs curvature
      into `a`. §33 saw this as slope ratios running ~8% high. I expect a bias of that
      order rather than exact agreement, and a bias that is consistent in sign across
      gamma is evidence for the prefactor rather than against the form.
  P3  d(Omega_x)/d(ln tau) = -1/c EXACTLY, with no k and no combinatorics in it. Fitted
      over a decade in tau it must match -1/c computed independently from the same
      T(N) sweep. This is the cleanest single number in the section because both sides
      come from different places.
  P4  The renewal model's own validity (std/mean ~ 1, §33's P2) must hold at every cell
      used, INCLUDING the small-Omega cells where the committed state is marginal. If
      it degrades there, those cells are inadmissible and the usable gamma range
      shrinks -- reported, not quietly dropped.
  P5  Cells whose crossover needs N = k*Omega_x beyond the exactly-solvable range are
      SKIPPED AND REPORTED. Larger gamma means a shallower barrier, smaller c, and a
      LARGER Omega_x, so the high-gamma end should fall out of reach first -- the
      prediction of which cells fail is itself part of the test.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time
from math import comb, ceil, log

import numpy as np

from crnl.cme import first_passage_moments
from crnl.deterministic import jacobian
from crnl.networks.am_reversible import am_reversible, delta_star, fixed_points


def committed_state(gamma: float, total: int) -> np.ndarray:
    att = [p for p in fixed_points(gamma) if p["kind"] == "attractor"][0]
    nb = int(round(att["b"] * total))
    rest = total - nb
    sep = int(round(delta_star(gamma) * total))
    if (rest - sep) % 2:
        sep -= 1
    return np.array([(rest + sep) // 2, (rest - sep) // 2, nb], dtype=np.int64)


def relax_time(gamma: float) -> float:
    att = [p for p in fixed_points(gamma) if p["kind"] == "attractor"][0]
    x = np.array([att["x"], att["y"], att["b"]], dtype=float)
    ev = np.linalg.eigvals(jacobian(am_reversible(gamma), x)).real
    ev = np.sort(ev[np.abs(ev) > 1e-9])
    return float(1.0 / abs(ev[-1]))


def lifetime(gamma: float, total: int, theta: float) -> dict:
    net = am_reversible(gamma)
    thr = max(2, int(round(theta * delta_star(gamma) * total)))
    r = first_passage_moments(net, total, float(total), committed_state(gamma, total),
                              lambda s: int(s[1]) - int(s[0]) >= thr)
    return {"mean": r["mean_time"], "std": r["std_time"], "valid": bool(r["valid"])}


def remerge(T: float, k: int, tau: float) -> float:
    m = ceil((k + 1) / 2)
    q = 1.0 - np.exp(-tau / T)
    pf = sum(comb(k, j) * q ** j * (1.0 - q) ** (k - j) for j in range(m, k + 1))
    return tau / pf if pf > 0 else float("inf")


def predicted_crossover(c: float, a: float, k: int, tau: float) -> float:
    m = ceil((k + 1) / 2)
    return ((m - 1) * (a - log(tau)) - log(comb(k, m))) / (c * (k - m))


def measured_crossover(life: dict, k: int, tau: float, omegas) -> float | None:
    """Continuous zero of ln(L_hold / L_remerge), never the largest integer winner."""
    xs, ys = [], []
    for om in omegas:
        s, h = life.get(om), life.get(k * om)
        if not (s and h and s["valid"] and h["valid"] and s["mean"] and h["mean"]):
            continue
        xs.append(float(om))
        ys.append(float(np.log(h["mean"] / remerge(s["mean"], k, tau))))
    if len(xs) < 2:
        return None
    xs, ys = np.array(xs), np.array(ys)
    sign = np.sign(ys)
    idx = np.where(np.diff(sign) != 0)[0]
    if not len(idx):
        return None
    i = idx[0]
    return float(xs[i] - ys[i] * (xs[i + 1] - xs[i]) / (ys[i + 1] - ys[i]))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.15, 0.20, 0.25, 0.30, 0.35])
    ap.add_argument("--ks", type=int, nargs="+", default=[3, 5, 7])
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--omegas", type=int, nargs="+",
                    default=list(range(6, 37, 2)))
    ap.add_argument("--tau-mults", type=float, nargs="+",
                    default=[0.5, 1.0, 2.0, 3.0, 5.0, 8.0])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/crossover_law.json"))
    args = ap.parse_args()

    t0 = time.time()
    print("Omega_x = [(m-1)(a - ln tau) - ln C(k,m)] / [c (k-m)],  "
          "c and a from HOLD data only")
    out = {}
    for g in args.gammas:
        tr = relax_time(g)
        need = sorted({om for om in args.omegas}
                      | {k * om for k in args.ks for om in args.omegas})
        life = {N: lifetime(g, N, args.theta) for N in need if N <= 260}
        ok = {N: v for N, v in life.items()
              if v["valid"] and v["mean"] and np.isfinite(v["mean"])}
        if len(ok) < 4:
            print(f"\ngamma={g}: only {len(ok)} trustworthy solves -- SKIPPED")
            continue
        N = np.array(sorted(ok)); T = np.array([ok[n]["mean"] for n in N])
        c, a = np.polyfit(N, np.log(T), 1)
        r2 = 1 - np.var(np.log(T) - (c * N + a)) / np.var(np.log(T))
        sm = np.array([ok[n]["std"] / ok[n]["mean"] for n in N])
        print(f"\n=== gamma = {g}   t_relax = {tr:.4f}")
        print(f"  hold-only fit over N = {N.min()}..{N.max()} ({len(N)} sizes): "
              f"c = {c:.6f}  a = {a:.6f}  R^2 = {r2:.5f}")
        print(f"  P4 renewal validity: std/mean over those sizes = "
              f"{sm.min():.4f}-{sm.max():.4f}"
              f"   -> {'ok' if abs(sm - 1).max() < 0.10 else 'DEGRADED'}")
        rows = []
        print(f"  {'k':>3}{'tau/trlx':>10}{'predicted':>12}{'measured':>11}"
              f"{'pred/meas':>11}")
        for k in args.ks:
            for mult in args.tau_mults:
                tau = mult * tr
                pred = predicted_crossover(c, a, k, tau)
                meas = measured_crossover(ok, k, tau, args.omegas)
                rows.append({"k": k, "tau_mult": mult, "tau": tau,
                             "predicted": pred, "measured": meas})
                ms = f"{meas:.2f}" if meas is not None else "unreached"
                rt = f"{pred/meas:.4f}" if meas else "-"
                print(f"  {k:>3}{mult:>10.1f}{pred:>12.2f}{ms:>11}{rt:>11}")
        out[str(g)] = {"c": float(c), "a": float(a), "r2": float(r2),
                       "t_relax": tr, "rows": rows,
                       "std_over_mean": [float(x) for x in sm]}

        got = [r for r in rows if r["measured"] is not None]
        if got:
            rr = np.array([r["predicted"] / r["measured"] for r in got])
            print(f"  P2 absolute: pred/meas over {len(got)} reachable cells = "
                  f"{rr.min():.4f}-{rr.max():.4f}, mean {rr.mean():.4f}")
        miss = len(rows) - len(got)
        if miss:
            print(f"  P5: {miss} of {len(rows)} cells never crossed in the Omega grid "
                  f"(reported, not dropped)")

        print(f"  P1 k-independence at tau = t_relax: ", end="")
        pk = [predicted_crossover(c, a, k, tr) for k in args.ks]
        print(" ".join(f"k={k}:{p:.2f}" for k, p in zip(args.ks, pk))
              + f"   spread {100*(max(pk)-min(pk))/np.mean(pk):.2f}%")

        print(f"  P3 d(Omega_x)/d(ln tau) must be -1/c = {-1/c:.4f}: ", end="")
        k0 = args.ks[0]
        lt = np.array([log(m * tr) for m in args.tau_mults])
        px = np.array([predicted_crossover(c, a, k0, m * tr) for m in args.tau_mults])
        s_pred = np.polyfit(lt, px, 1)[0]
        mx = [(log(m * tr), measured_crossover(ok, k0, m * tr, args.omegas))
              for m in args.tau_mults]
        mx = [(x, y) for x, y in mx if y is not None]
        if len(mx) >= 3:
            s_meas = np.polyfit([x for x, _ in mx], [y for _, y in mx], 1)[0]
            print(f"predicted-form {s_pred:.4f}, MEASURED {s_meas:.4f}"
                  f"   ratio {s_meas/(-1/c):.4f}")
            out[str(g)]["dOm_dlntau_measured"] = float(s_meas)
        else:
            print(f"predicted-form {s_pred:.4f}, too few reachable cells to measure")

    print(f"\n=== P2 across all gamma (the absolute test)")
    allr = []
    for g, v in out.items():
        rr = [r["predicted"] / r["measured"] for r in v["rows"] if r["measured"]]
        if rr:
            allr += rr
            print(f"  gamma={g}: {len(rr)} cells, pred/meas mean "
                  f"{np.mean(rr):.4f}  ({min(rr):.4f}-{max(rr):.4f})")
    if allr:
        allr = np.array(allr)
        print(f"  OVERALL: {len(allr)} cells, mean {allr.mean():.4f}, "
              f"sd {allr.std(ddof=1):.4f}, range {allr.min():.4f}-{allr.max():.4f}")
        print(f"  bias is {'CONSISTENT in sign' if (allr > 1).all() or (allr < 1).all() else 'NOT consistent in sign'}"
              f" -> {'supports the Kramers-prefactor reading (P2)' if (allr > 1).all() or (allr < 1).all() else 'not the predicted systematic'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
