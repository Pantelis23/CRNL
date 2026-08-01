"""T16: does concatenated AM show a threshold, or does the ceiling survive?

FINDINGS 1's wall and 12.1's depth ceiling say a single AM stage has a fidelity
ceiling. The threshold theorem of fault-tolerant computing says the opposite about a
restoring code: below a threshold physical error rate, CONCATENATION buys arbitrary
fidelity. Nothing in this project has ever concatenated -- every result here is level-0.

THE CONSTRUCTION IS A POOL MERGE, AND THAT IS THE WHOLE DESIGN (rule 10). The obvious
version -- run k tanks and take the majority of their answers in numpy -- inserts a
free, noiseless, perfectly reliable majority gate, the exact class of error that has
already cost this project three withdrawn results. Here, k independent AM tanks of size
Omega each run to commitment, and their CONTENTS ARE PHYSICALLY COMBINED into one tank
of size k*Omega which then runs AM itself. The merged tank's initial margin is the sum
of the k committed margins -- positive iff a majority answered correctly -- and the
combining stage carries its own noise because it IS an AM tank. No sign(), no free
comparison, nothing the chemistry could not do.

Everything is computed EXACTLY, with no Monte Carlo anywhere:

    p_1 = sum_j  C(k,j) * p_0^j * (1-p_0)^(k-j) * p_merge(j)

where p_0 is the exact CME error of one Omega tank and p_merge(j) is the exact CME
error of a k*Omega tank started from the merge of j wrong and (k-j) right outputs.

TWO QUESTIONS, AND ONLY THE SECOND IS INFORMATIVE.

  (a) Does voting suppress error at all -- is p_1 < p_0 with a log-log exponent of 2?
      **This is close to definitional** and is reported as a check on the machinery, not
      as a finding: if p_merge is ~0 for j <= 1 and ~1 for j >= 2, then p_1 = 3p_0^2 -
      2p_0^3 identically and the exponent is 2 by construction. §24.2 flagged its own
      near-definitional arm the same way.

  (b) **Is voting the best use of k*Omega molecules?** The control is a SINGLE tank of
      size k*Omega started at the same RELATIVE margin -- the same molecules, not
      divided. This is the comparison that decides whether concatenation buys anything,
      and it has an answer that can be predicted in advance from §27's collapse:

          p_0(Omega) ~ exp(-Omega*c)   =>   3*p_0^2 ~ 3*exp(-2*Omega*c)
          p_pool(k*Omega)              ~        exp(-3*Omega*c)     at k = 3

      **Voting SQUARES the error; pooling CUBES the exponent.** So pooling should win,
      and by a margin that grows exponentially in Omega. That is a sharp, quantitative,
      falsifiable prediction rather than a direction.

PREDICTIONS, written before running:

  P1  (near-definitional, stated for completeness) p_1 < p_0, and a log-log fit of p_1
      against p_0 across an eps sweep gives an exponent near 2.
  P2  THE RESULT. p_1 > p_pool at every cell: concatenation is WORSE than pooling the
      same molecules into one tank. If instead p_1 < p_pool anywhere, my reading of the
      exponents is wrong and voting buys something pooling does not -- which would be
      the more interesting outcome and is why the control is computed at every cell.
  P3  THE ABSOLUTE TEST (rule 16, and the reason this is not merely a fitted trend).
      ln(p_1 / p_pool) should be LINEAR IN Omega with slope c, where c is the SAME
      collapse rate §27/§28 measured independently for this gamma and eps. Fitting the
      slope and comparing it against `-2*V_exact` -- a number this project already has,
      computed from closed forms with no free parameter -- is a prediction in absolute
      terms, not a shape check. If the fitted slope disagrees with the known collapse
      rate, the exponent story is wrong however good the fit looks.
  P4  p_merge(j) is ~0 for j <= 1 and ~1 for j >= 2. **If p_merge(1) is NOT negligible,
      the merged stage has its own error floor**, p_1 saturates, and the ceiling
      survives for a second and independent reason -- which would be a stronger result
      than P2 and must not be missed by only looking at the totals.
  P5  RULE 13, the approximation's own parameter as a second axis: the readout
      convention (level-0 output taken at the ATTRACTOR versus at the THRESHOLD) must
      not change the exponent. If it does, the answer is a property of the convention.

WHY THE EXPECTED OUTCOME IS THE INTERESTING ONE. In QEC there is no pool operation: the
physical error rate is fixed and cannot be lowered by using more of the same qubit, so
concatenation is the only lever available. Here error falls exponentially in Omega, so
a bigger tank is a lever -- and a better one. **If P2 holds, the difference between
chemical and quantum restoration is not that one restores and the other does not, but
that chemistry has a cheaper knob and therefore never needs the code.**
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time
from math import comb

import numpy as np

from crnl.cme import splitting_probability
from crnl.networks.am_reversible import am_reversible, delta_star, fixed_points
from experiments.collapse_slope_absolute import V_exact


def attractor_counts(gamma: float, omega: int, readout: str, theta: float):
    """Level-0 output state (n_X, n_Y, n_B) for a tank that committed to X."""
    att = [p for p in fixed_points(gamma) if p["kind"] == "attractor"]
    b = att[0]["b"]
    nb = int(round(b * omega))
    rest = omega - nb
    d = delta_star(gamma) if readout == "attractor" else theta * delta_star(gamma)
    sep = int(round(d * omega))
    if (rest - sep) % 2:
        sep -= 1
    nx, ny = (rest + sep) // 2, (rest - sep) // 2
    return np.array([nx, ny, nb], dtype=np.int64)


def start_at_margin(gamma: float, omega: int, eps_frac: float):
    """A tank started at relative margin eps_frac * delta*, X ahead.

    Returns the state and the REALISED eps, which is not the requested one: the
    integer lattice plus the parity fix move it by up to ~10%, and §27 measured that
    that wobble alone bounces raw local slopes in Omega by 40%. Every fit in Omega
    below is eps-controlled with these values, as in §27/§28.2/§28.3.
    """
    nb = int(round(omega * gamma / (1.0 + gamma)))
    rest = omega - nb
    d0 = max(1, int(round(eps_frac * delta_star(gamma) * omega)))
    if (rest - d0) % 2:
        d0 -= 1
    nx, ny = (rest + d0) // 2, (rest - d0) // 2
    realised = d0 / (delta_star(gamma) * omega)
    return np.array([nx, ny, nb], dtype=np.int64), float(realised)


def p_error(gamma: float, total: int, start, theta: float) -> float:
    """Exact CME P(the tank commits to Y | start). X is the correct answer."""
    net = am_reversible(gamma)
    thr = max(2, int(round(theta * delta_star(gamma) * total)))
    res = splitting_probability(
        net, total, float(total), start,
        lambda s: abs(int(s[0]) - int(s[1])) >= thr,
        lambda s: int(s[0]) > int(s[1]))
    return (1.0 - res["split"]) if res["valid"] else float("nan")


def concatenated(gamma: float, omega: int, k: int, eps_frac: float,
                 theta: float, readout: str) -> dict:
    n0, eps0 = start_at_margin(gamma, omega, eps_frac)
    p0 = p_error(gamma, omega, n0, theta)
    win = attractor_counts(gamma, omega, readout, theta)
    lose = np.array([win[1], win[0], win[2]], dtype=np.int64)   # mirror image
    merged, p1 = [], 0.0
    for j in range(k + 1):
        st = (k - j) * win + j * lose
        pm = p_error(gamma, k * omega, st, theta)
        w = comb(k, j) * p0 ** j * (1.0 - p0) ** (k - j)
        merged.append({"j": j, "start": st.tolist(), "weight": float(w),
                       "p_merge": float(pm)})
        p1 += w * pm
    npool, eps_pool = start_at_margin(gamma, k * omega, eps_frac)
    pool = p_error(gamma, k * omega, npool, theta)
    return {"omega": omega, "eps_frac": eps_frac, "readout": readout,
            "start": n0.tolist(), "p0": float(p0), "p1": float(p1),
            "p_pool": float(pool), "merged": merged,
            "eps_realised": eps0, "eps_pool_realised": eps_pool}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.25)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--eps-fracs", type=float, nargs="+",
                    default=[0.10, 0.15, 0.20, 0.25, 0.30, 0.35])
    ap.add_argument("--omega-fixed", type=int, default=40)
    ap.add_argument("--omegas", type=int, nargs="+", default=[20, 28, 36, 44, 52, 60])
    ap.add_argument("--eps-fixed", type=float, default=0.25)
    ap.add_argument("--readouts", type=str, nargs="+",
                    default=["attractor", "threshold"])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/concatenation.json"))
    args = ap.parse_args()

    t0 = time.time()
    g, k = args.gamma, args.k
    print(f"gamma={g} k={k} theta={args.theta}   pool merge, everything exact (no MC)")

    all_rows = {}
    for readout in args.readouts:
        print(f"\n=== readout convention: {readout}  (rule 13 -- second axis)")
        print(f"\n--- eps sweep at Omega={args.omega_fixed}  "
              f"(k*Omega = {k*args.omega_fixed})")
        print(f"{'eps':>6}{'p0':>12}{'p1 (vote)':>13}{'p_pool':>13}"
              f"{'p1/p_pool':>12}{'p1/p0':>10}")
        eps_rows = []
        for e in args.eps_fracs:
            r = concatenated(g, args.omega_fixed, k, e, args.theta, readout)
            eps_rows.append(r)
            print(f"{e:>6.2f}{r['p0']:>12.4e}{r['p1']:>13.4e}{r['p_pool']:>13.4e}"
                  f"{r['p1']/r['p_pool']:>12.3e}{r['p1']/r['p0']:>10.3e}")

        print(f"\n--- Omega sweep at eps={args.eps_fixed}")
        print(f"{'Omega':>6}{'p0':>12}{'p1 (vote)':>13}{'p_pool':>13}"
              f"{'p1/p_pool':>12}{'ln(p1/pool)':>13}")
        om_rows = []
        for om in args.omegas:
            r = concatenated(g, om, k, args.eps_fixed, args.theta, readout)
            om_rows.append(r)
            ratio = r["p1"] / r["p_pool"] if r["p_pool"] > 0 else float("nan")
            print(f"{om:>6}{r['p0']:>12.4e}{r['p1']:>13.4e}{r['p_pool']:>13.4e}"
                  f"{ratio:>12.3e}{np.log(ratio) if ratio > 0 else float('nan'):>13.4f}")
        all_rows[readout] = {"eps": eps_rows, "omega": om_rows}

        print(f"\n  P1 (near-definitional): exponent of p1 vs p0")
        x = np.log([r["p0"] for r in eps_rows]); y = np.log([r["p1"] for r in eps_rows])
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() >= 3:
            c = np.polyfit(x[ok], y[ok], 1)
            res = y[ok] - np.polyval(c, x[ok])
            print(f"    p1 ~ p0^{c[0]:.4f}   R^2 = {1 - res.var()/y[ok].var():.5f}"
                  f"   (2 expected)")

        print(f"  P4: is the merged stage itself reliable?")
        mid = eps_rows[len(eps_rows) // 2]
        for m in mid["merged"]:
            print(f"    j={m['j']} (of {k} wrong)  start {str(m['start']):>16}  "
                  f"p_merge = {m['p_merge']:.6e}   weight {m['weight']:.3e}")

        print(f"  P2/P3: does voting beat pooling, and by the known collapse rate?")
        om = np.array([r["omega"] for r in om_rows], float)
        lr = np.array([np.log(r["p1"] / r["p_pool"]) if r["p_pool"] > 0 else np.nan
                       for r in om_rows])
        good = np.isfinite(lr)
        beats = [r["omega"] for r in om_rows if r["p1"] < r["p_pool"]]
        print(f"    cells where voting beats pooling: {len(beats)}/{len(om_rows)}"
              + (f" -> Omega {beats}" if beats else " (none)"))
        if good.sum() >= 4:
            net = am_reversible(g)
            predicted = -2.0 * V_exact(net, args.eps_fixed * delta_star(g))
            raw = np.polyfit(om[good], lr[good], 1)
            rres = lr[good] - np.polyval(raw, om[good])
            print(f"    RAW slope = {raw[0]:.6f}   R^2 = "
                  f"{1 - rres.var()/lr[good].var():.5f}   <- contaminated by the "
                  f"realised-eps wobble (§27)")
            # eps-controlled: both realised eps enter, with different sensitivities
            # because ln p1 ~ 2 ln p0 while ln p_pool ~ ln p_pool.
            e0 = np.array([r["eps_realised"] for r in om_rows])[good]
            ep = np.array([r["eps_pool_realised"] for r in om_rows])[good]
            A = np.vstack([om[good], e0 - e0.mean(), ep - ep.mean(),
                           np.ones(good.sum())]).T
            c, *_ = np.linalg.lstsq(A, lr[good], rcond=None)
            res = lr[good] - A @ c
            r2 = 1 - res.var() / lr[good].var()
            print(f"    eps-CONTROLLED slope in Omega = {c[0]:.6f}   R^2 = {r2:.5f}")
            print(f"    independently-known collapse rate -2*V_exact = "
                  f"{predicted:.6f}  -> |slope| / |rate| = "
                  f"{abs(c[0])/abs(predicted):.4f}")
            print(f"    realised eps ranged {e0.min():.4f}-{e0.max():.4f} at level 0 "
                  f"and {ep.min():.4f}-{ep.max():.4f} in the pool")
            all_rows[readout]["slope"] = float(c[0])
            all_rows[readout]["slope_raw"] = float(raw[0])
            all_rows[readout]["slope_r2"] = float(r2)
            all_rows[readout]["predicted_rate"] = float(predicted)

    if len(args.readouts) > 1:
        print(f"\n=== P5 (rule 13): does the readout convention change the answer?")
        for readout in args.readouts:
            s = all_rows[readout].get("slope")
            print(f"  {readout:>10}: ln(p1/p_pool) slope = "
                  + (f"{s:.6f}" if s is not None else "n/a"))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"gamma": g, "k": k, "theta": args.theta,
                                    "results": all_rows}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
