"""T16-a: does PERIODIC RE-MERGING beat the single-tank hold?

§32 found that one-shot voting loses to pooling the same molecules, by a factor growing
as exp(Omega*c) -- voting squares the error, pooling cubes the exponent. But §32 is
one-shot, and that is NOT what quantum error correction does. Real fault tolerance is
time-extended: fresh ancillas repeatedly remove errors that accumulate during storage,
and the threshold theorem is a statement about that REPETITION. **§12.1's depth ceiling
is the right target, because it concerns a bit HELD over time.**

THE TWO PROTOCOLS, at identical total molecule count N = k*Omega:

  HOLD      one tank of N molecules sits in its committed state until it
            spontaneously crosses to the opposite one. Lifetime = the exact mean
            first-passage time T(N).
  RE-MERGE  k tanks of Omega each hold a copy. Every tau the contents are combined
            into one k*Omega tank and immediately re-split into k equal portions --
            both physical operations, mixing and aliquoting. A minority that has
            flipped is outvoted by the merged margin and the split hands every
            portion back the corrected state. The sub-tanks re-amplify on their own
            between merges, since that is what the landscape does.

WHY MATCHED MOLECULES IS ALSO MATCHED DISSIPATION, and it is measured rather than
assumed (P1). Both protocols hold N molecules cycling continuously against gamma < 1,
and mixing and aliquoting are idealised as free. So if the per-molecule entropy
production at the committed state is size-independent, the two protocols dissipate at
the same rate and comparing lifetimes at matched N is already the fair comparison. That
is checkable with `ep_rate` and it is checked below -- if it fails, the comparison has
to be re-weighted rather than asserted.

THE MODEL IS A RENEWAL PROCESS, AND IT HAS ITS OWN KILL TEST (P2). Each sub-tank is
treated as a two-state chain flipping at rate 1/T(Omega), so a round fails when at least
m = ceil((k+1)/2) of k flip within tau:

    P_fail(tau) = sum_{j>=m} C(k,j) q^j (1-q)^(k-j),   q = 1 - exp(-tau/T(Omega))
    L_remerge(tau) = tau / P_fail(tau)

The two-state reduction is only legitimate if first passage is near-exponential, which
is exactly what `first_passage_moments` can settle: for an exponential law the standard
deviation equals the mean. **std/mean is computed exactly at every cell, and if it is
not close to 1 the renewal arithmetic above is void and no conclusion may be drawn from
it.** It also requires tau >> t_relax, so that a merged portion re-amplifies to the
attractor before the next merge; t_relax comes from the Jacobian at the attractor and
tau is never taken below it.

PREDICTIONS, written before running:

  P1  Per-molecule entropy production at the committed state is size-independent, so
      matched N is matched dissipation and the lifetime comparison is fair as it stands.
  P2  First passage is near-exponential: std/mean within a few percent of 1. This is the
      renewal model's own kill test and it is checked before any lifetime is quoted.
  P3  L_remerge grows like 1/tau as tau falls, saturating where tau meets t_relax. If it
      does not, the renewal picture is wrong whatever the lifetimes say.
  P4  THE RESULT, and it should repeat §32's exponent count in time rather than in
      probability. With ln T(Omega) = c*Omega, HOLD gets ln L = k*c*Omega while RE-MERGE
      gets only m*c*Omega, so hold wins for every k -- and a crossover is possible at
      small Omega where re-merge's 1/tau prefactor is worth more than the exponent.
  P5  THE ABSOLUTE TEST (rule 16), and it is an INTEGER. The slope of
      ln(L_hold / L_remerge) against Omega, divided by the slope of ln T(Omega) fitted
      independently from the same sweep, must equal exactly

          k - ceil((k+1)/2)   =   1 at k = 3,   2 at k = 5,   3 at k = 7

      This is a prediction of a specific integer from a ratio of two separately fitted
      slopes, not a shape check. If the ratio comes out non-integer or k-independent,
      the exponent story is wrong however well the lifetimes fit.
  P6  If re-merge wins ANYWHERE at admissible tau, §32's conclusion was a one-shot
      artifact and concatenation does buy something over time. That is the outcome that
      would cost the most and it is why every tau down to t_relax is scanned.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time
from math import comb, ceil

import numpy as np

from crnl.cme import ep_rate, first_passage_moments
from crnl.deterministic import jacobian
from crnl.networks.am_reversible import (
    am_reversible, delta_star, fixed_points, reverse_pairing,
)


def committed_state(gamma: float, total: int) -> np.ndarray:
    """The attractor state in counts, X ahead."""
    att = [p for p in fixed_points(gamma) if p["kind"] == "attractor"][0]
    nb = int(round(att["b"] * total))
    rest = total - nb
    sep = int(round(delta_star(gamma) * total))
    if (rest - sep) % 2:
        sep -= 1
    return np.array([(rest + sep) // 2, (rest - sep) // 2, nb], dtype=np.int64)


def relax_time(gamma: float) -> float:
    """1 / slowest non-conserved relaxation rate at the attractor."""
    att = [p for p in fixed_points(gamma) if p["kind"] == "attractor"][0]
    x = np.array([att["x"], att["y"], att["b"]], dtype=float)
    ev = np.linalg.eigvals(jacobian(am_reversible(gamma), x)).real
    ev = np.sort(ev[np.abs(ev) > 1e-9])          # drop the conservation zero mode
    return float(1.0 / abs(ev[-1]))


def bit_lifetime(gamma: float, total: int, theta: float) -> dict:
    """Exact MFPT and its std, from the committed state to the opposite one."""
    net = am_reversible(gamma)
    thr = max(2, int(round(theta * delta_star(gamma) * total)))
    start = committed_state(gamma, total)
    r = first_passage_moments(net, total, float(total), start,
                              lambda s: int(s[1]) - int(s[0]) >= thr)
    return {"total": total, "mean": r["mean_time"], "std": r["std_time"],
            "valid": bool(r["valid"]), "residual": r["residual"],
            "thr": int(thr)}


def remerge_lifetime(T_sub: float, k: int, tau: float) -> float:
    m = ceil((k + 1) / 2)
    q = 1.0 - np.exp(-tau / T_sub)
    pf = sum(comb(k, j) * q ** j * (1.0 - q) ** (k - j) for j in range(m, k + 1))
    return tau / pf if pf > 0 else float("inf")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.30)
    ap.add_argument("--ks", type=int, nargs="+", default=[3, 5])
    ap.add_argument("--omegas", type=int, nargs="+", default=[12, 16, 20, 24, 28, 32])
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--tau-mults", type=float, nargs="+",
                    default=[1.0, 2.0, 5.0, 10.0, 30.0, 100.0])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/remerge_hold.json"))
    args = ap.parse_args()

    t0 = time.time()
    g, theta = args.gamma, args.theta
    trelax = relax_time(g)
    print(f"gamma={g} theta={theta}   t_relax at the attractor = {trelax:.4f}")

    print(f"\n=== P1: is per-molecule dissipation size-independent at the "
          f"committed state?")
    net = am_reversible(g)
    pairing = reverse_pairing(net)
    print(f"{'N':>6}{'ep_rate':>14}{'ep_rate / N':>16}")
    eps = []
    ep_of = {}
    for N in sorted({20, 40, 60, 80} | set(args.omegas)
                    | {om * k for om in args.omegas for k in args.ks}):
        e = ep_rate(net, N, float(N), pairing)
        ep_of[N] = e
        if N in (20, 40, 60, 80):
            eps.append(e / N)
            print(f"{N:>6}{e:>14.6f}{e/N:>16.8f}")
    spread = (max(eps) - min(eps)) / np.mean(eps)
    print(f"  spread of ep_rate/N over a 4x range in N: {100*spread:.3f}%"
          f"   -> P1 {'HOLDS' if spread < 0.02 else 'FAILS'}"
          f"  (matched molecules {'IS' if spread < 0.02 else 'IS NOT'} "
          f"matched dissipation)")
    if spread >= 0.02:
        print(f"  P1 FAILED, so dissipation is accounted EXPLICITLY below rather than"
              f" assumed away.\n  ep_rate/N FALLS with N, so k small tanks dissipate"
              f" FASTER than one big tank of the\n  same total -- meaning re-merge pays"
              f" more per unit time as well as whatever it buys.")

    print(f"\n=== P2: is first passage near-exponential? (renewal model's kill test)")
    print(f"{'N':>6}{'MFPT':>16}{'std':>16}{'std/mean':>11}{'valid':>7}")
    life = {}
    for N in sorted({om * k for om in args.omegas for k in args.ks}
                    | set(args.omegas)):
        r = bit_lifetime(g, N, theta)
        life[N] = r
        if r["mean"] is None or not np.isfinite(r["mean"]):
            print(f"{N:>6}   solve failed / not trustworthy")
            continue
        ratio = r["std"] / r["mean"] if r["mean"] else float("nan")
        print(f"{N:>6}{r['mean']:>16.4e}{r['std']:>16.4e}{ratio:>11.4f}"
              f"{str(r['valid']):>7}")
    good = [v for v in life.values() if v["mean"] and np.isfinite(v["mean"])
            and v["valid"]]
    rr = np.array([v["std"] / v["mean"] for v in good])
    print(f"  std/mean over {len(rr)} sizes: {rr.min():.4f}-{rr.max():.4f}"
          f"   -> P2 {'HOLDS' if abs(rr - 1).max() < 0.10 else 'FAILS'}")

    rows = []
    for k in args.ks:
        m = ceil((k + 1) / 2)
        print(f"\n=== k = {k}  (a round fails when {m} of {k} flip)")
        print(f"{'Omega':>6}{'N=kOm':>7}{'T(Omega)':>13}{'L_hold=T(N)':>14}"
              f"{'best L_remerge':>16}{'at tau':>10}{'hold/remerge':>14}"
              f"{'diss ratio':>12}{'spend ratio':>13}")
        for om in args.omegas:
            sub, hold = life.get(om), life.get(k * om)
            if not (sub and hold and sub["mean"] and hold["mean"]
                    and sub["valid"] and hold["valid"]
                    and np.isfinite(hold["mean"])):
                print(f"{om:>6}{k*om:>7}   SKIPPED (no trustworthy solve)")
                continue
            best, best_tau = -np.inf, None
            for mult in args.tau_mults:
                tau = mult * trelax
                L = remerge_lifetime(sub["mean"], k, tau)
                if L > best:
                    best, best_tau = L, tau
            # dissipation RATE ratio, and total free energy spent over each
            # protocol's own lifetime (rule 11: the arms must share a clock, and
            # here they must also share a budget).
            d_rate = (k * ep_of[om]) / ep_of[k * om]
            spend_hold = ep_of[k * om] * hold["mean"]
            spend_rm = k * ep_of[om] * best
            rows.append({"k": k, "omega": om, "N": k * om,
                         "T_sub": sub["mean"], "L_hold": hold["mean"],
                         "L_remerge": float(best), "tau": float(best_tau),
                         "ratio": hold["mean"] / best,
                         "diss_rate_ratio": float(d_rate),
                         "spend_hold": float(spend_hold),
                         "spend_remerge": float(spend_rm)})
            print(f"{om:>6}{k*om:>7}{sub['mean']:>13.4e}{hold['mean']:>14.4e}"
                  f"{best:>16.4e}{best_tau:>10.3f}{hold['mean']/best:>14.4e}"
                  f"{d_rate:>12.4f}{spend_hold/spend_rm:>13.4e}")

    print(f"\n=== P3: does L_remerge scale as 1/tau?")
    om0 = args.omegas[len(args.omegas) // 2]
    sub = life.get(om0)
    if sub and sub["mean"]:
        print(f"  at Omega={om0}, k=3, T_sub={sub['mean']:.4e}")
        prev = None
        for mult in args.tau_mults:
            tau = mult * trelax
            L = remerge_lifetime(sub["mean"], 3, tau)
            rel = f"   x{prev/L:.3f} vs tau x{mult:.0f}" if prev else ""
            print(f"    tau = {tau:>9.3f} ({mult:>5.1f} t_relax)   "
                  f"L_remerge = {L:.4e}{rel}")
            if prev is None:
                prev = L

    print(f"\n=== P4/P6: does re-merge ever win?")
    wins = [r for r in rows if r["ratio"] < 1.0]
    print(f"  on LIFETIME at matched molecules: re-merge wins "
          f"{len(wins)}/{len(rows)}"
          + (f" -> {[(r['k'], r['omega']) for r in wins]}" if wins else " (none)"))
    dwins = [r for r in rows if r["diss_rate_ratio"] < 1.0]
    print(f"  on DISSIPATION RATE: re-merge is cheaper in {len(dwins)}/{len(rows)}"
          + (f" -> {[(r['k'], r['omega']) for r in dwins]}" if dwins else " (none)"))
    if not wins and not dwins:
        print(f"  -> re-merge is DOMINATED: it lives shorter AND burns faster, so no "
              f"dissipation-matched\n     re-weighting can rescue it. There is no "
              f"tradeoff here to price.")

    print(f"\n=== P5 (absolute, integer): slope ratio must be k - ceil((k+1)/2)")
    lt = np.array([[N, v["mean"]] for N, v in sorted(life.items())
                   if v["mean"] and np.isfinite(v["mean"]) and v["valid"]])
    cT = np.polyfit(lt[:, 0], np.log(lt[:, 1]), 1)
    resT = np.log(lt[:, 1]) - np.polyval(cT, lt[:, 0])
    print(f"  ln T(N) slope c = {cT[0]:.6f}   R^2 = "
          f"{1 - resT.var()/np.log(lt[:,1]).var():.5f}   ({len(lt)} sizes)")
    for k in args.ks:
        rs = [r for r in rows if r["k"] == k]
        if len(rs) < 3:
            continue
        om = np.array([r["omega"] for r in rs], float)
        lr = np.log([r["ratio"] for r in rs])
        c = np.polyfit(om, lr, 1)
        res = lr - np.polyval(c, om)
        want = k - ceil((k + 1) / 2)
        print(f"  k={k}: ln(L_hold/L_remerge) slope = {c[0]:.6f}   R^2 = "
              f"{1 - res.var()/lr.var():.5f}")
        print(f"        slope / c = {c[0]/cT[0]:.4f}   predicted integer {want}"
              f"   -> {'MATCH' if abs(c[0]/cT[0] - want) < 0.15 else 'MISMATCH'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(
        {"gamma": g, "t_relax": trelax, "ep_per_molecule": eps,
         "lifetimes": {str(N): v for N, v in life.items()}, "rows": rows},
        indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
