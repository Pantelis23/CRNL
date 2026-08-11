"""T-OPT-a: free per-class rates, and a calibration of the ruler everything rests on

§57 searched 696 networks on a ONE-PARAMETER rate slice (all forwards 1, all reverses
gamma) and found nothing beating AM by more than 9%. §56 says that slice is a
one-dimensional curve through a high-dimensional space, so the search may have been
sampling the wrong thing -- the same class of mistake as the screening bug §57 caught.
This optimises over the full rate space.

**AND IT CHECKS THE RULER.** §40 named its own leading suspect in advance: "our absorbing
set is TWO-SIDED, `|delta| >= thr`, where the standard statement is one-sided." The
first-passage TUR is a statement about a current reaching a threshold -- ONE threshold.
Every "how far is AM from the bound" number in this project, Q_min = 5.39 included, is
measured against a two-sided set that is not the TUR's setting. **Nothing has ever measured
how much that costs.**

Two free knobs the slice could not reach, and both are already named by other sections:

  * **gamma_dis, the drive on the disagreement channel, separately from gamma_rec.** §54
    proves `X+Y->2B` is self-mirror, so d = 0 and it contributes IDENTICALLY ZERO to P.
    A drive on it is entropy for no signal. §57's winner set it to detailed balance by
    accident of the enumeration; free rates can set it deliberately.
  * **rho = k_dis/k_rec**, §44's free lever, worth 43-50% on cost. The slice fixes rho = 1.

PREDICTIONS, written before running:

  P1  GATE, two parts. (a) Q is invariant under a UNIFORM rescale of every rate -- §44's
      P1a proved Sigma exactly invariant and the relative variance is dimensionless, so
      this must hold to solver precision, and it is what makes the search
      (2m-1)-dimensional rather than 2m. (b) The optimiser started at §57's slice point
      reproduces §57's Q = 5.4750 for AM before it moves.
  P2  **THE CALIBRATION.** One-sided absorption (at +thr only) against two-sided
      (|delta| >= thr), same network, same everything else. **I expect one-sided Q to be
      LARGER**, because waiting for the correct outcome alone admits wrong-way excursions
      that two-sided absorption terminates -- inflating both <T> and Var(T). If instead it
      is smaller, then §40's 5.39 OVERSTATES how close AM is to the bound and every
      "ribosome-grade" reading needs restating.
  P3  **THE OPTIMISATION.** Free rates beat §57's 5.0045. The margin over AM grows beyond
      9%.
  P4  **THE STRUCTURE OF THE OPTIMUM, predicted from three prior sections.** The free
      optimum should have (i) **gamma_dis = 1**, detailed balance on the self-mirror
      channel, because §54 says it carries no signal; (ii) **rho = k_dis/k_rec large**,
      because §44 measured that lever; (iii) gamma_rec small. If all three appear together
      out of an optimiser that was told none of them, three independent results converge.
  P5  How close to Q = 1 does chemistry get? Reported as a number, with §40's rule binding:
      **any Q < 1 is a suspected instrument failure, not a discovery** -- and P2 is exactly
      the instrument being suspected.
  P6  If free rates do NOT beat the slice materially, §57's headline strengthens from "on a
      slice" to "in the family", which is the more surprising outcome and is reported as
      such.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
import scipy.sparse.linalg as spla
from scipy.optimize import minimize

from crnl.cme import enumerate_states, first_passage_moments, generator
from crnl.networks.am_reversible import reverse_pairing
from crnl.reactions import Reaction, ReactionNetwork
from crnl.vectorized import compile_network
from experiments.cost_of_reliability import sigma_local
from experiments.optimal_element import SPECIES, _counts, landscape, symmetric_classes
from experiments.slaving_axis import slaved


def build_free(classes, rates):
    """Network with an independent forward and reverse rate per class.

    `rates` is [kf_0, kr_0, kf_1, kr_1, ...] in class order. Exchange symmetry is
    preserved because every reaction in a class shares its rate -- which is required for
    §43's divisibility and hence for P to exist at all.
    """
    rx, seen = [], set()
    for cls, kf in zip(classes, rates[0::2]):
        for l, r in cls:
            if (l, r) not in seen:
                seen.add((l, r))
                rx.append(Reaction(_counts(l), _counts(r), float(kf), name=f"f:{l}->{r}"))
    for cls, kr in zip(classes, rates[1::2]):
        for l, r in cls:
            if (r, l) not in seen:
                seen.add((r, l))
                rx.append(Reaction(_counts(r), _counts(l), float(kr), name=f"r:{r}->{l}"))
    return ReactionNetwork(species=list(SPECIES), reactions=rx, name="free")


def q_of(net, ds, omega, eps=0.35, theta=0.80, one_sided=False):
    """§40's Q. `one_sided` absorbs only at +thr, which is the TUR's actual setting."""
    st = slaved(net, eps * ds)
    if st is None or min(st) < 0:
        return None
    nb = int(round(st[2] * omega))
    rest = omega - nb
    d0 = max(1, int(round(eps * ds * omega)))
    if (rest - d0) % 2:
        d0 -= 1
    thr = max(2, int(round(theta * ds * omega)))
    if thr <= d0 or rest - d0 < 0 or nb < 0:
        return None
    n0 = np.array([(rest + d0) // 2, (rest - d0) // 2, nb], dtype=np.int64)
    try:
        pairing = reverse_pairing(net)
    except Exception:
        return None
    if (pairing < 0).any():
        return None

    if one_sided:
        absorbed = lambda s, t=thr: int(s[0]) - int(s[1]) >= t
    else:
        absorbed = lambda s, t=thr: abs(int(s[0]) - int(s[1])) >= t

    fp = first_passage_moments(net, int(omega), float(omega), n0, absorbed)
    if not fp["valid"] or not np.isfinite(fp["var_time"]) or fp["var_time"] <= 0:
        return None
    if not np.isfinite(fp["mean_time"]) or fp["mean_time"] <= 0:
        return None
    states, index = enumerate_states(3, int(omega))
    absorb = np.array([absorbed(s) for s in states])
    if absorb.all() or not absorb.any():
        return None
    Q = generator(net, int(omega), float(omega))
    tr = np.where(~absorb)[0]
    tmap = {int(i): r for r, i in enumerate(tr)}
    comp = compile_network(net, float(omega))
    sig = sigma_local(net, comp, states, pairing)[tr]
    try:
        Sig = float(spla.spsolve(Q[tr][:, tr].tocsr(), -sig)[tmap[index[tuple(n0)]]])
    except Exception:
        return None
    if not np.isfinite(Sig) or Sig <= 0:
        return None
    rel = fp["var_time"] / fp["mean_time"] ** 2
    return {"Q": float(rel * Sig / 2.0), "Sigma": Sig, "mean_T": fp["mean_time"],
            "rel_var": float(rel), "delta_star": ds, "thr": int(thr)}


def evaluate(classes, logr, omega, one_sided=False):
    rates = np.exp(np.asarray(logr, dtype=float))
    rates = rates / rates.max()                       # fix the scale (P1a)
    net = build_free(classes, rates)
    ds = landscape(net)
    if ds is None or ds < 0.05:
        return None, None
    return q_of(net, ds, omega, one_sided=one_sided), net


AM_IDX = None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=200)
    ap.add_argument("--restarts", type=int, default=4)
    ap.add_argument("--maxiter", type=int, default=120)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/free_rate_optimum.json"))
    args = ap.parse_args()

    t0 = time.time()
    cls = symmetric_classes()
    dis = next(i for i, c in enumerate(cls) if c[0] == (("X", "Y"), ("B", "B")))
    rec = next(i for i, c in enumerate(cls) if c[0] == (("B", "X"), ("X", "X")))
    am = [cls[dis], cls[rec]]

    print("=== P1 GATE (a): is Q invariant under a uniform rate rescale?")
    base = np.array([1.0, 0.05, 1.0, 0.05])
    ref, _ = evaluate(am, np.log(base), args.omega)
    print(f"{'lambda':>9}{'Q':>12}{'rel dev':>12}")
    worst = 0.0
    for lam in (0.01, 1.0, 30.0, 1000.0):
        r, _ = evaluate(am, np.log(base * lam), args.omega)
        dev = abs(r["Q"] - ref["Q"]) / ref["Q"]
        worst = max(worst, dev)
        print(f"{lam:>9.2f}{r['Q']:>12.6f}{dev:>12.2e}")
    print(f"  -> P1a {'HOLDS' if worst < 1e-9 else 'FAILS'}")

    print(f"\n=== P1 GATE (b): §57's AM slice point")
    print(f"  AM at gamma=0.05, two-sided: Q = {ref['Q']:.4f}   (§57 measured 5.4750)")

    print(f"\n=== P2 CALIBRATION: one-sided vs two-sided absorption, AM")
    print(f"{'gamma':>7}{'two-sided Q':>14}{'one-sided Q':>14}{'ratio':>9}"
          f"{'<T> 2s':>10}{'<T> 1s':>10}")
    calib = []
    for g in (0.05, 0.10, 0.20, 0.30):
        b = np.log(np.array([1.0, g, 1.0, g]))
        r2, _ = evaluate(am, b, args.omega, one_sided=False)
        r1, _ = evaluate(am, b, args.omega, one_sided=True)
        if r2 is None or r1 is None:
            continue
        calib.append({"gamma": g, "two": r2["Q"], "one": r1["Q"]})
        print(f"{g:>7.2f}{r2['Q']:>14.4f}{r1['Q']:>14.4f}{r1['Q']/r2['Q']:>9.3f}"
              f"{r2['mean_T']:>10.3f}{r1['mean_T']:>10.3f}")
    if calib:
        rr = np.array([c["one"] / c["two"] for c in calib])
        print(f"  one/two ratio: {rr.min():.3f}..{rr.max():.3f}")
        print(f"  -> {'one-sided is LARGER, as predicted' if rr.min() > 1 else 'one-sided is SMALLER -- §40 5.39 OVERSTATES how close AM is'}")

    print(f"\n=== P3/P4: free-rate optimisation (two-sided, comparable to §57)")
    cands = {
        "AM": [cls[dis], cls[rec]],
        "AM+revdis": [cls[dis], cls[rec],
                      next(c for c in cls if c[0] == (("B", "B"), ("X", "Y")))],
    }
    results = {}
    for name, chosen in cands.items():
        m = 2 * len(chosen)
        best = None
        for s in range(args.restarts):
            rng = np.random.default_rng(1000 + s)
            x0 = np.log(np.concatenate([[1.0, 0.05]] * len(chosen))) \
                if s == 0 else rng.uniform(-4, 2, m)

            def obj(z):
                r, _ = evaluate(chosen, z, args.omega)
                return r["Q"] if r is not None else 1e6

            res = minimize(obj, x0, method="Nelder-Mead",
                           options={"maxiter": args.maxiter, "xatol": 1e-3,
                                    "fatol": 1e-4})
            if res.fun < 1e5 and (best is None or res.fun < best[0]):
                best = (float(res.fun), res.x.copy())
        if best is None:
            print(f"  {name}: no valid optimum found")
            continue
        q, z = best
        rates = np.exp(z); rates /= rates.max()
        net = build_free(chosen, rates)
        results[name] = {"Q": q, "rates": rates.tolist(),
                         "reactions": [(r.name, r.k) for r in net.reactions]}
        print(f"\n  {name}: best Q = {q:.4f}")
        for r in net.reactions:
            print(f"    {r.name:<26} k = {r.k:.5f}")
        kf = rates[0::2]; kr = rates[1::2]
        print(f"    gamma per class (kr/kf): {[round(a/b, 4) for a, b in zip(kr, kf)]}")
        if len(chosen) >= 2:
            print(f"    rho = k_dis/k_rec = {kf[0]/kf[1]:.4f}")

    print(f"\n=== P3/P5/P6 verdict")
    if results:
        bq = min(r["Q"] for r in results.values())
        print(f"  best free-rate Q = {bq:.4f}")
        print(f"  §57 slice best   = 5.0045   AM on slice = 5.4750")
        print(f"  free vs §57: {5.0045/bq:.3f}x;  free vs AM: {5.4750/bq:.3f}x")
        if bq < 1.0:
            print(f"  ⚠ Q < 1: per §40's pre-registered rule this is a SUSPECTED INSTRUMENT")
            print(f"    FAILURE, and P2 above measures exactly the suspect.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"gate_worst": worst, "calibration": calib,
                                    "optima": results}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
