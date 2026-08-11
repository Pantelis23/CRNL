"""T-OPT-b: is AM's frontier the frontier? — and does the approach to the bound have a law?

§58 found AM sits ON the (time, Q) Pareto frontier at its fast end, and named its own limit:
**§57's 696-network enumeration was never re-run with free rates.** So "AM is
Pareto-optimal" is, so far, a statement about AM's own two classes. This tests it against
other networks, and asks whether the frontier has a universal shape.

**§54 MAKES THE SEARCH CHEAP.** A network with m classes has 2m rates, and a blind grid over
that is hopeless past m = 2. But §54 classifies every class by d_r = S_X - S_Y:

  * **d_r != 0: signal-carrying.** These are the only classes that enter P at all.
  * **d_r = 0: self-mirror or p=q.** These contribute IDENTICALLY ZERO to the drift (§51,
    §54, §55 -- one fact in three places), so they act only through the pool and the noise.

§58's optimum has exactly the shape that predicts: the non-signal channel fast (rho large,
§44's lever) and weakly driven, the signal channel's drive gamma_s tracing the frontier. So
the search collapses to **(rho_ns, gamma_ns, gamma_s)** whatever m is -- the theory buying
back the dimensions.

PREDICTIONS, written before running:

  P1  GATE. AM's frontier here reproduces §58's, which ran (4.04, 5.40) at the fast end to
      (783, 1.25) at the slow end. If the reduced parameterisation cannot recover §58's own
      points, it has thrown away the axis that mattered and nothing below counts.
  P2  **THE TEST. Does any network's frontier lie strictly BELOW AM's?** If none does,
      §58's headline strengthens from "AM is on its family's frontier" to **"AM is on the
      frontier of the enumerated family"** -- which is the strongest form the founding claim
      can take here. If one does, §58's headline is a two-class statement and the better
      network gets named.
  P3  **UNIVERSALITY.** Fit `Q - 1 = a * t^(-b)` per network. §58 noted its own points sit
      near b = 1/2. **Do the exponents agree across networks?** A shared b would be a
      quantitative time-cost law for approaching the thermodynamic bound: to get within eps
      of it costs t ~ eps^(-1/b).
  P4  **THE PRIOR IS THAT THEY WILL NOT AGREE**, and it is stated so the confirming outcome
      is not the flattering one. §39.2 found a 1/sep coefficient that did not transfer
      between axes; §46 found the SCALING itself did not. Two prior attempts at
      transferable exponents in this project have failed, so a shared b would be the
      surprise and a scattered one the expectation.
  P5  If some network's frontier crosses AM's -- better at one speed, worse at another --
      then there is no single best element and the answer is a frontier of frontiers.
      Reported as such rather than resolved by picking a speed.
  P6  §40's rule still binds: any Q < 1 is a suspected instrument failure, and §58 measured
      the suspect (the two-sided set differs from the one-sided one by e^0.63*Omega).
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from experiments.amplification_signature import mirror_pairs
from experiments.free_rate_optimum import build_free, evaluate
from experiments.optimal_element import symmetric_classes


def class_is_signal(cls_reactions, net_builder):
    """Does this class carry signal (some d_r != 0) under §54?"""
    net = net_builder([cls_reactions], [1.0, 0.3])
    try:
        pairs = mirror_pairs(net, "X", "Y")
    except Exception:
        return None
    if pairs is None:
        return None
    for d, idx in pairs:
        r = net.reactions[idx]
        if r.reactants.get("X", 0) != r.reactants.get("Y", 0) and d != 0:
            return True
    return False


def frontier(points):
    """Lower-left Pareto frontier of (time, Q)."""
    out, best = [], np.inf
    for t, q in sorted(points):
        if q < best - 1e-12:
            out.append((t, q))
            best = q
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=200)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/frontier_universality.json"))
    args = ap.parse_args()

    t0 = time.time()
    cls = symmetric_classes()
    idx_of = {c[0]: i for i, c in enumerate(cls)}
    DIS = idx_of[(("X", "Y"), ("B", "B"))]
    REC = idx_of[(("B", "X"), ("X", "X"))]
    REVDIS = idx_of[(("B", "B"), ("X", "Y"))]

    # classify every class as signal-carrying or not, using §54
    sig = {}
    for i, c in enumerate(cls):
        net = build_free([c], [1.0, 0.3])
        try:
            pairs = mirror_pairs(net, "X", "Y")
        except Exception:
            pairs = None
        s = False
        if pairs:
            for d, j in pairs:
                r = net.reactions[j]
                if r.reactants.get("X", 0) != r.reactants.get("Y", 0) and d != 0:
                    s = True
        sig[i] = s
    n_sig = sum(sig.values())
    print(f"§54 classification of the 16 classes: {n_sig} signal-carrying,"
          f" {16 - n_sig} not")
    print(f"  DIS (X+Y->2B) signal-carrying? {sig[DIS]}   "
          f"REC (B+X->2X)? {sig[REC]}   REVDIS (2B->X+Y)? {sig[REVDIS]}")

    cands = {
        "AM  {dis, rec}": [DIS, REC],
        "AM+revdis": [DIS, REC, REVDIS],
        "rec only": [REC],
    }
    # add a few more capable multi-class networks from the enumeration
    extra = []
    for i in range(len(cls)):
        if i in (DIS, REC, REVDIS) or not sig[i]:
            continue
        extra.append(i)
    for i in extra[:3]:
        cands[f"AM+cls{i}"] = [DIS, REC, i]
        cands[f"{{cls{i}, dis}}"] = [DIS, i]

    RHO = (1.0, 20.0, 400.0)
    GNS = (1e-4, 0.01, 0.2, 1.0)
    GS = (0.02, 0.05, 0.10, 0.20, 0.35, 0.44, 0.49)

    results = {}
    for name, ids in cands.items():
        chosen = [cls[i] for i in ids]
        pts = []
        for rho in RHO:
            for gns in GNS:
                for gs in GS:
                    rates = []
                    for i in ids:
                        if sig[i]:
                            rates += [1.0, gs]
                        else:
                            rates += [rho, gns * rho]
                    r, _ = evaluate(chosen, np.log(np.array(rates)), args.omega)
                    if r is not None and np.isfinite(r["Q"]) and r["Q"] > 0:
                        pts.append((r["mean_T"], r["Q"]))
        if not pts:
            print(f"\n{name}: no valid cells")
            continue
        fr = frontier(pts)
        results[name] = {"classes": ids, "points": pts, "frontier": fr}
        print(f"\n{name}: {len(pts)} cells, frontier of {len(fr)}")
        print("   " + "  ".join(f"({t:.1f},{q:.2f})" for t, q in fr[:9]))

    print(f"\n=== P1 GATE: does AM's frontier reproduce §58's?")
    am = results.get("AM  {dis, rec}")
    if am:
        f = am["frontier"]
        print(f"  fast end: ({f[0][0]:.2f}, {f[0][1]:.3f})   §58 had (4.04, 5.400)")
        print(f"  slow end: ({f[-1][0]:.2f}, {f[-1][1]:.3f})  §58 had (783.3, 1.253)")
        ok = f[0][1] < 6.5 and f[-1][1] < 1.6
        print(f"  -> P1 {'HOLDS' if ok else 'FAILS'}")

    print(f"\n=== P2: does any network's frontier lie BELOW AM's?")
    if am:
        at, aq = np.array([p[0] for p in am["frontier"]]), np.array([p[1] for p in am["frontier"]])
        print(f"{'network':>18}{'best Q':>9}{'at time':>10}{'beats AM at same t?':>22}")
        for name, r in results.items():
            if name.startswith("AM  "):
                continue
            fr = r["frontier"]
            better = 0
            for t, q in fr:
                amq = np.interp(t, at, aq)
                if q < amq - 1e-9:
                    better += 1
            bq = min(q for _, q in fr)
            bt = [t for t, q in fr if q == bq][0]
            print(f"{name:>18}{bq:>9.3f}{bt:>10.1f}{better:>15} of {len(fr)}")
        anybetter = any(
            any(q < np.interp(t, at, aq) - 1e-9 for t, q in r["frontier"])
            for n, r in results.items() if not n.startswith("AM  "))
        print(f"  -> P2 {'SOME network beats AM somewhere -- §58 was a two-class statement' if anybetter else 'NO network beats AM at any speed -- §58 STRENGTHENS'}")

    print(f"\n=== P3/P4: fit Q - 1 = a * t^(-b) per network")
    print(f"{'network':>18}{'b':>9}{'a':>10}{'points':>8}")
    exps = {}
    for name, r in results.items():
        fr = [(t, q) for t, q in r["frontier"] if q > 1.0 + 1e-6]
        if len(fr) < 4:
            continue
        x = np.log(np.array([t for t, _ in fr]))
        y = np.log(np.array([q - 1.0 for _, q in fr]))
        b, la = np.polyfit(x, y, 1)
        exps[name] = -b
        print(f"{name:>18}{-b:>9.3f}{np.exp(la):>10.3f}{len(fr):>8}")
    if len(exps) > 1:
        v = np.array(list(exps.values()))
        print(f"  exponents span {v.min():.3f}..{v.max():.3f}, spread"
              f" {100*(v.max()-v.min())/v.mean():.0f}%")
        print(f"  -> {'SHARED exponent -- a time-cost law for approaching the bound' if (v.max()-v.min())/v.mean() < 0.25 else 'NOT shared, as P4 predicted from §39.2 and §46'}")

    below = [(n, t, q) for n, r in results.items() for t, q in r["frontier"] if q < 1.0]
    print(f"\n=== P6: cells reading Q < 1: {len(below)}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"signal_classes": {str(k): v for k, v in sig.items()},
                                    "results": results, "exponents": exps},
                                   indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
