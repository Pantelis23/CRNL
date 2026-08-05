"""Is AM near the thermodynamic bound on decision-making? — the first external standard

Every cost number this project has produced -- §38's G, §39's closed form, §37's R -- is
a measurement with no external reference. Nothing says whether AM is a GOOD decision
element or merely a measurable one. A thermodynamic uncertainty relation supplies that
reference, because it is a BOUND rather than a fit.

The first-passage TUR (Gingrich & Horowitz) bounds the relative fluctuation of a
first-passage time by the entropy produced up to it:

    Var(T) / <T>^2   >=   2 / <Sigma>

so the dimensionless **uncertainty product Q = (Var(T)/<T>^2) * <Sigma> / 2** satisfies
Q >= 1, with Q = 1 at saturation. **Q is how far a decision sits from the thermodynamic
limit on decision-making.**

**BOTH SIDES ARE ALREADY EXACT AND ALREADY BUILT.** `first_passage_moments` solves
`Q_tt m2 = -2T` for the second moment and returns Var(T) -- added earlier this session
and cross-checked against the SSA in `tests/test_cme.py`. §37's solve gives <Sigma> from
`Q_tt Sigma = -sigma_local`. Same generator, same absorbing set, same on-manifold start
(§36), theta scaled with delta*. One solve each.

**WHAT A VIOLATION WOULD MEAN, stated first so it is not rationalised later.** The TUR in
this form is derived for first passage of a current-like observable to a threshold.
`delta = n_X - n_Y` is current-like (it changes by exactly +-1 per reaction, and only
four of the six reactions move it -- §30). But our absorbing set is TWO-SIDED,
`|delta| >= thr`, where the standard statement is one-sided. **If Q < 1 the leading
suspect is that the inequality does not apply in this form, NOT that thermodynamics is
wrong** -- and that must be reported before any physical reading is attempted.

PREDICTIONS, written before running:

  P1  GATE. Var(T) > 0 and finite at every cell (the solver returns invalid on negative
      variance rather than clipping), and <Sigma> > 0. Cells failing either are dropped
      and reported.
  P2  Q is Omega-INDEPENDENT. For a barrier crossing T becomes deterministic as Omega
      grows, so Var(T)/<T>^2 ~ 1/Omega, while <Sigma> ~ Omega. Their product is O(1).
      **This is a real structural check**: if Q drifts with Omega it is not a property of
      the chemistry and nothing below means anything.
  P3  THE TEST. Q >= 1 at every gamma. Reported with the caveat above if violated.
  P4  THE INTERESTING ONE. If Q >= 1, does its MINIMUM coincide with §38's optimal drive
      gamma ~ 0.20? The drive that minimises cost per e-fold of gain and the drive that
      comes closest to the thermodynamic bound are logically independent quantities.
      **If they coincide, that is a genuine unification** -- the optimum is not an
      artifact of the cost measure but the point where AM is closest to
      thermodynamically optimal. If they differ, they are two different optima and the
      design principle needs restating.
  P5  Q is O(1) rather than O(100). A decision element operating orders of magnitude from
      the bound would say the motif is thermodynamically mediocre and that its
      ubiquity in biology is about something else -- robustness, simplicity, speed.
      That is a real outcome and is reported as one.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.cme import first_passage_moments
from crnl.networks.am_reversible import am_reversible, delta_star
from experiments.cost_of_reliability import cell
from experiments.slaving_axis import slaved

THETA = 0.80


def tur_cell(gamma, omega, eps, theta=THETA):
    net = am_reversible(gamma)
    ds = delta_star(gamma)
    x0 = eps * ds
    st = slaved(net, x0)
    nb = int(round(st[2] * omega))
    rest = omega - nb
    d0 = max(1, int(round(x0 * omega)))
    if (rest - d0) % 2:
        d0 -= 1
    n0 = np.array([(rest + d0) // 2, (rest - d0) // 2, nb], dtype=np.int64)
    thr = max(2, int(round(theta * ds * omega)))
    fp = first_passage_moments(net, int(omega), float(omega), n0,
                               lambda s, t=thr: abs(int(s[0]) - int(s[1])) >= t)
    if not fp["valid"] or not np.isfinite(fp["var_time"]) or fp["var_time"] <= 0:
        return None
    c = cell(gamma, omega, eps, theta)
    if not np.isfinite(c["Sigma"]) or c["Sigma"] <= 0:
        return None
    rel_var = fp["var_time"] / fp["mean_time"] ** 2
    return {"gamma": gamma, "omega": omega, "eps": eps,
            "mean_T": fp["mean_time"], "var_T": fp["var_time"],
            "rel_var": float(rel_var), "Sigma": c["Sigma"],
            "Q": float(rel_var * c["Sigma"] / 2.0),
            "L": c["L"]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40])
    ap.add_argument("--omegas", type=int, nargs="+", default=[150, 300, 450, 600])
    ap.add_argument("--eps-frac", type=float, default=0.35)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/tur_bound.json"))
    args = ap.parse_args()

    t0 = time.time()
    print("Q = (Var(T)/<T>^2) * <Sigma> / 2     TUR says Q >= 1, saturation at 1")
    print(f"{'gamma':>7}{'Omega':>7}{'<T>':>11}{'rel var':>12}{'<Sigma>':>12}"
          f"{'Q':>10}")
    rows, dropped = [], []
    for g in args.gammas:
        for om in args.omegas:
            try:
                r = tur_cell(g, om, args.eps_frac)
            except Exception as e:
                dropped.append((g, om, type(e).__name__)); continue
            if r is None:
                dropped.append((g, om, "invalid solve")); continue
            rows.append(r)
            print(f"{g:>7.2f}{om:>7}{r['mean_T']:>11.4f}{r['rel_var']:>12.3e}"
                  f"{r['Sigma']:>12.3f}{r['Q']:>10.4f}")
        print()

    if dropped:
        print(f"  dropped {len(dropped)} cells (reported, not hidden): {dropped}\n")

    print("=== P2 gate: is Q Omega-independent?")
    print(f"{'gamma':>7}" + "".join(f"{o:>10}" for o in args.omegas) + f"{'spread':>9}")
    conv = {}
    for g in args.gammas:
        qs = [r["Q"] for o in args.omegas
              for r in rows if r["gamma"] == g and r["omega"] == o]
        if len(qs) < 2:
            continue
        q = np.array(qs)
        conv[g] = q[-1]
        print(f"{g:>7.2f}" + "".join(f"{x:>10.4f}" for x in q)
              + f"{100*(q.max()-q.min())/q.mean():>8.2f}%")
    if conv:
        allsp = [100 * (np.ptp([r["Q"] for r in rows if r["gamma"] == g])
                        / np.mean([r["Q"] for r in rows if r["gamma"] == g]))
                 for g in conv]
        print(f"  worst spread across Omega: {max(allsp):.2f}%"
              f"   -> P2 {'HOLDS' if max(allsp) < 10 else 'FAILS'}")

    print("\n=== P3/P5: does the bound hold, and how close is AM to it?")
    q = np.array([r["Q"] for r in rows])
    print(f"  Q over {len(q)} cells: {q.min():.4f} .. {q.max():.4f}")
    below = [r for r in rows if r["Q"] < 1.0]
    if below:
        print(f"  ⚠ {len(below)} cells BELOW 1 -> the inequality does not hold in this "
              f"form here. Leading suspect is the TWO-SIDED absorbing set, not physics.")
        print(f"    worst: Q = {min(r['Q'] for r in below):.4f}")
    else:
        print(f"  bound HOLDS at every cell. Closest approach Q = {q.min():.4f}"
              f"  ({100*(q.min()-1):.1f}% above saturation)")

    print("\n=== P4: does the TUR optimum coincide with §38's cost optimum (gamma ~ 0.20)?")
    if conv:
        gs = np.array(sorted(conv)); qs = np.array([conv[g] for g in gs])
        print(f"{'gamma':>7}{'Q':>10}")
        for gg, qq in zip(gs, qs):
            print(f"{gg:>7.2f}{qq:>10.4f}")
        i = int(np.argmin(qs))
        interior = 0 < i < len(qs) - 1
        print(f"\n  minimum Q = {qs[i]:.4f} at gamma = {gs[i]:.2f}"
              f"   ({'interior' if interior else 'at the edge'})")
        print(f"  §38's cost optimum: gamma ~ 0.20  -> "
              f"{'COINCIDE' if abs(gs[i]-0.20) < 0.06 else 'DIFFERENT optima'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
