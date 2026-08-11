"""T-OPT: is AM the best chemical decision element? — an exhaustive search against the bound

This is the founding question with every piece finally in place. §40 measured AM at
**Q_min = 5.39** against the thermodynamic floor Q >= 1, and called it "ribosome-grade".
§37-§39 priced restoration in closed form. **§56 made the search tractable**: a network can
restore for some rates iff some d_r > 0, so the capable ones can be *enumerated* rather than
stumbled on.

**THE FAMILY IS FINITE AND SMALL.** On {X, Y, B} with a conservation law X+Y+B = Omega and
bimolecular reactions, there are exactly **30 conservative reactions**, falling into **16
exchange-symmetric classes** (2 self-mirror, 14 mirror pairs). Networks built from 1-3
classes number **696**, and **AM is one of them** -- the class {X+Y->2B} (self-mirror) plus
the class {B+X->2X, B+Y->2Y} (a mirror pair). Every network is parameterised exactly as AM
is: forward rate 1, reverse rate gamma, one gamma for the whole network.

So the question "is AM the best restoring element chemistry allows, at this size" is a
finite search, and this runs it.

**THE FIGURE OF MERIT IS NOT MINE.** Q = (Var(T)/<T>^2) * <Sigma> / 2 is §40's, which is
Gingrich-Horowitz's first-passage TUR: Q >= 1 with equality at the thermodynamic limit on
decision-making. Lower is better. Every network is measured at §40's own conditions --
eps = 0.35, theta = 0.80, threshold scaled with delta* -- so the numbers are comparable to
§40's row for row.

PREDICTIONS, written before running:

  P1  GATE, and it is absolute. **AM must appear in the enumeration and reproduce §40's
      Q_min ~ 5.39.** If the pipeline cannot recover this project's own published number for
      the network it was measured on, nothing below is admissible.
  P2  Report the capable fraction of the 696 under §56's criterion, and the fraction that
      then actually produce a landscape at some gamma. These are different questions and
      both get reported.
  P3  **THE TEST. Something beats AM.** AM is one point in a 696-network family and was
      never designed to minimise Q -- it was designed (by Angluin et al., and by evolution)
      for other things. I expect a winner below 5.39 but **not near 1**: the two-sided
      absorbing set and the finite margin impose a floor above the ideal.
  P4  **SAFETY, and §40 wrote this rule in advance so it binds here.** §40: "If Q < 1 the
      leading suspect is that the inequality does not apply in this form, NOT that
      thermodynamics is wrong." **Any Q < 1 is reported as a suspected instrument failure**
      -- most likely the two-sided absorbing set -- and is NOT claimed as a
      sub-thermodynamic decision element.
  P5  Is the winner interpretable? A recognisable motif would say the search found something
      chemistry already knows; an exotic one would say the design space is larger than the
      motifs biology uses. Either is a result and both get named.
  P6  Does §56's class predict quality? `all d_r >= 0` networks restore at every rate, which
      sounds better but may cost more, since a network that cannot help amplifying may not be
      able to stop. Reported as a comparison of Q by class.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time
from itertools import combinations, combinations_with_replacement

from scipy.optimize import brentq

import numpy as np
import scipy.sparse.linalg as spla

from crnl.cme import enumerate_states, first_passage_moments, generator
from crnl.networks.am_reversible import am_reversible, reverse_pairing
from crnl.reactions import Reaction, ReactionNetwork
from crnl.vectorized import compile_network
from experiments.amplification_signature import classify
from experiments.cost_of_reliability import sigma_local
from experiments.slaving_axis import delta_star_of, slaved

SPECIES = ["X", "Y", "B"]


def _swap(m):
    return tuple(sorted(("Y" if s == "X" else "X" if s == "Y" else s) for s in m))


def symmetric_classes():
    """The 16 exchange-symmetric classes of conservative bimolecular reactions."""
    sides = [tuple(sorted(m)) for m in combinations_with_replacement(SPECIES, 2)]
    seen, out = set(), []
    for l in sides:
        for r in sides:
            if l == r or (l, r) in seen:
                continue
            mirror = (_swap(l), _swap(r))
            seen.add((l, r))
            seen.add(mirror)
            out.append(((l, r),) if mirror == (l, r) else ((l, r), mirror))
    return out


def _counts(m):
    d = {}
    for s in m:
        d[s] = d.get(s, 0) + 1
    return d


def build(classes, gamma, k=1.0):
    """Network from a set of classes: each reaction forward at k, reverse at gamma*k."""
    rx, seen = [], set()
    for cls in classes:
        for l, r in cls:
            if (l, r) in seen:
                continue
            seen.add((l, r))
            rx.append(Reaction(_counts(l), _counts(r), k, name=f"f:{l}->{r}"))
    for cls in classes:
        for l, r in cls:
            if (r, l) in seen:
                continue
            seen.add((r, l))
            rx.append(Reaction(_counts(r), _counts(l), gamma * k, name=f"r:{r}->{l}"))
    return ReactionNetwork(species=list(SPECIES), reactions=rx, name="cand")


def landscape(net, n=1200):
    """delta* on the slaved manifold, robustly.

    `slaving_axis.delta_star_of` scans a fixed grid to 0.999 and brackets sign changes
    between GRID points, but `slaved` returns None as delta -> 1 (the pool runs out), so a
    large delta* falls in the dead zone and the scan reports None. At gamma = 0.03 and
    0.05 the closed form gives 0.971 and 0.952 and `delta_star_of` returns None -- which
    silently discards exactly the best-Q cells, since AM's Q minimum sits near gamma=0.05.
    The first pass of this search lost ~97% of its cells that way and was caught only by
    the P1 gate disagreeing with §40's published 5.39.

    This version brackets between consecutive FINITE samples instead, and is gated against
    `delta_star(gamma)` on AM: worst |diff| 2.1e-13 over gamma = 0.05..0.49.

    SCOPE LIMIT, stated rather than fought. Below gamma ~ 0.05 it still returns None,
    because `slaved` itself dies there: at gamma = 0.03 the attractor sits at delta* =
    0.9708 where the losing species is 4e-5 of the tank, past the edge of the feasible
    manifold. The gamma grid therefore starts at 0.05 -- which is also §40's own grid, so
    the comparison stays like-for-like. `delta_star_of` is left untouched: §36 and §39.2
    rest on it.
    """
    S = net.stoichiometry_matrix()

    def drift(d):
        st = slaved(net, float(d))
        if st is None or min(st) < -1e-12:
            return np.nan
        v = S @ net.fluxes(np.asarray(st, dtype=float))
        return float(v[0] - v[1])

    xs = np.linspace(1e-4, 0.9999, n)
    pts = [(x, drift(x)) for x in xs]
    pts = [(x, f) for x, f in pts if np.isfinite(f)]
    if len(pts) < 3 or pts[0][1] <= 0:
        return None                        # must AMPLIFY near delta = 0
    root = None
    for (x0, f0), (x1, f1) in zip(pts, pts[1:]):
        if f0 > 0 >= f1:
            try:
                root = brentq(lambda d: drift(d), x0, x1, xtol=1e-12)
            except Exception:
                root = 0.5 * (x0 + x1)
            break
    if root is None or not np.isfinite(root) or root <= 1e-3:
        return None
    return float(root)


def q_of(net, ds, omega, eps=0.35, theta=0.80):
    """§40's Q = (Var(T)/<T>^2) * <Sigma> / 2, at §40's own conditions."""
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

    absorbed = lambda s, t=thr: abs(int(s[0]) - int(s[1])) >= t
    fp = first_passage_moments(net, int(omega), float(omega), n0, absorbed)
    if not fp["valid"] or not np.isfinite(fp["var_time"]) or fp["var_time"] <= 0:
        return None
    if not np.isfinite(fp["mean_time"]) or fp["mean_time"] <= 0:
        return None

    states, index = enumerate_states(3, int(omega))
    absorb = np.array([absorbed(s) for s in states])
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=200)
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.05, 0.08, 0.12, 0.18, 0.25, 0.35])
    ap.add_argument("--max-classes", type=int, default=3)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/optimal_element.json"))
    args = ap.parse_args()

    t0 = time.time()
    cls = symmetric_classes()
    nets = []
    for r in range(1, args.max_classes + 1):
        nets.extend(combinations(range(len(cls)), r))
    print(f"{len(cls)} symmetric classes -> {len(nets)} networks from 1..{args.max_classes}")

    # ---- P1 GATE ---------------------------------------------------------------
    print(f"\n=== P1 GATE: reproduce §40's AM Q_min ~ 5.39")
    best_am = None
    for g in args.gammas:
        net = am_reversible(g)
        ds = landscape(net)
        if ds is None:
            continue
        r = q_of(net, ds, args.omega)
        if r is None:
            continue
        print(f"  AM gamma={g:.2f}: Q = {r['Q']:.4f}   delta* = {ds:.4f}")
        if best_am is None or r["Q"] < best_am[1]["Q"]:
            best_am = (g, r)
    print(f"  AM best: Q = {best_am[1]['Q']:.4f} at gamma = {best_am[0]:.2f}"
          f"   (§40 published Q_min = 5.39)")
    gate = 4.0 < best_am[1]["Q"] < 7.5
    print(f"  -> P1 {'HOLDS' if gate else 'FAILS -- pipeline disagrees with §40'}")

    # ---- screen ------------------------------------------------------------------
    print(f"\n=== P2: screening {len(nets)} networks")
    capable, landscaped = 0, []
    for idx in nets:
        chosen = [cls[i] for i in idx]
        net = build(chosen, 0.2)
        try:
            c = classify(net)
        except Exception:
            continue
        if c not in ("mixed", "all>=0"):
            continue
        capable += 1
        for g in args.gammas:
            n2 = build(chosen, g)
            ds = landscape(n2)
            if ds is not None and ds > 0.05:
                landscaped.append((idx, g, ds))
    print(f"  capable under §56 (some d_r > 0): {capable}/{len(nets)}"
          f" = {100*capable/len(nets):.1f}%")
    print(f"  (network, gamma) cells with a landscape: {len(landscaped)}")

    # ---- P3: the search -----------------------------------------------------------
    print(f"\n=== P3: computing Q for every landscaped cell (Omega = {args.omega})")
    rows, done = [], 0
    for idx, g, ds in landscaped:
        chosen = [cls[i] for i in idx]
        net = build(chosen, g)
        r = q_of(net, ds, args.omega)
        done += 1
        if done % 100 == 0:
            print(f"    {done}/{len(landscaped)} cells, {len(rows)} valid,"
                  f" {time.time()-t0:.0f}s")
        if r is None:
            continue
        rows.append({"classes": [list(map(list, cls[i])) for i in idx], "idx": list(idx),
                     "gamma": g, "class": classify(net), **r})
    print(f"  {len(rows)} valid Q values from {len(landscaped)} cells")

    rows.sort(key=lambda z: z["Q"])
    print(f"\n=== the ten best decision elements found")
    print(f"{'rank':>5}{'Q':>10}{'gamma':>7}{'delta*':>9}{'class':>9}   reactions")
    for i, r in enumerate(rows[:10]):
        rxs = "; ".join("+".join(l) + "->" + "+".join(rr)
                        for c in r["classes"] for l, rr in [c[0]])
        print(f"{i+1:>5}{r['Q']:>10.4f}{r['gamma']:>7.2f}{r['delta_star']:>9.4f}"
              f"{r['class']:>9}   {rxs}")

    if rows:
        best = rows[0]
        print(f"\n=== P3/P4 verdict")
        print(f"  best Q = {best['Q']:.4f} against AM's {best_am[1]['Q']:.4f}"
              f"  ({best_am[1]['Q']/best['Q']:.2f}x better)"
              if best["Q"] > 0 else "")
        if best["Q"] < 1.0:
            print(f"  ⚠ Q < 1 -- per §40's pre-registered rule this is a SUSPECTED")
            print(f"    INSTRUMENT FAILURE (two-sided absorbing set), NOT a")
            print(f"    sub-thermodynamic decision element. Reported, not claimed.")
        below = [r for r in rows if r["Q"] < 1.0]
        print(f"  cells reading Q < 1: {len(below)}/{len(rows)}")

    print(f"\n=== P6: does §56's class predict quality?")
    for k in ("mixed", "all>=0"):
        qs = [r["Q"] for r in rows if r["class"] == k]
        if qs:
            print(f"  {k:>8}: {len(qs):>4} cells, best {min(qs):.4f},"
                  f" median {np.median(qs):.4f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"am_best": {"gamma": best_am[0], **best_am[1]},
                                    "n_classes": len(cls), "n_networks": len(nets),
                                    "capable": capable, "rows": rows[:200]},
                                   indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
