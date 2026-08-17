"""T-COST-n: does the escape action have a thermodynamic price? Asked against the RIGHT exponent

The founding question -- does reliability cost dissipation -- was pursued four ways and closed
each time (SYNTHESIS §3). **Every one of those closures priced a quantity that §80 has since
shown is not the one that governs anything.** §80/§81 established that for a chemically-coupled
cascade the per-stage error is escape, not misreading, so the exponent that sets reliability and
therefore depth is the escape action A. A has never been asked the thermodynamic question.

**The test is available exactly, because the drive of these networks is ONE number.** AM's
reversible pairs span a one-dimensional cycle space, so its entire non-equilibrium force is the
cycle affinity -3 ln(gamma) (§16, and `cycle_affinity` derives it generically rather than
assuming it). Schloegl's two reversible pairs likewise span one cycle, with affinity
ln(k1a k2r / k1r k2b). **So the force can be held EXACTLY fixed while the kinetics move**, which
is rule 9's opposite sweep in its strongest form: not a nuisance moved against a suspected cause,
but a suspected cause pinned to a constant by construction while everything else varies.

  If A is constant on a level set of the affinity, reliability has a thermodynamic price in the
  strict sense, and §3's four closures were all looking at the wrong observable.
  If A moves on that level set, dissipation does not price the exponent that governs depth, and
  the closure is complete rather than provisional.

**WHY THE SADDLE IS PINNED AT r2 = 1 in the Schloegl arm.** The affinity ln(e1 e2 / e3) is
invariant under rescaling all three roots (the lambda^3 cancels), while A = -int ln(mu/lam) dx
scales LINEARLY with the concentration scale (§75's collapse in lambda*Omega, §78's eta ~ lambda).
So a uniform rescaling moves A at fixed affinity for free -- and would be correctly dismissed as
a change of units. Fixing r2 = 1 removes that freedom, leaving a genuine one-parameter family of
DIFFERENT LANDSCAPES at one force and one scale.

**SCOPE, recorded before the measurement, not discovered after it.** The deep-barrier AM regime
(gamma^3 = 0.008, affinity 4.828) was scouted first and BOTH exact instruments failed there: the
stationary route dropped most cells on the engine's own trustworthiness guard, and the direct 2-D
first-passage solve returned NEGATIVE mean times with ln T saturating at ~35 in every surviving
cell -- the double-precision ceiling for that conditioning, not a barrier. AM is therefore
measured only in the shallow-barrier regime, and that is a limitation of the instrument, not a
choice of flattering data. The 1-D arm carries the section because it is exact to quadrature.

PREDICTIONS, written before running.

  P1  GATE. AM's local A = d ln T / d Omega must CONVERGE (rule 20: convergence, not a fixed
      tolerance), and any cell with a negative mean first-passage time is an instrument failure
      and not a small number (§80 P1). Non-converging and negative cells are EXCLUDED and
      COUNTED, never fitted -- §81.1 is one session old, where a saturated numerator produced a
      smooth, plausible series that inflated the headline 22x.
  P2  **THE TEST, 1-D and exact.** Schloegl roots swept along the level set e1 e2 / e3 = const at
      fixed saddle. **Predicted: A moves, and by a large factor** -- A is a functional of the
      rate functions over the whole barrier, and the affinity is a single ratio of rate
      constants that cannot see the shape between them.
  P3  **THE TEST, 2-D.** AM with per-reaction reverse ratios (g1, g2, g2) on the level set
      g1 g2^2 = const, which holds the affinity exactly fixed AND keeps X <-> Y exchange
      symmetry intact so the saddle stays on the diagonal. **Predicted: A moves here too.** A
      second substrate is required because rule 9 has failed three times on a single axis.
  P4  **THE OPPOSITE SWEEP, so the dissociation runs BOTH ways.** P2/P3 show affinity does not
      determine A. The converse -- that A does not determine affinity -- needs pairs with the
      SAME A and DIFFERENT affinity. **Predicted: they exist**, since A is one number extracted
      from a two-parameter family. One-way insensitivity would leave open that A still bounds
      the drive.
  P5  **WHAT DOES PRICE IT, AND IT IS A SUSPECT (rule 17), NOT A RESULT.** The obvious reading of
      P2/P3 is "A tracks the landscape". **The kill test is whether the landscape's headline
      number is enough: predicted NO** -- Delta alone must not determine A, because Delta is one
      scalar and A integrates a function over the whole barrier. Pairs with equal Delta and
      different A are searched for explicitly. If they do not exist, the suspect survives one
      test and is still a suspect.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import scipy.sparse.linalg as spla
from scipy.integrate import quad
from scipy.optimize import brentq, fsolve

from crnl.cme import enumerate_states, generator
from crnl.networks.am_reversible import cycle_affinity, reverse_pairing
from crnl.reactions import Reaction, ReactionNetwork
from experiments.cascade_schlogl import schlogl_consts

# ----------------------------------------------------------------- Schloegl, exact


def schlogl_affinity(r1, r3, r2=1.0):
    """The cycle affinity ln(k1a k2r / k1r k2b) in terms of the roots.

    With k1r = 1 and (k1a, k2b, k2r) = (e1, e3, e2), this is ln(e1 e2 / e3) -- a function of
    the roots alone, and invariant under rescaling them all (hence the pinned saddle).
    """
    e1, e2, e3 = r1 + r2 + r3, r1 * r2 + r1 * r3 + r2 * r3, r1 * r2 * r3
    return float(np.log(e1 * e2 / e3))


def schlogl_A(r1, r3, r2=1.0):
    """A = -int ln(mu/lam) dx, §81's exact quadrature."""
    k1a, k1r, k2b, k2r = schlogl_consts(r1, r2, r3)
    val, err = quad(lambda x: np.log((k1r * x ** 3 + k2r * x) / (k1a * x ** 2 + k2b)),
                    r2, r3, limit=200)
    return -float(val), float(err)


def r3_at_affinity(r1, target, r2=1.0, hi=400.0):
    """The partner root putting (r1, r2, r3) on the affinity level set."""
    return float(brentq(lambda t: schlogl_affinity(r1, t, r2) - target, r2 * 1.0001, hi))


# ----------------------------------------------------------------- AM, per-reaction gammas


def am_gammas(g1, g2, g3=None, k=1.0):
    """Reversible AM with the three reverse ratios set independently.

    g3 defaults to g2, which keeps X <-> Y exchange symmetry exact and the saddle on the
    diagonal. With g2 != g3 the element is asymmetric, the saddle leaves the diagonal, and
    (scouted) the bistable structure can collapse to a single fixed point -- a different
    experiment, and not this one.
    """
    if g3 is None:
        g3 = g2
    return ReactionNetwork(
        species=["X", "Y", "B"],
        reactions=[
            Reaction({"X": 1, "Y": 1}, {"B": 2}, k, name="f1:X+Y->2B"),
            Reaction({"B": 1, "X": 1}, {"X": 2}, k, name="f2:B+X->2X"),
            Reaction({"B": 1, "Y": 1}, {"Y": 2}, k, name="f3:B+Y->2Y"),
            Reaction({"B": 2}, {"X": 1, "Y": 1}, g1 * k, name="r1:2B->X+Y"),
            Reaction({"X": 2}, {"B": 1, "X": 1}, g2 * k, name="r2:2X->B+X"),
            Reaction({"Y": 2}, {"B": 1, "Y": 1}, g3 * k, name="r3:2Y->B+Y"),
        ],
        name=f"am-g({g1},{g2},{g3})",
    )


def am_affinity(g1, g2):
    net = am_gammas(g1, g2)
    return float(cycle_affinity(net, reverse_pairing(net)))


def am_rail(g1, g2):
    """The high-x attractor of the mass-action field, found numerically (g1 != g2 moves it)."""
    net = am_gammas(g1, g2)
    S = net.stoichiometry_matrix()
    best = None
    for d0 in np.linspace(0.05, 0.95, 25):
        x0 = np.clip(np.array([(0.8 + d0) / 2, (0.8 - d0) / 2, 0.2]), 1e-6, None)
        x0 /= x0.sum()
        s, _, ier, _ = fsolve(lambda z: np.append((S @ net.fluxes(z))[:2], z.sum() - 1.0),
                              x0, full_output=True)
        if ier == 1 and (s > -1e-9).all() and s[0] - s[1] > 1e-6:
            if best is None or s[0] - s[1] > best[0] - best[1]:
                best = s
    return best


def am_ln_mfpt(g1, g2, omega):
    """ln of the exact mean first passage from the high rail to the diagonal n_X = n_Y.

    Transient set is {n_X > n_Y}; everything at or past the diagonal has escaped. Qtt is a
    strictly diagonally dominant M-matrix, so the solve is far better conditioned than the
    stationary null-space problem -- but NOT unconditionally: at large A*Omega the solution
    spans e^(A Omega) and returns negative times (P1's guard, and §80's).
    """
    net = am_gammas(g1, g2)
    states, index = enumerate_states(3, int(omega))
    d = states[:, 0].astype(np.int64) - states[:, 1].astype(np.int64)
    trans = np.where(d > 0)[0]
    Qtt = generator(net, int(omega), float(omega))[trans][:, trans].tocsc()
    T = spla.spsolve(Qtt, -np.ones(len(trans)))
    x = am_rail(g1, g2)
    if x is None:
        return None
    nx = int(round(x[0] * omega))
    ny = int(round(x[1] * omega))
    if nx <= ny:
        return None
    r = int(np.where(trans == index[(nx, ny, int(omega) - nx - ny)])[0][0])
    t = float(T[r])
    if not np.isfinite(t) or t <= 0 or float(T.min()) <= 0:
        return None                      # P1: a negative mean time is a broken instrument
    return float(np.log(t))


def am_A(g1, g2, omegas=(120, 180, 240, 300), rel=0.02):
    """A from the convergence of d ln T / d Omega. Returns (A, series) or (None, series)."""
    pts = []
    for om in omegas:
        lt = am_ln_mfpt(g1, g2, om)
        if lt is not None:
            pts.append((om, lt))
    loc = [((pts[i + 1][1] - pts[i][1]) / (pts[i + 1][0] - pts[i][0]))
           for i in range(len(pts) - 1)]
    if len(loc) < 2 or abs(loc[-1] - loc[-2]) / abs(loc[-1]) > rel:
        return None, loc
    return loc[-1], loc


# ----------------------------------------------------------------- verdict rules
# Factored out and unit-tested on engineered data BEFORE this ran. Rule 19: for each, name the
# data that makes it print the other verdict, and check that data is the thing it means.


def spread_verdict(vals, thresh=2.0):
    """P2/P3: does A move on the affinity level set? Other verdict <- a constant A.

    An empty or singleton family is NOT evidence of constancy -- it is no measurement, and
    max()/min() on it either throws or returns 1.0, which reads as "constant" (§71 P2 printed
    HOLDS off a None). It returns None, and the caller must say so.
    """
    vals = [v for v in vals if v is not None and np.isfinite(v) and v > 0]
    if len(vals) < 3:
        return None, None
    s = max(vals) / min(vals)
    return s, bool(s > thresh)


def both_ways_verdict(affs, min_n=2, min_span=0.2):
    """P4: at equal A, does the affinity still vary? Other verdict <- affinity pinned by A.

    With fewer than two solved points the span is 0 by construction and would print the
    "A bounds the drive" verdict from no data at all.
    """
    if len(affs) < min_n:
        return None, None
    s = max(affs) - min(affs)
    return s, bool(s > min_span)


# ----------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/action_is_not_priced.json"))
    args = ap.parse_args()
    out = {}

    F0 = schlogl_affinity(0.1, 1.9)
    print("=== P2: Schloegl, EXACT. Roots on the affinity level set, saddle pinned at r2 = 1")
    print(f"    affinity held at ln(e1 e2 / e3) = {F0:.6f} for every row")
    print(f"{'r1':>8}{'r3':>10}{'affinity':>12}{'A':>12}{'quad err':>11}{'Delta':>9}")
    sch = []
    for r1 in (0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5):
        try:
            r3 = r3_at_affinity(r1, F0)
        except ValueError:
            print(f"{r1:>8}    no root on the level set")
            continue
        A, err = schlogl_A(r1, r3)
        sch.append({"r1": r1, "r3": r3, "aff": schlogl_affinity(r1, r3), "A": A,
                    "delta": (r3 - r1) / 2})
        print(f"{r1:>8}{r3:>10.4f}{schlogl_affinity(r1, r3):>12.6f}{A:>12.6f}"
              f"{err:>11.1e}{(r3 - r1) / 2:>9.4f}")
    spread, moved = spread_verdict([s["A"] for s in sch])
    affspread = max(s["aff"] for s in sch) - min(s["aff"] for s in sch)
    out["schlogl"] = sch
    print(f"  affinity varies by {affspread:.2e} nats across the family (held fixed by "
          f"construction)")
    if spread is None:
        print("  -> P2 UNDECIDED: fewer than three usable landscapes on the level set")
    else:
        print(f"  A spans {min(s['A'] for s in sch):.6f} .. {max(s['A'] for s in sch):.6f}"
              f"  -- a factor of {spread:.0f}")
        print(f"  -> P2 {'HOLDS: A is NOT priced by the affinity' if moved else 'FAILS: A is constant on the level set -- reliability HAS a thermodynamic price'}")

    print("\n=== P1/P3: AM, 2-D. Reverse ratios on the level set g1 g2^2 = const")
    C = 0.35 ** 3
    print(f"    affinity held at -ln(g1 g2^2) = {-np.log(C):.6f} for every row")
    print(f"{'g1':>8}{'g2':>7}{'affinity':>11}{'A':>10}{'Delta*':>9}  local A series")
    am, excluded = [], []
    for g2 in (0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60):
        g1 = C / g2 ** 2
        if not 0.02 < g1 <= 1.0:
            continue
        A, loc = am_A(g1, g2)
        ser = ", ".join(f"{v:.4f}" for v in loc)
        if A is None:
            excluded.append((round(g1, 4), g2, ser))
            print(f"{g1:>8.4f}{g2:>7.2f}{am_affinity(g1, g2):>11.6f}{'EXCLUDED':>10}"
                  f"{'':>9}  {ser}")
            continue
        x = am_rail(g1, g2)
        am.append({"g1": g1, "g2": g2, "aff": am_affinity(g1, g2), "A": A,
                   "delta": float(x[0] - x[1]), "series": loc})
        print(f"{g1:>8.4f}{g2:>7.2f}{am_affinity(g1, g2):>11.6f}{A:>10.5f}"
              f"{x[0] - x[1]:>9.4f}  {ser}")
    out["am"], out["am_excluded"] = am, excluded
    amspread, ammoved = spread_verdict([a["A"] for a in am])
    print(f"  cells excluded on P1 (non-convergent or negative mean time): {len(excluded)}")
    if amspread is None:
        print(f"  -> P3 UNDECIDED: only {len(am)} cells survived the P1 gate; AM cannot"
              " settle this and the 1-D arm carries the section")
    else:
        amaff = max(a["aff"] for a in am) - min(a["aff"] for a in am)
        print(f"  affinity varies by {amaff:.2e} nats across the family")
        print(f"  A spans {min(a['A'] for a in am):.5f} .. {max(a['A'] for a in am):.5f}"
              f"  -- a factor of {amspread:.1f}")
        print(f"  -> P3 {'HOLDS on a second substrate too' if ammoved else 'FAILS: A is constant on AMs level set'}")

    print("\n=== P4: the opposite sweep -- same A, DIFFERENT affinity?")
    target = sch[1]["A"]
    print(f"    target A = {target:.6f} (the reference element 0.1/1.0/1.9)")
    print(f"{'r1':>8}{'r3':>10}{'A':>12}{'affinity':>12}")
    pairs = []
    for r1 in (0.05, 0.08, 0.1, 0.12, 0.15):
        try:
            r3 = brentq(lambda t: schlogl_A(r1, t)[0] - target, 1.0001, 60.0)
        except ValueError:
            continue
        pairs.append({"r1": r1, "r3": r3, "A": schlogl_A(r1, r3)[0],
                      "aff": schlogl_affinity(r1, r3)})
        print(f"{r1:>8}{r3:>10.4f}{schlogl_A(r1, r3)[0]:>12.6f}"
              f"{schlogl_affinity(r1, r3):>12.6f}")
    out["same_A"] = pairs
    fspan, bothways = both_ways_verdict([p["aff"] for p in pairs])
    if fspan is None:
        print(f"  -> P4 UNDECIDED: only {len(pairs)} element(s) hit the target A")
    else:
        aspan = max(p["A"] for p in pairs) / min(p["A"] for p in pairs)
        print(f"  A is equal to {100*(aspan-1):.4f}% across these;"
              f" affinity spans {fspan:.4f} nats")
        print(f"  -> P4 {'HOLDS: the dissociation runs BOTH ways -- A does not bound the drive either' if bothways else 'FAILS: equal A forces equal affinity'}")

    print("\n=== P5 (SUSPECT, rule 17): does the landscape headline number Delta determine A?")
    print("    NOTE: this swept r1 over 0.05..0.35 first -- a 1.4x move in the saddle's")
    print("    relative position -- and printed 'Delta alone reproduces A' off a 3% span.")
    print("    Widened to the full admissible range 0.02 < r1 < 0.97 below (§82.1).")
    print(f"{'Delta':>8}{'r1':>7}{'r3':>9}{'saddle frac':>13}{'A':>12}{'affinity':>11}")
    eq, dspans = [], []
    for D in (0.9, 1.5164, 3.0):
        rows = []
        for r1 in (0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.9, 0.97):
            r3 = r1 + 2 * D
            if not r1 < 1.0 < r3:
                continue
            A, _ = schlogl_A(r1, r3)
            frac = (1.0 - r1) / (r3 - r1)
            rows.append({"delta": D, "r1": r1, "r3": r3, "frac": frac, "A": A,
                         "aff": schlogl_affinity(r1, r3)})
            print(f"{D:>8}{r1:>7}{r3:>9.4f}{frac:>13.4f}{A:>12.6f}"
                  f"{schlogl_affinity(r1, r3):>11.5f}")
        eq += rows
        As = [r["A"] for r in rows]
        fr = [r["frac"] for r in rows]
        ff = [r["aff"] for r in rows]
        dspans.append(max(As) / min(As))
        print(f"    -> at Delta = {D}: saddle position moved {max(fr)/min(fr):.1f}x,"
              f" A spans {max(As)/min(As):.3f}x, affinity spans {max(ff)-min(ff):.3f} nats")
    out["equal_delta"] = eq
    worst = max(dspans)
    print(f"  worst A span at IDENTICAL Delta: {worst:.3f}x"
          f"   (against {spread:.0f}x at identical AFFINITY)")
    if worst > 1.2:
        print("  -> P5 as predicted: Delta does NOT determine A either -- but it MISSES by a")
        print(f"     factor of {worst:.2f} where the affinity misses by {spread:.0f}. The suspect")
        print("     remains 'A is a functional of the whole field'; it is not yet confirmed.")
    else:
        print("  -> P5 Delta alone reproduces A across the full range -- the suspect survives")
        print("     this test and is STILL a suspect (rule 17), not a mechanism.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
