"""T-COST-e: does the cost optimum survive real Arrhenius kinetics, or is gamma the only job T has?

§38 and §39 locate an optimal drive gamma* ~ 0.20-0.24 by two computationally independent
routes. `crnl/cooling.py` maps gamma to a temperature, gamma = exp(-dE/T) -- and its own
docstring flags what it leaves out:

    "The forward rates have no temperature dependence, so this is the minimal change
     that lets the balance move, not a thermochemistry."

Real forward rates are Arrhenius too. **If temperature moves the forward rates as well as
the balance, the optimum in T is not simply the optimum in gamma re-labelled**, and the
"optimal temperature" reading of §38 -- which is the whole basis for asking which physical
substrate a restoring element should be built from -- is unearned.

WHAT T ACTUALLY CONTROLS ONCE EVERY RATE IS ARRHENIUS. Write k_r = A exp(-E_r/T). Then:

  * A UNIFORM activation energy is a pure clock change. Q -> lambda Q and sigma -> lambda
    sigma together, so Sigma = Q^-1 sigma is EXACTLY invariant. Temperature acting on all
    forward rates equally cannot move the cost optimum -- §5.1's uniform-order argument one
    level up. **So the Arrhenius upgrade has content only if the channels differ.**
  * The channels do differ, physically: AM's disagreement X+Y->2B and its recruitment
    B+X->2X need not share an activation energy, and in any real implementation they would
    not. Their ratio

        rho(T) = k_dis / k_rec = exp(-dEa/T),      dEa = E_dis - E_rec

    is a landscape change, not a clock change -- **and this rig has never varied it.**
    Every result in FINDINGS is at rho = 1.

So with genuine Arrhenius kinetics T has TWO jobs -- it sets the drive gamma AND the
kinetic asymmetry rho -- and the optimal temperature is a compromise between two optima
that need not agree. That is the claim to test.

THE GEOMETRY, derived by hand and checked before use. The off-symmetry fixed point still
needs b* = gamma/(1+gamma) (§30's identity is untouched: rho multiplies the disagreement
channel, which cancels in the difference), so with s = 1/(1+gamma),

    delta*(gamma, rho)^2 = (rho - gamma - 4 rho gamma^3) / [ (1+gamma)^2 (rho - gamma) ]

reducing to `delta_star(gamma)` at rho = 1. **It has a trap and the trap was hit.** For
rho < gamma both numerator and denominator are negative, so the expression returns a
plausible POSITIVE value -- 0.806 at (0.25, 0.05), and 6.37 and 11.4 near rho = gamma,
where delta* <= 1 is a hard bound -- in a region where the true landscape is absent
(checked against the ODE nullcline: 0.0). The landscape exists iff

    rho > rho_c(gamma) = gamma / (1 - 4 gamma^3),   and rho_c > gamma always,

so the single guard rho > rho_c covers it. Cells failing it are DROPPED AND REPORTED.

COST. §38's G = Sigma / (Omega ln(theta/eps)), k_B per molecule per e-fold of gain, with
theta and eps as fractions of delta* so the gain is constant across the sweep by
construction (the scaling whose absence killed §9.2). Directly comparable to §39's table.

PREDICTIONS, written before running:

  P1  GATE, two parts. (a) Multiplying every rate by lambda leaves Sigma unchanged to
      solver precision -- analytic, so a failure means the cost measure is clock-dependent
      and §37-§39 are all in question. (b) At rho = 1 the gamma-sweep reproduces §38's
      gamma* ~ 0.20 and §39's 0.240. Nothing below is admissible if either fails.
  P2  At dEa = 0, T is a pure relabelling of gamma, so T* = dE/ln(1/gamma*) identically,
      and with T_c = dE/ln 2 the death of the landscape, **T*/T_c = ln2/ln(1/gamma*) ~
      0.43-0.49: the optimum sits near 46% of the critical temperature.** A check on the
      mapping, not new physics -- but if it fails the mapping is wrong.
  P3  THE MECHANISM TEST, and it is the cheap one: sweep rho at FIXED gamma, no
      temperature anywhere. **Cost has an interior minimum in rho, at rho* < 1.** Reason:
      the disagreement reaction moves delta by EXACTLY ZERO -- that is §30's first
      cancellation, it consumes both species equally -- so it produces entropy and no
      signal, and raising rho should raise cost; but rho -> rho_c destroys the landscape.
      An interior optimum below 1 would mean **AM's own convention k_dis = k_rec is not
      the cheapest kinetics**, which no section has ever questioned. If rho* = 1 to
      resolution, the convention is optimal and that is the more surprising result.
  P4  THE OPPOSITE SWEEP (the confound-breaking rule). sign(dEa) decides whether warming
      raises or lowers rho. Run dEa > 0 and dEa < 0 at the same dE. **If T* moves in
      OPPOSITE directions from the dEa = 0 baseline, the rho channel is real and separable
      from the gamma channel.** If T* is unmoved under both, temperature acts only through
      gamma, the Arrhenius upgrade is cosmetic, and cooling.py's minimal model was
      sufficient after all -- reported as the real outcome it would be.
  P5  CONTROL, so a cost effect cannot be a geometry rescaling in disguise. delta* varies
      only weakly with rho (measured before writing this: 0.730 -> 0.773 over rho = 0.4 ->
      5 at gamma = 0.25, a 6% span against cost changes expected to be far larger).
      delta* is reported in every row so the reader can see the geometry barely moves.
  P6  ARITHMETIC, labelled as arithmetic and not as a measurement. §16 pins the cycle
      affinity at A = 3 ln(1/gamma) (in `verify_base`), so gamma* ~ 0.20-0.24 gives
      A* = 4.28-4.83 k_B T, i.e. an optimal fuel drop of ~0.11-0.12 eV at 300 K against
      ATP's ~20 k_B T. **This is NOT merged with §40's Q_min = 5.39.** Two numbers near 5
      arrived along axes chosen for other reasons, and rule 9 exists because that has
      already failed three times here.

=============================================================================
SECOND PASS. The first pass ran with P1-P6 as above; its verdicts are kept in FINDINGS
whatever this pass says. Two of them need re-analysis, and the criteria below are fixed
HERE, before re-running, so the filter cannot be chosen for the answer it gives.

P3 WAS REFUTED, AND IN THE OPPOSITE DIRECTION. G falls monotonically as rho rises, over
the whole grid rho = 0.35..6 at three gammas and two Omegas -- no interior optimum, the
minimum pinned at the upper edge. The reasoning behind P3 conflated FLUX with
DISSIPATION: the disagreement channel does move delta by zero, but driving a reaction
fast in BOTH directions pushes it toward local equilibrium, and a locally equilibrated
channel carries large flux at small NET entropy production. Raising rho makes the
non-signal channel cheap, not expensive.

P4's PRE-REGISTERED VERDICT WAS "FAILS", AND ITS ARGMIN IS NOT ADMISSIBLE. In the
dEa = +0.6 arm the reported T* = 1.114 sits at delta* = 0.227 with rho/rho_c = 1.04 --
the landscape is 4% from death. **§9.2 was withdrawn for exactly this**: G falls there
not because the decision is cheap but because it is SMALL, the threshold theta*delta**Omega
having collapsed with delta*. The dEa != 0 sweeps drive delta* over a 4x range, far wider
than any previous section, and G is not comparable across it.

  P7  DIAGNOSIS. Admissible cells are those with **delta* >= 0.40**, chosen because the
      rho = 1 gamma-sweep that P1b validates against §38 itself spans delta* = 0.40..0.97
      -- so that is the range over which this instrument is known to agree with the
      established result, and nowhere else. Under that filter the dEa = +0.6 arm should
      show its interior minimum near T ~ 0.57, BELOW the dEa = 0 baseline of 0.637.
  P8  With the filter, T* splits by sign(dEa) in OPPOSITE directions, so P4's substance
      holds even though its pre-registered form failed. **If the filtered +0.6 optimum is
      still ABOVE the baseline, P4 is refuted for real**, temperature acts through gamma
      alone, and cooling.py's minimal model was sufficient -- reported as such.
  P9  P3's REPLACEMENT HYPOTHESIS, with its own kill test (rule 17 -- the mechanism above
      is a suspect, not a result). If cheapness at large rho is local equilibration of the
      disagreement channel, then **G must ASYMPTOTE as rho -> infinity**, to the reduced
      model in which X+Y <-> 2B is equilibrated and recruitment carries the whole net
      current. Extend rho to 100. If G instead keeps falling without bound the mechanism
      is wrong again and the edge optimum needs a different account.
  P10 THE SKEPTICAL CHECK on P9, because G prices dissipation and not time, and §39.1
      found the cost residual is entirely a TIME residual. A higher rate constant is
      thermodynamically free -- catalysts do not change dG -- so if large rho also lowers
      the mean decision time, it is a free win on both axes with no trade-off, which is a
      strong claim and must be shown rather than assumed. Mean first-passage time is
      reported in every row. **If time RISES with rho, the "free win" reading is wrong
      and the cheapness is bought with slowness.**
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
import scipy.sparse.linalg as spla

from crnl.cme import enumerate_states, first_passage, generator
from crnl.networks.am_reversible import GAMMA_C, delta_star, reverse_pairing
from crnl.reactions import Reaction, ReactionNetwork
from crnl.vectorized import compile_network
from experiments.cost_of_reliability import sigma_local
from experiments.slaving_axis import slaved


# P7's admissibility floor, fixed before the second pass: the rho = 1 gamma-sweep that
# P1b validates against §38 spans delta* = 0.40..0.97, so that is the range over which
# this instrument is known to agree with the established result.
DSTAR_MIN = 0.40


def am_rho(gamma: float, rho: float, k: float = 1.0) -> ReactionNetwork:
    """AM with the disagreement channel at rate rho*k and recruitment at k.

    rho = k_dis/k_rec is the kinetic asymmetry between the two forward channels, which
    every previous section has held at 1. Reverse rates stay at gamma times their own
    forward, so the cycle affinity is 3 ln(1/gamma) exactly as in §16 -- rho changes the
    kinetics, not the drive.
    """
    kd, kr = rho * k, k
    return ReactionNetwork(
        species=["X", "Y", "B"],
        reactions=[
            Reaction({"X": 1, "Y": 1}, {"B": 2}, kd, name="f1:X+Y->2B"),
            Reaction({"B": 1, "X": 1}, {"X": 2}, kr, name="f2:B+X->2X"),
            Reaction({"B": 1, "Y": 1}, {"Y": 2}, kr, name="f3:B+Y->2Y"),
            Reaction({"B": 2}, {"X": 1, "Y": 1}, gamma * kd, name="r1:2B->X+Y"),
            Reaction({"X": 2}, {"B": 1, "X": 1}, gamma * kr, name="r2:2X->B+X"),
            Reaction({"Y": 2}, {"B": 1, "Y": 1}, gamma * kr, name="r3:2Y->B+Y"),
        ],
        name=f"am-rho-g{gamma}-r{rho}",
    )


def rho_critical(gamma: float) -> float:
    """Below this rho the landscape is absent. Always > gamma, which is why one guard
    suffices for the double sign flip documented in the module docstring."""
    if gamma <= 0.0:
        return 0.0
    denom = 1.0 - 4.0 * gamma ** 3
    return float("inf") if denom <= 0.0 else gamma / denom


def delta_star_rho(gamma: float, rho: float) -> float:
    """Attractor separation at kinetic asymmetry rho; 0 where there is no landscape.

    Verified against the ODE nullcline at 30 (gamma, rho) pairs to 1e-9, and against
    `delta_star(gamma)` at rho = 1 to machine precision, BEFORE this was used anywhere.
    """
    if gamma >= GAMMA_C and rho <= 1.0:
        return 0.0
    if rho <= rho_critical(gamma):
        return 0.0
    v = (rho - gamma - 4.0 * rho * gamma ** 3) / ((1.0 + gamma) ** 2 * (rho - gamma))
    if not (v > 0.0):
        return 0.0
    ds = float(np.sqrt(v))
    # delta <= s = 1/(1+gamma) is a hard bound; a violation means the guard is wrong
    return ds if ds < 1.0 / (1.0 + gamma) else 0.0


def cost_cell(gamma, rho, omega, eps, theta):
    """§38's G = Sigma/(Omega ln(theta/eps)), on the rho-generalised network.

    Mirrors `cost_of_reliability.cell` exactly except for the network and delta*, so the
    rho = 1 column is comparable to §37-§40 row for row.
    """
    ds = delta_star_rho(gamma, rho)
    if ds <= 0.0:
        return None
    net = am_rho(gamma, rho)
    x0 = eps * ds
    st = slaved(net, x0)
    if st is None:
        return None
    nb = int(round(st[2] * omega))
    rest = omega - nb
    d0 = max(1, int(round(x0 * omega)))
    if (rest - d0) % 2:
        d0 -= 1
    if d0 < 1 or rest - d0 < 0:
        return None
    n0 = np.array([(rest + d0) // 2, (rest - d0) // 2, nb], dtype=np.int64)
    thr = max(2, int(round(theta * ds * omega)))
    if thr <= d0:
        return None

    states, index = enumerate_states(3, int(omega))
    absorb = np.array([abs(int(s[0]) - int(s[1])) >= thr for s in states])
    Q = generator(net, int(omega), float(omega))
    tr = np.where(~absorb)[0]
    tmap = {int(i): r for r, i in enumerate(tr)}
    Qtt = Q[tr][:, tr].tocsr()
    si = tmap[index[tuple(n0)]]

    fav = np.array([int(s[1]) > int(s[0]) for s in states])[absorb].astype(float)
    b = -(Q[tr][:, np.where(absorb)[0]].tocsr() @ fav)
    p_err = float(spla.spsolve(Qtt, b)[si])

    comp = compile_network(net, float(omega))
    sig = sigma_local(net, comp, states, reverse_pairing(net))[tr]
    Sig = float(spla.spsolve(Qtt, -sig)[si])
    if not np.isfinite(Sig) or Sig <= 0:
        return None

    fp = first_passage(net, int(omega), float(omega), n0,
                       lambda s, t=thr: abs(int(s[0]) - int(s[1])) >= t)
    mt = float(fp["mean_time"]) if fp.get("valid", True) else float("nan")

    gain = np.log(theta / eps)
    return {"gamma": gamma, "rho": rho, "omega": omega, "delta_star": ds,
            "p_err": p_err, "L": float(-np.log(max(p_err, 1e-300))), "Sigma": Sig,
            "G": Sig / (omega * gain), "thr": int(thr), "mean_time": mt,
            "admissible": bool(ds >= DSTAR_MIN)}


def argmin_parabolic(xs, ys):
    """Vertex of a parabola through the three points around the discrete minimum."""
    i = int(np.argmin(ys))
    if i == 0 or i == len(xs) - 1:
        return float(xs[i]), False
    x0, x1, x2 = xs[i - 1], xs[i], xs[i + 1]
    y0, y1, y2 = ys[i - 1], ys[i], ys[i + 1]
    d = (x0 - x1) * (x0 - x2) * (x1 - x2)
    if abs(d) < 1e-300:
        return float(x1), True
    a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / d
    b = (x2 ** 2 * (y0 - y1) + x1 ** 2 * (y2 - y0) + x0 ** 2 * (y1 - y2)) / d
    return (float(-b / (2 * a)) if abs(a) > 1e-300 else float(x1)), True


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omegas", type=int, nargs="+", default=[200, 300])
    ap.add_argument("--eps", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--dE", type=float, default=1.0, help="fuel drop per reaction, sets gamma(T)")
    ap.add_argument("--dEas", type=float, nargs="+", default=[0.0, 0.6, -0.6])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/arrhenius_optimum.json"))
    args = ap.parse_args()

    t0 = time.time()
    eps, th, dE = args.eps, args.theta, args.dE
    out = {"eps": eps, "theta": th, "dE": dE}

    # ---- P1a: clock invariance, analytic and exact -------------------------------
    print("=== P1a GATE: is Sigma invariant under a uniform rate rescale?")
    base = cost_cell(0.25, 1.0, 200, eps, th)
    print(f"{'lambda':>9}{'Sigma':>14}{'rel dev':>12}")
    worst = 0.0
    for lam in (0.1, 1.0, 7.0, 100.0):
        net = am_rho(0.25, 1.0, k=lam)
        ds = delta_star_rho(0.25, 1.0)
        st = slaved(net, eps * ds)
        nb = int(round(st[2] * 200)); rest = 200 - nb
        d0 = max(1, int(round(eps * ds * 200)))
        if (rest - d0) % 2:
            d0 -= 1
        n0 = np.array([(rest + d0) // 2, (rest - d0) // 2, nb], dtype=np.int64)
        thr = max(2, int(round(th * ds * 200)))
        states, index = enumerate_states(3, 200)
        absorb = np.array([abs(int(s[0]) - int(s[1])) >= thr for s in states])
        Q = generator(net, 200, 200.0)
        tr = np.where(~absorb)[0]
        tmap = {int(i): r for r, i in enumerate(tr)}
        comp = compile_network(net, 200.0)
        sig = sigma_local(net, comp, states, reverse_pairing(net))[tr]
        Sg = float(spla.spsolve(Q[tr][:, tr].tocsr(), -sig)[tmap[index[tuple(n0)]]])
        dev = abs(Sg - base["Sigma"]) / base["Sigma"]
        worst = max(worst, dev)
        print(f"{lam:>9.1f}{Sg:>14.6f}{dev:>12.2e}")
    p1a = worst < 1e-9
    print(f"  -> P1a {'HOLDS' if p1a else 'FAILS -- the cost measure is clock-dependent'}"
          f"  (worst {worst:.2e})")
    out["p1a_worst"] = worst

    # ---- P1b: reproduce §38/§39's gamma* at rho = 1 -------------------------------
    print(f"\n=== P1b GATE: gamma-sweep at rho = 1 must reproduce gamma* ~ 0.20-0.24")
    gammas = [0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.32, 0.36, 0.40]
    print(f"{'gamma':>7}" + "".join(f"{f'G(Om={o})':>13}" for o in args.omegas)
          + f"{'delta*':>10}")
    grows, gstars = [], {}
    for om in args.omegas:
        gstars[om] = []
    for g in gammas:
        line, ds = f"{g:>7.2f}", delta_star_rho(g, 1.0)
        for om in args.omegas:
            c = cost_cell(g, 1.0, om, eps, th)
            line += f"{c['G']:>13.4f}" if c else f"{'--':>13}"
            gstars[om].append(c["G"] if c else np.nan)
            if c:
                grows.append(c)
        print(line + f"{ds:>10.4f}")
    p1b = {}
    for om in args.omegas:
        ys = np.array(gstars[om])
        ok = np.isfinite(ys)
        gs, _ = argmin_parabolic(np.array(gammas)[ok], ys[ok])
        p1b[om] = gs
        print(f"  Omega={om}: gamma* = {gs:.4f}")
    p1b_ok = all(0.14 < v < 0.32 for v in p1b.values())
    print(f"  §38 CME 0.20, §39 closed form 0.240  -> P1b "
          f"{'HOLDS' if p1b_ok else 'FAILS'}")
    out["gamma_star_rho1"] = {str(k): v for k, v in p1b.items()}

    # ---- P3/P5: the rho axis at fixed gamma, no temperature -----------------------
    print(f"\n=== P3/P5/P9/P10: sweep rho at FIXED gamma -- is k_dis = k_rec optimal?")
    rhos = [0.35, 0.5, 0.7, 0.85, 1.0, 1.3, 1.8, 2.5, 4.0, 6.0,
            10.0, 20.0, 50.0, 100.0]      # P9 extends the first pass's 0.35..6
    rho_star = {}
    rrows = []
    for g in (0.16, 0.24, 0.32):
        print(f"\n  gamma = {g}   (rho_c = {rho_critical(g):.4f})")
        print(f"{'rho':>8}" + "".join(f"{f'G(Om={o})':>12}" for o in args.omegas)
              + f"{'L(Om=200)':>11}{'time(200)':>11}{'delta*':>9}")
        cols = {om: [] for om in args.omegas}
        xs, dropped = [], []
        for r in rhos:
            ds = delta_star_rho(g, r)
            if ds <= 0:
                dropped.append(r)
                continue
            line, ok, first = f"{r:>8.2f}", True, None
            for om in args.omegas:
                c = cost_cell(g, r, om, eps, th)
                if c is None:
                    ok = False
                    line += f"{'--':>12}"
                else:
                    cols[om].append(c["G"])
                    rrows.append(c)
                    line += f"{c['G']:>12.4f}"
                    if first is None:
                        first = c
            if ok:
                xs.append(r)
            extra = (f"{first['L']:>11.2f}{first['mean_time']:>11.3f}"
                     if first else f"{'--':>11}{'--':>11}")
            print(line + extra + f"{ds:>9.4f}")
        if dropped:
            print(f"    dropped (no landscape, reported not hidden): rho = {dropped}")
        for om in args.omegas:
            if len(cols[om]) == len(xs) and len(xs) >= 3:
                rs, interior = argmin_parabolic(np.array(xs), np.array(cols[om]))
                rho_star[(g, om)] = (rs, interior)
                y = np.array(cols[om])
                tail = 100.0 * (y[-2] - y[-1]) / y[-1] if len(y) > 2 else float("nan")
                print(f"    Omega={om}: rho* = {rs:.4f}"
                      f"  ({'interior' if interior else 'AT THE EDGE'})"
                      f"   last step G falls {tail:.2f}%  ->"
                      f" {'ASYMPTOTING' if tail < 2.0 else 'still falling'}")
    out["rho_star"] = {f"{g}_{om}": v for (g, om), v in rho_star.items()}

    interior_below = [v for v, i in rho_star.values() if i and v < 1.0]
    interior_any = [v for v, i in rho_star.values() if i]
    print(f"\n  P3: {len(interior_any)}/{len(rho_star)} cells have an INTERIOR rho*;"
          f" {len(interior_below)} of those are below 1")
    if interior_any:
        print(f"      rho* range {min(interior_any):.3f} .. {max(interior_any):.3f}"
              f"   -> AM's k_dis = k_rec is "
              f"{'NOT optimal' if interior_below else 'optimal or beaten from above'}")

    # ---- P2/P4: temperature sweeps, with the opposite sweep ----------------------
    print(f"\n=== P2/P4: sweep T, with gamma(T) = exp(-{dE}/T) and rho(T) = exp(-dEa/T)")
    Tc = dE / np.log(2.0)
    print(f"  T_c = dE/ln2 = {Tc:.4f} (landscape death at gamma = 1/2)")
    Ts = [round(x, 4) for x in np.linspace(0.30, 1.25, 15)]
    tstar, tstar_filt = {}, {}
    trows = []
    for dEa in args.dEas:
        print(f"\n  dEa = {dEa:+.2f}   ({'no rho channel' if dEa == 0 else 'rho moves with T'})")
        print(f"{'T':>7}{'T/Tc':>7}{'gamma':>8}{'rho':>8}"
              + "".join(f"{f'G(Om={o})':>12}" for o in args.omegas)
              + f"{'L(200)':>9}{'delta*':>9}{'adm?':>6}")
        cols = {om: [] for om in args.omegas}
        xs, adm, dropped = [], [], []
        for T in Ts:
            g = float(np.exp(-dE / T))
            r = float(np.exp(-dEa / T))
            ds = delta_star_rho(g, r)
            if ds <= 0:
                dropped.append(round(T, 4))
                continue
            line, ok, first = f"{T:>7.3f}{T/Tc:>7.3f}{g:>8.4f}{r:>8.4f}", True, None
            for om in args.omegas:
                c = cost_cell(g, r, om, eps, th)
                if c is None:
                    ok = False
                    line += f"{'--':>12}"
                else:
                    cols[om].append(c["G"])
                    c2 = dict(c); c2["T"] = T; c2["dEa"] = dEa
                    trows.append(c2)
                    line += f"{c['G']:>12.4f}"
                    if first is None:
                        first = c
            if ok:
                xs.append(T)
                adm.append(ds >= DSTAR_MIN)
            print(line + (f"{first['L']:>9.2f}" if first else f"{'--':>9}")
                  + f"{ds:>9.4f}{'yes' if ds >= DSTAR_MIN else 'NO':>6}")
        if dropped:
            print(f"    dropped (no landscape): T = {dropped}")
        for om in args.omegas:
            if len(cols[om]) == len(xs) and len(xs) >= 3:
                X, Y, A = np.array(xs), np.array(cols[om]), np.array(adm)
                ts, interior = argmin_parabolic(X, Y)
                tstar[(dEa, om)] = (ts, interior)
                msg = (f"    Omega={om}: pre-registered T* = {ts:.4f}"
                       f" (T/Tc {ts/Tc:.3f}, {'interior' if interior else 'AT THE EDGE'})")
                if A.sum() >= 3:
                    tf, fi = argmin_parabolic(X[A], Y[A])
                    tstar_filt[(dEa, om)] = (tf, fi)
                    msg += (f"   |  delta*>={DSTAR_MIN} T* = {tf:.4f}"
                            f" (T/Tc {tf/Tc:.3f}, {'interior' if fi else 'at edge'})")
                else:
                    msg += f"   |  too few admissible cells ({int(A.sum())})"
                print(msg)
    out["T_star"] = {f"{d}_{om}": v for (d, om), v in tstar.items()}
    out["T_star_filtered"] = {f"{d}_{om}": v for (d, om), v in tstar_filt.items()}
    out["Tc"] = Tc
    out["dstar_min"] = DSTAR_MIN

    print(f"\n=== P2: does the dEa = 0 optimum sit where gamma* says it should?")
    for om in args.omegas:
        if (0.0, om) in tstar and om in p1b:
            ts = tstar[(0.0, om)][0]
            pred = dE / np.log(1.0 / p1b[om])
            print(f"  Omega={om}: T* = {ts:.4f}   predicted dE/ln(1/gamma*) = {pred:.4f}"
                  f"   ratio {ts/pred:.4f}   T*/T_c = {ts/Tc:.4f}")

    print(f"\n=== P4/P8: the opposite sweep -- does sign(dEa) move T* in opposite directions?")
    print("  BOTH criteria reported. P4 is the pre-registered one; P8 restricts to the")
    print(f"  delta* >= {DSTAR_MIN} band where P1b validates this instrument against §38.")
    out["verdicts"] = {}
    for tag, table in (("P4 (pre-registered, all cells)", tstar),
                       (f"P8 (delta* >= {DSTAR_MIN})", tstar_filt)):
        print(f"\n  --- {tag}")
        verdict = []
        for om in args.omegas:
            if not all((d, om) in table for d in (0.0, 0.6, -0.6)):
                print(f"    Omega={om}: incomplete, no verdict")
                continue
            b = table[(0.0, om)][0]
            up, dn = table[(0.6, om)][0], table[(-0.6, om)][0]
            print(f"    Omega={om}: T*(0) = {b:.4f};  +0.6 -> {up:.4f} "
                  f"({100*(up-b)/b:+.1f}%);  -0.6 -> {dn:.4f} ({100*(dn-b)/b:+.1f}%)")
            verdict.append((up - b) * (dn - b) < 0
                           and max(abs(up - b), abs(dn - b)) / b > 0.02)
        ok = bool(verdict) and all(verdict)
        out["verdicts"][tag] = ok
        if not verdict:
            continue
        if ok:
            print("    -> HOLDS: T* splits by sign(dEa). The rho channel is real and")
            print("       separable, so T does NOT act through gamma alone and cooling.py's")
            print("       minimal model is not sufficient for an optimal-temperature claim.")
        else:
            print("    -> FAILS: T* does not split by sign(dEa) under this criterion.")

    print(f"\n=== P6 (ARITHMETIC on §16's affinity + §38, not a measurement)")
    for om in args.omegas:
        if om in p1b:
            A = 3.0 * np.log(1.0 / p1b[om])
            print(f"  Omega={om}: gamma* = {p1b[om]:.4f} -> A* = 3ln(1/gamma*) ="
                  f" {A:.3f} k_B T  -> {A*0.025852:.4f} eV at 300 K")
    print("  ATP hydrolysis ~20 k_B T. NOT merged with §40's Q_min = 5.39 (rule 9).")

    out["rows"] = {"gamma_sweep": grows, "rho_sweep": rrows, "T_sweep": trows}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
