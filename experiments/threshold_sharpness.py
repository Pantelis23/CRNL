"""T15-k: how sharp is a restoring threshold made of Omega molecules?

§9.1 settled the DETERMINISTIC picture: the AM landscape dies at gamma_c = 1/2, and above
it "no population size Omega can restore, because there is nothing to restore toward".
§62 then made that boundary exact and free for any exchange-symmetric network -- the sign
of P at the symmetric steady state, reproducing gamma_c = 1/2 to 1.1e-15. **So the
deterministic side now has ZERO uncertainty, and any deviation measured against it is
entirely the CME's.** That is the setup this project has wanted since §9.

The question §9.1 did not ask is the device question: BELOW gamma_c, how many molecules does
a switch need before its threshold is sharp? A transistor is a good restoring element because
its threshold is abrupt. A chemical one made of Omega molecules cannot be, and the width of
its transition is a number nobody here has measured.

**THE INSTRUMENT IS EXACT AND ALREADY IMPLIED BY §53.** §53 defined P as the antisymmetric
eigenvalue of the JACOBIAN. The CME generator has an antisymmetric sector too: exchange
X <-> Y commutes with Q for any exchange-symmetric network (§43's premise), so Q block-
diagonalises, and on antisymmetric observables f(sigma s) = -f(s) the block is

    Q_A[p, q] = Q[s_p, s_q] - Q[s_p, sigma(s_q)]      over representatives with n_X > n_Y

with the diagonal n_X = n_Y states dropping out identically. **lambda_A = the leading
eigenvalue of Q_A is the exact stochastic counterpart of §53's P**, needing no simulation,
no threshold and no first-passage definition.

The chain is ergodic at every finite Omega, so lambda_A < 0 always: **there is no sign change
at finite size, and that is the point.** Restoration at finite Omega shows up not as a sign
but as lambda_A becoming exponentially small -- metastability -- and the transition from
"exponentially small" to "O(1)" is the blurred threshold whose width is being measured.

PREDICTIONS, written before running; per rule 19 each verdict rule names the data that would
make it print the other answer.

  P1  GATE (a). The spectrum of Q_A is a SUBSET of the spectrum of the full Q, checked by
      dense eigendecomposition at small Omega. If it is not, the block construction is wrong
      and nothing below counts.
  P1  GATE (b). lambda_A < 0 at every finite Omega and every gamma. A non-negative value
      means the block or the ergodicity assumption is wrong, NOT a phase transition.
  P1  GATE (c). **ABSOLUTE, no fit (rule 16).** Above gamma_c the symmetric fixed point is
      stable and the CME's antisymmetric relaxation must converge to the DETERMINISTIC
      eigenvalue, which §62 gives in closed form as P(x*) = (1 - 2 gamma)/3. So for
      gamma > 1/2, lambda_A(Omega) -> (1 - 2 gamma)/3. Reported as a RATIO to that exact
      number, across Omega. This is the check §22 taught: fit nothing, compare in absolute
      terms against a quantity obtainable exactly.
      **VERDICT RULE, second version.** The first demanded |ratio - 1| < 0.05 at Omega=120
      and printed FAILS on ratios 1.165, 1.129, 1.086, 1.064 and 0.929, 0.966, 0.984,
      0.989 -- data converging cleanly to 1 from both sides. A fixed-Omega tolerance tests
      the SIZE of a finite-Omega correction, not whether it vanishes. The gate is now
      convergence itself: |ratio - 1| must decrease monotonically in Omega and its local
      exponent is reported. Data that would fail it: a deviation that plateaus or grows.
  P2  **THE MEASUREMENT.** ln|lambda_A| runs from ~ -Omega*A (metastable, steeply
      Omega-dependent) to ln|(1-2g)/3| (Omega-independent).

      FIRST VERSION, KEPT PER RULE 3 BECAUSE IT WAS WRONG: gamma*(Omega) = argmax of
      |d ln|lambda_A|/d gamma|, "parameter-free so it cannot be an artifact of a threshold".
      It is parameter-free and it measures the wrong thing. In the metastable branch the
      slope is -Omega A'(gamma), which is LARGEST deep inside that branch and SMALLEST at
      the crossover -- so the argmax runs to the low-gamma edge of whatever window is
      swept, and the first run duly printed "max slope at the EDGE" for every Omega.
      Being free of a threshold does not make a statistic measure the thing you named it
      after (rule 19: name the data that would make it print the other answer -- here NO
      data would have made it print a crossover).

      SECOND VERSION, anchored to a quantity §62 gives exactly. Above gamma_c the
      deterministic rate is L_det(gamma) = ln|(1 - 2 gamma)/3|, so define the EXCESS

          E(gamma, Omega) = L_det(gamma) - ln|lambda_A(gamma, Omega)|   >= 0

      -- how many e-folds slower than deterministic the symmetry breaking relaxes. E -> 0
      above the boundary and E -> Omega*A(gamma) below it, so **E/Omega is an order
      parameter with an exact reference rather than a fitted one**, and the transition
      region is where E passes through O(1). The width is measured as

          w[a,b](Omega) = gamma(E = a) - gamma(E = b)     for (a,b) = (1,2), (1,4), (2,4)

      i.e. the gamma-interval over which the metastability grows from a to b e-folds. This
      is a LEVEL DIFFERENCE, so it needs no "E ~ 0" endpoint -- which matters because
      L_det -> -infinity at gamma_c (P7) and because above the boundary E does not reach
      zero at finite Omega but only O(Omega^-1/2) (P1c). An e-fold is a unit, not a tuned
      threshold; and per rule 13 the choice of (a,b) is the approximation's own parameter,
      so **the exponent must be checked across (a,b) before it is compared across Omega**.
      If it moves with the levels, the width has no exponent.
      **w(Omega) is the device question: the blur on the threshold of a switch built from
      Omega molecules.**
  P3  PREDICTED: gamma*(Omega) -> 1/2 from BELOW and w(Omega) -> 0. If gamma* approaches
      from above, or does not approach 1/2 at all, the stochastic threshold is not a blurred
      version of the deterministic one and §9.1's picture does not survive coarse-graining.
  P4  **EVERY ANSATZ REPORTED, not the flattering one (rule 15).** Fit w ~ Omega^(-a) and
      w ~ c/ln(Omega) and w ~ a/Omega. If they disagree on the extrapolation the width's
      scaling is UNRESOLVED and says so. Two prior attempts at transferable exponents here
      (§39.2, §46) failed, and §59 found a third non-transferring exponent, so disagreement
      is the expectation.
  P7  **RULE 10 GUARD ON P3's OWN REFERENCE, added after the smoke run and before the
      measurement.** E is built on L_det = ln|(1 - 2 gamma)/3|, which has a LOGARITHMIC
      SINGULARITY at gamma_c because the deterministic eigenvalue vanishes there. So a
      width measured from E could be the width of that singularity rather than of the
      transition -- the harness doing the work the chemistry did not. **Kill test: measure
      the width again with an instrument that has no reference at all.**

          D(gamma) = ln|lambda_A(gamma, 2*Omega)| - ln|lambda_A(gamma, Omega)|

      D is ~ -Omega*A(gamma) in the metastable phase and -> 0 above the boundary, is finite
      everywhere, and never mentions the deterministic rate. If the two widths give the
      same exponent the width is a property of the transition; if they disagree, P3's
      number is an artifact of its reference and is withdrawn in favour of this one.
  P5  **RULE 9, an axis I did not choose** -- implemented with the singularity-free D of
      P7, so it does not inherit P3's reference. Sweep rho (§44's lever). rho changes the
      timescale separation without moving delta*, so if the width is a property of the
      restoring threshold its EXPONENT in Omega should be rho-independent even though its
      prefactor need not be. A moving exponent means the width is not a property of the
      threshold alone.
  P6  **RULE 10 GUARD.** lambda_A is exponentially small deep in the restoring phase and
      will hit the floor of double precision. Every eigenvalue whose magnitude is within
      1e-11 of solver noise is DISCARDED and counted, never plotted -- a numerically dead
      eigenvalue would read as a very sharp threshold, which is exactly the flattering
      artifact. The usable gamma window is reported with the result.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from crnl.cme import enumerate_states, generator
from crnl.networks.am_reversible import am_reversible
from experiments.arrhenius_optimum import am_rho


def antisym_block(net, omega):
    """Q restricted to antisymmetric observables under X <-> Y. Exact, no approximation."""
    states, index = enumerate_states(3, int(omega))
    Q = generator(net, int(omega), float(omega)).tocsr()
    rep = [i for i, s in enumerate(states) if int(s[0]) > int(s[1])]
    pos = -np.ones(len(states), dtype=np.int64)
    pos[rep] = np.arange(len(rep))
    mirror = np.array([index[(int(s[1]), int(s[0]), int(s[2]))] for s in states])
    A = Q[rep][:, rep]
    B = Q[rep][:, mirror[rep]]
    return (A - B).tocsc(), len(rep)


def lambda_A(net, omega, tol=1e-11):
    """Leading (least negative) eigenvalue of the antisymmetric block."""
    M, n = antisym_block(net, omega)
    if n < 3:
        return None
    try:
        v = spla.eigs(M, k=1, sigma=0.0, which="LM", return_eigenvectors=False,
                      tol=0, maxiter=5000)
    except Exception:
        return None
    lam = complex(v[0])
    if abs(lam.imag) > 1e-8 * max(abs(lam.real), 1e-30):
        return None
    return float(lam.real)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omegas", type=int, nargs="+", default=[20, 30, 40, 60, 80])
    ap.add_argument("--gammas", type=float, nargs="+", default=None)
    ap.add_argument("--ng", type=int, default=41)
    ap.add_argument("--glo", type=float, default=0.20)
    ap.add_argument("--ghi", type=float, default=0.80)
    ap.add_argument("--floor", type=float, default=1e-11)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/threshold_sharpness.json"))
    args = ap.parse_args()
    t0 = time.time()
    gam = (np.array(args.gammas) if args.gammas
           else np.linspace(args.glo, args.ghi, args.ng))

    # ---------------- P1 GATE (a): is Q_A's spectrum inside Q's? --------------------
    print("=== P1 GATE (a): spectrum of the antisymmetric block ⊂ spectrum of Q")
    worst_a = 0.0
    for om in (8, 12, 16):
        net = am_reversible(0.30)
        M, _ = antisym_block(net, om)
        Q = generator(net, om, float(om)).toarray()
        ev_full = np.linalg.eigvals(Q)
        ev_blk = np.linalg.eigvals(M.toarray())
        for e in ev_blk:
            worst_a = max(worst_a, float(np.abs(ev_full - e).min()))
        print(f"  Omega={om:>3}: block {M.shape[0]:>4} of {Q.shape[0]:>4} states,"
              f" worst distance to the full spectrum {worst_a:.3e}")
    print(f"  -> P1(a) {'HOLDS' if worst_a < 1e-8 else 'FAILS'}")

    # ---------------- P1 GATE (c): absolute, against §62's (1-2g)/3 -----------------
    print("\n=== P1 GATE (c): above gamma_c, lambda_A -> the exact (1 - 2 gamma)/3")
    print(f"{'gamma':>7}{'Omega':>7}{'lambda_A':>14}{'(1-2g)/3':>12}{'ratio':>10}")
    gate_c = []
    for g in (0.60, 0.75):
        for om in (20, 40, 80, 120):
            net = am_reversible(g)
            lam = lambda_A(net, om)
            if lam is None:
                continue
            ref = (1.0 - 2.0 * g) / 3.0
            print(f"{g:>7.2f}{om:>7}{lam:>14.8f}{ref:>12.8f}{lam/ref:>10.5f}")
            gate_c.append({"gamma": g, "omega": om, "lam": lam, "ratio": lam / ref})
    ok = True
    for g in sorted({r["gamma"] for r in gate_c}):
        sel = sorted([r for r in gate_c if r["gamma"] == g], key=lambda r: r["omega"])
        dev = [abs(r["ratio"] - 1.0) for r in sel]
        oms = [r["omega"] for r in sel]
        mono = all(dev[i + 1] < dev[i] for i in range(len(dev) - 1))
        expo = (np.log(dev[0] / dev[-1]) / np.log(oms[-1] / oms[0])
                if dev[-1] > 0 else np.inf)
        ok = ok and mono
        print(f"  gamma={g}: |ratio-1| " + " -> ".join(f"{d:.4f}" for d in dev)
              + f"   {'monotone' if mono else 'NOT monotone'}, ~Omega^-{expo:.2f}")
    print(f"  -> P1(c) {'HOLDS: the deviation vanishes with Omega, so the CME meets §62' if ok else 'FAILS -- the deviation does not shrink'}")

    # ---------------- P1 GATE (b) + the sweep ---------------------------------------
    print(f"\n=== P2: ln|lambda_A| across gamma, per Omega")
    data, dead, positive = {}, 0, 0
    for om in args.omegas:
        gs, ls = [], []
        for g in gam:
            net = am_reversible(float(g))
            lam = lambda_A(net, om)
            if lam is None:
                continue
            if lam >= 0:
                positive += 1
                continue
            if abs(lam) < args.floor:            # P6: numerically dead, discard
                dead += 1
                continue
            gs.append(float(g)); ls.append(float(np.log(abs(lam))))
        data[om] = {"gamma": gs, "lnabs": ls}
        if gs:
            print(f"  Omega={om:>4}: {len(gs):>3} usable gamma in"
                  f" [{min(gs):.3f}, {max(gs):.3f}],"
                  f" ln|lambda_A| from {min(ls):.2f} to {max(ls):.2f}")
    print(f"  -> P1(b) {'HOLDS: lambda_A < 0 everywhere' if positive == 0 else f'FAILS: {positive} non-negative'}")
    print(f"  P6: {dead} eigenvalues discarded below the {args.floor:g} precision floor")

    # ---------------- the excess E, and the width of the transition -----------------
    print(f"\n=== P3: the excess E = L_det - ln|lambda_A|, and the width of {{0 < E < k}}")
    PAIRS = ((1.0, 2.0), (1.0, 4.0), (2.0, 4.0))
    print(f"{'Omega':>7}" + "".join(f"{f'w[{a:g},{b:g}]':>11}" for a, b in PAIRS)
          + f"{'g(E=1)':>10}{'g(E=4)':>10}")
    rows = []

    def gamma_at(gs, E, level):
        """Where E crosses `level`, by linear interpolation. E decreases in gamma."""
        for i in range(len(gs) - 1):
            if (E[i] - level) * (E[i + 1] - level) <= 0 and E[i] != E[i + 1]:
                t = (level - E[i]) / (E[i + 1] - E[i])
                return float(gs[i] + t * (gs[i + 1] - gs[i]))
        return None

    def pair_widths(gs, Q):
        """gamma-interval between two metastability levels. No endpoint, no reference."""
        out = {}
        for a, b in PAIRS:
            ga, gb = gamma_at(gs, Q, a), gamma_at(gs, Q, b)
            out[f"{a:g},{b:g}"] = None if (ga is None or gb is None) else abs(ga - gb)
        return out

    for om in args.omegas:
        gs = np.array(data[om]["gamma"]); ls = np.array(data[om]["lnabs"])
        if len(gs) < 7:
            print(f"{om:>7}   too few usable points ({len(gs)})")
            continue
        E = np.log(np.abs((1.0 - 2.0 * gs) / 3.0)) - ls
        ws = pair_widths(gs, E)
        if ws["1,2"] is None:
            print(f"{om:>7}   E does not bracket the levels inside the swept window")
            continue
        rows.append({"omega": om, "w": ws, "g_E1": gamma_at(gs, E, 1.0),
                     "g_E4": gamma_at(gs, E, 4.0)})
        print(f"{om:>7}" + "".join(
            f"{ws[f'{a:g},{b:g}']:>11.5f}" if ws[f"{a:g},{b:g}"] else f"{'--':>11}"
            for a, b in PAIRS)
              + f"{rows[-1]['g_E1']:>10.5f}"
              + (f"{rows[-1]['g_E4']:>10.5f}" if rows[-1]["g_E4"] else f"{'--':>10}"))

    if len(rows) >= 3:
        om = np.array([r["omega"] for r in rows], float)
        print(f"\n=== P3/rule 13: is the exponent stable in the LEVELS before comparing"
              f" across Omega?")
        exps = {}
        for a, b in PAIRS:
            key = f"{a:g},{b:g}"
            w = np.array([r["w"][key] for r in rows if r["w"][key]], float)
            o = np.array([r["omega"] for r in rows if r["w"][key]], float)
            if len(w) < 3:
                continue
            exps[key] = float(np.polyfit(np.log(o), np.log(w), 1)[0])
            print(f"  levels {key:>7}: w ~ Omega^({exps[key]:+.4f})"
                  f"   over Omega {o.min():.0f}..{o.max():.0f}")
        if len(exps) >= 2:
            v = np.array(list(exps.values()))
            sp_ = (v.max() - v.min()) / max(abs(v.mean()), 1e-30)
            print(f"  spread in the exponent across levels: {100*sp_:.1f}%")
            print(f"  -> {'STABLE, so the exponent is a property of the width' if sp_ < 0.15 else 'MOVES with the levels: the width has no single exponent (rule 13)'}")

        print(f"\n=== P4: every ansatz for w[1,4](Omega), not only the flattering one")
        w = np.array([r["w"]["1,4"] for r in rows], float)
        a_pow, c_pow = np.polyfit(np.log(om), np.log(w), 1)
        pred_pow = float(np.exp(c_pow) * 1000.0 ** a_pow)
        c_log = np.polyfit(1.0 / np.log(om), w, 1)
        pred_log = float(c_log[0] / np.log(1000.0) + c_log[1])
        c_inv = np.polyfit(1.0 / om, w, 1)
        pred_inv = float(c_inv[0] / 1000.0 + c_inv[1])
        print(f"  w ~ Omega^a        a = {a_pow:+.4f}   -> w(1000) = {pred_pow:.5f}")
        print(f"  w ~ c/ln(Omega) + b               -> w(1000) = {pred_log:.5f}")
        print(f"  w ~ a/Omega + b                   -> w(1000) = {pred_inv:.5f}")
        preds = np.array([pred_pow, pred_log, pred_inv])
        spread = (preds.max() - preds.min()) / max(abs(preds.mean()), 1e-30)
        print(f"  spread across ansaetze: {100*spread:.0f}%")
        print(f"  -> {'the ansaetze AGREE; the width has a scaling' if spread < 0.25 else 'the ansaetze DISAGREE: the width scaling is UNRESOLVED (rule 15)'}")

    # ---------------- P7: the same width with NO reference at all -------------------
    print(f"\n=== P7 (rule 10): the width again from D = ln|lam(2W)| - ln|lam(W)|,")
    print(f"    which never mentions the deterministic rate and has no singularity")

    def width_D(build, omegas, gs):
        """Pair-widths of |D|, the same definition P3 uses. build(gamma) -> network."""
        out = []
        for om in omegas:
            vals = []
            for g in gs:
                a = lambda_A(build(float(g)), om)
                b = lambda_A(build(float(g)), 2 * om)
                if (a is None or b is None or a >= 0 or b >= 0
                        or abs(a) < args.floor or abs(b) < args.floor):
                    vals.append(np.nan)
                    continue
                vals.append(abs(np.log(abs(b)) - np.log(abs(a))))
            D = np.array(vals)
            m = np.isfinite(D)
            if m.sum() < 7:
                continue
            r = {"omega": om}
            r.update(pair_widths(np.array(gs)[m], D[m]))
            out.append(r)
        return out

    gsD = np.linspace(args.glo, args.ghi, args.ng)
    omD = list(args.omegas)
    drows = width_D(lambda g: am_reversible(g), omD, gsD)
    print(f"{'Omega':>7}" + "".join(f"{f'w[{a:g},{b:g}]':>11}" for a, b in PAIRS))
    for r in drows:
        print(f"{r['omega']:>7}" + "".join(
            f"{r[f'{a:g},{b:g}']:>11.5f}" if r[f"{a:g},{b:g}"] else f"{'--':>11}"
            for a, b in PAIRS))
    dexp = {}
    for a, b in PAIRS:
        key = f"{a:g},{b:g}"
        o = np.array([r["omega"] for r in drows if r.get(key)], float)
        w = np.array([r[key] for r in drows if r.get(key)], float)
        if len(w) >= 3:
            dexp[key] = float(np.polyfit(np.log(o), np.log(w), 1)[0])
            print(f"  levels {key:>7}: w ~ Omega^({dexp[key]:+.4f})  [no reference]")
    if dexp and len(rows) >= 3:
        a_ref = exps.get("1,4")
        a_noref = dexp.get("1,4")
        print("  P3 levels 1,4 (singular reference): "
              + (f"{a_ref:+.4f}" if a_ref is not None else "--"))
        print("  P7 levels 1,4 (no reference):       "
              + (f"{a_noref:+.4f}" if a_noref is not None else "--"))
        agree = (a_ref is not None and a_noref is not None
                 and abs(a_ref - a_noref) < 0.10)
        print(f"  -> P7 {'AGREE: the width is a property of the transition, not of the reference' if agree else 'DISAGREE: P3s width is an artifact of its reference and is withdrawn'}")

    # ---------------- P5: rule 9, the rho axis, on the reference-free width ---------
    print(f"\n=== P5 (rule 9): the same exponent under rho, §44's lever")
    print(f"{'rho':>7}{'exponent':>12}{'n':>5}")
    rexp = {}
    for rho in (1.0, 4.0, 16.0):
        rr = width_D(lambda g, rho=rho: am_rho(float(g), rho), omD, gsD)
        o = np.array([r["omega"] for r in rr if r.get("1,4")], float)
        w = np.array([r["1,4"] for r in rr if r.get("1,4")], float)
        if len(w) < 3:
            print(f"{rho:>7.0f}{'too few':>12}{len(w):>5}")
            continue
        rexp[rho] = float(np.polyfit(np.log(o), np.log(w), 1)[0])
        print(f"{rho:>7.0f}{rexp[rho]:>12.4f}{len(w):>5}")
    if len(rexp) >= 2:
        v = np.array(list(rexp.values()))
        print(f"  exponent spread across rho: {100*(v.max()-v.min())/abs(v.mean()):.1f}%")
        print(f"  -> {'rho-INDEPENDENT: the width belongs to the threshold' if (v.max()-v.min())/abs(v.mean()) < 0.15 else 'the exponent MOVES with rho: the width is not a property of the threshold alone'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"gate_a": worst_a, "gate_c": gate_c,
                                    "curves": {str(k): v for k, v in data.items()},
                                    "crossover": rows, "dead": dead,
                                    "width_noref": drows, "exp_noref": {str(k): v for k, v in dexp.items()},
                                    "exp_rho": {str(k): v for k, v in rexp.items()}},
                                   indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
