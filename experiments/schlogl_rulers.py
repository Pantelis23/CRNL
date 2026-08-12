"""T-THM-b: run this project's rulers on a restoring element with NO exchange symmetry

Every quantitative result in CRNL -- G ~ 2 k_B per molecule per e-fold (§38), the affinity
floor A > 3 ln 2 (§9.1), gamma_c = 1/2 (§9.1/§62), Q = 5.475 and the (time, Q) frontier
exponent 0.58 (§57-§59) -- is measured on ONE network, Approximate Majority, and on its
exchange-symmetric relatives. §65's theorem is conditioned on that symmetry outright, and
§42 measured residuals of 1.9e1 and 2.9e2 when it is broken. **Meanwhile the founding object
is asymmetric**: an inverter drives toward one rail, not toward whichever it started nearest.

So the answer to the founding question currently has n = 1, and this project has found three
times that exponents and coefficients do NOT transfer between axes (§39.2's 1/sep coefficient,
§46's scaling, §59's frontier exponent). **The prior must therefore be that G, Q and the
frontier exponent are AM facts, not restoration facts.**

THE SUBSTRATE. Schloegl's model, the textbook bistable CRN, with A and B chemostatted and
folded into the rate constants so X is the only dynamic species:

    2X <-> 3X     forward k1a = k1*[A], reverse k1r
     0 <-> X      forward k2b = k2*[B], reverse k2r

No exchange symmetry exists to break: there is one species. It is 1-D, so the birth-death
chain is EXACT -- splitting probability in closed form (§61's method), MFPT and its variance
and the mean entropy production to absorption all tridiagonal solves. **This is a better
instrument than AM's, so a failure to transfer cannot be blamed on numerics.**

PARAMETERISATION, chosen so the landscape is the knob. Placing the three fixed points at
x0 - m, x0, x0 + m makes the cubic (x-x0)^3 - m^2 (x-x0), so

    k1a/k1r = 3*x0,   k2r/k1r = 3*x0^2 - m^2,   k2b/k1r = x0*(x0^2 - m^2)

with x0 the operating point and m the half-width of the landscape. m -> 0 is where bistability
dies -- Schloegl's analogue of gamma -> gamma_c.

**A STRUCTURAL DIFFERENCE FROM AM, WORTH STATING BEFORE ANY NUMBER.** In AM, gamma sets the
landscape and the affinity independently enough that §9.1 could speak of a floor. Here the
affinity is DETERMINED by the landscape: A = ln[k1a*k2r/(k1r*k2b)] = ln[3(3x0^2 - m^2)/(x0^2 - m^2)].
There is no free drive knob at fixed landscape. That is a fact about the substrate, not a
limitation of the measurement.

PREDICTIONS. P3 was DERIVED analytically before running and then confirmed against the
engine's independent `cycle_affinity`, which knows nothing of the formula (it takes the null
space of the per-pair stoichiometry); that ordering is the point of rule 16.

  P1  GATE. The deterministic cubic has roots exactly x0 - m, x0, x0 + m, and the
      birth-death chain's drift (lambda - mu)/Omega agrees with the ODE field.
      **SECOND VERSION: the first demanded agreement to 1e-12, which is FALSE and could
      never have passed** -- the chain uses falling factorials n(n-1) and n(n-1)(n-2) where
      mass action uses x^2 and x^3, so the two differ by O(1/Omega) BY CONSTRUCTION. That is
      the discreteness this whole project measures, not an error. The gate is now that
      Omega * |drift - field| is bounded and Omega-independent, i.e. the disagreement is
      exactly first order. Caught by the pre-run test (§66's convention) before any cell ran;
      a gate that demands something false is the same defect class as a verdict that cannot
      fail.
  P2  GATE. `cycle_affinity` on the built network equals ln[3(3x0^2 - m^2)/(x0^2 - m^2)] to
      1e-12, at several (x0, m). Checked at six cells before this file was written: it does,
      to 10 decimals.
  P3  **THE PRE-REGISTERED ABSOLUTE TEST. As m -> 0 the affinity tends to ln 9 = 2 ln 3 =
      2.1972245773, INDEPENDENT OF x0.** So Schloegl has an affinity floor at the death of
      bistability, exactly as AM does, and the two numbers are
          **AM: 3 ln 2 = 2.0794415417        Schloegl: 2 ln 3 = 2.1972245773**
      differing by 5.7%. **Predicted: they are DIFFERENT, so the floor is not universal** --
      but both are ln(small integer) and within 6% of each other with no shared chemistry,
      which is the striking part either way. §9.1 read its own floor as "3 reactions x ln 2";
      this substrate has 2 reversible pairs and gives 2 x ln 3, so that reading is testable
      here rather than decorative.
  P4  **G, THE COST PER E-FOLD.** §38 measured ~2 k_B per molecule per e-fold of
      amplification for AM. Measured here as (mean entropy production to absorption) /
      (Omega * e-folds of spread growth). **PREDICTED: NOT 2**, on the strength of three
      transfer failures. If it lands near 2 anyway, that is a substrate-independent price of
      restoration and by far the bigger result -- which is why the prediction is written
      against it.
  P5  **Q, §40's first-passage TUR ratio**, two-sided as §57-§58 used, floor 1. AM's best is
      5.475 and its frontier runs to 1.115 at t = 3747. **PREDICTED: Schloegl's Q differs
      from AM's by more than the 9% that separated 696 AM-family networks in §57.**
  P6  **RULE 10 GUARD.** The chain is unbounded above and must be capped. A cap acts as a
      REFLECTING WALL, which would push probability back toward the upper rail and
      manufacture restoration the chemistry did not do. Every reported quantity is checked
      for cap-independence and the cap is raised until it moves nothing at 1e-9.
  SECOND PASS, after the first run and before P4/P5 are read as a result. The raw numbers
  came out G = 8.3-27.1 against AM's ~2, and Q = 388-7101 against AM's 5.475. **Those are not
  a like-for-like comparison and reporting them as "does not transfer" would have been a
  rule-11 error** -- the control must share a clock with its arm. AM is CLOSED and
  conservative: all of its entropy production belongs to the decision. Schloegl is
  CHEMOSTATTED: the 0 <-> X channel runs at O(Omega) forever, so the system dissipates
  steadily while sitting still, and Sigma-to-absorption is dominated by a HOUSEKEEPING term
  proportional to the time taken rather than to the decision. Sigma = 14293 at Omega = 800
  against a cycle affinity of 2.26 is that term, not restoration.

  FIRST REPAIR ATTEMPT, AND IT FAILED -- kept per rule 3. Subtracting a housekeeping term
  Sigma_ex = Sigma - sigma_local(x0) * <T> gave a NEGATIVE excess in all 15 cells
  (-85 to -5878), because sigma_local at the operating point exceeds its average along the
  path, so the subtraction overshoots. A negative excess is not a small number, it is an
  invalid definition, and it is reported rather than tuned away.

  WHY NO SUBTRACTION WORKS. A 1-D birth-death chain has ZERO stationary probability current,
  so in the adiabatic/non-adiabatic decomposition essentially all of sigma_local is
  housekeeping and the non-adiabatic part is the system term ln[pi(n0)/pi(n_f)], which is
  NEGATIVE for a trajectory running from the unstable point down to a rail. **So §38's
  "entropy produced per e-fold of amplification" has no counterpart here at all** -- it is a
  quantity defined by AM's closed conservative bookkeeping, not by restoration.

  THE COMPARABLE QUANTITY, which needs no split. Report **Sigma / (Omega * A)**: dissipation
  per molecule measured in units of the network's own cycle affinity, i.e. **how many driven
  cycles the element turns per molecule to make one decision**. It is dimensionless, defined
  for closed and chemostatted networks alike, and can be computed for AM from its published
  numbers. That is the cross-substrate price of restoration, if there is one.
  P7  **VERDICT RULES, unit-tested on engineered data before this runs** (§66's convention,
      tests/test_schlogl_rulers.py). "Transfers" and "does not transfer" must both be
      reachable, and the comparison must be against AM's PUBLISHED numbers with their own
      stated uncertainties, not against a fresh AM run whose noise differs.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
from scipy.linalg import solve_banded

from crnl.networks.am_reversible import cycle_affinity, reverse_pairing
from crnl.reactions import Reaction, ReactionNetwork

AM_FLOOR = 3.0 * np.log(2.0)          # §9.1, published
AM_G = 2.0                            # §38, published, ~2 k_B/molecule/e-fold
AM_Q = 5.475                          # §57, published (the §40 value was 5.39)


def schlogl(x0, m, k1r=1.0):
    """Fixed points at x0 - m, x0, x0 + m. Chemostats folded into the constants."""
    k1a = 3.0 * x0 * k1r
    k2r = (3.0 * x0 ** 2 - m ** 2) * k1r
    k2b = x0 * (x0 ** 2 - m ** 2) * k1r
    rx = [Reaction({"X": 2}, {"X": 3}, k1a, name="f1:2X->3X"),
          Reaction({"X": 3}, {"X": 2}, k1r, name="r1:3X->2X"),
          Reaction({}, {"X": 1}, k2b, name="f2:->X"),
          Reaction({"X": 1}, {}, k2r, name="r2:X->")]
    return ReactionNetwork(species=["X"], reactions=rx, name="schlogl")


def _consts(x0, m, k1r=1.0):
    return (3.0 * x0 * k1r, k1r, x0 * (x0 ** 2 - m ** 2) * k1r,
            (3.0 * x0 ** 2 - m ** 2) * k1r)


def channels(n, omega, x0, m):
    """The four propensities at count n. Exact, no continuum approximation."""
    k1a, k1r, k2b, k2r = _consts(x0, m)
    n = np.asarray(n, dtype=float)
    a1 = k1a * n * (n - 1.0) / omega                     # 2X -> 3X
    a2 = k1r * n * (n - 1.0) * (n - 2.0) / omega ** 2    # 3X -> 2X
    a3 = np.full_like(n, k2b * omega)                    # 0 -> X
    a4 = k2r * n                                         # X -> 0
    return np.maximum(a1, 0.0), np.maximum(a2, 0.0), a3, np.maximum(a4, 0.0)


def birth_death(omega, x0, m, cap):
    ns = np.arange(0, cap + 1, dtype=float)
    a1, a2, a3, a4 = channels(ns, omega, x0, m)
    return a1 + a3, a2 + a4, ns          # lambda(n), mu(n)


def sigma_local(omega, x0, m, cap):
    """Exact local entropy-production rate: each channel against its own reverse."""
    ns = np.arange(0, cap + 1, dtype=float)
    a1, a2, a3, a4 = channels(ns, omega, x0, m)
    up1, up2 = channels(ns + 1, omega, x0, m)[1], channels(ns + 1, omega, x0, m)[3]
    dn1, dn2 = channels(ns - 1, omega, x0, m)[0], channels(ns - 1, omega, x0, m)[2]
    s = np.zeros_like(ns)
    for fwd, rev in ((a1, up1), (a3, up2), (a2, dn1), (a4, dn2)):
        ok = (fwd > 0) & (rev > 0)
        s[ok] += fwd[ok] * np.log(fwd[ok] / rev[ok])
    return s


def _tridiag(lo, hi, lam, mu, source):
    """Solve L f = -source on (lo, hi) with f(lo) = f(hi) = 0. L is the generator."""
    idx = np.arange(lo + 1, hi)
    k = len(idx)
    ab = np.zeros((3, k))
    ab[0, 1:] = lam[idx[:-1]]
    ab[1, :] = -(lam[idx] + mu[idx])
    ab[2, :-1] = mu[idx[1:]]
    return idx, solve_banded((1, 1), ab, -source[idx])


def cell(omega, x0, m, eps, theta, cap_mult=3.0):
    """All rulers at once on the exact chain. Start eps below the unstable point."""
    cap = int(np.ceil(cap_mult * (x0 + m) * omega))
    lam, mu, ns = birth_death(omega, x0, m, cap)
    lo = int(round((x0 - theta * m) * omega))
    hi = int(round((x0 + theta * m) * omega))
    n0 = int(round((x0 - eps * m) * omega))
    if not (lo < n0 < hi) or hi >= cap:
        return None

    # exact splitting probability in logs (§61's construction, no cancellation)
    lp, acc = np.zeros(hi - lo + 1), 0.0
    for k in range(lo + 1, hi + 1):
        if lam[k] <= 0 or mu[k] <= 0:
            return None
        acc += np.log(mu[k]) - np.log(lam[k])
        lp[k - lo] = acc

    def lse(v):
        M = v.max()
        return M + np.log(np.exp(v - M).sum())

    i = n0 - lo
    p_down = float(np.exp(lse(lp[i:hi - lo]) - lse(lp[0:hi - lo])))

    one = np.ones(cap + 1)
    idx, T1 = _tridiag(lo, hi, lam, mu, one)             # mean first passage
    j = int(np.where(idx == n0)[0][0])
    src2 = np.zeros(cap + 1)
    src2[idx] = 2.0 * T1
    _, T2 = _tridiag(lo, hi, lam, mu, src2)              # second moment
    sig = sigma_local(omega, x0, m, cap)
    _, SS = _tridiag(lo, hi, lam, mu, sig)               # mean EP to absorption

    mT, m2, Sig = float(T1[j]), float(T2[j]), float(SS[j])
    var = m2 - mT ** 2
    efolds = np.log(theta / eps)
    # rule 11: AM is closed, Schloegl is chemostatted and burns entropy sitting still.
    # The comparable quantity is the EXCESS over that housekeeping rate.
    hk = float(sig[int(round(x0 * omega))])
    Sig_ex = Sig - hk * mT          # refuted definition, kept and reported (rule 3)
    Q = (var / mT ** 2) * Sig / 2.0 if mT > 0 else np.nan
    return {"omega": omega, "x0": x0, "m": m, "eps": eps, "theta": theta, "cap": cap,
            "lo": lo, "hi": hi, "n0": n0, "p_down": p_down, "mean_T": mT,
            "var_T": var, "Sigma": Sig, "sigma_hk": hk, "Sigma_ex": float(Sig_ex),
            "Q": float(Q), "efolds": float(efolds),
            "G": float(Sig / (omega * efolds)),
            "G_ex": float(Sig_ex / (omega * efolds))}


def verdict_transfer(measured, published, tol_frac, label):
    """Reachable both ways; compares against a PUBLISHED number, not a fresh run."""
    rel = abs(measured - published) / abs(published)
    if rel <= tol_frac:
        return "transfers", (f"{label}: {measured:.4f} against AM's published "
                             f"{published:.4f} -- {100*rel:.1f}%, inside {100*tol_frac:.0f}%. "
                             f"TRANSFERS, which the prior said it would not.")
    return "does-not", (f"{label}: {measured:.4f} against AM's published {published:.4f} "
                        f"-- {100*rel:.1f}%, outside {100*tol_frac:.0f}%. Does NOT transfer, "
                        f"as three prior failures predicted.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--x0", type=float, default=1.0)
    ap.add_argument("--ms", type=float, nargs="+", default=[0.3, 0.4, 0.5, 0.6, 0.7])
    ap.add_argument("--omegas", type=int, nargs="+", default=[200, 400, 800])
    ap.add_argument("--eps", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/schlogl_rulers.json"))
    args = ap.parse_args()
    t0 = time.time()

    print("=== P1 GATE: roots exact, and the chain's drift equals the ODE field")
    worst_r = worst_d = 0.0
    for x0, m in ((args.x0, 0.5), (2.0, 0.4), (0.7, 0.3)):
        net = schlogl(x0, m)
        k1a, k1r, k2b, k2r = _consts(x0, m)
        roots = np.sort(np.roots([-k1r, k1a, -k2r, k2b]).real)
        worst_r = max(worst_r, float(np.abs(roots - np.array([x0 - m, x0, x0 + m])).max()))
        scaled = []
        for om in (200, 400, 800):
            lam, mu, ns = birth_death(om, x0, m, cap=int(3 * (x0 + m) * om))
            n = int(1.2 * x0 * om)
            drift = (lam[n] - mu[n]) / om
            field = float(net.rhs(np.array([n / om]))[0])
            scaled.append(om * abs(drift - field))
        worst_d = max(worst_d, float(np.ptp(scaled) / max(np.mean(scaled), 1e-30)))
    print(f"  worst root error {worst_r:.3e}")
    print(f"  Omega*|drift - field| spread across Omega = {worst_d:.3e}"
          f"  (0 would mean exactly first order)")
    print(f"  -> P1 {'HOLDS: roots exact, discreteness exactly O(1/Omega)' if worst_r < 1e-10 and worst_d < 0.05 else 'FAILS'}")

    print("\n=== P2/P3: the affinity, and its floor as the landscape dies")
    print(f"{'x0':>6}{'m':>8}{'A (engine)':>14}{'A (formula)':>14}{'diff':>11}")
    worst_a = 0.0
    for x0, m in ((args.x0, 0.5), (args.x0, 0.9), (2.5, 0.5), (0.7, 0.3),
                  (args.x0, 1e-6), (2.5, 1e-6), (0.4, 1e-6)):
        net = schlogl(x0, m)
        A = cycle_affinity(net, reverse_pairing(net))
        pred = float(np.log(3 * (3 * x0 ** 2 - m ** 2) / (x0 ** 2 - m ** 2)))
        worst_a = max(worst_a, abs(A - pred))
        print(f"{x0:>6.2f}{m:>8.6g}{A:>14.10f}{pred:>14.10f}{abs(A-pred):>11.2e}")
    print(f"  -> P2 {'HOLDS' if worst_a < 1e-12 else 'FAILS'}")
    floor = float(np.log(9.0))
    print(f"\n  **Schloegl floor  = 2 ln 3 = {floor:.10f}**")
    print(f"  **AM floor (§9.1) = 3 ln 2 = {AM_FLOOR:.10f}**")
    print(f"  difference {100*abs(floor-AM_FLOOR)/AM_FLOOR:.2f}%   ratio {floor/AM_FLOOR:.6f}")
    print(f"  -> P3: the floor is NOT universal, but both are ln(small integer)"
          f" and agree to {100*abs(floor-AM_FLOOR)/AM_FLOOR:.1f}%")

    print("\n=== P6 GUARD: is anything moved by the cap (a reflecting wall)?")
    base = cell(args.omegas[0], args.x0, args.ms[len(args.ms) // 2], args.eps, args.theta,
                cap_mult=3.0)
    worst_cap = 0.0
    for cm in (4.0, 6.0):
        alt = cell(args.omegas[0], args.x0, args.ms[len(args.ms) // 2], args.eps,
                   args.theta, cap_mult=cm)
        for k in ("p_down", "mean_T", "Sigma", "Q"):
            worst_cap = max(worst_cap, abs(alt[k] - base[k]) / max(abs(base[k]), 1e-30))
    print(f"  worst relative change from cap 3x -> 6x: {worst_cap:.3e}")
    print(f"  -> P6 {'HOLDS: the cap is not doing work' if worst_cap < 1e-9 else 'FAILS: the cap is acting as a wall'}")

    print(f"\n=== P4/P5: G and Q on the exact chain")
    print("  raw Sigma is NOT comparable to AM's (rule 11: AM closed, Schloegl chemostatted).")
    print("  Sig_ex is the REFUTED housekeeping subtraction, shown because it came out"
          " negative (rule 3).")
    print("  The comparable column is cyc/mol = Sigma/(Omega*A): driven cycles per molecule"
          " per decision.")
    print(f"{'m':>6}{'Om':>6}{'A':>9}{'mean_T':>10}{'Sigma':>11}{'Sig_ex':>11}"
          f"{'cyc/mol':>10}{'Q':>11}")
    rows = []
    for m in args.ms:
        A = float(np.log(3 * (3 * args.x0 ** 2 - m ** 2) / (args.x0 ** 2 - m ** 2)))
        for om in args.omegas:
            r = cell(om, args.x0, m, args.eps, args.theta)
            if r is None:
                continue
            r["A"] = A
            rows.append(r)
            r["cyc_per_mol"] = r["Sigma"] / (om * A)
            print(f"{m:>6.2f}{om:>6}{A:>9.4f}{r['mean_T']:>10.4f}{r['Sigma']:>11.2f}"
                  f"{r['Sigma_ex']:>11.1f}{r['cyc_per_mol']:>10.4f}{r['Q']:>11.2f}")

    if rows:
        big = [r for r in rows if r["omega"] == args.omegas[-1]]
        exs = np.array([r["Sigma_ex"] for r in big])
        print(f"\n  P4 FIRST ATTEMPT REFUTED: the housekeeping-subtracted excess is NEGATIVE "
              f"in {int((exs < 0).sum())}/{len(exs)} cells ({exs.min():.0f}..{exs.max():.0f})."
              f" That definition is invalid, not small. §38's G has no counterpart here.")
        cyc = np.array([r["cyc_per_mol"] for r in big])
        # AM's comparable number, from PUBLISHED values only (rule: no fresh-run baseline)
        am_cyc = AM_G * np.log(0.80 / 0.35) / (-3.0 * np.log(0.20))
        print(f"\n  P4 COMPARABLE: driven cycles per molecule per decision")
        print(f"    Schloegl  {cyc.min():.3f}..{cyc.max():.3f} at Omega={args.omegas[-1]}")
        print(f"    AM        {am_cyc:.3f}   (from published G~{AM_G} k_B/molecule/e-fold,"
              f" {np.log(0.80/0.35):.3f} e-folds, A = -3 ln 0.2 = {-3*np.log(0.2):.3f})")
        cg, mg = verdict_transfer(float(np.median(cyc)), float(am_cyc), 0.25,
                                  "P4  cycles/molecule")
        print(f"  -> {mg}")
        # ...but is it even intensive? If it moves with Omega it inherits the same defect.
        mmid = args.ms[len(args.ms) // 2]
        ser = [(r["omega"], r["Sigma"] / (r["omega"] * r["A"]))
               for r in rows if r["m"] == mmid]
        print(f"  CAVEAT, and it disqualifies this measure too: cycles/molecule at m={mmid}"
              f" runs " + ", ".join(f"{v:.3f}@W={o}" for o, v in ser))
        print(f"  It GROWS with Omega, because Sigma is dominated by housekeeping ~ Omega*<T>"
              f" and <T> itself grows. So it is not intensive and not a substrate-independent"
              f" price either.")
        print(f"  -> **NO dissipation-based cost measure ported.** The chemostatted element's"
              f" dissipation is dominated by a term proportional to the decision TIME, which"
              f" a closed element has no analogue for.")
        Qs = np.array([r["Q"] for r in big])
        print(f"\n  P5  Q (raw Sigma, so an UPPER bound and not comparable to AM's):"
              f" {Qs.min():.1f}..{Qs.max():.1f} against AM's published {AM_Q}")
        print(f"  §40's floor Q >= 1: {int((Qs < 1.0).sum())} cells below it")
        print(f"  -> P5 NOT DECIDED: Q inherits the same non-comparability as Sigma;"
              f" reported, not compared.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"rows": rows, "floor_schlogl": floor,
                                    "floor_am": AM_FLOOR, "gate_affinity": worst_a,
                                    "gate_cap": worst_cap}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
