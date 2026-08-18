"""T15-n.1: §84's closed form sits BELOW the exact action. Two suspects die; the third is testable

§84 derived the escape action in closed form with no fitted parameter and validated it in
absolute terms -- the ratio to the exact first-passage action runs 0.8630 -> 0.9776 over
gamma = 0.30..0.46. **It approaches 1 from BELOW, and T15-n.1 asks why.** Three things could
produce a systematic deficit, and they are separable:

  S1  TRUNCATION of the slow manifold. §84 expanded b(u) = b0 + b2 u^2 and dropped u^4. The fix
      is to solve ds/dt = 0 for b exactly at each u, with no expansion.
  S2  THE DIFFUSION APPROXIMATION. §84 used the Fokker-Planck quasipotential 2 int F/D du. The
      exact action for a JUMP process is int ln(lam/mu) du, and these differ at third order:
      with r = lam/mu, ln r - 2(r-1)/(r+1) ~ (r-1)^3/12 > 0, so **the diffusion form must
      underestimate** -- which is the observed sign.
  S3  THE ADIABATIC ELIMINATION ITSELF -- the fast variable's fluctuations, which §84's P2 named
      as the one approximation §83 says to distrust.

**PREDICTIONS AS THEY WERE MADE, and two of them were wrong.** S2 was the motivating suspicion
for this section: it has the right sign and an elementary mechanism. Both S1 and S2 were
predicted to close most of the gap. They were scouted first, and **both are refuted**; the
predictions stay with what killed them (rule 3), and the section became the elimination of
candidates rather than the confirmation of one.

**Note on the sign, which rules out the easy version of S3.** The measured action is LARGER than
the reduced one. Fast-variable fluctuations feeding the lead would add to its effective
diffusion, and more noise means a SMALLER barrier. So whatever S3 is, it is not simply "the lead
is noisier than the reduction thinks", and a correction of that form would move the wrong way.

PREDICTIONS, written before running.

  P1  GATE. The reduced chain must reproduce the exact dynamics it claims to reduce: the slow
      manifold must satisfy b(u*) = b* = gamma/(1+gamma) to machine precision (it is the same
      fixed point), and the reduced drift lam_u - mu_u must equal the exactly-factored lead drift
      (1+gamma) u (b - b*) identically in u, not just at the fixed points.
  P2  **S1, and it was predicted to help.** Using the exact slow manifold instead of §84's
      quadratic truncation. **Predicted: the ratio moves closer to 1.**
  P3  **S2, the motivating suspect, with the right sign.** int ln(lam_u/mu_u) du against
      2 int F/D du on the same manifold. **Predicted: the WKB form closes most of the deficit.**
  P4  **S3, and the kill test is cheap and decisive.** The pair X+Y -> 2B and 2B -> X+Y has
      Delta u = 0: it moves the fast variable and does NOT touch the lead's rates. Multiplying
      both by M therefore sharpens the timescale separation, leaves the cycle affinity exactly
      unchanged (the ratio of the pair is untouched), and makes the adiabatic elimination
      asymptotically exact. **Predicted: the residual falls monotonically toward 0 as M grows.**
      If it does not, S3 is refuted too and the reduced chain is built wrongly rather than
      approximately.
  P5  **RULE 15, AND IT IS AWKWARD FOR §84.** Report that the TRUNCATED closed form agrees better
      than the EXACT reduction. **Predicted: that is a cancellation of errors, not an
      improvement** -- and the test is whether A_nf/A_exact-manifold is non-monotone in gamma, a
      systematic improvement being necessarily one-signed. §84's nu = 2 is untouched either way,
      since it is the leading behaviour and both reductions share it; what would be luck is the
      numerical closeness.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from crnl.cme import enumerate_states, generator
from crnl.networks.am_reversible import delta_star
from crnl.reactions import Reaction, ReactionNetwork
from experiments.nu_is_two import A_nf

# §84's measured actions, quoted as stored numbers (rule 16) -- not recomputed.
MEASURED = {0.30: 0.125373, 0.35: 0.071577, 0.38: 0.046444, 0.40: 0.032622,
            0.42: 0.021153, 0.44: 0.012080, 0.45: 0.008458, 0.46: 0.005460}


def slow_manifold(u, g, M=1.0):
    """Solve ds/dt = 0 for b at fixed lead u -- exactly, with no expansion in u.

    M scales the u-neutral pair (X+Y <-> 2B), which is the only part of the field that moves
    the fast variable without touching the lead's own rates.
    """
    def gg(b):
        s = 1.0 - b
        return (M * (-(s * s - u * u) / 2.0 + 2 * g * b * b)
                + b * s - g * (s * s + u * u) / 2.0)
    return float(brentq(gg, 1e-12, 1.0 - 1e-12))


def rates_u(u, g, M=1.0):
    """Birth (Delta u = +1) and death (Delta u = -1) rates of the lead on the slow manifold.

    Only f2/r2/f3/r3 change the lead, so these do NOT depend on M except through b(u).
    """
    b = slow_manifold(u, g, M)
    s = 1.0 - b
    x, y = (s + u) / 2.0, (s - u) / 2.0
    return b * x + g * y * y, g * x * x + b * y


def u_star(g, M=1.0):
    """The attractor: where lam_u = mu_u, i.e. where the slow manifold crosses b = b*."""
    if M == 1.0:
        return float(delta_star(g))
    bs = g / (1.0 + g)

    def h(u):
        return slow_manifold(u, g, M) - bs
    lo, hi = 1e-6, 0.999
    if h(lo) * h(hi) > 0:
        return None
    return float(brentq(h, lo, hi))


def A_fp(g, M=1.0):
    """Fokker-Planck quasipotential 2 int F/D du -- the form §84's closed formula approximates."""
    us = u_star(g, M)
    if us is None:
        return None

    def integrand(u):
        lam, mu = rates_u(u, g, M)
        return 2.0 * (lam - mu) / (lam + mu)
    return float(quad(integrand, 0.0, us, limit=200)[0])


def A_wkb(g, M=1.0):
    """The exact jump-process action on the same reduced chain: int ln(lam_u/mu_u) du."""
    us = u_star(g, M)
    if us is None:
        return None

    def integrand(u):
        lam, mu = rates_u(u, g, M)
        return float(np.log(lam / mu))
    return float(quad(integrand, 0.0, us, limit=200)[0])


def am_fast(g, M, k=1.0):
    """Reversible AM with the u-neutral pair X+Y <-> 2B scaled by M.

    Delta u = 0 for both, so the lead's own rates are untouched and the cycle affinity is
    unchanged (M cancels in the ratio of the pair).
    """
    return ReactionNetwork(
        species=["X", "Y", "B"],
        reactions=[
            Reaction({"X": 1, "Y": 1}, {"B": 2}, M * k, name="f1"),
            Reaction({"B": 1, "X": 1}, {"X": 2}, k, name="f2"),
            Reaction({"B": 1, "Y": 1}, {"Y": 2}, k, name="f3"),
            Reaction({"B": 2}, {"X": 1, "Y": 1}, M * g * k, name="r1"),
            Reaction({"X": 2}, {"B": 1, "X": 1}, g * k, name="r2"),
            Reaction({"Y": 2}, {"B": 1, "Y": 1}, g * k, name="r3"),
        ], name=f"am-fast-g{g}-M{M}")


def ln_mfpt_fast(g, M, omega):
    import scipy.sparse.linalg as spla
    net = am_fast(g, M)
    states, index = enumerate_states(3, int(omega))
    d = states[:, 0].astype(np.int64) - states[:, 1].astype(np.int64)
    trans = np.where(d > 0)[0]
    Qtt = generator(net, int(omega), float(omega))[trans][:, trans].tocsc()
    T = spla.spsolve(Qtt, -np.ones(len(trans)))
    us = u_star(g, M)
    if us is None:
        return None
    b = slow_manifold(us, g, M)
    s = 1.0 - b
    nx, ny = int(round((s + us) / 2 * omega)), int(round((s - us) / 2 * omega))
    if nx <= ny or nx + ny > omega:
        return None
    t = float(T[int(np.where(trans == index[(nx, ny, int(omega) - nx - ny)])[0][0])])
    if not np.isfinite(t) or t <= 0 or float(T.min()) <= 0:
        return None
    return float(np.log(t))


def A_measured_fast(g, M, band=(6.0, 20.0), rel=0.02):
    """Same instrument and the same Omega band as §84, re-targeted with the reduced A."""
    est = A_fp(g, M) or A_nf(g)
    # dedupe: at large M the action grows, the targets round onto the same Omega, and a
    # repeated Omega puts a zero in the local-slope denominator.
    oms = sorted({max(50, int(round(t / est / 50.0) * 50)) for t in np.linspace(*band, 4)})
    if len(oms) < 3:
        oms = sorted({max(50, int(round(t / est / 25.0) * 25))
                      for t in np.linspace(*band, 5)})
    pts = []
    for om in oms:
        v = ln_mfpt_fast(g, M, om)
        if v is not None and 0.0 < v < 28.0:
            pts.append((om, v))
    loc = [(pts[i + 1][1] - pts[i][1]) / (pts[i + 1][0] - pts[i][0])
           for i in range(len(pts) - 1)]
    if len(loc) < 2 or abs(loc[-1] - loc[-2]) / abs(loc[-1]) > rel:
        return None, loc, oms
    return loc[-1], loc, oms


def resid_exponent(gs, ratios):
    d = np.log(0.5 - np.asarray(gs, dtype=float))
    return float(np.polyfit(d, np.log(1.0 - np.asarray(ratios)), 1)[0])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/where_the_deficit_lives.json"))
    args = ap.parse_args()
    out = {}
    gs = sorted(MEASURED)

    print("=== P1 GATE: does the reduced chain reproduce the dynamics it claims to reduce?")
    worst_b, worst_f = 0.0, 0.0
    for g in (0.30, 0.40, 0.46):
        bs = g / (1.0 + g)
        worst_b = max(worst_b, abs(slow_manifold(float(delta_star(g)), g) - bs))
        for u in np.linspace(0.01, float(delta_star(g)) * 0.99, 12):
            lam, mu = rates_u(u, g)
            exact = (1.0 + g) * u * (slow_manifold(u, g) - bs)
            worst_f = max(worst_f, abs((lam - mu) - exact))
    print(f"  worst |b(u*) - b*| = {worst_b:.2e}   (the attractor is the same fixed point)")
    print(f"  worst |(lam_u - mu_u) - (1+g)u(b-b*)| over 36 points = {worst_f:.2e}")
    ok1 = worst_b < 1e-12 and worst_f < 1e-12
    print(f"  -> P1 {'HOLDS: the reduction is faithful; any deficit is the reduction, not a bug' if ok1 else 'FAILS'}")
    assert ok1

    print("\n=== P2/P3/P5: S1 (truncation) and S2 (diffusion approximation)")
    print(f"{'gamma':>7}{'measured':>11}{'§84 trunc':>11}{'ratio':>8}"
          f"{'exact mfld FP':>15}{'ratio':>8}{'exact mfld WKB':>16}{'ratio':>8}")
    rows = []
    for g in gs:
        m = MEASURED[g]
        a_t, a_f, a_w = A_nf(g), A_fp(g), A_wkb(g)
        rows.append({"gamma": g, "meas": m, "trunc": a_t, "fp": a_f, "wkb": a_w,
                     "r_trunc": a_t / m, "r_fp": a_f / m, "r_wkb": a_w / m})
        print(f"{g:>7}{m:>11.6f}{a_t:>11.6f}{a_t/m:>8.4f}{a_f:>15.6f}{a_f/m:>8.4f}"
              f"{a_w:>16.6f}{a_w/m:>8.4f}")
    out["reductions"] = rows
    e_t = resid_exponent(gs, [r["r_trunc"] for r in rows])
    e_f = resid_exponent(gs, [r["r_fp"] for r in rows])
    e_w = resid_exponent(gs, [r["r_wkb"] for r in rows])
    print(f"  residual ~ (gc-g)^p :  truncated {e_t:.3f},  exact manifold FP {e_f:.3f},"
          f"  exact manifold WKB {e_w:.3f}")
    s1_helps = e_f > e_t
    print(f"  -> P2 (S1) {'HOLDS: the exact manifold helps' if s1_helps else 'REFUTED: the EXACT manifold is WORSE than §84s truncation -- truncation is not the deficit'}")
    gap_fp_wkb = max(abs(r["wkb"] / r["fp"] - 1) for r in rows)
    print(f"  worst |WKB/FP - 1| on the same manifold = {100*gap_fp_wkb:.2f}%")
    s2_helps = gap_fp_wkb > 0.20
    print(f"  -> P3 (S2) {'HOLDS: the jump-process form closes the gap' if s2_helps else 'REFUTED: WKB and FP agree to well under the deficit. lam_u/mu_u stays near 1 across the barrier, so the cubic term is negligible -- the diffusion approximation is NOT the deficit'}")
    ratios = [r["trunc"] / r["fp"] for r in rows]
    nonmono = not (all(np.diff(ratios) > 0) or all(np.diff(ratios) < 0))
    print(f"  A_trunc/A_exact-manifold: " + ", ".join(f"{v:.4f}" for v in ratios))
    print(f"  -> P5 {'as predicted: NON-MONOTONE, so §84s closer agreement is a CANCELLATION of errors and not an improvement. Its nu = 2 is untouched; its numerical closeness is luck' if nonmono else 'monotone -- the truncation is systematically better, which needs explaining'}")

    print("\n=== P4: S3, the kill test. Speed up the fast variable and the elimination must fix")
    print("    (X+Y <-> 2B has Delta u = 0: it moves b, not the lead, and M cancels in the")
    print("     affinity, so this sharpens the timescale separation and changes nothing else)")
    print(f"{'gamma':>7}{'M':>5}{'u*':>9}{'A reduced':>12}{'A measured':>12}{'ratio':>8}"
          f"{'1-ratio':>10}")
    p4, excluded = [], 0
    for g in (0.40, 0.44):
        for M in (1, 2, 4, 8, 16):
            ared = A_fp(g, float(M))
            am, loc, oms = A_measured_fast(g, float(M))
            us = u_star(g, float(M))
            if ared is None or am is None:
                excluded += 1
                print(f"{g:>7}{M:>5}{(us if us else float('nan')):>9.4f}"
                      f"{(ared if ared else float('nan')):>12.6f}{'EXCLUDED':>12}")
                continue
            p4.append({"gamma": g, "M": M, "ustar": us, "A_red": ared, "A_meas": am,
                       "ratio": ared / am})
            print(f"{g:>7}{M:>5}{us:>9.4f}{ared:>12.6f}{am:>12.6f}{ared/am:>8.4f}"
                  f"{1-ared/am:>10.4f}")
    out["p4"], out["p4_excluded"] = p4, excluded
    print(f"  cells excluded on the convergence gate: {excluded}")
    verdicts = []
    for g in (0.40, 0.44):
        ser = [r for r in p4 if r["gamma"] == g]
        if len(ser) < 3:
            continue
        res = [1 - r["ratio"] for r in ser]
        falls = all(np.diff(res) < 0)
        verdicts.append(falls)
        print(f"  gamma={g}: 1-ratio over M = " + ", ".join(f"{v:.4f}" for v in res)
              + f"   ({'falls monotonically' if falls else 'NOT monotone'})")
    if not verdicts:
        print("  -> P4 UNDECIDED: too few M values survived the gate")
    elif all(verdicts):
        print("  -> P4 HOLDS: the deficit is the ADIABATIC ELIMINATION. Sharpening the")
        print("     timescale separation removes it, and S1/S2 are excluded above.")
    else:
        print("  -> P4 REFUTED: the deficit does NOT vanish with timescale separation, so it is")
        print("     not the elimination either -- and all three suspects are then dead.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
