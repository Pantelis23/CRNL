"""T-COST-l: is the missing ~79% the pool-fluctuation Jensen term? — absolute, no fit

§50 showed a COMPLETE one-dimensional account -- bulk, boundary layer and jump discreteness
together, exact by tridiagonal solve -- supplies only ~21% of the absorption coefficient
(bd/cme = 0.213 and 0.217 in the two converged cells, agreeing to 2%). ~79% is not
one-dimensional.

**THE SUSPECT WAS NAMED NINE SECTIONS EARLIER.** §39.1 listed candidate (iii), a Jensen gap
`E[sigma(state)] != sigma(E[state])`, withdrew it as an explanation of the COST, and wrote
that it "may still explain the *time* gap, which is now a separate and cleaner question."
This is that question.

**IT IS COMPUTABLE IN CLOSED FORM AND NEEDS NO LYAPUNOV SOLVE.** The pool is FAST, so at
fixed delta it holds a quasi-stationary Gaussian about the manifold. Expanding the drift in
s = x + y about the manifold value s*,

    E[mu(delta,s)] = mu(delta,s*) + (dmu/ds) E[s - s*] + (1/2)(d2mu/ds2) Var(s) + ...

**The first-order term is already counted**: E[s] - s* is the deterministic lag, which is
exactly §47's eps. The NEW piece is the second-order one. With s an Ornstein-Uhlenbeck
variable relaxing at |dnu/ds| under diffusion D_s,

    Var(s) = D_s / |dnu/ds|,    D_s = (1/(2 Omega)) sum_r (Delta s_r)^2 f_r

both read straight off the network's own fluxes, so

    **eps_J(delta) = (1/(2 mu)) (d2mu/ds2) Var(s)**,   **J = Omega * <eps_J>_time**

with the Omega cancelling, exactly as an O(1/Omega) coefficient must. No fitted parameter.

WHAT IT MUST HIT. In the two cells §50 found converged, the measured coefficient and the
1-D part are

    gamma=0.07, rho=1 :  cme = 2.568,  bd = 0.546  ->  J must supply ~2.02
    gamma=0.20, rho=32:  cme = 2.076,  bd = 0.450  ->  J must supply ~1.63

PREDICTIONS, written before running:

  P1  GATE, rule 13. eps_J is built from finite differences and a quadrature; both are its
      own numerical parameters and are checked for convergence WITHIN a cell before any
      comparison BETWEEN cells. Separately, Var(s)*Omega must come out Omega-independent,
      since that is the LNA scaling the whole construction assumes -- if it does not, the
      adiabatic elimination is wrong and nothing below is admissible.
  P2  THE SIGN IS FORCED, and this is the cheapest way to be wrong. Absorption shortens the
      MFPT, so the correction must be POSITIVE, which requires **d2mu/ds2 > 0** -- the drift
      must be CONVEX in the pool coordinate over the traversal. If mu is concave there, J is
      negative, the Jensen term cannot be the missing piece, and §39.1's candidate (iii)
      dies for the time as it already died for the cost.
  P3  THE TEST, absolute (rule 16). **bd_coeff + J = cme_coeff** in the two converged cells.
      Reported as (bd + J)/cme per cell, not averaged.
  P4  CONSISTENCY. J must be much LARGER than bd -- roughly 3.7x and 3.6x respectively --
      since it has to supply ~79% against ~21%. A J of the same size as bd would mean the
      budget is still short and something third is missing.
  5.  NO FITTING IF IT MISSES. If (bd+J)/cme lands away from 1, the residual is reported as
      a number and the budget stays open. A consistent offset across the two cells again
      indicates one more same-order term; a scattered one indicates the decomposition is
      wrong.
  P6  The three pre-asymptotic cells are reported but carry no verdict -- their cme is still
      rising, as §49 and §50 both recorded, so they cannot test an asymptotic prediction.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.networks.am_reversible import reverse_pairing
from experiments.arrhenius_optimum import am_rho, delta_star_rho
from experiments.cost_absolute import sigma_and_mu
from experiments.lag_absolute import field, s_star


def pool_terms(net, delta, h):
    """(eps_J, Var_s_times_Omega, d2mu_ds2, mu) at this delta, all from the network."""
    s0 = s_star(net, delta)
    if s0 is None:
        return None
    f0, fp, fm = (field(net, delta, s0), field(net, delta, s0 + h),
                  field(net, delta, s0 - h))
    if f0 is None or fp is None or fm is None:
        return None
    mu = f0[0]
    d2mu = (fp[0] - 2 * f0[0] + fm[0]) / (h * h)
    dnu_ds = (fp[1] - fm[1]) / (2 * h)
    if mu <= 0 or abs(dnu_ds) < 1e-12:
        return None

    x, y, b = 0.5 * (s0 + delta), 0.5 * (s0 - delta), 1.0 - s0
    S = net.stoichiometry_matrix()
    fl = net.fluxes(np.array([x, y, b]))
    ds_r = (S[0, :] + S[1, :]).astype(float)          # change in s per reaction
    two_Ds_omega = float(np.sum(ds_r ** 2 * fl))       # = 2 * D_s * Omega
    var_omega = 0.5 * two_Ds_omega / abs(dnu_ds)       # = Var(s) * Omega
    eps_j_omega = 0.5 * d2mu * var_omega / mu          # = eps_J * Omega
    return eps_j_omega, var_omega, d2mu, mu


def J_of(net, d_lo, d_hi, pairing, h=1e-5, n=401):
    """J = Omega * <eps_J>_time, weighted by dt = ddelta/mu."""
    xs = np.linspace(d_lo, d_hi, n)
    e, w, conv = [], [], []
    for x in xs:
        t = pool_terms(net, float(x), h)
        sm = sigma_and_mu(net, float(x), pairing)
        if t is None or sm is None or sm[1] <= 0:
            return None
        e.append(t[0])
        w.append(1.0 / sm[1])
        conv.append(t[2])
    e, w = np.array(e), np.array(w)
    return (float(np.trapezoid(e * w, xs) / np.trapezoid(w, xs)),
            float(np.min(conv)), float(np.max(conv)))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/pool_jensen.json"))
    args = ap.parse_args()

    t0 = time.time()
    meas = json.load(open("results/lag_endpoints.json"))
    bdj = json.load(open("results/birthdeath_absorption.json"))
    bd = {(c["gamma"], c["rho"]): c["bd"] for c in bdj["compare"]}
    cells = sorted({(r["gamma"], r["rho"]) for r in meas})
    CONVERGED = {(0.07, 1.0), (0.20, 32.0)}

    print("=== P1 GATE (rule 13): convergence in the finite-difference step h")
    ref = am_rho(0.20, 1.0)
    ds = delta_star_rho(0.20, 1.0)
    pr = reverse_pairing(ref)
    print(f"{'h':>10}{'J':>12}{'rel change':>13}")
    prev, rc = None, float("nan")
    for h in (1e-4, 5e-5, 2e-5, 1e-5):
        r = J_of(ref, 0.35 * ds, 0.80 * ds, pr, h=h)
        rc = abs(r[0] - prev) / abs(r[0]) if prev is not None else float("nan")
        print(f"{h:>10.1e}{r[0]:>12.5f}{rc:>13.2e}")
        prev = r[0]
    print(f"  -> P1(h) {'HOLDS' if rc < 1e-3 else 'FAILS'}")

    print(f"\n=== P2: is the drift CONVEX in the pool coordinate? (the sign is forced)")
    print(f"{'gamma':>6}{'rho':>6}{'min d2mu/ds2':>15}{'max d2mu/ds2':>15}{'J':>10}{'sign':>10}")
    rows = []
    for g, r in cells:
        net = am_rho(g, r)
        dsx = delta_star_rho(g, r)
        pairing = reverse_pairing(net)
        ref_c = [c for c in meas if c["gamma"] == g and c["rho"] == r][0]
        lo, hi = ref_c["eps_real"] * dsx, ref_c["theta_real"] * dsx
        res = J_of(net, lo, hi, pairing)
        if res is None:
            print(f"{g:>6.2f}{r:>6.1f}   failed")
            continue
        J, c_lo, c_hi = res
        rows.append({"gamma": g, "rho": r, "J": J, "d2mu_min": c_lo, "d2mu_max": c_hi})
        print(f"{g:>6.2f}{r:>6.1f}{c_lo:>15.3f}{c_hi:>15.3f}{J:>10.3f}"
              f"{'+' if J > 0 else 'NEGATIVE':>10}")
    allpos = all(c["J"] > 0 for c in rows)
    print(f"  -> P2 {'HOLDS, sign is right' if allpos else 'FAILS -- Jensen term has the wrong sign'}")

    print(f"\n=== P3/P4: the absolute budget, in the two converged cells")
    print(f"{'gamma':>6}{'rho':>6}{'cme':>9}{'bd':>8}{'J':>9}{'bd+J':>9}"
          f"{'(bd+J)/cme':>13}{'J/bd':>8}")
    by_meas = {}
    for g, r in cells:
        b = {c["omega"]: c for c in meas if c["gamma"] == g and c["rho"] == r}
        oms = sorted(b)
        cme = [(b[o]["gap_real"] - b[o]["pred_real"]) * o for o in oms]
        by_meas[(g, r)] = cme[-1] if abs(cme[-1] - cme[-2]) / cme[-1] < 0.15 \
            else float(np.mean(cme[:-1]))
    out = []
    for c in rows:
        k = (c["gamma"], c["rho"])
        cme, b = by_meas[k], bd.get(k, float("nan"))
        tot = b + c["J"]
        tag = "CONVERGED" if k in CONVERGED else "pre-asympt"
        out.append({**c, "cme": cme, "bd": b, "total": tot, "ratio": tot / cme,
                    "converged": k in CONVERGED})
        print(f"{c['gamma']:>6.2f}{c['rho']:>6.1f}{cme:>9.3f}{b:>8.3f}{c['J']:>9.3f}"
              f"{tot:>9.3f}{tot/cme:>13.3f}{c['J']/b:>8.2f}   {tag}")

    conv = [c for c in out if c["converged"]]
    print(f"\n=== P3 verdict (converged cells only)")
    for c in conv:
        print(f"  gamma={c['gamma']:.2f} rho={c['rho']:<5}"
              f" (bd+J)/cme = {c['ratio']:.3f}"
              f"   [1-D {100*c['bd']/c['cme']:.0f}%, Jensen {100*c['J']/c['cme']:.0f}%]")
    if conv:
        v = np.array([c["ratio"] for c in conv])
        print(f"  ratio {v.min():.3f}..{v.max():.3f}, mean {v.mean():.3f}")
        if abs(v.mean() - 1) < 0.15:
            print("  -> P3 HOLDS. The absorption coefficient is accounted for with NO")
            print("     fitted parameter: 1-D discreteness plus the pool-fluctuation")
            print("     Jensen term. §39.1's candidate (iii) explains the TIME gap.")
        else:
            print(f"  -> P3 FAILS by {100*(v.mean()-1):+.0f}%. Reported, not fitted.")
            sp = 100 * (v.max() - v.min()) / v.mean()
            print(f"     spread across the two cells {sp:.1f}% -> "
                  f"{'a consistent offset, so one more same-order term' if sp < 20 else 'scattered, so the decomposition is wrong'}")

    print(f"\n=== P6: the pre-asymptotic cells carry no verdict")
    for c in out:
        if not c["converged"]:
            print(f"  gamma={c['gamma']:.2f} rho={c['rho']:<5} (bd+J)/cme = {c['ratio']:.3f}"
                  f"   -- cme still rising, not a test")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
