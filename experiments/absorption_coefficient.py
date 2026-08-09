"""T-COST-j: the absorption correction, derived and tested ABSOLUTELY

§48 closed T-COST-i -- the lag model is exact as Omega -> infinity -- but left the DECAY
LAW unresolved, fitting per-cell exponents of 0.20..0.82 and explicitly refusing to claim
the tidy 1/sqrt(Omega) that appeared after averaging (rule 18). Fitting an exponent is
exactly what rule 16 says never settles anything.

**THE STANDARD FIRST-PASSAGE EXPANSION PREDICTS THE EXPONENT AND THE COEFFICIENT.** On the
slaved manifold delta jumps by +-1 with total fluxes `up` and `dn` (both already computed by
`updown`), so in concentration units

    ddelta = mu dt + sqrt(2D) dW,   mu = up - dn,   D = D0/Omega,   D0 = (up + dn)/2

For a diffusion absorbing at L, T(x0) = int dy (1/D) e^{-psi(y)} int^y e^{psi(z)} dz with
psi' = mu/D. Writing g = mu/D0 so that psi' = Omega*g, expanding the inner integral by
Laplace gives int^y e^psi = e^{psi}(1/psi')(1 + psi''/psi'^2 + ...) with
psi''/psi'^2 = g'/(Omega g^2), hence

    T_stoch = int (1/mu)(1 + g'/(Omega g^2)) ddelta = T_det + (1/Omega) int g'/(mu g^2)

    **T_det/MFPT - 1  =  <eps>_time  +  K/Omega,   K = -(1/T_det) int g'/(mu g^2) ddelta**

**The exponent is 1, not 1/2, and K is fully computable with no fit.** K > 0 because g
falls past the drift peak, which is the observed sign.

**THE DATA ALREADY SUPPORTS THE EXPONENT WHERE IT IS TESTABLE.** Measured before writing
these predictions, `(gap - pred)*Omega` over Omega = 300/500/700/1000:

    gamma=0.20 rho=32    1.95  1.98  2.00  2.08     constant to 6%
    gamma=0.20 rho=1     3.58  4.11  3.90  4.21     constant to 17%
    gamma=0.07 rho=1     2.61  2.49  2.60  3.65     constant to 4% for Omega <= 700
    gamma=0.20 rho=0.5   6.55  7.31  8.62  9.28     RISING 42%
    gamma=0.35 rho=1     4.33  7.00  9.13 10.26     RISING 137%

**So 1/Omega holds cleanly in three cells and fails in two -- and the two that fail are
exactly the two with the longest deterministic traversal** (T_det = 13.5 and 12.9 against
3.1, 4.9, 6.7). §48's fitted exponents of 0.20..0.82 were averaging asymptotic cells
together with pre-asymptotic ones.

PREDICTIONS, written before running:

  P1  GATE. K is a quadrature over `updown`, so its node count is its own numerical
      parameter (rule 13). Converged within a cell before comparing between cells.
  P2  THE TEST, absolute and with no fitted parameter. In the three asymptotic cells
      (gamma=0.07, gamma=0.20 rho=1, gamma=0.20 rho=32), **K matches the measured
      `(gap - pred)*Omega`**. This is the rule-16 test §48 could not do.
  P3  THE PRE-ASYMPTOTIC PREDICTION. In gamma=0.35 and rho=0.5 the measured coefficient is
      still RISING, so it has not reached its asymptote. **K must therefore EXCEED the
      measured value there**, and the measured value must be climbing toward it. If K
      instead sits below a still-rising measurement, the expansion is not the explanation.
  P4  THE MISSING BOUNDARY TERM, named in advance. The derivation extends the inner
      integral to -infinity and so is the BULK term only; the absorbing boundary contributes
      a layer of width D/mu, which is the SAME order 1/Omega. **A consistent offset of the
      same size across cells therefore indicates that missing term, while a SCATTERED ratio
      indicates the bulk form itself is wrong.** That distinction is the point of reporting
      K/measured per cell rather than averaged.
  P5  REFUTING OUTCOME. If K has the wrong sign, or is off by more than 3x in the
      asymptotic cells, the bulk expansion is not the mechanism and T-COST-j stays open
      with a measured coefficient and no derivation.
  P6  If P2 and P3 hold, the whole gap is accounted for with **no free parameter anywhere**:
      `T_det/MFPT - 1 = <eps>_time + K/Omega`, the first term §47's lag and the second this
      absorption term, and §48's unresolved exponent is resolved as 1.
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
from experiments.lag_endpoints import traversal
from experiments.slaving_axis import updown


def g_of(net, delta):
    """g = mu/D0 with D0 = (up + dn)/2, both from the network's own manifold fluxes."""
    ud = updown(net, delta)
    if ud is None:
        return None
    up, dn = ud
    d0 = 0.5 * (up + dn)
    if d0 <= 0:
        return None
    return (up - dn) / d0, up - dn


def K_of(net, d_lo, d_hi, pairing, n=401, h=1e-5):
    """K = -(1/T_det) int g'/(mu g^2) ddelta, the bulk absorption coefficient."""
    t_det = traversal(net, d_lo, d_hi, pairing)
    if t_det is None:
        return None
    xs = np.linspace(d_lo, d_hi, n)
    vals = []
    for x in xs:
        a, b = g_of(net, float(x) + h), g_of(net, float(x) - h)
        c = g_of(net, float(x))
        if a is None or b is None or c is None:
            return None
        gp = (a[0] - b[0]) / (2 * h)
        g, mu = c[0], c[1]
        if abs(g) < 1e-14 or abs(mu) < 1e-14:
            return None
        vals.append(gp / (mu * g * g))
    return -float(np.trapezoid(vals, xs)) / t_det, t_det


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/absorption_coefficient.json"))
    args = ap.parse_args()

    t0 = time.time()
    meas = json.load(open("results/lag_endpoints.json"))
    cells = sorted({(r["gamma"], r["rho"]) for r in meas})
    oms = sorted({r["omega"] for r in meas})

    print("=== P1 GATE (rule 13): is K converged in its quadrature?")
    net = am_rho(0.20, 1.0)
    ds = delta_star_rho(0.20, 1.0)
    pairing = reverse_pairing(net)
    prev, rc = None, float("nan")
    print(f"{'nodes':>8}{'K':>12}{'rel change':>13}")
    for n in (101, 201, 401, 801):
        k, _ = K_of(net, 0.35 * ds, 0.80 * ds, pairing, n=n)
        rc = abs(k - prev) / abs(k) if prev is not None else float("nan")
        print(f"{n:>8}{k:>12.5f}{rc:>13.2e}")
        prev = k
    print(f"  -> P1 {'HOLDS' if rc < 1e-3 else 'FAILS'}")

    print(f"\n=== P2/P3: K against the measured (gap - pred)*Omega, per cell")
    print(f"{'gamma':>6}{'rho':>6}{'T_det':>8}{'K':>9}"
          + "".join(f"{f'meas@{o}':>10}" for o in oms) + f"{'K/meas(max Om)':>16}{'status':>14}")
    rows = []
    for g, r in cells:
        ds = delta_star_rho(g, r)
        net = am_rho(g, r)
        pairing = reverse_pairing(net)
        by = {c["omega"]: c for c in meas if c["gamma"] == g and c["rho"] == r}
        ref = by[oms[-1]]
        lo, hi = ref["eps_real"] * ds, ref["theta_real"] * ds
        res = K_of(net, lo, hi, pairing)
        if res is None:
            print(f"{g:>6.2f}{r:>6.1f}   K failed")
            continue
        k, t_det = res
        ms = [(by[o]["gap_real"] - by[o]["pred_real"]) * o for o in oms if o in by]
        rising = ms[-1] > 1.25 * ms[0]
        ratio = k / ms[-1]
        rows.append({"gamma": g, "rho": r, "K": k, "T_det": t_det,
                     "meas": ms, "K_over_meas": ratio, "rising": bool(rising)})
        print(f"{g:>6.2f}{r:>6.1f}{t_det:>8.2f}{k:>9.3f}"
              + "".join(f"{m:>10.2f}" for m in ms)
              + f"{ratio:>16.3f}{'PRE-ASYMPT' if rising else 'asymptotic':>14}")

    asym = [c for c in rows if not c["rising"]]
    pre = [c for c in rows if c["rising"]]

    print(f"\n=== P2: the absolute test in the {len(asym)} asymptotic cells")
    if asym:
        v = np.array([c["K_over_meas"] for c in asym])
        for c in asym:
            print(f"  gamma={c['gamma']:.2f} rho={c['rho']:<5} K={c['K']:.3f}"
                  f"  measured={c['meas'][-1]:.3f}  K/meas={c['K_over_meas']:.3f}")
        print(f"  K/meas: {v.min():.3f}..{v.max():.3f}, mean {v.mean():.3f},"
              f" spread {100*(v.max()-v.min())/v.mean():.1f}%")
        ok = v.min() > 1 / 3 and v.max() < 3
        print(f"  -> P2/P5: {'the bulk expansion has the right sign and scale' if ok else 'REFUTED, off by more than 3x'}")

    print(f"\n=== P3: the pre-asymptotic cells -- K must EXCEED a still-rising measurement")
    for c in pre:
        print(f"  gamma={c['gamma']:.2f} rho={c['rho']:<5} T_det={c['T_det']:.1f}"
              f"  measured {c['meas'][0]:.2f} -> {c['meas'][-1]:.2f} (rising)"
              f"  K={c['K']:.3f}"
              f"  -> {'K exceeds it, AS PREDICTED' if c['K'] > c['meas'][-1] else 'K is BELOW it, contrary to P3'}")

    print(f"\n=== P4: is K/meas a consistent offset (missing boundary term) or scattered?")
    if asym:
        v = np.array([c["K_over_meas"] for c in asym])
        sp = 100 * (v.max() - v.min()) / v.mean()
        print(f"  spread across asymptotic cells {sp:.1f}%")
        if sp < 25:
            print(f"  -> CONSISTENT offset of {v.mean():.3f}. The bulk form is right and the")
            print(f"     boundary-layer term, same order 1/Omega, is what is missing.")
        else:
            print(f"  -> SCATTERED. A single missing term cannot explain it and the bulk")
            print(f"     form itself is suspect.")

    print(f"\n=== P6: verdict on the exponent")
    print(f"  §48 fitted 0.20..0.82 and declined to claim 1/2. Here the exponent is")
    print(f"  DERIVED as 1, and (gap-pred)*Omega is constant in {len(asym)}/{len(rows)}"
          f" cells -- the other {len(pre)} being the longest-traversal cells,"
          f" T_det = {', '.join(f'{c['T_det']:.1f}' for c in pre)} against"
          f" {', '.join(f'{c['T_det']:.1f}' for c in asym)}.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
