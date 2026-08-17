"""T-DEPTH-d: the last free number -- nats of reliability per molecule, and does IT transfer?

§76 reduced everything this project set out to measure to a single element quantity: the
per-stage error probability eps. Reliability is eps; depth is c*/eps; there is no third
thing. **So the whole founding question now rests on one coefficient: how fast does ln(1/eps)
grow with molecule count?** Define

    eta  =  d ln(1/eps) / d Omega        [nats of reliability per molecule]

**This is the transfer test that matters, and the earlier ones did not.** §67 priced
dissipation per e-fold and found no counterpart on a driven element; §68 priced the affinity
floor and found it varies fourfold inside one family. Both were pricing quantities that §73
later showed are invisible to composition. eta is the only quantity left, so if anything about
a restoring element is substrate-independent, it is this -- and if eta is substrate-specific
then nothing is, and the chemistry's entire contribution is one number.

DEFINITION, stated because it is a choice. eps is the Gaussian readout of the element's own
intrinsic noise: eps = Phi(-Delta/sigma) with Delta the rail separation and sigma the rail's
quasi-stationary width, both measured exactly from the CME as in §75. That is precisely the
quantity the chemically-coupled cascade uses (§75), so it is the right eps, but it is a
readout convention and not the only conceivable one.

PREDICTIONS, written before running.

  P1  GATE, and per rule 20 it is a CONVERGENCE test, not a tolerance. ln(1/eps) must become
      LINEAR in Omega at fixed landscape, i.e. the local slope eta must settle. Report the
      slope between successive Omega and require it to stop moving. If it does not settle,
      eta does not exist and nothing below counts.
  P2  **AM's eta looks already-linear in §75's own numbers** -- Delta/sigma = 14.740, 20.915,
      29.622 at Omega = 60, 120, 240 gives ratios 1.4189 and 1.4163 against sqrt(2) = 1.4142,
      so (Delta/sigma)^2 ~ Omega to 0.3% and eta = 1.83 at gamma = 0.20, 9.88 at gamma = 0.05.
      **Predicted: eta is constant in Omega and strongly dependent on gamma.**
  P3  **Schloegl's looks NOT already-linear** -- Delta/sigma = 3.0037 at lambda*Omega = 800 and
      6.8718 at 3200, a ratio of 2.288 against sqrt(4) = 2. **Predicted: that is a small-Omega
      contamination of the rail width, not a different physics**, because §75 established the
      exact collapse in lambda*Omega and a collapse cannot coexist with two different
      exponents. It should straighten at larger Omega. If it does not, §75's collapse and this
      section disagree and one of them is wrong.
  P4  **THE TEST. Does eta transfer?** Compare AM's eta at matched conditions against
      Schloegl's. **Predicted: NO** -- eta should depend on the landscape shape, since it is
      the one place the chemistry can still show up. A transferring eta would be the single
      most surprising result in this project and would mean restoration has a universal
      exchange rate after all.
  P5  **RULE 9.** Sweep the landscape shape within each substrate too -- gamma for AM, the
      root spacing for Schloegl -- so that "depends on the landscape" is measured rather than
      asserted.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
from scipy.stats import norm

from crnl.networks.am_reversible import delta_star
from experiments.chemical_channel_noise import am_rail_width, rail_width


def ln_inv_eps(delta, sigma):
    """ln(1/eps) with eps = Phi(-Delta/sigma), computed in logs so the tail does not vanish."""
    z = delta / sigma
    return -float(norm.logcdf(-z))


def schlogl_eta(lam, omega, base=(0.1, 1.0, 1.9)):
    r = rail_width(omega, lam * base[0], lam * base[1], lam * base[2])
    if r is None:
        return None
    delta = lam * (base[2] - base[0]) / 2.0
    return ln_inv_eps(delta, r["sd_exact"]), r["sd_exact"], delta


def am_eta(gamma, omega):
    """None when the 2-D stationary solve is not trustworthy -- the engine guards this
    itself, and the first version of this file swallowed the failure in a bare `except`
    and skipped the cell silently. Failures are counted and printed instead (rule 10)."""
    try:
        sd, _ = am_rail_width(gamma, omega)
    except RuntimeError:
        return None
    ds = float(delta_star(gamma))
    return ln_inv_eps(ds, sd), sd, ds


def local_slopes(xs, ys):
    return [(ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]) for i in range(len(xs) - 1)]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/nats_per_molecule.json"))
    args = ap.parse_args()
    out = {"am": [], "schlogl": []}

    skipped = []
    print("=== P1/P2: AM -- does ln(1/eps) become linear in Omega, and what is eta?")
    print(f"{'gamma':>7}{'Omega':>7}{'delta*':>9}{'sigma':>10}{'ln(1/eps)':>12}"
          f"{'local eta':>11}")
    for g in (0.30, 0.20, 0.05):
        oms, ys = [], []
        for om in (60, 120, 240, 400):
            r = am_eta(g, om)
            if r is None:
                skipped.append((g, om))
                print(f"{g:>7.2f}{om:>7}   EXCLUDED: the stationary solve is not trustworthy")
                continue
            y, sd, ds = r
            oms.append(om); ys.append(y)
            out["am"].append({"gamma": g, "omega": om, "ln_inv_eps": y, "sigma": sd})
            sl = local_slopes(oms, ys)
            print(f"{g:>7.2f}{om:>7}{ds:>9.5f}{sd:>10.6f}{y:>12.3f}"
                  + (f"{sl[-1]:>11.4f}" if sl else f"{'--':>11}"))
        if len(ys) >= 3:
            sl = local_slopes(oms, ys)
            drift = abs(sl[-1] - sl[-2]) / abs(sl[-1])
            print(f"         -> eta = {sl[-1]:.4f} nats/molecule,"
                  f" last-step drift {100*drift:.2f}%")

    print("\n=== P1/P3: Schloegl -- does it straighten at larger Omega?")
    print(f"{'lambda':>8}{'Omega':>8}{'sigma':>10}{'ln(1/eps)':>12}{'local eta':>11}")
    for lam in (1.0, 4.0):
        oms, ys = [], []
        for om in (400, 1600, 6400, 25600):
            r = schlogl_eta(lam, om)
            if r is None:
                continue
            y, sd, delta = r
            oms.append(om); ys.append(y)
            out["schlogl"].append({"lam": lam, "omega": om, "ln_inv_eps": y, "sigma": sd})
            sl = local_slopes(oms, ys)
            print(f"{lam:>8.1f}{om:>8}{sd:>10.6f}{y:>12.3f}"
                  + (f"{sl[-1]:>11.6f}" if sl else f"{'--':>11}"))
        if len(ys) >= 3:
            sl = local_slopes(oms, ys)
            drift = abs(sl[-1] - sl[-2]) / abs(sl[-1])
            print(f"          -> eta = {sl[-1]:.6f} nats/molecule,"
                  f" last-step drift {100*drift:.2f}%"
                  f"   {'SETTLED' if drift < 0.10 else 'still moving'}")

    print("\n=== P4/P5: does eta transfer, and what does it depend on?")
    print(f"{'element':>22}{'landscape':>14}{'eta':>12}")
    tab = []
    for g in (0.30, 0.20, 0.05):
        pts = [(om, am_eta(g, om)) for om in (120, 240, 400)]
        pts = [(om, r[0]) for om, r in pts if r is not None]
        if len(pts) < 2:
            print(f"{'AM':>22}{f'gamma={g}':>14}{'  no usable pair':>12}")
            continue
        e = (pts[-1][1] - pts[-2][1]) / (pts[-1][0] - pts[-2][0])
        tab.append(("AM", f"gamma={g}", e))
        print(f"{'AM':>22}{f'gamma={g}':>14}{e:>12.4f}"
              + f"   (Omega {pts[-2][0]}->{pts[-1][0]})")
    for lam in (1.0, 4.0):
        ys = [schlogl_eta(lam, om)[0] for om in (6400, 25600)]
        e = (ys[1] - ys[0]) / (25600 - 6400)
        tab.append(("Schlogl", f"lambda={lam}", e))
        print(f"{'Schloegl':>22}{f'lambda={lam}':>14}{e:>12.6f}")
    for spread in (0.6, 0.9):
        r1, r3 = 1.0 - spread, 1.0 + spread
        ys = []
        for om in (6400, 25600):
            rr = rail_width(om, r1, 1.0, r3)
            ys.append(ln_inv_eps((r3 - r1) / 2.0, rr["sd_exact"]))
        e = (ys[1] - ys[0]) / (25600 - 6400)
        tab.append(("Schlogl", f"spread={spread}", e))
        print(f"{'Schloegl':>22}{f'spread={spread}':>14}{e:>12.6f}")

    vals = [t[2] for t in tab]
    print(f"\n  eta spans {min(vals):.6f} .. {max(vals):.4f}"
          f"  -- a factor of {max(vals)/min(vals):.0f}")
    print(f"  -> P4 {'eta does NOT transfer: it is the one place the chemistry survives' if max(vals)/min(vals) > 2 else 'eta TRANSFERS -- restoration has a universal exchange rate'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    print(f"\n  cells excluded on an untrustworthy stationary solve: {len(skipped)}"
          + (f"  {skipped}" if skipped else ""))
    args.out.write_text(json.dumps({"rows": out, "skipped": skipped,
                                    "eta": [[a, b, c] for a, b, c in tab]},
                                   indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
