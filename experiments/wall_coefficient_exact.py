"""The wall coefficient kappa, measured exactly instead of expanded.

FINDINGS 2 derives `c(eps) = kappa * eps^2` for irreversible AM with
`kappa = 3/2`, and FINDINGS 12 carries it to gamma > 0 as

    kappa_12(gamma) = (9/2) * lambda_antisym(gamma) = (3/2) * (1 - 2 gamma).

That is `lambda / D0` with the RESTORING GAIN taken at gamma and the DIFFUSION
taken at gamma = 0. But the reverse reactions are extra jumps along the decision
mode, so they add noise. At the symmetric point, with v = (1, -1, 0),

    D0(gamma) = sum_r (v . S_r)^2 a_r = (2/9) * (1 + gamma),

(the two recruitments contribute 1/9 each and the two homodimer reverses
gamma/9 each; the disagreement reaction and its reverse are orthogonal to v),
so the coefficient should be

    kappa(gamma) = lambda / D0 = (3/2) * (1 - 2 gamma) / (1 + gamma).

Both agree at gamma = 0, so FINDINGS 2 is untouched. They differ by 1/(1+gamma),
which is 31% at gamma = 0.45.

THE TEST. `quasipotential.landscape` gives W = -(1/Omega) ln P_ss exactly, and
the ridge curvature at the saddle IS kappa -- no expansion, no fit to an escape
probability. Two limits have to be taken and both are taken here rather than
assumed:

  * eps -> 0, because the quartic terms make the measured curvature drift with
    the fit window (FINDINGS 2's own table shows this drift, 1.586 -> 1.809).
    kappa(w) is fit as kappa_0 + c*w^2 and extrapolated.
  * Omega -> infinity, because W is only Omega-independent to leading WKB order.
    kappa_0(Omega) is fit as kappa_inf + c/Omega and extrapolated.

Prediction, written before running: the extrapolated ratio to
(3/2)(1-2g)/(1+g) is 1.000, and to FINDINGS 12's kappa it is 1/(1+gamma).

Note this does NOT explain FINDINGS 12's fitted slopes (0.795 / 0.626 / 0.419 /
0.497 at gamma = 0.05 / 0.15 / 0.30 / 0.45), which are non-monotone and so
cannot be a smooth 1/(1+gamma). Whether the corrected kappa improves 12's
collapse is a separate measurement -- see `--recollapse`.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from crnl.networks.am_reversible import am_reversible, delta_star, lambda_antisym
from crnl.quasipotential import P_FLOOR, barriers, landscape, ridge_profile

WINDOWS = (0.10, 0.15, 0.20, 0.25, 0.30, 0.35)


def kappa_predicted(gamma: float) -> float:
    """lambda / D0 with the diffusion evaluated at the same gamma."""
    return lambda_antisym(gamma) / ((2.0 / 9.0) * (1.0 + gamma))


def kappa_findings12(gamma: float) -> float:
    """(9/2) * lambda: the diffusion frozen at its gamma = 0 value."""
    return 1.5 * (1.0 - 2.0 * gamma)


#: Two decades of headroom above the probability floor. Probabilities a few x
#: the floor are already round-off-contaminated even though they are nonzero,
#: and a curvature fit through them looks perfectly well behaved -- it produced
#: the gamma=0.35, Omega=400 outlier (kappa 0.343 -> 0.309) in the first run.
SAFETY = 100.0

#: Minimum lattice sites across the NARROWEST fit window. eps is quantised at
#: 1/Omega, so a narrow window at modest Omega is a parabola through a handful
#: of points; that produced the gamma=0.45, Omega=150 outlier. 21 sites is +-10.
MIN_SITES = 21


def usable_omega(gamma: float, omega: float, barrier: float) -> tuple[bool, str]:
    """Whether Omega falls in this instrument's window, and why not if not.

    The two guards pull in OPPOSITE directions, which is the awkward and
    important fact about this route: discreteness needs large Omega, the
    probability floor needs small Omega. The window can be empty -- it is for
    gamma <= 0.2, where the barrier is too deep to resolve before the lattice is
    fine enough -- and an empty window is a real limit, not a tuning failure.
    """
    lo = MIN_SITES / (2.0 * min(WINDOWS) * delta_star(gamma))
    hi = (-np.log(P_FLOOR) - np.log(SAFETY)) / max(barrier, 1e-12)
    if omega < lo:
        return False, f"lattice: need Omega >= {lo:.0f}"
    if omega > hi:
        return False, f"floor: need Omega <= {hi:.0f}"
    return True, ""


def curvature(land: dict, gamma: float) -> tuple[float, dict]:
    """Ridge curvature extrapolated to a zero-width fit window.

    kappa = -c0 of the quadratic fit, which is independent of where the ridge
    peak actually sits -- so this stays correct for a tilted landscape whose
    peak is not at eps = 0.
    """
    ds = delta_star(gamma)
    eps, W = ridge_profile(land)
    ok = np.isfinite(W)
    ceiling = land["dW_max"] - np.log(SAFETY) / land["omega"]
    per_window = {}
    for w in WINDOWS:
        inside = np.abs(eps) < w * ds
        m = ok & inside
        if m.sum() < MIN_SITES or not np.all(ok[inside]):
            continue
        if W[m].max() > ceiling:          # too close to the round-off floor
            continue
        per_window[w] = -float(np.polyfit(eps[m], W[m], 2)[0])
    if len(per_window) < 3:
        return float("nan"), per_window
    ws = np.array(sorted(per_window))
    ks = np.array([per_window[w] for w in ws])
    k0 = float(np.polyfit(ws ** 2, ks, 1)[1])      # extrapolate w -> 0
    return k0, per_window


def run(gammas, omegas, out: pathlib.Path | None) -> list[dict]:
    rows = []
    for gamma in gammas:
        pred, p12 = kappa_predicted(gamma), kappa_findings12(gamma)
        print(f"\ngamma = {gamma}   lambda/D0 = {pred:.6f}   "
              f"FINDINGS-12 kappa = {p12:.6f}   ratio = {pred / p12:.6f}")
        print(f"{'Omega':>7} {'unres':>8} {'dWmax':>7} {'kappa(w->0)':>12} "
              f"{'/pred':>8} {'/k12':>8}")
        got = []
        # locate the usable window from a cheap barrier estimate, then walk it
        probe = 100
        try:
            est = barriers(am_reversible(gamma), probe, float(probe))["dW_x"]
        except (ValueError, RuntimeError):
            est = float("nan")
        if not np.isfinite(est):
            print("  no barrier estimate at Omega=100; instrument window is empty")
            continue

        oms = [o for o in omegas if usable_omega(gamma, o, est)[0]]
        if not oms:
            lo = MIN_SITES / (2.0 * min(WINDOWS) * delta_star(gamma))
            hi = (-np.log(P_FLOOR) - np.log(SAFETY)) / est
            print(f"  instrument window is EMPTY: lattice needs Omega >= {lo:.0f}, "
                  f"floor needs Omega <= {hi:.0f} (barrier ~ {est:.4f})")
            continue
        for om in oms:
            try:
                land = landscape(am_reversible(gamma), int(om), float(om))
            except RuntimeError as e:
                print(f"{om:>7}   solve rejected: {str(e)[:50]}")
                continue
            k0, per_w = curvature(land, gamma)
            if not np.isfinite(k0):
                print(f"{om:>7} {land['n_unresolved']:>8} {land['dW_max']:>7.4f}"
                      "   ridge not resolved across enough windows")
                continue
            print(f"{om:>7} {land['n_unresolved']:>8} {land['dW_max']:>7.4f} "
                  f"{k0:>12.6f} {k0 / pred:>8.4f} {k0 / p12:>8.4f}")
            got.append((float(om), k0))
            rows.append({"gamma": gamma, "omega": float(om), "kappa": k0,
                         "kappa_pred": pred, "kappa_findings12": p12,
                         "per_window": {str(k): v for k, v in per_w.items()},
                         "n_unresolved": land["n_unresolved"]})
        if len(got) >= 3:
            oms = np.array([g[0] for g in got])
            ks = np.array([g[1] for g in got])
            kinf = float(np.polyfit(1.0 / oms, ks, 1)[1])
            print(f"  Omega -> inf:  kappa = {kinf:.6f}   "
                  f"/pred = {kinf / pred:.4f}   /k12 = {kinf / p12:.4f}   "
                  f"(1/(1+g) = {1 / (1 + gamma):.4f})")
            rows.append({"gamma": gamma, "omega": "inf", "kappa": kinf,
                         "kappa_pred": pred, "kappa_findings12": p12})
    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(rows, indent=2))
        print(f"\nwrote {out}")
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.20, 0.30, 0.35, 0.40, 0.45])
    ap.add_argument("--omegas", type=int, nargs="+",
                    default=[100, 150, 200, 300, 400, 600])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/wall_coefficient_exact.json"))
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    if args.quick:
        args.gammas, args.omegas = [0.35, 0.45], [100, 150, 200]
    run(args.gammas, args.omegas, args.out)


if __name__ == "__main__":
    main()
