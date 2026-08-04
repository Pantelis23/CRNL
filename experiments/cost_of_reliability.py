"""What does a nat of reliability COST? — the founding question, measured exactly

This project exists to ask what restoration costs. It has a minimum affinity
(§9.1: A_c = 3 ln 2), an error exponent (§1: P ~ exp(-c Omega)), and a fuel lifetime
(§20), and **no relation connecting them**. §9.2 once claimed a dissipation minimum near
gamma ~ 0.3 and it is in THEORIES §4 as a withdrawn result -- killed because the decision
threshold was held fixed while delta*(gamma) shrank, so above gamma ~ 0.42 the threshold
sat OUTSIDE the landscape and "deciding" meant fluctuating past the attractor.

THREE THINGS ARE AVAILABLE NOW THAT WERE NOT THEN, and all three bear on exactly the
failure that killed it:

  * theta is scaled with delta*, which is the fix §9.2's withdrawal named.
  * The start is placed ON the slaved manifold (§36), which is the fourteen-section
    misattribution corrected an hour ago -- an off-manifold start charges the decision
    for a relaxation transient that is not part of it.
  * P(error) is solved for DIRECTLY (§35), so reliability is exact at any depth rather
    than floored at 1e-12 by cancellation.

THE QUANTITY. Both sides of the trade are exact linear solves on the same generator:

    L(Omega)  = -ln P(error)                    [nats of reliability]
    Sigma(Omega) = E[ entropy produced before absorption ]   [k_B]

`first_passage` solves `Q_tt T = -1` for the mean time; the expected entropy production
solves the SAME system with the local entropy rate as its source,

    Q_tt Sigma = -sigma_local,   sigma_local(n) = sum_j a_j(n) ln[ a_j(n) / a_rev(n+S_j) ]

so the cost of a decision is exact, not a rate times a time. The ratio

    **R = Sigma / L   [k_B per nat of reliability]**

is the thing the founding question is about, and it should be Omega-INDEPENDENT: both
numerator and denominator are extensive, so R is a property of the CHEMISTRY, not of how
many molecules were thrown at it.

PREDICTIONS, written before running:

  P1  GATE, and the whole framing rests on it. R is Omega-independent once Omega is
      large enough for the exponential regime -- Sigma and L both grow linearly and
      their ratio flattens. If R keeps drifting with Omega, R is not a property of the
      chemistry and nothing below means anything.
  P2  R has a MINIMUM in gamma. This is §9.2's claim, redone with the instrument that
      killed it plus §36's start-state fix. A minimum would be a genuine design
      principle: **the drive at which restoration is cheapest per nat**, and a candidate
      answer to what a near-ideal restoring switch should cost.
  P3  IF a minimum exists, its VALUE is a candidate floor on the cost of reliability in
      this network class. It is compared against recognisable scales -- ln 2, 1, 2, and
      the cycle affinity floor A_c = 3 ln 2 = 2.079 -- **as a comparison, not as a fit**.
      Landing on one of them would be a real result; landing between them is the more
      likely outcome and gets reported as a number, not decorated.
  P4  IF R is monotone in gamma with no interior minimum, §9.2's withdrawal is confirmed
      with better instruments and the honest answer is "harder drive buys cheaper
      reliability, without limit until the affinity floor". That is the outcome that
      costs a hoped-for result and it is reported the same way.
  P5  The cost per nat must EXCEED some positive floor at every gamma -- restoration
      cannot be free. R > 0 is trivially true; what is informative is whether R stays
      bounded away from 0 as the drive grows, or whether it can be made arbitrarily
      cheap by driving harder. The second would be the more surprising result and would
      say the affinity floor, not the dissipation, is the binding constraint.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
import scipy.sparse.linalg as spla

from crnl.approximations import propensities_batch
from crnl.cme import enumerate_states, generator
from crnl.networks.am_reversible import am_reversible, delta_star, reverse_pairing
from crnl.vectorized import compile_network
from experiments.slaving_axis import slaved


def sigma_local(net, comp, states, pairing):
    """sum_j a_j(n) ln[a_j(n)/a_rev(n+S_j)], vectorised over all states."""
    S = net.stoichiometry_matrix().astype(np.int64)
    A = propensities_batch(comp, states.astype(float))
    out = np.zeros(len(states))
    for j in range(net.n_reactions):
        rev = int(pairing[j])
        n2 = states + S[:, j]
        ok = (n2 >= 0).all(axis=1) & (A[:, j] > 0)
        if not ok.any():
            continue
        Arev = propensities_batch(comp, n2[ok].astype(float))[:, rev]
        good = Arev > 0
        idx = np.where(ok)[0][good]
        out[idx] += A[idx, j] * np.log(A[idx, j] / Arev[good])
    return out


def cell(gamma, omega, eps, theta):
    net = am_reversible(gamma)
    ds = delta_star(gamma)
    x0 = eps * ds
    st = slaved(net, x0)
    nb = int(round(st[2] * omega))
    rest = omega - nb
    d0 = max(1, int(round(x0 * omega)))
    if (rest - d0) % 2:
        d0 -= 1
    n0 = np.array([(rest + d0) // 2, (rest - d0) // 2, nb], dtype=np.int64)
    thr = max(2, int(round(theta * ds * omega)))

    states, index = enumerate_states(3, int(omega))
    absorb = np.array([abs(int(s[0]) - int(s[1])) >= thr for s in states])
    Q = generator(net, int(omega), float(omega))
    tr = np.where(~absorb)[0]
    tmap = {int(i): r for r, i in enumerate(tr)}
    Qtt = Q[tr][:, tr].tocsr()
    si = tmap[index[tuple(n0)]]

    # reliability, solved for directly (no 1 - split)
    fav = np.array([int(s[1]) > int(s[0]) for s in states])[absorb].astype(float)
    b = -(Q[tr][:, np.where(absorb)[0]].tocsr() @ fav)
    p_err = float(spla.spsolve(Qtt, b)[si])

    # expected entropy produced before absorption
    comp = compile_network(net, float(omega))
    sig = sigma_local(net, comp, states, reverse_pairing(net))[tr]
    Sig = float(spla.spsolve(Qtt, -sig)[si])

    return {"omega": omega, "p_err": p_err, "L": -np.log(p_err), "Sigma": Sig,
            "R": Sig / (-np.log(p_err)), "thr": int(thr),
            "eps_realised": float(d0 / (ds * omega))}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45])
    ap.add_argument("--omegas", type=int, nargs="+", default=[100, 200, 300, 400])
    ap.add_argument("--eps-frac", type=float, default=0.35)
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/cost_of_reliability.json"))
    args = ap.parse_args()

    t0 = time.time()
    print("R = Sigma / L  =  k_B of entropy produced per NAT of reliability bought")
    print("start ON the slaved manifold (§36); theta scaled with delta* (§9.2's fix)\n")
    print(f"{'gamma':>7}{'Omega':>7}{'P(error)':>13}{'L (nats)':>11}"
          f"{'Sigma (k_B)':>13}{'R':>10}")
    out = {}
    for g in args.gammas:
        rows = []
        for om in args.omegas:
            try:
                r = cell(g, om, args.eps_frac, args.theta)
            except Exception as e:
                print(f"{g:>7.2f}{om:>7}   SKIPPED ({type(e).__name__})")
                continue
            if not (np.isfinite(r["R"]) and r["L"] > 0):
                print(f"{g:>7.2f}{om:>7}   SKIPPED (degenerate)")
                continue
            rows.append(r)
            print(f"{g:>7.2f}{om:>7}{r['p_err']:>13.4e}{r['L']:>11.4f}"
                  f"{r['Sigma']:>13.4f}{r['R']:>10.4f}")
        if len(rows) >= 2:
            Rs = np.array([r["R"] for r in rows])
            drift = 100 * (Rs[-1] - Rs[-2]) / Rs[-1]
            print(f"{'':>7}{'':>7}{'':>13}{'':>11}{'R drift last step':>13}"
                  f"{drift:>9.2f}%")
        out[str(g)] = rows
        print()

    print("=== P1 gate: is R Omega-independent?")
    print(f"{'gamma':>7}{'R at largest Omega':>21}{'drift over last step':>23}")
    conv = {}
    for g, rows in out.items():
        if len(rows) < 2:
            continue
        Rs = [r["R"] for r in rows]
        d = 100 * (Rs[-1] - Rs[-2]) / Rs[-1]
        conv[float(g)] = Rs[-1]
        print(f"{float(g):>7.2f}{Rs[-1]:>21.4f}{d:>22.2f}%")

    print("\n=== P2/P3/P4: the cost of reliability against the drive")
    if conv:
        gs = np.array(sorted(conv))
        Rs = np.array([conv[g] for g in gs])
        A = -3.0 * np.log(gs)
        print(f"{'gamma':>7}{'affinity A':>12}{'R (k_B/nat)':>14}")
        for g, a, r in zip(gs, A, Rs):
            print(f"{g:>7.2f}{a:>12.4f}{r:>14.4f}")
        i = int(np.argmin(Rs))
        interior = 0 < i < len(Rs) - 1
        print(f"\n  minimum R = {Rs[i]:.4f} k_B/nat at gamma = {gs[i]:.2f} "
              f"(A = {A[i]:.4f})")
        print(f"  interior minimum? {'YES -- P2' if interior else 'NO -- P4, monotone'}")
        for name, val in (("ln 2", np.log(2)), ("1", 1.0), ("2", 2.0),
                          ("A_c = 3 ln 2", 3 * np.log(2))):
            print(f"    against {name:>13}: R_min / {name} = {Rs[i]/val:.4f}")
        print(f"  monotone in gamma? "
              f"{'yes' if np.all(np.diff(Rs) > 0) or np.all(np.diff(Rs) < 0) else 'NO'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
