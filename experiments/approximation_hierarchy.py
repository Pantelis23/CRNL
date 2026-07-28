"""Where each level of coarse-graining breaks on a restoration observable.

The observable is the one with a known exact answer: P(error) from a start biased
toward X by eps, with "decided" at |n_X - n_Y| >= theta * delta* * Omega. The exact
CME gives it as a splitting probability with no sampling error at all, so every
approximate level is measured against a reference rather than against each other.

    ODE        deterministic mass action
    CLE        Gaussian noise of variance = propensity (crnl/approximations)
    tau-leap   Poisson firings over a fixed window
    SSA        exact jumps -- the anchor, must agree with CME to sampling error
    CME        exact reference

PREDICTIONS, written before running:

  P1  The ODE fails CATEGORICALLY, not quantitatively. It reports p = 0 at every
      Omega and every eps > 0 where the truth is finite, and there is no
      refinement parameter that improves it. Not an inaccurate number -- the wrong
      kind of number.
  P2  The SSA agrees with the CME within sampling error at every cell. This is the
      anchor; if it fails, the harness is wrong and nothing else here means
      anything.
  P3  The CLE gets the EXPONENT wrong, and in a signed direction. It is the
      quadratic truncation of the jump Hamiltonian, so it agrees with the exact
      rate function only near the saddle -- which is exactly the limit FINDINGS 15
      verified to 0.1%. Away from it the exact barrier is larger (`ln r` versus
      `2(r-1)/(r+1)`), so **the CLE should overestimate the failure probability,
      by a factor that grows with the barrier** i.e. with eps and with Omega.
  P4  Tau-leaping interpolates: indistinguishable from the SSA below some tau, and
      systematically wrong above it. It should break at the BOUNDARY (the minority
      species approaching zero) before it breaks at the barrier.
  P5  Cost: the exactness is bought at O(Omega) events per unit time for the SSA,
      or O(Omega^2) memory for the CME. The cheap levels buy their speed
      specifically by getting the tail wrong -- which is the only part that matters
      for restoration.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.approximations import cle_run, tau_leap_run
from crnl.cme import first_passage
from crnl.deterministic import integrate
from crnl.networks.am_reversible import am_reversible, delta_star, reverse_pairing
from crnl.thermo import gillespie_instrumented
from crnl.vectorized import compile_network


def _setup(gamma: float, omega: int, eps_frac: float, theta: float):
    ds = delta_star(gamma)
    d0 = max(1, int(round(eps_frac * ds * omega)))
    nb = int(round(omega * gamma / (1.0 + gamma)))
    rest = omega - nb
    if (rest - d0) % 2:
        d0 -= 1
    nx, ny = (rest + d0) // 2, (rest - d0) // 2
    thr = max(2, int(round(theta * ds * omega)))
    return np.array([nx, ny, nb], dtype=np.int64), thr, ds


def p_cme(gamma, omega, eps_frac, theta):
    n0, thr, _ = _setup(gamma, omega, eps_frac, theta)
    net = am_reversible(gamma)
    fp = first_passage(net, int(omega), float(omega), n0,
                       lambda s, thr=thr: abs(int(s[0]) - int(s[1])) >= thr,
                       reverse_pairing(net))
    return (1.0 - float(fp["split"])) if fp["valid"] else float("nan")


def p_ode(gamma, omega, eps_frac, theta):
    n0, thr, _ = _setup(gamma, omega, eps_frac, theta)
    traj = integrate(am_reversible(gamma), n0 / omega, t_span=(0.0, 400.0))
    x, y = traj.final()[0], traj.final()[1]
    return 0.0 if x > y else 1.0


def _sample(runner, gamma, omega, eps_frac, theta, trials, seed, **kw):
    n0, thr, _ = _setup(gamma, omega, eps_frac, theta)
    comp = compile_network(am_reversible(gamma), float(omega))
    rng = np.random.default_rng(seed)
    stop = lambda n: abs(n[0] - n[1]) >= thr
    wrong = ok = 0
    retries = budget = 0
    for _ in range(trials):
        r = runner(comp, n0, rng, stop=stop, t_max=4000.0, **kw)
        retries += r.retries
        if r.hit_budget:
            budget += 1
            continue
        ok += 1
        if r.n_final[0] <= r.n_final[1]:
            wrong += 1
    if ok == 0:
        return float("nan"), 0, retries, budget
    return wrong / ok, ok, retries, budget


def p_ssa(gamma, omega, eps_frac, theta, trials, seed):
    """Exact SSA, stopped at FIRST PASSAGE -- the same observable as the CME.

    Must use a stop predicate. Running to a large t_max and reading the final
    state is a different quantity: by then the chain has relaxed into an attractor
    and may have crossed back, so it measures an occupancy rather than a splitting
    probability, and it would break the P2 anchor for a reason that has nothing to
    do with coarse-graining.
    """
    n0, thr, _ = _setup(gamma, omega, eps_frac, theta)
    net = am_reversible(gamma)
    comp = compile_network(net, float(omega))
    pair = reverse_pairing(net)
    rng = np.random.default_rng(seed)
    wrong = 0
    for _ in range(trials):
        r = gillespie_instrumented(
            comp, n0, rng, pair,
            stop=lambda n: abs(int(n[0]) - int(n[1])) >= thr,
            max_steps=20_000_000)
        wrong += int(r.n_final[0] <= r.n_final[1])
    return wrong / trials


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.30)
    ap.add_argument("--omegas", type=int, nargs="+", default=[40, 60, 80, 100])
    ap.add_argument("--eps-fracs", type=float, nargs="+", default=[0.20, 0.35])
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--dts", type=float, nargs="+", default=[0.02])
    ap.add_argument("--taus", type=float, nargs="+", default=[0.02, 0.10, 0.40])
    ap.add_argument("--trials", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=20260728)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/approximation_hierarchy.json"))
    args = ap.parse_args()

    t0 = time.time()
    rows = []
    for ef in args.eps_fracs:
        print(f"\n=== eps/delta* = {ef}, gamma = {args.gamma}, theta = {args.theta}")
        hdr = f"{'Omega':>6} {'CME (exact)':>13} {'ODE':>7} {'SSA':>12} {'SSA/CME':>8}"
        hdr += f" | {'CLE':>12} {'CLE/CME':>8}"
        for tau in args.taus:
            hdr += f" | {'tau=' + str(tau):>11} {'ratio':>7}"
        print(hdr)
        for om in args.omegas:
            pc = p_cme(args.gamma, om, ef, args.theta)
            po = p_ode(args.gamma, om, ef, args.theta)
            ps = p_ssa(args.gamma, om, ef, args.theta, args.trials, args.seed + om)
            line = (f"{om:>6} {pc:>13.3e} {po:>7.1f} {ps:>12.3e} "
                    f"{ps/pc if pc > 0 else float('nan'):>8.3f}")
            rec = {"gamma": args.gamma, "omega": om, "eps_frac": ef,
                   "p_cme": pc, "p_ode": po, "p_ssa": ps}
            pl, nok, rt, bd = _sample(cle_run, args.gamma, om, ef, args.theta,
                                      args.trials, args.seed + 1000 + om,
                                      dt=args.dts[0])
            line += f" | {pl:>12.3e} {pl/pc if pc > 0 else float('nan'):>8.3f}"
            rec.update({"p_cle": pl, "cle_retries": rt, "cle_budget": bd})
            for tau in args.taus:
                pt, nok2, rt2, bd2 = _sample(
                    tau_leap_run, args.gamma, om, ef, args.theta, args.trials,
                    args.seed + 2000 + om, tau=tau)
                line += f" | {pt:>11.3e} {pt/pc if pc > 0 else float('nan'):>7.3f}"
                rec[f"p_tau_{tau}"] = pt
                rec[f"tau_{tau}_retries"] = rt2
            print(line)
            rows.append(rec)

    # the exponent is the physics: fit -ln p against Omega for each level
    print("\n=== error exponent c = d(-ln p)/dOmega, by level")
    print(f"{'eps/d*':>7} {'CME':>9} {'SSA':>9} {'CLE':>9} " +
          " ".join(f"{'tau=' + str(t):>9}" for t in args.taus))
    for ef in args.eps_fracs:
        sel = [r for r in rows if r["eps_frac"] == ef]
        om = np.array([r["omega"] for r in sel], dtype=float)
        cells = []
        for key in (["p_cme", "p_ssa", "p_cle"]
                    + [f"p_tau_{t}" for t in args.taus]):
            p = np.array([r[key] for r in sel])
            good = np.isfinite(p) & (p > 0)
            cells.append(np.polyfit(om[good], -np.log(p[good]), 1)[0]
                         if good.sum() >= 3 else float("nan"))
        print(f"{ef:>7.2f} " + " ".join(f"{c:>9.5f}" for c in cells))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
