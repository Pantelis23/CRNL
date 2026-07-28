"""T11a's kill test: does FINDINGS 21's cliff survive a bigger alphabet?

§21 measured the coarse-graining ladder on AM (n = 2) and found a cliff: every
level keeping ANY noise recovers the restoration error exponent to 2-12%, and the
ODE, keeping none, is categorically wrong. THEORIES T11a names the obvious way for
that to be an AM property rather than a restoration property -- run it at n >= 3.

WHAT HAD TO CHANGE, and one of them was a real bug waiting to happen:

  * `cme.first_passage` hardcodes its favoured set as `n[0] > n[1]`, which is
    exactly right at n = 2 and silently WRONG above it -- in a 3-winner race a
    state where X3 has won can still satisfy n_X1 > n_X2 and would be scored as a
    success. `cme.splitting_probability` takes the predicate instead. It
    reproduces `first_passage` to 0.00e+00 on AM, where both are defined.
  * There is no `delta_star` for n > 2, so the landscape width is measured: the
    ODE is integrated from a biased start and the champion-minus-best-rival
    separation at the attractor IS the scale. Bias and threshold are then the same
    fractions of it that §21 used, so the two sections are asking the same question
    of their respective networks.

PREDICTION, written before running: if the cliff is a property of restoration, the
CLE's exponent error at n = 3 stays within a few percent as it did at n = 2, and
the ODE still reports exactly zero. If the cliff is an AM property, the CLE error
grows -- the breaking mode at n = 3 is (0.816, -0.408, -0.408, 0) rather than
(0.707, -0.707, 0), so the noise the CLE has to get right is spread over more
directions and there is no reason its quadratic truncation should hold up equally.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.approximations import cle_run, tau_leap_run
from crnl.cme import splitting_probability
from crnl.deterministic import integrate
from crnl.networks.am_reversible import reverse_pairing
from crnl.networks.n_winner_reversible import (
    gamma_critical, n_winner_reversible, symmetric_state,
)
from crnl.thermo import gillespie_instrumented
from crnl.vectorized import compile_network


def landscape_width(n: int, gamma: float) -> float:
    """Champion minus best rival at the attractor: the n > 2 stand-in for delta*."""
    x, b = symmetric_state(n, gamma)
    x0 = np.full(n + 1, x)
    x0[n] = b
    x0[0] += 0.02
    x0[1:n] -= 0.02 / (n - 1)
    traj = integrate(n_winner_reversible(n, gamma), x0, t_span=(0.0, 800.0))
    f = traj.final()[:n]
    return float(f.max() - np.sort(f)[-2])


def setup(n, gamma, omega, eps_frac, theta, width):
    """Start state whose champion-minus-best-rival margin is EXACTLY the target.

    The obvious construction -- rivals at `(rest-m)//n`, champion taking the
    remainder -- lets the realised margin overshoot by up to n-1 counts, because
    the floor division dumps the remainder onto the champion. At Omega = 60, n = 3
    that is a 15% error in eps, and since it varies erratically with Omega it
    corrupts the exponent fit, which requires eps held FIXED. It showed up as a
    non-monotone CME error probability (0.115 -> 0.064 -> 0.092 -> 0.026), which
    is not a thing a first-passage probability can do.

    Here the rival maximum is pinned first and the remainder is taken off the
    rivals instead, so `max(rivals) = champion - m` holds exactly; the assert is
    the guard, not a comment.
    """
    _, b = symmetric_state(n, gamma)
    nb = int(round(b * omega))
    rest = omega - nb
    m = max(1, int(round(eps_frac * width * omega)))
    for dm in (0, 1, -1, 2, -2):
        mm = m + dm
        if mm < 1:
            continue
        R = -((-(rest - mm)) // n)              # ceil((rest-mm)/n)
        c = R + mm
        rivals = [R] * (n - 1)
        excess = sum(rivals) + c - rest
        if excess < 0 or excess > n - 2:        # must leave one rival at R
            continue
        for i in range(excess):
            rivals[-(i + 1)] -= 1
        n0 = np.array([c] + rivals + [nb], dtype=np.int64)
        if n0.sum() == omega and (n0 >= 0).all() and c - max(rivals) == mm:
            thr = max(2, int(round(theta * width * omega)))
            return n0, thr, float(mm) / omega
    raise ValueError(f"no exact-margin start at Omega={omega}, n={n}")


def _absorbing(n_committed, thr):
    def f(s):
        v = np.sort(np.asarray(s[:n_committed]))
        return (v[-1] - v[-2]) >= thr
    return f


def _favoured(n_committed):
    def f(s):
        v = np.asarray(s[:n_committed])
        return int(np.argmax(v)) == 0
    return f


def p_cme(n, gamma, omega, n0, thr):
    r = splitting_probability(n_winner_reversible(n, gamma), int(omega),
                              float(omega), n0, _absorbing(n, thr), _favoured(n))
    return (1.0 - float(r["split"])) if r["valid"] else float("nan")


def p_ode(n, gamma, omega, n0):
    traj = integrate(n_winner_reversible(n, gamma), n0 / omega, t_span=(0.0, 800.0))
    return 0.0 if int(np.argmax(traj.final()[:n])) == 0 else 1.0


def _sampled(runner, n, gamma, omega, n0, thr, trials, seed, **kw):
    comp = compile_network(n_winner_reversible(n, gamma), float(omega))
    rng = np.random.default_rng(seed)
    stop = _absorbing(n, thr)
    wrong = ok = 0
    for _ in range(trials):
        r = runner(comp, n0, rng, stop=stop, t_max=6000.0, **kw)
        if r.hit_budget:
            continue
        ok += 1
        wrong += int(np.argmax(r.n_final[:n]) != 0)
    return (wrong / ok if ok else float("nan")), ok


def p_ssa(n, gamma, omega, n0, thr, trials, seed):
    net = n_winner_reversible(n, gamma)
    comp = compile_network(net, float(omega))
    pair = reverse_pairing(net)   # derived, not assumed; n-winner pairs resolve
    rng = np.random.default_rng(seed)
    stop = _absorbing(n, thr)
    wrong = 0
    for _ in range(trials):
        r = gillespie_instrumented(comp, n0, rng, pair, stop=stop,
                                   max_steps=20_000_000)
        wrong += int(np.argmax(r.n_final[:n]) != 0)
    return wrong / trials


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--gamma-frac", type=float, default=0.60,
                    help="gamma as a fraction of gamma_c(n); 0.60 matches §21's AM run")
    ap.add_argument("--omegas", type=int, nargs="+", default=[30, 45, 60, 80])
    ap.add_argument("--eps-fracs", type=float, nargs="+", default=[0.25, 0.40])
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--taus", type=float, nargs="+", default=[0.05, 0.30])
    ap.add_argument("--trials", type=int, default=6000)
    ap.add_argument("--seed", type=int, default=20260728)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/approximation_hierarchy_nwinner.json"))
    args = ap.parse_args()

    n = args.n
    gamma = args.gamma_frac * gamma_critical(n)
    width = landscape_width(n, gamma)
    t0 = time.time()
    print(f"n={n}  gamma={gamma:.5f} (= {args.gamma_frac} x gamma_c={gamma_critical(n):.5f})"
          f"  landscape width = {width:.4f}")
    rows = []
    for ef in args.eps_fracs:
        print(f"\n=== eps/width = {ef}")
        hdr = (f"{'Omega':>6} {'eps_real':>7} {'CME (exact)':>13} {'ODE':>6} {'SSA':>12} {'/CME':>7}"
               f" | {'CLE':>12} {'/CME':>7}")
        for tau in args.taus:
            hdr += f" | {'tau=' + str(tau):>11} {'/CME':>7}"
        print(hdr)
        for om in args.omegas:
            n0, thr, realised = setup(n, gamma, om, ef, args.theta, width)
            pc = p_cme(n, gamma, om, n0, thr)
            po = p_ode(n, gamma, om, n0)
            ps = p_ssa(n, gamma, om, n0, thr, args.trials, args.seed + om)
            line = (f"{om:>6} {realised:>7.4f} {pc:>13.3e} {po:>6.1f} {ps:>12.3e} "
                    f"{ps/pc if pc > 0 else float('nan'):>7.3f}")
            rec = {"n": n, "gamma": gamma, "omega": om, "eps_frac": ef,
                   "realised_margin": realised, "p_cme": pc, "p_ode": po,
                   "p_ssa": ps}
            pl, _ = _sampled(cle_run, n, gamma, om, n0, thr, args.trials,
                             args.seed + 1000 + om, dt=args.dt)
            line += f" | {pl:>12.3e} {pl/pc if pc > 0 else float('nan'):>7.3f}"
            rec["p_cle"] = pl
            for tau in args.taus:
                pt, _ = _sampled(tau_leap_run, n, gamma, om, n0, thr, args.trials,
                                 args.seed + 2000 + om, tau=tau)
                line += f" | {pt:>11.3e} {pt/pc if pc > 0 else float('nan'):>7.3f}"
                rec[f"p_tau_{tau}"] = pt
            print(line)
            rows.append(rec)

    # Fit against eps^2 * Omega, not Omega. The lattice cannot place the same
    # margin at every Omega (at Omega = 30 the realised eps is 22% off target),
    # and since c ~ kappa eps^2 a 5% drift in eps moves the exponent by 10% --
    # comparable to the effect being measured. Regressing on eps^2*Omega absorbs
    # it exactly and the slope IS kappa.
    print(f"\n=== kappa = d(-ln p)/d(eps^2 Omega), n = {n}")
    print(f"{'eps':>6} {'CME':>9} {'SSA':>9} {'CLE':>9} " +
          " ".join(f"{'tau=' + str(t):>9}" for t in args.taus))
    for ef in args.eps_fracs:
        sel = [r for r in rows if r["eps_frac"] == ef]
        om = np.array([r["omega"] * r["realised_margin"] ** 2 for r in sel],
                      dtype=float)
        cells = []
        for key in ["p_cme", "p_ssa", "p_cle"] + [f"p_tau_{t}" for t in args.taus]:
            p = np.array([r[key] for r in sel])
            g = np.isfinite(p) & (p > 0)
            cells.append(np.polyfit(om[g], -np.log(p[g]), 1)[0]
                         if g.sum() >= 3 else float("nan"))
        print(f"{ef:>6.2f} " + " ".join(f"{c:>9.5f}" for c in cells))
        base = cells[0]
        if np.isfinite(base) and base != 0:
            print(f"       ratios to exact: " + " ".join(
                f"{c/base:>9.4f}" for c in cells))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
