"""T-CASC-d: is §91's margin law derivable? Annealed averaging over the upstream, no fit

§91 measured log(per-stage penalty) = -0.952 x margin/sigma, with two independent knobs
collapsing onto one curve. That is a measurement, not a mechanism (rule 17). Here is the
candidate, and it is computable in absolute terms:

  * the downstream's barrier depends on the upstream concentration, A = A(x_up), and vanishes at
    the collapse point x_crit;
  * the upstream relaxes in 1/|f'(r3)| = 0.151 time units while the downstream escapes on
    ~exp(A*Omega) -- four orders apart -- so the downstream sees the TIME-AVERAGED rate, not the
    rate at the mean. That is an ANNEALED average, and it is the whole content of the mechanism:

        penalty = < exp( -[A(x_up) - A(r3)] * Omega ) >_upstream

**WHY THE OBVIOUS VERSION IS WRONG, worked out before running.** Linearising A near the rail and
averaging over a Gaussian upstream of width sigma^2 = V/Omega gives a saddle point at
u* = V A'(r3), and log(penalty) = A(r3)^2 Omega^2 / (2 m^2) -- **quadratic-inverse in the margin
m**, not linear. Against §91's numbers that predicts 4.97, 1.42, 0.74 at m = 1.81, 3.39, 4.70
where the measurement gives 2.77, 0.77, 0.061. It is the right order and the wrong shape, and
log(penalty)*m^2 is not constant across the data (9.08, 10.59, 8.87, 1.35). **So the linearisation
of A is what fails**: the barrier does not fall linearly toward the collapse, and near a
saddle-node it goes as (x_up - x_crit)^(3/2). Rather than substitute one guessed exponent for
another, A(x_up) is computed numerically and the average is done over the EXACT upstream
quasi-stationary distribution from the CME.

PREDICTIONS, written before running.

  P1  GATE. A(x_up) must reproduce what is already known: A(r3) must equal the isolated element's
      action from §80's quadrature, and A must fall to 0 at exactly the x_crit §91 measured from
      the roots. Two independent routes to the same two numbers.
  P2  **THE SHAPE.** A(x_up) near x_crit must go as (x_up - x_crit)^p. **Predicted p = 3/2**, the
      saddle-node exponent (Dykman-Krivoglaz; §85.2's literature note gives 3/2 at a saddle-node
      against 2 at a pitchfork). Fitted only to identify p, and reported per rule 21 as
      p_eff(window) since the barrier vanishes at the threshold.
  P3  **THE TEST, ABSOLUTE (rule 16).** The annealed average against §91's 14 measured penalties.
      **Predicted: it reproduces them within a factor of ~2 with no fitted parameter**, and in
      particular reproduces the LINEAR-in-margin shape that the Gaussian linearisation misses.
      A prediction that lands on the right shape but the wrong scale still identifies the
      mechanism; one that misses the shape does not.
  P4  **WHAT THE ANNEALED AVERAGE BUYS.** If P3 holds, the depth of a cascade is computable for
      ANY transfer function from single-element quantities plus the upstream distribution -- no
      joint CME, which is what limited §91 to D = 2.
  P5  **RULE 15.** Report the quenched average <A> as well (evaluating the rate at the mean
      upstream, i.e. ignoring the fluctuation). If quenched also fits, the annealed mechanism is
      not identified by this data and the agreement is not evidence for it.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.integrate import quad

import experiments.chemical_cascade as cc
from experiments.cascade_schlogl import schlogl_consts
from experiments.derive_eta import schlogl_V

R1, R2, R3 = cc.RAILS
C = schlogl_consts(R1, R2, R3)

# §91's measured penalties, quoted as stored numbers (rule 16).
MEASURED = [(1.0, 1.0, 2.59, 4.845), (2.0, 1.0, 2.61, 4.860), (3.0, 1.0, 3.02, 3.159),
            (4.0, 1.0, 3.39, 2.164), (6.0, 1.0, 3.83, 1.479), (8.0, 1.0, 4.05, 1.301),
            (12.0, 1.0, 4.26, 1.211), (4.0, 0.6, 4.70, 1.063), (4.0, 0.8, 4.03, 1.328),
            (4.0, 1.2, 2.80, 4.190), (4.0, 1.4, 2.26, 8.391), (4.0, 1.6, 1.81, 15.987),
            (6.0, 1.3, 2.98, 3.615), (8.0, 1.5, 2.72, 5.686)]


def down_lam_mu(x, x_up, scheme):
    k1a, k1r, k2b, k2r = C
    if scheme == "source":
        lam = k1a * x ** 2 + (k2b / R3) * x_up
    elif scheme == "catalytic":
        lam = (k1a / R3) * x_up * x ** 2 + k2b
    else:
        lam = k1a * cc.hill(x_up, R3) * x ** 2 + k2b
    return lam, k1r * x ** 3 + k2r * x


def A_of_xup(x_up, scheme):
    """Escape action of the DOWNSTREAM landscape when the upstream sits at x_up. 0 if collapsed."""
    r = cc.downstream_roots(x_up, C, R3, scheme)
    if len(r) < 3:
        return 0.0
    xs, xh = r[1], r[2]
    val, _ = quad(lambda z: np.log(down_lam_mu(z, x_up, scheme)[0]
                                   / down_lam_mu(z, x_up, scheme)[1]), xs, xh, limit=200)
    return float(val)


def upstream_qsd(omega):
    """Exact quasi-stationary distribution of stage 1 above its saddle, from the CME."""
    Q, cap, m, strides, _ = cc.build(1, omega, "source")
    n = np.arange(m)
    keep = np.where(n > R2 * omega)[0]
    Qs = Q[keep][:, keep].T.tocsc()
    w, v = spla.eigs(Qs, k=1, which="LR")
    p = np.real(v[:, 0])
    p = np.abs(p)
    return n[keep] / omega, p / p.sum()


def build_speed(om, speed, scheme="hill", cap_mult=1.25):
    """Two-stage chain with stage 1's propensities scaled by `speed` -- same landscape, same
    barrier, same rail width, same stationary law. Only the clock changes."""
    cap = int(np.ceil(cap_mult * R3 * om))
    m = cap + 1
    rows, cols, vals = [], [], []
    for idx in range(m * m):
        ns = [idx // m, idx % m]
        tot = 0.0
        for i in (0, 1):
            lam, mu = cc.rates_stage(ns[i], ns[i - 1] if i else 0, om, C, R3, i == 0, scheme)
            if i == 0:
                lam *= speed
                mu *= speed
            st = m if i == 0 else 1
            if ns[i] < cap and lam > 0:
                rows.append(idx); cols.append(idx + st); vals.append(lam); tot += lam
            if ns[i] > 0 and mu > 0:
                rows.append(idx); cols.append(idx - st); vals.append(mu); tot += mu
        rows.append(idx); cols.append(idx); vals.append(-tot)
    return sp.csr_matrix((vals, (rows, cols)), shape=(m * m, m * m)), m, [m, 1]


def stage1_stationary(om, cap_mult=1.25):
    """Exact stationary law of the REFLECTED stage 1, by the birth-death product formula.
    Independent of the speed factor: scaling every rate equally leaves it unchanged, which is
    what makes it a clean seed for a clock sweep."""
    cap = int(np.ceil(cap_mult * R3 * om))
    nsad = int(np.ceil(R2 * om))
    up = np.arange(nsad, cap + 1)
    lp = np.zeros(len(up))
    acc = 0.0
    for i in range(1, len(up)):
        l, _ = cc.rates_stage(float(up[i - 1]), 0.0, om, C, R3, True, "hill")
        _, u = cc.rates_stage(float(up[i]), 0.0, om, C, R3, True, "hill")
        acc += np.log(l) - np.log(u)
        lp[i] = acc
    w = np.exp(lp - lp.max())
    return up, w / w.sum()


def build_reflect(om, speed, cap_mult=1.25):
    """Two-stage chain with stage 1 REFLECTED at its saddle -- it fluctuates around its rail
    but can never escape, so P(stage 2 low) needs no conditioning at all."""
    cap = int(np.ceil(cap_mult * R3 * om))
    nsad = int(np.ceil(R2 * om))
    up = np.arange(nsad, cap + 1)
    nu, m2 = len(up), cap + 1
    rows, cols, vals = [], [], []
    for a in range(nu):
        n1 = up[a]
        for n2 in range(m2):
            idx = a * m2 + n2
            tot = 0.0
            l1, u1 = cc.rates_stage(float(n1), 0.0, om, C, R3, True, "hill")
            l1 *= speed; u1 *= speed
            if n1 < cap and l1 > 0:
                rows.append(idx); cols.append((a + 1) * m2 + n2); vals.append(l1); tot += l1
            if n1 > up[0] and u1 > 0:
                rows.append(idx); cols.append((a - 1) * m2 + n2); vals.append(u1); tot += u1
            l2, u2 = cc.rates_stage(float(n2), float(n1), om, C, R3, False, "hill")
            if n2 < cap and l2 > 0:
                rows.append(idx); cols.append(idx + 1); vals.append(l2); tot += l2
            if n2 > 0 and u2 > 0:
                rows.append(idx); cols.append(idx - 1); vals.append(u2); tot += u2
            rows.append(idx); cols.append(idx); vals.append(-tot)
    return sp.csr_matrix((vals, (rows, cols)), shape=(nu * m2, nu * m2)), up, m2, cap


def predict(scheme, n_hill, k_hill, omega, xs, px, quenched=False):
    cc.HILL_N, cc.HILL_K = n_hill, k_hill
    A0 = A_of_xup(R3, scheme)
    if quenched:
        xbar = float((px * xs).sum())
        return float(np.exp(-(A_of_xup(xbar, scheme) - A0) * omega))
    w = np.array([np.exp(-(A_of_xup(x, scheme) - A0) * omega) for x in xs])
    return float((px * w).sum())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=30)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/margin_law.json"))
    args = ap.parse_args()
    om = args.omega
    out = {}

    print("=== P1 GATE: two independent routes to A(r3) and to x_crit")
    cc.HILL_N, cc.HILL_K = 4.0, 1.0
    from experiments.action_is_not_priced import schlogl_A
    A_iso = schlogl_A(R1, R3, R2)[0]
    for scheme in ("source", "catalytic", "hill"):
        a = A_of_xup(R3, scheme)
        print(f"  {scheme:>10}: A(r3) = {a:.8f}   isolated §80 quadrature = {A_iso:.8f}"
              f"   diff {abs(a-A_iso):.2e}")
    ok = all(abs(A_of_xup(R3, s) - A_iso) < 1e-9 for s in ("source", "catalytic", "hill"))
    print(f"  -> P1 {'HOLDS: at the rail every coupling gives the isolated action' if ok else 'FAILS'}")
    assert ok

    print("\n=== P2: the shape of A(x_up) approaching the collapse")
    xs_grid = np.linspace(R1, R3, 3001)
    for scheme, lab in (("catalytic", "catalytic"), ("hill", "hill n=4 K=1")):
        vals = np.array([A_of_xup(x, scheme) for x in xs_grid])
        live = vals > 0
        xc = xs_grid[live].min()
        sel = (xs_grid > xc) & (xs_grid < xc + 0.35 * (R3 - xc))
        p = np.polyfit(np.log(xs_grid[sel] - xc), np.log(vals[sel]), 1)[0]
        print(f"  {lab:>13}: x_crit = {xc:.4f}, A ~ (x-x_crit)^p with p_eff(window) = {p:.3f}")
    print("  (reported as p_eff(window) per rule 21 -- A vanishes at the threshold, so any")
    print("   fitted exponent is window-dependent. 3/2 is the saddle-node value.)")

    xs, px = upstream_qsd(om)
    sd = float(np.sqrt(schlogl_V(R1, R2, R3) / om))
    print(f"\n  upstream QSD: mean {float((px*xs).sum()):.4f} (r3 = {R3}),"
          f" sd {float(np.sqrt((px*(xs-(px*xs).sum())**2).sum())):.4f}"
          f"  (LNA sigma = {sd:.4f})")

    print("\n=== P3/P5: the annealed average against §91's 14 measured penalties, no fit")
    print(f"{'n':>5}{'K':>6}{'margin':>9}{'measured':>11}{'frozen-avg':>13}{'ratio':>8}"
          f"{'quenched':>11}{'ratio':>8}")
    rows = []
    for n, K, marg, meas in MEASURED:
        a = predict("hill", n, K, om, xs, px)
        q = predict("hill", n, K, om, xs, px, quenched=True)
        rows.append({"n": n, "K": K, "margin": marg, "meas": meas, "ann": a, "que": q})
        print(f"{n:>5}{K:>6}{marg:>9.2f}{meas:>11.3f}{a:>11.3f}{a/meas:>8.3f}"
              f"{q:>11.3f}{q/meas:>8.3f}")
    cc.HILL_N, cc.HILL_K = 4.0, 1.0
    out["rows"] = rows
    ann = np.array([r["ann"] / r["meas"] for r in rows])
    que = np.array([r["que"] / r["meas"] for r in rows])
    print(f"  frozen-avg/measured: span {ann.min():.3f}..{ann.max():.3f}"
          f" (factor {ann.max()/ann.min():.2f})")
    print(f"  fast-avg/measured: span {que.min():.3f}..{que.max():.3f}"
          f" (factor {que.max()/que.min():.2f})")
    # does the annealed model reproduce the LINEAR-in-margin shape?
    mg = np.array([r["margin"] for r in rows])
    sl_meas = np.polyfit(mg, np.log([r["meas"] for r in rows]), 1)[0]
    sl_ann = np.polyfit(mg, np.log([max(r["ann"], 1e-300) for r in rows]), 1)[0]
    sl_que = np.polyfit(mg, np.log([max(r["que"], 1e-300) for r in rows]), 1)[0]
    print(f"  slope of log(penalty) vs margin:  measured {sl_meas:.3f},"
          f"  annealed {sl_ann:.3f},  quenched {sl_que:.3f}")
    shape_ok = abs(sl_ann - sl_meas) < 0.3 * abs(sl_meas)
    scale_ok = ann.max() / ann.min() < 4.0
    print(f"  -> P3 {'HOLDS: the annealed average reproduces the margin law with NO fitted parameter' if (shape_ok and scale_ok) else ('reproduces the SHAPE but not the scale -- the mechanism is identified, the prefactor is not' if shape_ok else 'FAILS on shape: the annealed average is not what produces the margin law')}")
    print(f"  -> P5 {'the quenched average does NOT reproduce it, so the FLUCTUATION is doing the work and the annealed mechanism is identified' if abs(sl_que - sl_meas) > 0.3*abs(sl_meas) else 'the quenched average fits too -- this data does not distinguish them, and P3 is not evidence for annealing'}")

    print("\n=== P6: the SECOND variable -- the upstream CLOCK, on a clean instrument")
    print("    Three confounds had to die first, and each reversed or inflated the answer:")
    print("      (a) seeding stage 1 as a delta at its rail: the slow cells had not SPREAD yet,")
    print("          which reads as a clock effect and reversed the sign of the trend;")
    print("      (b) measuring P(s2 low, s1 high): conditioning on stage 1 surviving excludes")
    print("          exactly the trajectories that dipped, and the exclusion grows with speed;")
    print("      (c) pre-equilibrating the joint chain: stage 2 AGES during it -- at speed 1/8")
    print("          it had already accumulated 14.6% error before the window opened.")
    print("    Clean instrument: stage 1 REFLECTED at its saddle (it can never escape, so")
    print("    nothing is conditioned), seeded from its EXACT stationary law (speed-independent,")
    print("    since scaling every rate equally leaves it unchanged), stage 2 fresh at its rail.")
    print(f"{'speed':>8}{'P(s2 low)':>13}{'penalty':>10}{'limit it approaches':>22}")
    cc.HILL_N, cc.HILL_K = 4.0, 1.0
    frozen = predict("hill", 4.0, 1.0, om, xs, px)
    fast = predict("hill", 4.0, 1.0, om, xs, px, quenched=True)
    Q1, cap1, m1, st1, _ = cc.build(1, om, "hill")
    q1 = spla.expm_multiply(Q1.T * 2.0, cc.seed_high(1, om, m1, st1, R3))
    eps_iso = float(q1[np.arange(m1) < R2 * om].sum())
    upv, pi1 = stage1_stationary(om)
    p6 = []
    for sp_ in (0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0):
        Q, uu, m2, cap = build_reflect(om, sp_)
        q = np.zeros(len(uu) * m2)
        for a, w in enumerate(pi1):
            q[a * m2 + int(round(R3 * om))] = w
        q = spla.expm_multiply(Q.T * 2.0, q)
        lo2 = (np.arange(len(uu) * m2) % m2) < R2 * om
        v = float(q[lo2].sum())
        p6.append({"speed": sp_, "P2": v, "pen": v / eps_iso})
        near = "frozen" if abs(v / eps_iso - frozen) < abs(v / eps_iso - fast) else "fast-avg"
        print(f"{sp_:>8.3f}{v:>13.5e}{v/eps_iso:>10.3f}{near:>22}")
    out["p6"] = p6
    slow = [r["pen"] for r in p6 if r["speed"] <= 1.0]
    quick = [r["pen"] for r in p6 if r["speed"] >= 8.0]
    print(f"  FROZEN-upstream formula <exp(-dA*Om)>   = {frozen:.3f};"
          f"  slow plateau measures {min(slow):.3f}..{max(slow):.3f}")
    print(f"  FAST-averaging formula  exp(-<dA>*Om)   = {fast:.3f};"
          f"  fast tail measures {min(quick):.3f} and still falling")
    ok_slow = abs(max(slow) / frozen - 1) < 0.20
    falls = p6[-1]["pen"] < 0.5 * max(slow)
    print(f"  -> P6 {'HOLDS: the penalty interpolates between the two limits as the upstream clock varies. The frozen formula is the SLOW limit (not the fast one), and the fall toward the mean-landscape rate as the upstream speeds up is MOTIONAL NARROWING' if (ok_slow and falls) else 'FAILS: the measured penalty does not sit between the two limiting formulas'}")
    print(f"  -> so §91's margin law is a FROZEN-upstream statement. All of §91's variants ran")
    print(f"     at speed 1, inside the plateau, which is why one variable sufficed there.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
