"""T-COST-o: two networks with IDENTICAL mass-action ODEs have different reliability

§82 left one suspect standing: A is a functional of the whole field and of no scalar summary of
it. The kill test as opened was to match Delta and the saddle placement and vary the curvature.
**There is a far sharper version, and it kills more than the suspect.**

A = -int ln(mu/lam) dx depends on lam and mu SEPARATELY. The mass-action ODE depends only on
their DIFFERENCE, f = lam - mu. So any reaction pair adding the same function to both leaves the
deterministic dynamics **exactly** unchanged -- every fixed point, Delta, the saddle, the whole
vector field -- while moving A. Mass action supplies such pairs directly: a birth and a death
whose REACTANT COMPLEXES are identical have identical propensities.

    Schloegl:   X -> 2X   and   X -> 0     at equal rate c.  Both propensities are c*n.
    AM:         X + Y -> 2X  and  X + Y -> 2Y  at equal rate c.  Both are c*n_X*n_Y/Omega.

Each pair is a futile cycle: it injects noise and does no deterministic work. In the limit of
large c, mu/lam -> 1 and **the barrier vanishes while the ODE never moves at all.**

**WHY THIS MATTERS BEYOND T-COST-o.** This project's entire method is the gap between the ODE
and the CME. If two networks share an ODE exactly and differ in A, then **the ODE does not
determine reliability, and does not determine the maximum composition depth** -- the quantity
§72-§81 spent nine sections on. §78/§81's "no master equation needed" survives (A is still
quadrature), but its gloss "a property of the deterministic field" was too loose: the Lyapunov
route needs D = S diag(a) S^T, which is the propensities and not the drift, and so does A.

**SCOPE, stated before running.** Both neutral pairs are IRREVERSIBLE as written, so their
dissipation is formally infinite. That is not a defect of the test -- it sharpens §82: here the
element pays an UNBOUNDED thermodynamic cost and gets WORSE reliability for it. Drive and
reliability move in opposite directions in the same network.

PREDICTIONS, written before running.

  P1  GATE, in two parts, and the section is void without both.
      (a) The ODE must be unchanged to machine precision -- f(x) identical across c at many x,
          and the three roots identical. If the pair perturbs the field at all, everything below
          is measuring a landscape change and nothing else.
      (b) RULE 16, ABSOLUTE. The quadrature A must match the EXACT mean first-passage time of
          the modified chain, for each c. A quadrature that drifts with c while the true barrier
          does not would produce this entire result out of nothing.
  P2  **THE TEST. A falls with c, monotonically, toward 0.** Because mu/lam -> 1 pointwise as the
      common term dominates. If A were instead constant in c, A would be a functional of the
      drift alone and the ODE would determine reliability after all.
  P3  **THE CONSEQUENCE, in the project's own units.** Report ln D_max = A*Omega - ln(t/c*)
      (§81 P5) across c at fixed Omega. **Predicted: orders of magnitude**, between elements
      that no ODE measurement can tell apart.
  P4  **DOES THE OTHER EXPONENT AGREE?** eta = Delta^2/(2V) with V = (lam+mu)/(2|f'|) at the
      rail. The common term raises lam+mu and leaves f' alone, so **eta must fall too**. If A
      fell and eta rose, one of the two is being computed wrongly and P2 would not be safe.
  P5  **RULE 9, SECOND SUBSTRATE.** AM with its own neutral pair, measured by the shallow-barrier
      first-passage instrument §82 established. **Predicted: the same collapse.** A conservative
      two-species element and an open one-species element share no chemistry, and rule 9 has
      failed three times on a single axis.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import scipy.sparse.linalg as spla
from scipy.integrate import quad
from scipy.special import logsumexp

from crnl.cme import enumerate_states, generator
from crnl.networks.am_reversible import am_reversible, delta_star
from crnl.reactions import Reaction, ReactionNetwork
from experiments.cascade_schlogl import schlogl_consts
from experiments.depth_is_error import c_star

# ------------------------------------------------------- Schloegl with a neutral futile pair


def lam_mu(x, r1, r2, r3, c):
    """Concentration-scaled birth and death rates. The c terms are IDENTICAL in both."""
    k1a, k1r, k2b, k2r = schlogl_consts(r1, r2, r3)
    return k1a * x ** 2 + k2b + c * x, k1r * x ** 3 + k2r * x + c * x


def drift(x, r1, r2, r3, c):
    lam, mu = lam_mu(x, r1, r2, r3, c)
    return lam - mu


def A_quad(r1, r2, r3, c):
    val, err = quad(lambda x: np.log(lam_mu(x, r1, r2, r3, c)[1]
                                     / lam_mu(x, r1, r2, r3, c)[0]),
                    r2, r3, limit=200)
    return -float(val), float(err)


def chain_rates(omega, r1, r2, r3, c, cap):
    """Stochastic propensities. X -> 2X and X -> 0 both contribute c*n, exactly."""
    k1a, k1r, k2b, k2r = schlogl_consts(r1, r2, r3)
    n = np.arange(0, cap + 1, dtype=float)
    lam = k1a * n * (n - 1.0) / omega + k2b * omega + c * n
    mu = k1r * n * (n - 1.0) * (n - 2.0) / omega ** 2 + k2r * n + c * n
    return np.maximum(lam, 0.0), np.maximum(mu, 0.0)


def ln_mfpt(omega, r1, r2, r3, c, cap_mult=1.8):
    """§80's exact log-space first passage, rail -> saddle, on the MODIFIED chain."""
    cap = int(np.ceil(cap_mult * r3 * omega))
    lam, mu = chain_rates(omega, r1, r2, r3, c, cap)
    a, n0 = int(round(r2 * omega)), int(round(r3 * omega))
    lp = np.full(cap + 1, -np.inf)
    lp[a] = 0.0
    acc = 0.0
    for k in range(a + 1, cap + 1):
        if lam[k - 1] <= 0 or mu[k] <= 0:
            break
        acc += np.log(lam[k - 1]) - np.log(mu[k])
        lp[k] = acc
    terms = []
    for k in range(a, n0):
        if lam[k] <= 0 or not np.isfinite(lp[k]):
            continue
        tail = lp[k + 1:cap + 1]
        tail = tail[np.isfinite(tail)]
        if tail.size:
            terms.append(-np.log(lam[k]) - lp[k] + logsumexp(tail))
    return float(logsumexp(terms)) if terms else np.nan


def A_mfpt(r1, r2, r3, c, oms=(3200, 6400)):
    a, b = ln_mfpt(oms[0], r1, r2, r3, c), ln_mfpt(oms[1], r1, r2, r3, c)
    return (b - a) / (oms[1] - oms[0])


def eta_lna(r1, r2, r3, c):
    """eta = Delta^2/(2V), V = (lam+mu)/(2|f'|) at the high rail. f' is c-INDEPENDENT."""
    k1a, k1r, k2b, k2r = schlogl_consts(r1, r2, r3)
    fp = 2 * k1a * r3 - 3 * k1r * r3 ** 2 - k2r          # the c terms cancel in f = lam - mu
    lam, mu = lam_mu(r3, r1, r2, r3, c)
    V = (lam + mu) / (2 * abs(fp))
    return ((r3 - r1) / 2) ** 2 / (2 * V), V


# ------------------------------------------------------- AM with its own neutral pair


def am_neutral(gamma, c, k=1.0):
    """Reversible AM plus X+Y->2X and X+Y->2Y at equal rate c.

    Both have reactant complex X+Y, so both have propensity c*n_X*n_Y/Omega; one raises the
    lead by 2 and the other lowers it by 2, so the drift contribution cancels EXACTLY while the
    lead's diffusion gains 8*c*x*y. Conservation X+Y+B is untouched.
    """
    net = am_reversible(gamma, k)
    rs = list(net.reactions)
    if c > 0:
        rs += [Reaction({"X": 1, "Y": 1}, {"X": 2}, c, name="n+:X+Y->2X"),
               Reaction({"X": 1, "Y": 1}, {"Y": 2}, c, name="n-:X+Y->2Y")]
    return ReactionNetwork(species=["X", "Y", "B"], reactions=rs,
                           name=f"am-neutral-g{gamma}-c{c}")


def am_ln_mfpt(gamma, c, omega):
    """§82's instrument: exact mean first passage from the rail to the diagonal n_X = n_Y."""
    net = am_neutral(gamma, c)
    states, index = enumerate_states(3, int(omega))
    d = states[:, 0].astype(np.int64) - states[:, 1].astype(np.int64)
    trans = np.where(d > 0)[0]
    Qtt = generator(net, int(omega), float(omega))[trans][:, trans].tocsc()
    T = spla.spsolve(Qtt, -np.ones(len(trans)))
    ds, b = float(delta_star(gamma)), gamma / (1.0 + gamma)
    nx = int(round((1 - b + ds) / 2 * omega))
    ny = int(round((1 - b - ds) / 2 * omega))
    if nx <= ny:
        return None
    r = int(np.where(trans == index[(nx, ny, int(omega) - nx - ny)])[0][0])
    t = float(T[r])
    if not np.isfinite(t) or t <= 0 or float(T.min()) <= 0:
        return None                       # §80/§82: a negative mean time is a broken instrument
    return float(np.log(t))


def am_A(gamma, c, omegas=(120, 180, 240, 300), rel=0.02):
    pts = [(om, am_ln_mfpt(gamma, c, om)) for om in omegas]
    pts = [(om, v) for om, v in pts if v is not None]
    loc = [(pts[i + 1][1] - pts[i][1]) / (pts[i + 1][0] - pts[i][0])
           for i in range(len(pts) - 1)]
    if len(loc) < 2 or abs(loc[-1] - loc[-2]) / abs(loc[-1]) > rel:
        return None, loc
    return loc[-1], loc


# ------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rails", type=float, nargs=3, default=[0.1, 1.0, 1.9])
    ap.add_argument("--cs", type=float, nargs="+", default=[0.0, 0.5, 1.0, 2.0, 5.0, 20.0])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/ode_does_not_determine_it.json"))
    args = ap.parse_args()
    r1, r2, r3 = args.rails
    out = {}

    print("=== P1(a) GATE: is the mass-action ODE really untouched?")
    xs = np.linspace(0.02, 2.5, 41)
    f0 = drift(xs, r1, r2, r3, 0.0)
    worst = 0.0
    for c in args.cs:
        worst = max(worst, float(np.abs(drift(xs, r1, r2, r3, c) - f0).max()))
    scale = float(np.abs(f0).max())
    print(f"  max |f_c(x) - f_0(x)| over 41 points and {len(args.cs)} values of c: {worst:.3e}")
    print(f"  (field scale max|f_0| = {scale:.4f}, so the relative perturbation is"
          f" {worst/scale:.1e})")
    ok_a = worst < 1e-12 * max(scale, 1.0)

    # ... and the SAME gate for AM, which was argued analytically above and must not be
    # taken on that argument: the AM half of the section means nothing without it.
    rng = np.random.default_rng(0)
    S0 = am_neutral(0.35, 0.0).stoichiometry_matrix()
    am_worst, am_scale = 0.0, 0.0
    for _ in range(200):
        x = rng.dirichlet([1.0, 1.0, 1.0])
        f0am = S0 @ am_neutral(0.35, 0.0).fluxes(x)
        am_scale = max(am_scale, float(np.abs(f0am).max()))
        for c in (0.25, 0.5, 1.0):
            net = am_neutral(0.35, c)
            fc = net.stoichiometry_matrix() @ net.fluxes(x)
            am_worst = max(am_worst, float(np.abs(fc - f0am).max()))
    print(f"  AM, 200 random simplex points x 3 values of c: max |f_c - f_0| = {am_worst:.3e}"
          f"  (field scale {am_scale:.4f})")
    ok_a = ok_a and am_worst < 1e-12 * max(am_scale, 1.0)
    print(f"  -> P1(a) {'HOLDS on BOTH substrates: the deterministic dynamics is identical' if ok_a else 'FAILS -- a pair is not neutral and nothing below means anything'}")
    assert ok_a, "a neutral pair perturbed the ODE; the experiment is void"

    print("\n=== P1(b) GATE (rule 16, ABSOLUTE): quadrature against the exact first passage")
    print(f"{'c':>7}{'A quadrature':>15}{'A exact MFPT':>15}{'ratio':>9}")
    rows = []
    for c in args.cs:
        Aq, err = A_quad(r1, r2, r3, c)
        Am = A_mfpt(r1, r2, r3, c)
        eta, V = eta_lna(r1, r2, r3, c)
        rows.append({"c": c, "A_quad": Aq, "A_mfpt": Am, "ratio": Aq / Am,
                     "quad_err": err, "eta": eta, "V": V})
        print(f"{c:>7}{Aq:>15.6f}{Am:>15.6f}{Aq/Am:>9.4f}")
    wr = max(abs(r["ratio"] - 1) for r in rows)
    print(f"  worst |quadrature/exact - 1| = {100*wr:.3f}%")
    ok_b = wr < 0.02
    print(f"  -> P1(b) {'HOLDS: the quadrature tracks the true barrier at every c' if ok_b else 'FAILS'}")
    out["schlogl"] = rows

    print("\n=== P2/P4: A and eta against the futile rate, at a FIXED deterministic field")
    print(f"{'c':>7}{'A':>12}{'A/A(0)':>10}{'eta':>12}{'eta/eta(0)':>13}")
    A0, e0 = rows[0]["A_quad"], rows[0]["eta"]
    for r in rows:
        print(f"{r['c']:>7}{r['A_quad']:>12.6f}{r['A_quad']/A0:>10.4f}"
              f"{r['eta']:>12.6f}{r['eta']/e0:>13.4f}")
    As = [r["A_quad"] for r in rows]
    etas = [r["eta"] for r in rows]
    mono_A = all(As[i + 1] < As[i] for i in range(len(As) - 1))
    mono_e = all(etas[i + 1] < etas[i] for i in range(len(etas) - 1))
    print(f"  A falls by a factor of {A0/As[-1]:.1f} across the sweep;"
          f" eta by {e0/etas[-1]:.1f}")
    print(f"  -> P2 {'HOLDS: A is NOT a functional of the drift -- the ODE does not determine reliability' if mono_A and A0/As[-1] > 1.5 else 'FAILS: A is constant at fixed ODE'}")
    print(f"  -> P4 {'the two exponents agree in direction, as they must' if mono_e else 'the exponents DISAGREE in direction -- one of them is computed wrongly and P2 is not safe'}")

    print("\n=== P3: the consequence, in depth. ln D_max = A*Omega - ln(t/c*)")
    cs, t = c_star(), 2.0
    print(f"{'c':>7}{'Omega':>8}{'ln D_max':>12}{'D_max':>14}")
    for r in rows:
        for om in (6400,):
            lnD = r["A_quad"] * om - np.log(t / cs)
            print(f"{r['c']:>7}{om:>8}{lnD:>12.1f}{np.exp(min(lnD, 700)):>14.3e}")
    span = (rows[0]["A_quad"] - rows[-1]["A_quad"]) * 6400
    print(f"  -> at Omega = 6400 the depth ceiling spans e^{span:.0f} between elements that"
          f" NO ODE measurement can tell apart")

    print("\n=== P5 (rule 9): AM, a conservative two-species element, with its own neutral pair")
    print(f"{'gamma':>7}{'c':>7}{'A':>10}{'A/A(0)':>10}  local A series")
    am, base = [], {}
    for g in (0.35, 0.40):
        for c in (0.0, 0.25, 0.5, 1.0):
            A, loc = am_A(g, c)
            ser = ", ".join(f"{v:.4f}" for v in loc)
            if A is None:
                print(f"{g:>7}{c:>7}{'EXCLUDED':>10}{'':>10}  {ser}")
                continue
            base.setdefault(g, A)
            am.append({"gamma": g, "c": c, "A": A, "series": loc})
            print(f"{g:>7}{c:>7}{A:>10.5f}{A/base[g]:>10.4f}  {ser}")
    out["am"] = am
    ok5 = None
    for g in (0.35, 0.40):
        row = [r for r in am if r["gamma"] == g]
        if len(row) >= 3:
            fell = all(row[i + 1]["A"] < row[i]["A"] for i in range(len(row) - 1))
            ok5 = fell if ok5 is None else (ok5 and fell)
            print(f"  gamma={g}: A falls by a factor of {row[0]['A']/row[-1]['A']:.2f}"
                  f"  ({'monotone' if fell else 'NOT monotone'})")
    if ok5 is None:
        print("  -> P5 UNDECIDED: too few AM cells survived the convergence gate")
    else:
        print(f"  -> P5 {'HOLDS on a second substrate: a conservative element collapses the same way' if ok5 else 'FAILS: AM does not follow'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
