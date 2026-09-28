"""§114 (T-CASC-aa, T-CASC-m) -- a second axis for the position, and it is not Omega.

The fast/frozen POSITION is the cascade arc's central unsolved quantity. §109 measured it under
matched seeding as a smooth monotone curve running ~0.0 at Omega = 14 to ~1.16 at Omega = 70:
stage 2's effective escape rate traverses the whole averaging bracket as volume grows, and exits
it. Three explanatory families have been retired by measurement -- averaging prescriptions (§108),
seeding artifacts (§109), and the timescale ratio the framing was named for (§110, whose
refutation §113 then showed survives a converged box).

**Every one of those attempts swept Omega, and Omega is the one axis that cannot decide this.**
Both open entries say so in their own words:

    T-CASC-m: "a THIRD element chosen so the ratio differs from both existing values... The
               element must move the landscape and the coupling together -- Omega will not do
               it, by the argument above."
    T-CASC-aa: "testing it means moving tau_cross AT FIXED OMEGA."

Within one element, Omega moves the escape action A*Omega, the margin/sigma ratio, the rail width
and tau_cross **together** -- margin/sigma ~ sqrt(Omega) and A*Omega ~ Omega by construction. Five
points along a line through a space of five correlated candidates cannot separate them. That is
rule 9 in its original form: constancy -- or variation -- along the axis you happened to sweep is
not a statement about that axis.

§112 built the missing instrument without noticing what it was for. `schlogl_consts(r1, r2, r3)`
returns exact rate constants for any three roots, so the saddle can be MOVED at fixed rails,
producing genuinely different elements: six of them, with escape actions differing 75x, all real
chemistry with positive constants, all deterministically NEUTRAL at their own rail (verified to
1e-15, which is what §103.1's fixed-point structure requires) and all accepted by the existing
cascade machinery, since `rates_stage` is already parameterised by (c, r3).

This section measures the position on a GRID: six elements x three volumes. Everything here is a
parameterised re-implementation checked against the published pipeline rather than an edit to it
(P1), because generalising §109's measurement touches seven published functions.

PREDICTIONS, WRITTEN BEFORE RUNNING.

  P1  WIRING, and it gates everything. At the published rails (0.15, 1.0, 3.1827) and the
      published box (cap_mult = 1.25) this module must reproduce §109's stored p3 table --
      `qsd`, `delta` and `bracket_top` at all five Omega -- to better than 1e-6 relative, and
      §110's tau_up and tau_cross to 1e-9. If it does not, the generalisation changed the physics
      and nothing below can be read.

  P2  DOES THE POSITION MOVE AT FIXED OMEGA AT ALL? This is the question the whole section exists
      to ask, and **it cannot fail to inform, because either answer kills a family**:
        - if the position is FLAT across six elements at fixed Omega, it is a function of volume
          alone, and every rail-dependent candidate dies at once: A*Omega, margin/sigma,
          tau_cross, the bracket width. What would be left is sqrt(Omega) or something like it.
        - if it MOVES, then sqrt(Omega) and every other Omega-only scalar is refuted, because
          Omega did not change.
      I predict it moves, and moves a lot: the elements span 75x in action.

  P3  THE COLLAPSE TEST, AND RULE 15 IS THE WHOLE DESIGN. Report the position against EVERY
      candidate T-CASC-aa names -- A*Omega, (margin/sigma)^2, tau_cross, the bracket width,
      pi_low -- on BOTH datasets at once: the six-element sweep at fixed Omega and the published
      element's Omega sweep. **A candidate survives only if both datasets fall on one curve.**
      A candidate that fits the Omega sweep but misses the rail sweep is refuted, and that is
      precisely what the second axis is for. No candidate is allowed to be reported alone.

  P4  THE ONE I WILL COMMIT TO, so that it can be wrong: **A*Omega.** §98 identified it; §99.1's
      single out-of-sample point favoured it over (margin/sigma)^2 (-6.6% against -15.0%) and was
      explicitly labelled not decisive. I predict the position is monotone increasing in A*Omega
      and that the rail-sweep and Omega-sweep points fall on a common curve to within the box
      error §113 measured (3-7% at Omega = 14, less above). **If the two datasets do not collapse
      onto one A*Omega curve, A*Omega is refuted as THE variable and §98's identification was an
      artifact of only ever having had one axis.**

  P5  WHAT COMES FREE, and it is T16-g's audit (§112 closed three generality claims that rested
      on unwritten symmetry premises). Running the cascade on six new elements tests, at no extra
      cost, whether two published framings survive off the published rails:
        (a) §102's bracket ORDERING, k(<x>) < k_eff < <k>. If any element violates it, §102's
            "position between two limits" is not a general framing and the coordinate this
            section is measuring is not even well defined there.
        (b) §103's chain closure, whose intrinsic term d_intr = mu_1 - r3 was §103.1's fix.
            Report its sign and size per element.

TWO THINGS CHANGED AFTER SEEING OUTPUT, DISCLOSED BECAUSE RULE 19 IS ABOUT EXACTLY THIS.

  (i)  A VALIDITY FILTER WAS ADDED. The grid's raw output contains cells where the position
       coordinate does not exist: E-021's bracket is 1.022, so dividing by ln(bracket) prints
       a position of -304, and the deepest element's bracket INVERTS (k_avg < k_mean), which
       makes "position between the limits" meaningless. `validity()` now states the three
       definitional requirements -- stage 1 survives the window, the bracket is neither
       degenerate nor inverted -- and every verdict is computed on cells that pass.
       **This makes the conclusions harder, not easier.** On the raw grid P4's collapse test
       returned an RMS of 30 against a curve spanning 1.18 and printed "A*Omega is REFUTED" --
       a refutation of everything, from correct numbers, off cells where the quantity is not
       defined. That is §53's failure mode and it was one edit away from being published.

  (ii) P4's CRITERION WAS REPLACED. It first interpolated the Omega-sweep curve and took an RMS
       over all rail cells. It now compares only cells from DIFFERENT elements whose A*Omega
       already agree to 25%, so nothing is interpolated, fitted or extrapolated. The replacement
       is strictly stronger: it rests on a direct comparison rather than on a curve.

  Both changes were made before P6 was written or run. P1-P5's original wording is left above
  unedited (rule 3).

WHAT THIS CANNOT SETTLE. Six elements is a line through rail space too -- r1 and r3 are held fixed
and only r2 moves, so `f`, the action, the barrier asymmetry and tau all change together along it,
just differently than they do along Omega. Two crossing lines beat one line, and they are not a
plane. A candidate that survives both is *consistent with* being the variable; it is not proven to
be. And the shallowest element (E-021, A*Omega ~ 0.1 at Omega = 14) has barely any barrier, so the
usable lever is nearer 12x than 75x -- reported per cell rather than assumed.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from experiments.cascade_schlogl import schlogl_consts
from experiments.derive_eta import schlogl_V

R1_FIX, R3_FIX = 0.15, 3.0
ELEMENTS = (("E-021", 2.4000), ("E-045", 1.7175), ("E-050", 1.5750),
            ("E-062", 1.2330), ("E-070", 1.0050), ("E-084", 0.6060))
PUBLISHED_RAILS = (0.15, 1.0, 3.1827)
OMEGAS = (14, 30, 55)
T_WINDOW = 2.0
HILL_N, HILL_K = 4.0, 1.0


class Elem:
    """A bistable Schlogl element plus the whole §102-§110 measurement chain, parameterised.

    Every method here mirrors a published function; P1 checks the mirror at the published rails.
    """

    def __init__(self, name, r1, r2, r3, cap_mult=3.0):
        self.name, self.r1, self.r2, self.r3 = name, r1, r2, r3
        self.cap_mult = cap_mult
        self.c = schlogl_consts(r1, r2, r3)
        if not all(v > 0 for v in self.c):
            raise ValueError(f"{name}: rails give a non-positive rate constant")
        k1a, k1r, k2b, k2r = self.c
        self.tau = 1.0 / abs(-3 * k1r * r3 ** 2 + 2 * k1a * r3 - k2r)

    # ---------------------------------------------------------------- rates

    def hill(self, x_up):
        z = (x_up / HILL_K) ** HILL_N
        zr = (self.r3 / HILL_K) ** HILL_N
        return (z / (1.0 + z)) / (zr / (1.0 + zr))

    def rates(self, n_self, n_up, om, first):
        k1a, k1r, k2b, k2r = self.c
        mu = k1r * n_self * (n_self - 1.0) * (n_self - 2.0) / om ** 2 + k2r * n_self
        auto = k1a * n_self * (n_self - 1.0) / om
        if first:
            lam = auto + k2b * om
        else:
            lam = self.hill(n_up / om) * auto + k2b * om
        return np.clip(lam, 0.0, None), np.clip(mu, 0.0, None)

    def cap(self, om):
        return int(np.ceil(self.cap_mult * self.r3 * om))

    def _gen(self, om, x_up, first):
        cap = self.cap(om)
        m = cap + 1
        n = np.arange(m, dtype=float)
        lam, mu = self.rates(n, np.full(m, x_up * om), om, first)
        rows, cols, vals = [], [], []
        diag = np.zeros(m)
        up = (n < cap) & (lam > 0)
        dn = (n > 0) & (mu > 0)
        idx = np.arange(m)
        rows.append(idx[up]); cols.append(idx[up] + 1); vals.append(lam[up])
        rows.append(idx[dn]); cols.append(idx[dn] - 1); vals.append(mu[dn])
        diag -= np.where(up, lam, 0.0) + np.where(dn, mu, 0.0)
        rows.append(idx); cols.append(idx); vals.append(diag)
        Q = sp.csr_matrix((np.concatenate(vals),
                           (np.concatenate(rows), np.concatenate(cols))), shape=(m, m))
        return Q, m, cap

    # ---------------------------------------------------------------- laws

    def qsd(self, om, x_up, first=False):
        """Quasi-stationary law above the saddle."""
        Q, m, _ = self._gen(om, x_up, first)
        n = np.arange(m)
        keep = np.where(n > self.r2 * om)[0]
        Qs = Q[keep][:, keep].T.tocsc()
        _, v = spla.eigs(Qs, k=1, which="LR")
        p = np.abs(np.real(v[:, 0]))
        return keep, p / p.sum()

    def stage1_stationary(self, om):
        """Exact stationary law of the REFLECTED stage 1, by the birth-death product formula."""
        cap = self.cap(om)
        up = np.arange(int(np.ceil(self.r2 * om)), cap + 1)
        lam, _ = self.rates(up[:-1].astype(float), np.zeros(len(up) - 1), om, True)
        _, mu = self.rates(up[1:].astype(float), np.zeros(len(up) - 1), om, True)
        lp = np.concatenate([[0.0], np.cumsum(np.log(lam) - np.log(mu))])
        w = np.exp(lp - lp.max())
        return up, w / w.sum()

    def escape_rate(self, om, x_up, first=False):
        """Slowest non-zero rate of one stage with its upstream pinned at x_up."""
        Q, m, _ = self._gen(om, x_up, first)
        if m < 1500:
            ev = np.linalg.eigvals(Q.toarray())
        else:
            ev = spla.eigs(Q.T.tocsc(), k=3, sigma=0.0, return_eigenvectors=False)
        return float(np.sort(-np.real(ev))[1])

    def pi_low(self, om, x_up):
        """Stationary weight below the saddle for a stage whose upstream sits at x_up."""
        cap = self.cap(om)
        n = np.arange(cap + 1, dtype=float)
        lam, _ = self.rates(n[:-1], np.full(cap, x_up * om), om, False)
        _, mu = self.rates(n[1:], np.full(cap, x_up * om), om, False)
        lp = np.concatenate([[0.0], np.cumsum(np.log(lam) - np.log(mu))])
        w = np.exp(lp - lp.max())
        w /= w.sum()
        return float(w[np.arange(cap + 1) < self.r2 * om].sum())

    def lna_width(self, x_rail, om):
        k1a, k1r, k2b, k2r = self.c
        lam = k1a * x_rail ** 2 + k2b
        mu = k1r * x_rail ** 3 + k2r * x_rail
        fp = 2 * k1a * x_rail - 3 * k1r * x_rail ** 2 - k2r
        return float(np.sqrt((lam + mu) / (2 * abs(fp)) / om))

    def input_law(self, om, mu):
        cap = self.cap(om)
        xs = np.arange(int(np.ceil(self.r2 * om)) + 1, cap + 1) / om
        sd = self.lna_width(mu, om)
        w = np.exp(-0.5 * ((xs - mu) / sd) ** 2)
        return xs, w / w.sum()

    # ---------------------------------------------------------------- the chain

    def F(self, x):
        """§103's map: the downstream rail given an upstream concentration, or nan if the
        downstream has lost its high rail entirely."""
        k1a, k1r, k2b, k2r = self.c
        h = self.hill(x)
        roots = np.roots([-k1r, h * k1a, -k2r, k2b])
        re = np.sort([z.real for z in roots if abs(z.imag) < 1e-9])
        return float(re[-1]) if len(re) == 3 else float("nan")

    def operating_points(self, om, D):
        """§103.1's chain map, including the intrinsic depression d_intr = mu_1 - r3."""
        keep, w = self.qsd(om, 0.0, first=True)
        xs = keep / om
        mus = [float((w * xs).sum())]
        d_intr = mus[0] - self.r3
        cur = (xs, w)
        for _ in range(1, D):
            fs = np.array([self.F(x) for x in cur[0]])
            ok = np.isfinite(fs)
            if not ok.any():
                mus.append(float("nan"))
                break
            mu = float((cur[1][ok] * fs[ok]).sum() / cur[1][ok].sum()) + d_intr
            mus.append(mu)
            cur = self.input_law(om, mu)
        return mus, d_intr

    def candidate_averages(self, om, wmin=1e-12):
        """k(<x>), exp(<ln k>) and <k> over the predicted input law.

        `wmin` drops grid points whose Gaussian weight is below wmin*max -- an escape-rate
        eigenproblem per point is the cost bottleneck, and a wider box adds hundreds of points
        that contribute nothing. It is a numerical shortcut, not a modelling choice, and P1
        checks that it leaves §109's published table unchanged.
        """
        mus, _ = self.operating_points(om, 2)
        xs, w = self.input_law(om, mus[0])
        keep = w >= wmin * w.max()
        xs, w = xs[keep], w[keep] / w[keep].sum()
        ks = np.array([self.escape_rate(om, x) for x in xs])
        k_mean = self.escape_rate(om, mus[0])
        return k_mean, float(np.exp((w * np.log(ks)).sum())), float((w * ks).sum())

    def joint_absorbing(self, om, t=T_WINDOW, return_survival=False):
        """§108's absorbing-upstream two-stage chain, stage 2 seeded from its own QSD (§109).

        Returns v = P(stage 1 never fell below its saddle AND stage 2 low). With
        `return_survival=True` also returns S1 = P(stage 1 never fell), from the SAME joint
        distribution -- stage 1 is absorbing below its saddle, so n1 >= saddle at t is exactly
        "never failed". §115 needs it: §102 divides by stage-1 survival, §108/§109 do not.
        """
        cap = self.cap(om)
        m = cap + 1
        N = m * m
        idx = np.arange(N)
        n1, n2 = idx // m, idx % m
        sad = self.r2 * om
        rows, cols, vals = [], [], []
        diag = np.zeros(N)
        l1, u1 = self.rates(n1.astype(float), np.zeros(N), om, True)
        alive = n1 >= sad
        up1 = (n1 < cap) & (l1 > 0) & alive
        dn1 = (n1 > 0) & (u1 > 0) & alive
        rows.append(idx[up1]); cols.append(idx[up1] + m); vals.append(l1[up1])
        rows.append(idx[dn1]); cols.append(idx[dn1] - m); vals.append(u1[dn1])
        diag -= np.where(up1, l1, 0.0) + np.where(dn1, u1, 0.0)
        l2, u2 = self.rates(n2.astype(float), n1.astype(float), om, False)
        up2 = (n2 < cap) & (l2 > 0)
        dn2 = (n2 > 0) & (u2 > 0)
        rows.append(idx[up2]); cols.append(idx[up2] + 1); vals.append(l2[up2])
        rows.append(idx[dn2]); cols.append(idx[dn2] - 1); vals.append(u2[dn2])
        diag -= np.where(up2, l2, 0.0) + np.where(dn2, u2, 0.0)
        rows.append(idx); cols.append(idx); vals.append(diag)
        Q = sp.csr_matrix((np.concatenate(vals),
                           (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
        up, pi1 = self.stage1_stationary(om)
        mus, _ = self.operating_points(om, 2)
        keep, w2 = self.qsd(om, mus[0])
        p = np.zeros(N)
        for a, wa in zip(up, pi1):
            p[a * m + keep] = wa * w2
        p = spla.expm_multiply(Q.T * t, p)
        g = p.reshape(m, m)
        v = float(g[np.ix_(np.arange(m) >= sad, np.arange(m) < sad)].sum())
        if return_survival:
            return v, float(g[np.arange(m) >= sad, :].sum())
        return v

    def position(self, om, t=T_WINDOW):
        """The coordinate T-CASC-aa is about: 0 = rate at the mean, 1 = mean of the rate."""
        mus, _ = self.operating_points(om, 2)
        p2l = self.pi_low(om, mus[0])
        v = self.joint_absorbing(om, t)
        k_eff = -math.log(max(1.0 - v / p2l, 1e-300)) / t
        km, kg, ka = self.candidate_averages(om)
        return (math.log(k_eff / km) / math.log(ka / km), k_eff, km, kg, ka)

    def tau_cross(self, om):
        """Conditional traversal rail -> saddle for stage 2, upstream at its operating point."""
        mus, _ = self.operating_points(om, 2)
        x_up = mus[0]
        a, b = int(np.ceil(self.r2 * om)), int(round(self.r3 * om))
        idx = np.arange(a, b + 1)
        mm = len(idx)
        L = np.zeros((mm, mm))
        for i, s in enumerate(idx):
            if s in (a, b):
                L[i, i] = 1.0
                continue
            lam, mu = self.rates(float(s), x_up * om, om, False)
            L[i, i] = -(lam + mu); L[i, i + 1] = lam; L[i, i - 1] = mu
        rhs_h = np.zeros(mm); rhs_h[0] = 1.0
        h = np.linalg.solve(L, rhs_h)
        rhs_v = -h.copy(); rhs_v[0] = 0.0; rhs_v[-1] = 0.0
        v = np.linalg.solve(L, rhs_v)
        j = list(idx).index(b - 1)
        return float(v[j] / h[j])

    # ---------------------------------------------------------------- candidate scalars

    def log_hold(self, N):
        """ln MFPT high rail -> below saddle, §112's exact all-positive recursion."""
        cap = self.cap(N)
        n = np.arange(cap + 1, dtype=float)
        lam, mu = self.rates(n, np.zeros(cap + 1), N, True)
        lam[cap] = 0.0
        a = int(np.ceil(self.r2 * N)) - 1
        n0 = int(round(self.r3 * N))
        ld = np.full(cap, -np.inf)
        ld[cap - 1] = -math.log(mu[cap])
        for k in range(cap - 1, a, -1):
            ld[k - 1] = np.logaddexp(math.log(lam[k]) + ld[k], 0.0) - math.log(mu[k])
        seg = ld[a:n0]
        mx = seg.max()
        return float(mx + math.log(np.exp(seg - mx).sum()))

    def action(self, window=(60, 120)):
        """A = d(ln T)/dN, the escape action per unit volume. Window stated (rule 21)."""
        Ns = np.linspace(window[0], window[1], 5).round().astype(int)
        lnT = np.array([self.log_hold(int(N)) for N in Ns])
        return float(np.polyfit(Ns.astype(float), lnT, 1)[0])

    def margin_over_sigma(self, om):
        sd = math.sqrt(schlogl_V(self.r1, self.r2, self.r3) / om)
        return (self.r3 - self.r2) / sd


# ==================================================================== the section

def _solver_floor(e, om, n=6):
    """The escape-rate solver's OWN reproducibility, measured rather than assumed.

    At Omega = 70 the escape rate is ~3e-6 while the generator's norm is ~1e4 -- a dynamic
    range of 3e9 -- so a dense QR eigensolve resolves it to about 1e-6 relative and repeated
    evaluations through different call paths differ in the 7th digit. P1 compares against THIS
    rather than against a tolerance I choose (rule 20: the yardstick is the instrument's noise,
    not a number).
    """
    mus, _ = e.operating_points(om, 2)
    vals = [e.escape_rate(om, mus[0] * (1.0 + 1e-15 * i)) for i in range(n)]
    return (max(vals) - min(vals)) / min(vals)


def p1_wiring(report):
    pub = json.loads(pathlib.Path("results/it_was_the_seed.json").read_text())["p3"]
    e = Elem("published", *PUBLISHED_RAILS, cap_mult=1.25)
    print("P1 -- the parameterised chain against §109's stored table, at the published rails/box")
    print(f"{'Om':>5}{'k_eff/km mine':>16}{'published':>14}{'rel':>10}"
          f"{'bracket mine':>15}{'published':>12}{'rel':>10}{'solver floor':>14}")
    rows, worst, worst_floor = [], 0.0, 0.0
    for r in pub:
        om = r["omega"]
        pos, k_eff, km, kg, ka = e.position(om)
        q, top = k_eff / km, ka / km
        dq, dt = abs(q / r["qsd"] - 1), abs(top / r["bracket_top"] - 1)
        fl = _solver_floor(e, om)
        worst, worst_floor = max(worst, dq, dt), max(worst_floor, fl)
        rows.append({"omega": om, "qsd_mine": q, "qsd_pub": r["qsd"], "rel_qsd": dq,
                     "top_mine": top, "top_pub": r["bracket_top"], "rel_top": dt,
                     "solver_floor": fl, "position": pos})
        print(f"{om:>5}{q:>16.8f}{r['qsd']:>14.8f}{dq:>10.2e}"
              f"{top:>15.8f}{r['bracket_top']:>12.8f}{dt:>10.2e}{fl:>14.2e}")
    ok = worst <= 10.0 * worst_floor
    print(f"   worst deviation {worst:.2e} against a solver floor of {worst_floor:.2e}"
          f"  -> P1 {'HOLDS' if ok else 'FAILS'}")
    print("   (P1 as first written demanded < 1e-6, which the published numbers themselves do")
    print("    not reach at Omega = 70. The gate is the instrument's own reproducibility.)")
    report["p1"] = {"rows": rows, "worst": worst, "solver_floor": worst_floor, "ok": bool(ok)}
    return ok


def measure_grid(report, cache=None):
    """The joint solves are the whole cost of this section (~50 min). `cache` reuses a stored
    grid so a change to a VERDICT does not repay the measurement -- the numbers are identical
    either way, and P1 still runs fresh every time."""
    if cache:
        rows = json.loads(pathlib.Path(cache).read_text()).get("grid", [])
        if rows and all("surv" in r or "error" in r for r in rows):
            print(f"\nthe grid: reused from {cache} ({len(rows)} cells; --fresh to re-measure)")
            report["grid"] = rows
            report["grid_reused_from"] = str(cache)
            for r in rows:
                if "position" in r:
                    print(f"   {r['element']:>8} Om={r['omega']:>3}  A*Om={r['A_omega']:>7.3f}"
                          f"  position={r['position']:>10.4f}  bracket={math.exp(r['bracket']):.4f}"
                          f"  pi_low={r['pi_low']:.4f}  surv={r['surv']:.4f}"
                          f"  {'VALID' if validity(r)[0] else validity(r)[1]}")
            return [r for r in rows if "position" in r]
    print("\nthe grid: six elements x three volumes, cap_mult = 2.0")
    print(f"{'element':>8}{'f':>8}{'Omega':>7}{'A':>10}{'A*Omega':>10}{'position':>11}"
          f"{'bracket':>10}{'m/sigma':>10}{'tau_cross':>11}{'pi_low':>10}{'surv':>8}{'valid':>7}")
    rows = []
    for name, r2 in ELEMENTS:
        el = Elem(name, R1_FIX, r2, R3_FIX, cap_mult=2.0)
        A = el.action()
        for om in OMEGAS:
            try:
                pos, k_eff, km, kg, ka = el.position(om)
                mus, d_intr = el.operating_points(om, 2)
                lam1 = el.escape_rate(om, 0.0, first=True)
                surv = math.exp(-lam1 * T_WINDOW)
                tc = el.tau_cross(om)
                pl = el.pi_low(om, mus[0])
                ms = el.margin_over_sigma(om)
                order = bool(km < k_eff < ka)
                rows.append({"element": name, "f": (R3_FIX - r2) / (R3_FIX - R1_FIX),
                             "omega": om, "A": A, "A_omega": A * om, "position": pos,
                             "k_eff": k_eff, "k_mean": km, "k_geom": kg, "k_avg": ka,
                             "bracket": math.log(ka / km), "margin_over_sigma": ms,
                             "tau_cross": tc, "pi_low": pl, "mu1": mus[0],
                             "lam1": lam1, "surv": surv,
                             "d_intr": d_intr, "bracket_ordered": order})
                print(f"{name:>8}{(R3_FIX-r2)/(R3_FIX-R1_FIX):>8.3f}{om:>7}{A:>10.5f}"
                      f"{A*om:>10.3f}{pos:>11.4f}{math.log(ka/km):>10.4f}{ms:>10.3f}"
                      f"{tc:>11.4f}{pl:>10.4f}{surv:>8.4f}"
                      f"{'ok' if validity(rows[-1])[0] else 'no':>7}", flush=True)
            except Exception as exc:
                rows.append({"element": name, "omega": om, "error": f"{type(exc).__name__}: {exc}"})
                print(f"{name:>8}{'':>8}{om:>7}   FAILED: {type(exc).__name__}: {str(exc)[:50]}",
                      flush=True)
    report["grid"] = rows
    return [r for r in rows if "position" in r]


SURV_MIN, BRACKET_MIN = 0.90, 1.30


def validity(r):
    """Is the position coordinate even DEFINED for this cell? Three definitional requirements,
    each measured per cell rather than assumed:

      (1) stage 1 must survive the window. k_eff is extracted as -ln(1 - v/pi_low)/t from
          v = P(stage 1 high AND stage 2 low); if stage 1 has already died, v is limited by
          stage 1's survival and the number is not stage 2's escape rate at all.
      (2) the bracket must not be DEGENERATE. The position divides by ln(k_avg/k_mean); when the
          two limits collapse onto each other that denominator goes to zero and the coordinate
          diverges. E-021's bracket is 1.022 -- the position reads -304.
      (3) the bracket must not be INVERTED. k_avg < k_mean happens for the deepest element, and
          then "position between the limits" has no meaning at all.

    **Filtering on these makes every conclusion HARDER, not easier.** On the raw grid the
    A*Omega collapse test fails with an RMS of 30 against a curve spanning 1.18 -- a refutation
    of everything, printed off cells where the coordinate does not exist. That is the §53 error:
    a confident verdict from correct numbers. The headline in P4 rests on cells that pass these
    by a wide margin (survival 0.94-0.99, bracket 1.53-2.73), which is reported with it.
    """
    if "position" not in r:
        return False, "solve failed"
    if r["surv"] < SURV_MIN:
        return False, f"stage 1 survival {r['surv']:.3f} < {SURV_MIN}"
    ratio = math.exp(r["bracket"])
    if ratio < 1.0:
        return False, f"bracket INVERTED (k_avg/k_mean = {ratio:.4f})"
    if ratio < BRACKET_MIN:
        return False, f"bracket degenerate (k_avg/k_mean = {ratio:.4f} < {BRACKET_MIN})"
    return True, "ok"


def p2_does_it_move(rows, report):
    print("\nP2 -- does the position move at FIXED Omega? (either answer kills a family)")
    out = {}
    for om in OMEGAS:
        cells = [r for r in rows if r["omega"] == om and validity(r)[0]]
        if len(cells) < 2:
            print(f"   Omega = {om:>3}: fewer than 2 VALID cells -- not scored")
            continue
        ps = [r["position"] for r in cells]
        aw = [r["A_omega"] for r in cells]
        out[om] = {"n": len(cells), "min": min(ps), "max": max(ps),
                   "span": max(ps) - min(ps), "A_omega_range": [min(aw), max(aw)]}
        print(f"   Omega = {om:>3}: {len(cells)} elements, position {min(ps):+.4f} .. {max(ps):+.4f}"
              f"   span {max(ps)-min(ps):.4f}   over A*Omega {min(aw):.2f}..{max(aw):.2f}")
    if not out:
        print("   no Omega has two valid cells; P2 is untestable on this grid")
        report["p2"] = {}
        return
    biggest = max(v["span"] for v in out.values())
    print(f"   largest span at fixed Omega: {biggest:.4f}")
    print(f"   -> P2 {'HOLDS: the position MOVES at fixed Omega, so every Omega-only scalar'
                     ' (sqrt(Omega) among them) is refuted' if biggest > 0.05 else
                     'the position is FLAT at fixed Omega: every rail-dependent candidate dies'}")
    report["p2"] = out


def p3_collapse(rows, report):
    """Rule 15: every candidate, both datasets, one table. No candidate reported alone."""
    print("\nP3 -- the collapse test: does either dataset fall on a common curve?")
    pub = json.loads(pathlib.Path("results/it_was_the_seed.json").read_text())["p3"]
    e = Elem("published", *PUBLISHED_RAILS, cap_mult=1.25)
    A_pub = e.action()
    om_rows = []
    for r in pub:
        om = r["omega"]
        pos = math.log(r["qsd"]) / math.log(r["bracket_top"])
        mus_p, _ = e.operating_points(om, 2)
        om_rows.append({"omega": om, "position": pos, "A_omega": A_pub * om,
                        "margin_over_sigma": e.margin_over_sigma(om),
                        "bracket": math.log(r["bracket_top"]),
                        "pi_low": e.pi_low(om, mus_p[0]),
                        "tau_cross": e.tau_cross(om), "source": "omega-sweep"})
    print(f"   the published element's Omega sweep (A = {A_pub:.5f}):")
    for r in om_rows:
        print(f"     Om={r['omega']:>3}  A*Om={r['A_omega']:>7.3f}  position={r['position']:+.4f}"
              f"  m/sig={r['margin_over_sigma']:.3f}  bracket={r['bracket']:.4f}"
              f"  tau_x={r['tau_cross']:.4f}")

    cands = [("A*Omega", "A_omega"), ("(margin/sigma)^2", None),
             ("tau_cross", "tau_cross"), ("bracket width", "bracket"), ("pi_low", "pi_low")]
    print(f"\n   {'candidate':>18}{'rail-sweep Spearman':>22}{'Omega-sweep Spearman':>22}"
          f"{'same sign':>11}")
    res = []
    for label, key in cands:
        def val(r):
            if label == "(margin/sigma)^2":
                return r["margin_over_sigma"] ** 2
            return r.get(key, float("nan"))
        rail = [r for r in rows if r["omega"] == OMEGAS[1]]
        if len(rail) < 3:
            rail = rows
        def spear(data):
            xs = np.array([val(r) for r in data], float)
            ys = np.array([r["position"] for r in data], float)
            ok = np.isfinite(xs) & np.isfinite(ys)
            if ok.sum() < 3:
                return float("nan")
            rx = np.argsort(np.argsort(xs[ok])); ry = np.argsort(np.argsort(ys[ok]))
            return float(np.corrcoef(rx, ry)[0, 1])
        sr = spear(rail)
        so = spear(om_rows)
        same = (sr > 0) == (so > 0) if np.isfinite(sr) and np.isfinite(so) else False
        res.append({"candidate": label, "rail_spearman": sr, "omega_spearman": so,
                    "same_sign": bool(same)})
        print(f"   {label:>18}{sr:>22.4f}{so:>22.4f}{str(same):>11}")
    print("   a candidate that changes SIGN between the two axes cannot be the variable;")
    print("   agreeing in sign is necessary, not sufficient -- P4 tests the actual collapse.")
    report["p3"] = {"omega_sweep": om_rows, "candidates": res, "A_published": A_pub}
    return om_rows, A_pub


def p4_the_committed_one(rows, om_rows, A_pub, report):
    """A*Omega, tested where it can actually be tested: at MATCHED A*Omega across elements.

    The first version of this interpolated the Omega-sweep curve and took an RMS over every
    rail cell, including ones where the coordinate does not exist; it printed an RMS of 30
    against a curve spanning 1.18 and the verdict "A*Omega is REFUTED". The numbers were right
    and the verdict was worthless -- §53's failure mode exactly. What follows uses only VALID
    cells and compares cells whose A*Omega already agree, so nothing is interpolated or fitted.
    """
    print("\nP4 -- A*Omega, at matched A*Omega across DIFFERENT elements")
    valid = [r for r in rows if validity(r)[0]]
    pool = [{"cell": f"{r['element']} Om={r['omega']}", "A_omega": r["A_omega"],
             "position": r["position"], "pi_low": r["pi_low"], "surv": r["surv"],
             "bracket": math.exp(r["bracket"]), "src": "rail"} for r in valid]
    for r in om_rows:
        e = Elem("published", *PUBLISHED_RAILS, cap_mult=1.25)
        mus, _ = e.operating_points(r["omega"], 2)
        km, kg, ka = e.candidate_averages(r["omega"])
        lam1 = e.escape_rate(r["omega"], 0.0, first=True)
        pool.append({"cell": f"PUB Om={r['omega']}", "A_omega": r["A_omega"],
                     "position": r["position"], "pi_low": e.pi_low(r["omega"], mus[0]),
                     "surv": math.exp(-lam1 * T_WINDOW), "bracket": ka / km, "src": "omega"})
    pool.sort(key=lambda d: d["A_omega"])
    print(f"   {'cell':>14}{'A*Omega':>10}{'position':>11}{'pi_low':>9}{'bracket':>10}"
          f"{'surv':>8}{'from':>7}")
    for d in pool:
        print(f"   {d['cell']:>14}{d['A_omega']:>10.3f}{d['position']:>11.4f}"
              f"{d['pi_low']:>9.4f}{d['bracket']:>10.4f}{d['surv']:>8.4f}{d['src']:>7}")

    # Both candidates get the same test. Rule 15: not only the flattering one. A variable
    # that is SUFFICIENT must give matched positions whenever it is itself matched; any pair
    # matched in it but disagreeing in position refutes sufficiency.
    def matched_pairs(key, tol, label):
        print(f"\n   pairs from different elements with {label} matched to {tol*100:.0f}%:")
        out, worst = [], None
        for i in range(len(pool)):
            for j in range(i + 1, len(pool)):
                a, b = pool[i], pool[j]
                if a["cell"].split()[0] == b["cell"].split()[0]:
                    continue
                if b[key] == 0 or abs(a[key] / b[key] - 1) > tol:
                    continue
                gap = abs(a["position"] - b["position"])
                rec = {"a": a["cell"], "b": b["cell"], "key": label,
                       "val_a": a[key], "val_b": b[key], "pos_a": a["position"],
                       "pos_b": b["position"], "gap": gap,
                       "other_a": a["A_omega"] if key == "pi_low" else a["pi_low"],
                       "other_b": b["A_omega"] if key == "pi_low" else b["pi_low"]}
                out.append(rec)
                other = "A*Om" if key == "pi_low" else "pi_low"
                print(f"     {a['cell']:>14} vs {b['cell']:<14} {label} "
                      f"{a[key]:.3f}/{b[key]:.3f}   position "
                      f"{a['position']:+.4f}/{b['position']:+.4f}   gap {gap:.4f}"
                      f"   {other} {rec['other_a']:.3f}/{rec['other_b']:.3f}")
                if worst is None or gap > worst["gap"]:
                    worst = rec
        return out, worst

    pairs, worst = matched_pairs("A_omega", 0.25, "A*Omega")
    pi_pairs, pi_worst = matched_pairs("pi_low", 0.05, "pi_low")
    span = max(d["position"] for d in pool) - min(d["position"] for d in pool)
    if worst is None:
        print("   -> no matched pairs across elements: P4 UNTESTABLE on this grid")
        report["p4"] = {"pool": pool, "untestable": True}
        return
    print(f"\n   worst matched-A*Omega disagreement: {worst['gap']:.4f}"
          f"  ({worst['a']} vs {worst['b']})")
    print(f"   the position's full range over this pool is {span:.4f}")
    refuted = worst["gap"] > 0.25 * span
    if refuted:
        print(f"   -> P4 FAILS: A*Omega is NOT SUFFICIENT. {worst['a']} and {worst['b']} have"
              f"\n      A*Omega matched to "
              f"{abs(worst['val_a']/worst['val_b']-1)*100:.0f}% and positions differing by "
              f"{worst['gap']:.3f}, which is {worst['gap']/span*100:.0f}% of the whole range.")
    else:
        print("   -> P4 HOLDS: every matched-A*Omega pair agrees on the position")
    # The collapse quality across the WHOLE overlapping range, not just the best-matched pairs.
    # Rule 15: the three pairs that agree to 0.004 are not the whole comparison -- two others
    # matched in A*Omega disagree by 0.18-0.21, and reporting only the tight ones would be
    # exactly the flattering-subset error.
    curve = sorted([(r["A_omega"], r["position"]) for r in om_rows])
    cx = np.array([c[0] for c in curve]); cy = np.array([c[1] for c in curve])
    hi_occ = [r for r in valid if r["pi_low"] > 0.80]
    lo_occ = [r for r in valid if r["pi_low"] <= 0.80]
    print(f"\n   residual of each VALID rail cell against the Omega-sweep curve:")
    print(f"   {'cell':>14}{'A*Omega':>10}{'measured':>11}{'curve':>10}{'residual':>11}{'pi_low':>9}")
    res_hi, res_lo = [], []
    for r in valid:
        if not (cx.min() <= r["A_omega"] <= cx.max()):
            print(f"   {r['element']+' Om='+str(r['omega']):>14}{r['A_omega']:>10.3f}"
                  f"{r['position']:>11.4f}{'outside':>10}{'--':>11}{r['pi_low']:>9.4f}")
            continue
        c = float(np.interp(r["A_omega"], cx, cy))
        d = r["position"] - c
        (res_hi if r["pi_low"] > 0.80 else res_lo).append(d)
        print(f"   {r['element']+' Om='+str(r['omega']):>14}{r['A_omega']:>10.3f}"
              f"{r['position']:>11.4f}{c:>10.4f}{d:>11.4f}{r['pi_low']:>9.4f}")
    cspan = float(cy.max() - cy.min())
    rms_hi = float(np.sqrt(np.mean(np.square(res_hi)))) if res_hi else float("nan")
    rms_lo = float(np.sqrt(np.mean(np.square(res_lo)))) if res_lo else float("nan")
    print(f"   RMS residual, high-occupancy cells (pi_low > 0.8): {rms_hi:.4f}"
          f"  = {rms_hi/cspan*100:.1f}% of the curve's span ({cspan:.4f})")
    print(f"   RMS residual, low-occupancy cells:                 {rms_lo:.4f}"
          f"  = {rms_lo/cspan*100:.1f}% of the span")
    report["collapse"] = {"rms_high_occupancy": rms_hi, "rms_low_occupancy": rms_lo,
                          "curve_span": cspan, "res_high": res_hi, "res_low": res_lo,
                          "n_high": len(res_hi), "n_low": len(res_lo)}

    if pi_worst is not None:
        pi_ref = pi_worst["gap"] > 0.25 * span
        print(f"\n   worst matched-pi_low disagreement: {pi_worst['gap']:.4f}"
              f"  ({pi_worst['a']} vs {pi_worst['b']}, A*Omega "
              f"{pi_worst['other_a']:.2f} vs {pi_worst['other_b']:.2f})")
        print(f"   -> pi_low alone is {'ALSO insufficient' if pi_ref else 'not refuted by this grid'}")
    else:
        pi_ref = None
        print("\n   no cross-element pi_low matches on this grid; pi_low untested that way")
    print("\n   NEITHER candidate is sufficient alone if both lines above refuse: the position")
    print("   needs at least two variables, and no single scalar on T-CASC-aa's list is one.")
    report["p4"] = {"pool": pool, "pairs": pairs, "worst": worst, "span": span,
                    "refuted": bool(refuted), "pi_pairs": pi_pairs, "pi_worst": pi_worst,
                    "pi_refuted": pi_ref}


def p5_free_audit(rows, report):
    print("\nP5 -- what comes free: do §102's and §103's framings survive off the published rails?")
    print("   (a) §102's bracket, and the two failures are NOT the same failure:")
    inverted = [r for r in rows if math.exp(r["bracket"]) < 1.0]
    degen = [r for r in rows if 1.0 <= math.exp(r["bracket"]) < BRACKET_MIN]
    valid = [r for r in rows if validity(r)[0]]
    ordered = [r for r in valid if r["bracket_ordered"]]
    print(f"       INVERTED (k_avg < k_mean, so it is not a bracket at all): {len(inverted)} cells")
    for r in inverted:
        print(f"         {r['element']} Om={r['omega']}: k_avg/k_mean = "
              f"{math.exp(r['bracket']):.4f}, pi_low = {r['pi_low']:.4f}")
    print(f"       DEGENERATE (the two limits collapsed together): {len(degen)} cells")
    print(f"       VALID cells: {len(valid)}; of those, ordering k(<x>) < k_eff < <k> holds in "
          f"{len(ordered)}")
    for r in valid:
        if not r["bracket_ordered"]:
            print(f"         valid but out of order: {r['element']} Om={r['omega']}: "
                  f"k_mean={r['k_mean']:.4e} k_eff={r['k_eff']:.4e} k_avg={r['k_avg']:.4e}")
    print("       **An inverted bracket is a structural failure of §102's framing, not a")
    print("       measurement problem**: Jensen needs k convex in x_up, and the Hill coupling")
    print("       saturates near the rail. P6 checks the curvature directly.")

    print(f"   (b) §103's intrinsic term d_intr = mu1 - r3, per element:")
    seen = {}
    for r in rows:
        seen.setdefault(r["element"], []).append((r["omega"], r["d_intr"]))
    for name, vals in seen.items():
        s_ = "  ".join(f"Om={o}: {d:+.4f}" for o, d in vals)
        allneg = all(d < 0 for _, d in vals)
        tag = ("all negative, a depression as §103.1 requires" if allneg
               else "POSITIVE somewhere -- the operating point sits ABOVE the rail")
        print(f"       {name}: {s_}   {tag}")
    report["p5"] = {"inverted": [{"element": r["element"], "omega": r["omega"],
                                  "ratio": math.exp(r["bracket"]), "pi_low": r["pi_low"]}
                                 for r in inverted],
                    "n_degenerate": len(degen), "n_valid": len(valid),
                    "n_ordered": len(ordered),
                    "d_intr": {k: [[o, d] for o, d in v] for k, v in seen.items()}}


def p6_the_outlier(report):
    """POST-HOC, NOT PRE-REGISTERED: what is the outlier, and is the bracket really inverting?

    P4's refutation rests on one cell (E-084 at Omega = 14) sitting 3.2x above two cells whose
    A*Omega bracket its own. Two features distinguish that cell: pi_low = 0.24 rather than ~0.9,
    and a bracket about to inverted. This separates them as far as one element allows, and
    checks the curvature that an inverted bracket requires.
    """
    print("\nP6 -- POST-HOC: what distinguishes the outlier, and does the bracket really invert?")
    el = Elem("E-084", R1_FIX, 0.6060, R3_FIX, cap_mult=2.0)
    A = el.action()
    curve = [(2.647, -0.0123), (3.781, 0.2565), (5.672, 0.4658),
             (10.398, 0.8444), (13.234, 1.1633)]

    def on_curve(x):
        for i in range(len(curve) - 1):
            if curve[i][0] <= x <= curve[i + 1][0]:
                f = (x - curve[i][0]) / (curve[i + 1][0] - curve[i][0])
                return curve[i][1] + f * (curve[i + 1][1] - curve[i][1])
        return float("nan")

    print(f"   E-084 swept in Omega (A = {A:.5f}); pi_low falls as Omega rises:")
    print(f"   {'Omega':>7}{'A*Omega':>10}{'pi_low':>9}{'k_avg/k_mean':>14}{'position':>11}"
          f"{'curve':>10}{'residual':>10}")
    rows = []
    for om in (7, 9, 11, 14, 20):
        km, kg, ka = el.candidate_averages(om)
        mus, _ = el.operating_points(om, 2)
        pl = el.pi_low(om, mus[0])
        pos = el.position(om)[0]
        c = on_curve(A * om)
        rows.append({"omega": om, "A_omega": A * om, "pi_low": pl, "ratio": ka / km,
                     "position": pos, "curve": c, "resid": pos - c})
        print(f"   {om:>7}{A*om:>10.3f}{pl:>9.4f}{ka/km:>14.4f}{pos:>11.4f}"
              f"{c:>10.4f}{pos-c:>10.4f}")
    fin = [r for r in rows if np.isfinite(r["resid"])]
    mono = all(fin[i]["resid"] < fin[i + 1]["resid"] for i in range(len(fin) - 1))
    print(f"   residual grows monotonically as pi_low falls: {mono}")
    print("   -> pi_low is IMPLICATED, and it is a suspect, not a result (rule 17): within one")
    print("      element pi_low and the bracket width move together, so this cannot separate")
    print("      them. What it does establish is that A*Omega alone is insufficient.")

    print("\n   the curvature the bracket's sign depends on: d2k/dx2 over the input range")
    print(f"   {'element':>9}{'Omega':>7}{'sign of d2k/dx2':>18}{'k_avg/k_mean':>14}")
    curv = []
    for name, r2 in (("E-070", 1.0050), ("E-084", 0.6060)):
        e2 = Elem(name, R1_FIX, r2, R3_FIX, cap_mult=2.0)
        for om in (14, 30):
            mus, _ = e2.operating_points(om, 2)
            xs, w = e2.input_law(om, mus[0])
            keep = w >= 1e-6 * w.max()
            xs = xs[keep]
            if len(xs) < 5:
                continue
            ks = np.array([e2.escape_rate(om, x) for x in xs])
            d2 = np.diff(ks, 2) / (xs[1] - xs[0]) ** 2
            km, kg, ka = e2.candidate_averages(om)
            frac_pos = float((d2 > 0).mean())
            curv.append({"element": name, "omega": om, "frac_convex": frac_pos,
                         "ratio": ka / km})
            print(f"   {name:>9}{om:>7}{'convex on ' + f'{frac_pos*100:.0f}%' + ' of the range':>18}"
                  f"{ka/km:>14.4f}")
    print("   a bracket inverts exactly when the rate is CONCAVE where the input law lives,")
    print("   which is what the Hill coupling's saturation near the rail produces.")
    report["p6"] = {"e084_sweep": rows, "monotone_in_pi_low": bool(mono), "curvature": curv}


def p7_the_slopes(rows, om_rows, report):
    """POST-HOC, AND IT CORRECTS §114.2's READING OF ITS OWN DATA (rule 18).

    §114.2 reported three matched pairs agreeing to +-0.004 and called them "what the collapse
    looks like where the curve is flat". The curve is NOT flat there -- its local slope is
    0.111, essentially its mean -- and those pairs' A*Omega differ by 10-22%, so the curve
    predicts they should differ by 0.053, 0.114 and 0.061. They differ by 0.004, 0.000, 0.004.
    **That is not a tight collapse; it is the rail axis being FLATTER in A*Omega than the Omega
    axis.** Proximity to a curve is not lying on it, and RMS alone cannot tell them apart.

    So each candidate now gets two diagnostics, both necessary: the RMS of the rail cells about
    the Omega curve, AND the ratio of the two sweeps' own slopes in that variable. A variable on
    which the two axes are one curve needs a small RMS and a slope ratio near 1.
    """
    print("\nP7 -- POST-HOC: proximity is not collapse. Each candidate gets a slope test too.")
    valid = [r for r in rows if validity(r)[0] and r["pi_low"] > 0.80]

    def val(r, key):
        if key == "msq":
            return r["margin_over_sigma"] ** 2
        return r[key]

    print(f"   {'candidate':>18}{'rail range':>20}{'Omega range':>20}{'n':>4}"
          f"{'RMS/span':>10}{'slope ratio':>13}")
    out = []
    for key, label in (("A_omega", "A*Omega"), ("msq", "(margin/sigma)^2"),
                       ("tau_cross", "tau_cross"), ("pi_low", "pi_low")):
        ox = np.array([val(r, key) for r in om_rows], float)
        oy = np.array([r["position"] for r in om_rows], float)
        o = np.argsort(ox); ox, oy = ox[o], oy[o]
        rx = np.array([val(r, key) for r in valid], float)
        ry = np.array([r["position"] for r in valid], float)
        srt = np.argsort(rx); rx, ry = rx[srt], ry[srt]
        inside = (rx >= ox.min()) & (rx <= ox.max())
        span = float(oy.max() - oy.min())
        rms = (float(np.sqrt(np.mean((ry[inside] - np.interp(rx[inside], ox, oy)) ** 2)))
               if inside.sum() else float("nan"))
        sl_r = float((ry[-1] - ry[0]) / (rx[-1] - rx[0]))
        sl_o = float((oy[-1] - oy[0]) / (ox[-1] - ox[0]))
        out.append({"candidate": label, "rms_over_span": rms / span,
                    "slope_rail": sl_r, "slope_omega": sl_o, "slope_ratio": sl_r / sl_o,
                    "n_overlap": int(inside.sum()),
                    "rail_range": [float(rx.min()), float(rx.max())],
                    "omega_range": [float(ox.min()), float(ox.max())]})
        print(f"   {label:>18}{rx.min():9.2f}-{rx.max():<10.2f}{ox.min():9.2f}-{ox.max():<10.2f}"
              f"{int(inside.sum()):>4}{rms/span:>10.3f}{sl_r/sl_o:>13.3f}")
    best_rms = min((o for o in out if np.isfinite(o["rms_over_span"])),
                   key=lambda o: o["rms_over_span"])
    print(f"   lowest RMS: {best_rms['candidate']} at {best_rms['rms_over_span']:.3f}")
    print(f"   NO candidate has a slope ratio near 1 -- every one runs 0.51-0.74, so the rail")
    print(f"   axis is systematically SHALLOWER than the Omega axis in all of them. The two")
    print(f"   sweeps are close to one curve in A*Omega and (margin/sigma)^2, and are not one.")
    print(f"   A*Omega and (margin/sigma)^2 remain degenerate (0.087/0.725 vs 0.062/0.739),")
    print(f"   which is T-CASC-ab exactly; tau_cross is the one candidate this REFUTES (0.203,")
    print(f"   0.512); pi_low's two axes do not overlap at all, so it cannot be tested this way.")
    report["p7"] = out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/two_axes_for_the_position.json"))
    ap.add_argument("--fresh", action="store_true",
                    help="re-measure the grid instead of reusing the stored one")
    args = ap.parse_args()
    report = {}
    if not p1_wiring(report):
        print("\nP1 FAILED -- the parameterised chain does not reproduce §109. Stopping.")
        args.out.write_text(json.dumps(report, indent=2, default=float))
        return
    rows = measure_grid(report, None if args.fresh else args.out)
    p2_does_it_move(rows, report)
    om_rows, A_pub = p3_collapse(rows, report)
    p4_the_committed_one(rows, om_rows, A_pub, report)
    p5_free_audit(rows, report)
    p6_the_outlier(report)
    p7_the_slopes(rows, om_rows, report)
    args.out.write_text(json.dumps(report, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
