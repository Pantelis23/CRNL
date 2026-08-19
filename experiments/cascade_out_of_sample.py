"""T-CASC-j: run the closed chain FORWARD on an element it was never built from

§91-§96 assembled a prediction chain for a chemically-coupled cascade, every link of which uses
only SINGLE-element quantities:

    1. the mean    mu_2 = <F(x_up)> over stage 1's own law, plus the element's intrinsic shift
                          (§96 -- computed exactly, not truncated at second order)
    2. the width   sigma_2 = the LNA width at that operating point           (§95)
    3. the penalty log(penalty) = -0.952 x margin/sigma                      (§91)

Each was calibrated on ONE Schloegl landscape (roots 0.15/1/3.1827) with ONE Hill coupling
(n = 4, K = 1) at Omega = 30. **A chain of three fitted-or-calibrated steps checked only on the
data it was built from is exactly what rule 16 exists to stop.** So it is run forward here on a
different element and a different coupling, with NOTHING measured on the cascade, against one
exact joint solve.

**The coefficient -0.952 is the exposed part.** The mean and width steps are derivations; the
penalty step is an empirical slope fitted on the old element in §91, and this is the first thing
that could show it to be element-specific.

PREDICTIONS, written before running.

  P1  GATE. The new element must be genuinely different and still a cascade: report its roots,
      escape action, collapse point and margin, and confirm the coupling TRANSMITS (§91 P1(b) --
      a downstream that keeps its high rail when the upstream is low carries no signal, and a
      neutrality gate alone would not catch it).
  P2  **THE MEAN.** Predicted from the exact static-transfer average. **Predicted: within ~0.5%**,
      since §96 got 0.12% on the calibration element and the only new error is that the residual
      (the finite correlation time) may differ.
  P3  **THE WIDTH.** LNA at the predicted operating point. **Predicted: within ~3%**, the LNA's own
      accuracy on the calibration element.
  P4  **THE PENALTY, and this is the real test.** §91's slope applied to the new element's margin.
      **Predicted: within §91's own 18% scatter IF the slope transfers.** If the mean and width
      land but the penalty does not, **-0.952 is element-specific and §91's law is a fit, not a
      law** -- which would be the most useful outcome here.
  P5  **RULE 15.** Report all three residuals whatever they are, and do not re-fit the slope to
      the new element. A coefficient re-fitted per element is not a law.
  P6  **IF THE SLOPE FAILS, THE FORMULA IT APPROXIMATES SHOULD NOT.** §92's frozen-upstream
      average <exp(-[A(x_up) - A(r3)] Omega)> is a DERIVATION, and it carries Omega and the
      barrier explicitly where §91's slope carries only the margin. This element's barrier is
      A*Omega = 13.2 against the calibration's 5.7, so a margin-only law cannot be right across
      both. **Predicted: the frozen average transfers where the fitted slope does not** -- and if
      it does, §91's law is exposed as a one-element parameterisation of it.
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

# ---- the NEW element and the NEW coupling. Nothing here was used to calibrate anything.
RAILS_NEW = (0.20, 1.0, 4.3505)      # §82's affinity level set, a different point on it
HILL_N_NEW, HILL_K_NEW = 6.0, 1.3
SLOPE_91 = -0.952                    # §91's fitted slope, quoted and NOT re-fitted (rule 16)


def consts():
    return schlogl_consts(*RAILS_NEW)


def hill_new(x_up):
    r3 = RAILS_NEW[2]
    f = lambda z: np.power(max(z, 0.0), HILL_N_NEW) / (
        HILL_K_NEW ** HILL_N_NEW + np.power(max(z, 0.0), HILL_N_NEW))
    return f(x_up) / f(r3)


def down_roots(x_up):
    k1a, k1r, k2b, k2r = consts()
    r = np.roots([-k1r, k1a * hill_new(x_up), -k2r, k2b])
    return np.sort([z.real for z in r if abs(z.imag) < 1e-9 and z.real > 1e-12])


def F_new(x_up):
    r = down_roots(x_up)
    return float(r[-1]) if len(r) >= 3 else np.nan


def A_new():
    r1, r2, r3 = RAILS_NEW
    k1a, k1r, k2b, k2r = consts()
    v, _ = quad(lambda x: np.log((k1r * x ** 3 + k2r * x) / (k1a * x ** 2 + k2b)), r2, r3,
                limit=200)
    return -float(v)


def rates(n, n_up, om, first):
    k1a, k1r, k2b, k2r = consts()
    auto = k1a * n * (n - 1.0) / om
    mu = k1r * n * (n - 1.0) * (n - 2.0) / om ** 2 + k2r * n
    lam = (auto + k2b * om) if first else (hill_new(n_up / om) * auto + k2b * om)
    return max(lam, 0.0), max(mu, 0.0)


def stage1_law(om, cap_mult=1.25):
    r1, r2, r3 = RAILS_NEW
    cap = int(np.ceil(cap_mult * r3 * om))
    up = np.arange(int(np.ceil(r2 * om)), cap + 1)
    lp, acc = np.zeros(len(up)), 0.0
    for i in range(1, len(up)):
        l, _ = rates(float(up[i - 1]), 0.0, om, True)
        _, u = rates(float(up[i]), 0.0, om, True)
        acc += np.log(l) - np.log(u)
        lp[i] = acc
    w = np.exp(lp - lp.max())
    return up, w / w.sum(), cap


def lna_width(x_rail, om):
    k1a, k1r, k2b, k2r = consts()
    lam = k1a * x_rail ** 2 + k2b
    mu = k1r * x_rail ** 3 + k2r * x_rail
    fp = 2 * k1a * x_rail - 3 * k1r * x_rail ** 2 - k2r
    return float(np.sqrt((lam + mu) / (2 * abs(fp)) / om))


def build(om, free_last, cap_mult=1.25):
    r1, r2, r3 = RAILS_NEW
    cap = int(np.ceil(cap_mult * r3 * om))
    ref = np.arange(int(np.ceil(r2 * om)), cap + 1)
    nr = len(ref)
    m2 = cap + 1 if free_last else nr
    rows, cols, vals = [], [], []
    diag = np.zeros(nr * m2)
    idx = np.arange(nr * m2)
    a, b = idx // m2, idx % m2
    c1 = ref[a].astype(float)
    c2 = (b if free_last else ref[b]).astype(float)
    l1 = np.array([rates(x, 0.0, om, True)[0] for x in c1])
    u1 = np.array([rates(x, 0.0, om, True)[1] for x in c1])
    l2 = np.array([rates(y, x, om, False)[0] for y, x in zip(c2, c1)])
    u2 = np.array([rates(y, x, om, False)[1] for y, x in zip(c2, c1)])
    lo2 = 0.0 if free_last else float(ref[0])
    for cnt, lam, mu, st, lo in ((c1, l1, u1, m2, float(ref[0])), (c2, l2, u2, 1, lo2)):
        up_ok = (cnt < cap) & (lam > 0)
        dn_ok = (cnt > lo) & (mu > 0)
        rows.append(idx[up_ok]); cols.append(idx[up_ok] + st); vals.append(lam[up_ok])
        rows.append(idx[dn_ok]); cols.append(idx[dn_ok] - st); vals.append(mu[dn_ok])
        diag -= np.where(up_ok, lam, 0.0) + np.where(dn_ok, mu, 0.0)
    rows.append(idx); cols.append(idx); vals.append(diag)
    Q = sp.csr_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(nr * m2, nr * m2))
    return Q, ref, nr, m2, cap


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=30)
    ap.add_argument("--t", type=float, default=2.0)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/cascade_out_of_sample.json"))
    args = ap.parse_args()
    om, t = args.omega, args.t
    r1, r2, r3 = RAILS_NEW
    out = {"rails": list(RAILS_NEW), "hill": [HILL_N_NEW, HILL_K_NEW], "omega": om}

    print("=== P1 GATE: is this a different element, and does the coupling TRANSMIT?")
    xs = np.linspace(r1, r3, 4001)
    xc = next((float(x) for x in xs[::-1] if len(down_roots(x)) < 3), None)
    print(f"  rails {RAILS_NEW}  (calibration element was 0.15/1.0/3.1827)")
    print(f"  Hill n={HILL_N_NEW}, K={HILL_K_NEW}  (calibration was n=4, K=1)")
    print(f"  escape action A = {A_new():.6f}   (calibration element: 0.190241)")
    print(f"  F(r3) = {F_new(r3):.10f} vs r3 = {r3}   -> neutral at the rail")
    print(f"  downstream roots at r1: {np.round(down_roots(r1), 4)}")
    trans = len(down_roots(r1)) == 1
    print(f"  collapse point x_crit = {xc}")
    ok1 = abs(F_new(r3) - r3) < 1e-9 and trans and xc is not None
    print(f"  -> P1 {'HOLDS: neutral at the rail AND it transmits' if ok1 else 'FAILS: not a cascade of this element'}")
    assert ok1

    up, pi1, cap = stage1_law(om)
    x1 = up / om
    mu1 = float((pi1 * x1).sum())
    sd1 = float(np.sqrt((pi1 * (x1 - mu1) ** 2).sum()))
    d_intr = mu1 - r3
    margin1 = (r3 - xc) / sd1
    print(f"\n  stage 1 (exact, 1-D): mean {mu1:.5f}  sd {sd1:.5f}"
          f"  intrinsic shift {d_intr:+.5f}  margin {margin1:.3f} sd")

    print("\n=== the PREDICTIONS, made before the joint solve, from single-element quantities")
    mu2_pred = float(sum(w * F_new(x) for w, x in zip(pi1, x1) if np.isfinite(F_new(x)))) + d_intr
    sd2_pred = lna_width(mu2_pred, om)
    margin2 = (r3 - xc) / sd2_pred
    pen_pred = float(np.exp(SLOPE_91 * (margin2 - margin1)))
    print(f"  mu_2      predicted = {mu2_pred:.5f}")
    print(f"  sigma_2   predicted = {sd2_pred:.5f}   (margin {margin2:.3f} sd)")
    print(f"  penalty   predicted = {pen_pred:.4f}   using §91's slope {SLOPE_91}, NOT re-fitted")

    print("\n=== the exact joint solve")
    Qr, ref, nr, m2, cap = build(om, free_last=False)
    p = np.zeros(nr * m2)
    pos = list(ref).index(int(round(r3 * om)))
    for i, w in enumerate(pi1):
        p[i * m2 + pos] = w
    p = spla.expm_multiply(Qr.T * t, p)
    idx = np.arange(nr * m2)
    c2 = ref[idx % m2].astype(float)
    w = p / p.sum()
    mu2 = float((w * c2).sum()) / om
    sd2 = float(np.sqrt((w * (c2 - mu2 * om) ** 2).sum())) / om
    Qf, reff, nrf, m2f, capf = build(om, free_last=True)
    pf = np.zeros(nrf * m2f)
    for i, wv in enumerate(pi1):
        pf[i * m2f + int(round(r3 * om))] = wv
    pf = spla.expm_multiply(Qf.T * t, pf)
    lo = (np.arange(nrf * m2f) % m2f) < r2 * om
    p_chain = float(pf[lo].sum())
    # the isolated reference: stage 2 with the upstream pinned at r3
    rows, cols, vals = [], [], []
    diag = np.zeros(cap + 1)
    for n in range(cap + 1):
        l, u = rates(float(n), r3 * om, om, False)
        if n < cap and l > 0:
            rows.append(n); cols.append(n + 1); vals.append(l); diag[n] -= l
        if n > 0 and u > 0:
            rows.append(n); cols.append(n - 1); vals.append(u); diag[n] -= u
    rows += list(range(cap + 1)); cols += list(range(cap + 1)); vals += list(diag)
    Qi = sp.csr_matrix((vals, (rows, cols)), shape=(cap + 1, cap + 1))
    q = np.zeros(cap + 1); q[int(round(r3 * om))] = 1.0
    q = spla.expm_multiply(Qi.T * t, q)
    eps_iso = float(q[np.arange(cap + 1) < r2 * om].sum())
    pen = p_chain / eps_iso
    print(f"  mu_2    measured = {mu2:.5f}")
    print(f"  sigma_2 measured = {sd2:.5f}")
    print(f"  penalty measured = {pen:.4f}   (eps_iso = {eps_iso:.4e})")

    print("\n=== P2/P3/P4/P5: the three residuals, reported whatever they are")
    print(f"{'quantity':>12}{'predicted':>12}{'measured':>12}{'ratio':>9}{'error':>9}")
    rows_out = [("mean", mu2_pred, mu2), ("width", sd2_pred, sd2), ("penalty", pen_pred, pen)]
    for name, pr, me in rows_out:
        print(f"{name:>12}{pr:>12.5f}{me:>12.5f}{pr/me:>9.4f}{100*(pr/me-1):>8.2f}%")
    out["pred"] = {"mu2": mu2_pred, "sd2": sd2_pred, "pen": pen_pred}
    out["meas"] = {"mu2": mu2, "sd2": sd2, "pen": pen, "eps_iso": eps_iso}
    e_mu = abs(mu2_pred / mu2 - 1)
    e_sd = abs(sd2_pred / sd2 - 1)
    e_pen = abs(pen_pred / pen - 1)
    print(f"  -> P2 {'HOLDS' if e_mu < 0.005 else 'FAILS'}: the mean is {100*e_mu:.2f}% off")
    print(f"  -> P3 {'HOLDS' if e_sd < 0.03 else 'FAILS'}: the width is {100*e_sd:.2f}% off")
    print(f"  -> P4 {'HOLDS: §91s slope TRANSFERS to an element it was not fitted on' if e_pen < 0.18 else 'FAILS: the slope -0.952 does NOT transfer -- §91s law is a fit calibrated on one element, not a law'}"
          f": the penalty is {100*e_pen:.1f}% off")
    print("  -> P5 all three reported; the slope was NOT re-fitted to this element.")

    print("\n=== P6: the DERIVATION behind §91's fit -- §92's frozen-upstream average")
    A0 = A_new()

    def A_at(x_up):
        r = down_roots(x_up)
        if len(r) < 3:
            return 0.0
        k1a, k1r, k2b, k2r = consts()
        v, _ = quad(lambda z: np.log((k1r * z ** 3 + k2r * z)
                                     / (k1a * hill_new(x_up) * z ** 2 + k2b)),
                    r[1], r[2], limit=200)
        return -float(v)

    w_frozen = np.array([np.exp(-(A_at(x) - A0) * om) for x in x1])
    frozen = float((pi1 * w_frozen).sum())
    fast = float(np.exp(-(A_at(mu1) - A0) * om))
    print(f"  A(r3) = {A0:.6f};  A*Omega = {A0*om:.1f}"
          f"  (calibration element: {0.190241*om:.1f})")
    print(f"  frozen-upstream average <exp(-dA*Om)> = {frozen:.4f}")
    print(f"  fast limit exp(-<dA>*Om)              = {fast:.4f}")
    print(f"  measured penalty                      = {pen:.4f}")
    out["p6"] = {"A": A0, "frozen": frozen, "fast": fast}
    lo_, hi_ = min(frozen, fast), max(frozen, fast)
    inside = lo_ * 0.9 <= pen <= hi_ * 1.1
    print(f"  frozen/measured = {frozen/pen:.3f};  fast/measured = {fast/pen:.3f}")
    print(f"  -> P6 {'HOLDS: the derivation brackets the measurement on an element it was never calibrated on, so §91s slope is a ONE-ELEMENT PARAMETERISATION of it and the margin alone is not the controlling variable -- the barrier depth A*Omega enters too' if inside else 'FAILS: even the derivation misses this element, so the frozen-average account does not transfer either'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
