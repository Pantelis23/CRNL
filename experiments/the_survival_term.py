"""§115 -- the survival term: §102 divides by stage-1 survival, §108/§109/§114 do not.

The fast/frozen POSITION of §108-§114 is `ln(k_eff/k_mean) / ln(k_avg/k_mean)`, and every
number in it is exact except one step: k_eff is inverted from a joint quantity,

    v = P(stage 1 never fell below its saddle  AND  stage 2 below its saddle)   at t = 2,
    k_eff = -ln(1 - v / pi_low) / t.

`pi_low` removes stage 2's own two-state occupancy (§106.2). Nothing removes stage 1's death.
Stage 1 is absorbing below its saddle, so every trajectory in which it failed is simply missing
from `v`, and for small lambda*t the inversion returns roughly S1 * k_true, where S1 = P(stage 1
never fell) over the window. §102's `rate_limits` does divide by a survival (`frac = pure / surv`);
§108's `effective_rate`, §109's P3, and §114's `Elem.position` -- a faithful mirror of §109 -- do
not. The two halves of the arc invert the same kind of quantity two different ways.

**Why this is a domino, and rule 9 is the reason.** S1 is not a constant. It RISES with Omega --
§114's diagnostic has the published element's stage-1 survival at 0.864 (Omega = 14), 0.992 (30),
1.000 (70) -- and Omega is exactly the axis along which the position "drifts" from ~0 to ~1.16.
A bias of ln(S1) in ln(k_eff) would pull the low-Omega end of that curve DOWN and steepen it,
without any physics. Three explanatory families have been retired for this drift (§108-§110).
None of those sections looked at the instrument's own survival term.

The corrected inversion uses the SAME joint distribution and needs no new model:

    k_eff_corr = -ln(1 - (v / S1) / pi_low) / t,    S1 = P(n1 >= saddle at t)

-- the conditional probability that stage 2 is low given its input never failed, which is what
"stage 2's escape rate while its input fluctuates in the high basin" means.

PREDICTIONS, WRITTEN BEFORE RUNNING.

  P1  WIRING. With the survival term ignored, this reproduces §109's stored `qsd` column and
      §114's stored positions to the solver's floor (~3e-6). The joint solve's `v` is unchanged
      by the new return path.

  P2  THE SIZE OF S1. The stage-1 seed is the REFLECTED stationary law, which §109 measured
      escaping 1.24-1.28x faster than the free QSD, so S1 should sit a little BELOW
      exp(-lambda1 * t): predict S1 ~ 0.80-0.86 for the published element at Omega = 14, above
      0.98 by Omega = 30, and indistinguishable from 1 by Omega = 55.

  P3  THE SHIFT. For small lambda*t the position moves up by about -ln(S1)/ln(k_avg/k_mean).
      With ln(bracket) = 0.873 at Omega = 14, predict the published curve's first point rises by
      0.15-0.25 and its last point by < 0.01, so **the drift's span shrinks by roughly a fifth**.
      It cannot remove the drift -- S1 is ~1 over most of the curve -- and it is not offered as
      the explanation of it. If the shift is < 0.05 the survival term is negligible and this is
      a null result, reported as one.

  P4  THE KILL TEST FOR THE CORRECTION ITSELF, which has to be independent of the position.
      A rate should not depend on the window it was measured in. Uncorrected, k_eff inherits
      S1(t), which decays with t; corrected, it should not. Measure both at t = 2 and t = 4:
      predict the UNcorrected ratio k(4)/k(2) departs from 1 by several percent at Omega = 14 and
      the corrected one by clearly less. If the corrected extraction is NOT more window-stable,
      dividing by S1 is not the right fix and P3's numbers mean nothing.

  P5  T-CASC-af, the domino, and it can print either verdict. §114.9 found the rail axis
      systematically shallower than the Omega axis -- slope ratios 0.725 (A*Omega) and 0.739
      ((margin/sigma)^2). The Omega sweep's steepest contribution is its low-Omega point, where
      S1 is smallest; the rail cells mostly sit at S1 > 0.98. **Predict the correction moves both
      slope ratios toward 1.** If they stay at ~0.73 the survival bias is not what makes the rail
      axis shallower, and T-CASC-af stays exactly as open as it was.

  P6  §114.3 SURVIVES, and rule 14 is why it is tested. The correction can only RAISE positions
      (S1 <= 1). E-084 at Omega = 14 is anomalously HIGH (1.47 against a bracketing pair at
      0.466), so the correction cannot rescue A*Omega's sufficiency unless it raises the pair far
      more than E-084 -- and the pair has the higher survival. Predict the matched triple's
      verdict is unchanged. A refutation is a claim too, and this is the first thing that could
      have moved it.

  P7  WHAT ELSE INHERITS. List every published number whose extraction omits S1 -- §108's
      true/k(<x>) table, §109's P3 `qsd` column, §110's position column, §114's grid -- with the
      corrected value printed BESIDE the original (rule 7). None is edited.

WHAT THIS CANNOT SETTLE. Dividing by S1 conditions on stage 1's survival exactly, but the
inversion through `pi_low` is still §106.2's two-state model, which assumes stage 2's own
relaxation is single-exponential. This removes one bias, identified and measured; it does not
certify the inversion.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib

import numpy as np

from experiments.two_axes_for_the_position import (
    ELEMENTS, PUBLISHED_RAILS, R1_FIX, R3_FIX, T_WINDOW, Elem,
)

OMEGAS_PUB = (14, 20, 30, 55, 70)
RAIL_CELLS = (("E-062", 30), ("E-062", 55), ("E-070", 30), ("E-070", 55), ("E-084", 14))
CAP = 2.0


def invert(v, s1, pl, t):
    """Both inversions from one joint solve. `s1=1` gives §108/§109's."""
    return -math.log(max(1.0 - (v / s1) / pl, 1e-300)) / t


def cell(el, om, t=T_WINDOW):
    mus, _ = el.operating_points(om, 2)
    pl = el.pi_low(om, mus[0])
    v, s1 = el.joint_absorbing(om, t, return_survival=True)
    km, kg, ka = el.candidate_averages(om)
    k0, k1 = invert(v, 1.0, pl, t), invert(v, s1, pl, t)
    lb = math.log(ka / km)
    return {"omega": om, "v": v, "S1": s1, "pi_low": pl, "k_mean": km, "k_avg": ka,
            "k_uncorr": k0, "k_corr": k1,
            "pos_uncorr": math.log(k0 / km) / lb, "pos_corr": math.log(k1 / km) / lb,
            "A_omega": el.action() * om, "msq": el.margin_over_sigma(om) ** 2}


def slope_ratio(rail, omg, xkey, ykey):
    ox = np.array([r[xkey] for r in omg]); oy = np.array([r[ykey] for r in omg])
    rx = np.array([r[xkey] for r in rail]); ry = np.array([r[ykey] for r in rail])
    o, s = np.argsort(ox), np.argsort(rx)
    ox, oy, rx, ry = ox[o], oy[o], rx[s], ry[s]
    sl_o = (oy[-1] - oy[0]) / (ox[-1] - ox[0])
    sl_r = (ry[-1] - ry[0]) / (rx[-1] - rx[0])
    inside = (rx >= ox.min()) & (rx <= ox.max())
    rms = float(np.sqrt(np.mean((ry[inside] - np.interp(rx[inside], ox, oy)) ** 2)))
    return float(sl_r / sl_o), rms / float(oy.max() - oy.min())


def pinned_control(el, om, t):
    """Stage 2 ALONE, upstream pinned at stage 1's operating point, QSD-seeded: the same
    two-state inversion with no stage 1 in the problem, so no death term can exist."""
    import scipy.sparse.linalg as spla
    mus, _ = el.operating_points(om, 2)
    pl = el.pi_low(om, mus[0])
    Q, m, _ = el._gen(om, mus[0], False)
    keep, w = el.qsd(om, mus[0])
    p = np.zeros(m); p[keep] = w
    p = spla.expm_multiply(Q.T * t, p)
    low = float(p[np.arange(m) < el.r2 * om].sum())
    return -math.log(max(1.0 - low / pl, 1e-300)) / t


def p4b(rep):
    """POST-HOC, DISCLOSED (rule 19). P4 as written -- "the corrected ratio is closer to 1 in
    every cell" -- printed FAILS. Read against the numbers it looks broken rather than the
    correction: at Omega = 30, where S1 ~ 1, BOTH extractions drift ~+6% between windows, so
    the two-state inversion has a window error of its own, and at Omega = 20 the uncorrected
    0.976 could be that +5% and a -7% survival loss cancelling by coincidence. "Closer to 1"
    cannot tell a correct extraction from two cancelling biases.

    That reading is itself a claim and gets an independent test rather than an argument: a
    control with the same stage-2 physics and NO stage 1. Criterion, written before the control
    was run: the corrected ratio must lie closer to the control's than the uncorrected ratio
    does, in every cell. The data that would print the opposite: control ratios near the
    uncorrected ones -- which is exactly the case in which dividing by S1 is wrong.
    """
    print("\nP4b -- POST-HOC control: stage 2 alone, upstream pinned, no stage 1 to die")
    e = Elem("published", *PUBLISHED_RAILS, cap_mult=1.25)
    rows = []
    for om in (14, 20, 30, 55):
        a, b = cell(e, om, 2.0), cell(e, om, 4.0)
        ctl = pinned_control(e, om, 4.0) / pinned_control(e, om, 2.0)
        ru, rc = b["k_uncorr"] / a["k_uncorr"], b["k_corr"] / a["k_corr"]
        closer = abs(rc - ctl) < abs(ru - ctl)
        rows.append({"omega": om, "control": ctl, "uncorr": ru, "corr": rc,
                     "S1_t2": a["S1"], "closer": bool(closer)})
        print(f"   Om={om:>3}  k(4)/k(2): control {ctl:.4f}   uncorrected {ru:.4f}"
              f"   corrected {rc:.4f}   S1(t=2) {a['S1']:.4f}   "
              f"{'corrected closer' if closer else 'UNCORRECTED closer'}", flush=True)
    ok = all(r["closer"] for r in rows if r["S1_t2"] < 0.999)
    print(f"   -> P4b {'HOLDS: wherever S1 < 1, dividing by it brings the rate to the no-death control' if ok else 'FAILS: the correction does not reproduce the no-death control'}")
    rep["p4b"] = {"rows": rows, "holds": bool(ok)}


def p8_window(rep):
    """POST-HOC PREDICTIONS, WRITTEN BEFORE THIS FUNCTION WAS RUN. P4b found that at Omega = 55,
    where S1 = 0.9999 and the survival correction is inert, the joint k_eff still drifts +10%
    between t = 2 and t = 4, against +2% for the no-stage-1 control. So k_eff at a fixed window
    is not a window-independent rate, and every position since §108 is built from k_eff at ONE
    window, T_WINDOW = 2.0, which no section swept. That is rule 13 exactly: an approximation's
    own numerical parameter is a second axis.

      P8a  the corrected position at Omega = 14 is near window-stable (survival removed, control
           drift ~+2%): predict it moves < 0.1 over t = 1..8.
      P8b  at Omega = 55 and 70 the position RISES with the window: +10% in k_eff per doubling
           is ~+0.13 in position per doubling, so predict +0.3 or more over t = 1..8.
      P8c  therefore the drift's span is window-dependent: predict it grows by > 0.2 from t = 1
           to t = 8. If instead the span is stable to < 0.1, the window is harmless for the
           drift's SIZE and only its absolute placement moves.
    """
    print("\nP8 -- the window was never swept: position vs t, survival-corrected, published element")
    e = Elem("published", *PUBLISHED_RAILS, cap_mult=1.25)
    ts = (1.0, 2.0, 4.0, 8.0)
    print(f"{'Om':>5}" + "".join(f"{'t='+str(t):>11}" for t in ts) + f"{'move':>9}")
    table = {}
    for om in OMEGAS_PUB:
        mus, _ = e.operating_points(om, 2)
        pl = e.pi_low(om, mus[0])
        km, kg, ka = e.candidate_averages(om)
        lb = math.log(ka / km)
        row = []
        for t in ts:
            v, s1 = e.joint_absorbing(om, t, return_survival=True)
            k = invert(v, s1, pl, t)
            row.append(math.log(k / km) / lb)
        table[om] = row
        print(f"{om:>5}" + "".join(f"{x:>11.4f}" for x in row) + f"{row[-1]-row[0]:>+9.4f}",
              flush=True)
    spans = [table[OMEGAS_PUB[-1]][i] - table[OMEGAS_PUB[0]][i] for i in range(len(ts))]
    print("   drift span (Om=70 minus Om=14) by window: "
          + "  ".join(f"t={t}: {sp:.4f}" for t, sp in zip(ts, spans)))
    rep["p8"] = {"windows": list(ts), "positions": {str(k): v for k, v in table.items()},
                 "spans": spans}


def _joint_alive(el, om):
    """The absorbing-upstream joint generator restricted to stage-1-ALIVE states (n1 >= saddle).
    Transitions out of the alive set become leakage on the diagonal: a proper sub-generator."""
    import scipy.sparse as sp
    cap = el.cap(om)
    m = cap + 1
    N = m * m
    idx = np.arange(N)
    n1, n2 = idx // m, idx % m
    sad = el.r2 * om
    rows, cols, vals = [], [], []
    diag = np.zeros(N)
    l1, u1 = el.rates(n1.astype(float), np.zeros(N), om, True)
    alive = n1 >= sad
    up1 = (n1 < cap) & (l1 > 0) & alive
    dn1 = (n1 > 0) & (u1 > 0) & alive
    rows.append(idx[up1]); cols.append(idx[up1] + m); vals.append(l1[up1])
    rows.append(idx[dn1]); cols.append(idx[dn1] - m); vals.append(u1[dn1])
    diag -= np.where(up1, l1, 0.0) + np.where(dn1, u1, 0.0)
    l2, u2 = el.rates(n2.astype(float), n1.astype(float), om, False)
    up2 = (n2 < cap) & (l2 > 0)
    dn2 = (n2 > 0) & (u2 > 0)
    rows.append(idx[up2]); cols.append(idx[up2] + 1); vals.append(l2[up2])
    rows.append(idx[dn2]); cols.append(idx[dn2] - 1); vals.append(u2[dn2])
    diag -= np.where(up2, l2, 0.0) + np.where(dn2, u2, 0.0)
    rows.append(idx); cols.append(idx); vals.append(diag)
    Q = sp.csr_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
    keep = np.where(alive)[0]
    return Q[keep][:, keep].tocsc()


def conditioned_gap(el, om, k=4):
    """lambda1 - lambda0 of the alive-restricted joint generator: stage 2's relaxation rate
    CONDITIONED on stage 1 never failing (the Q-process gap). No seed, no window, no S1."""
    import scipy.sparse.linalg as spla
    Qa = _joint_alive(el, om)
    ev = spla.eigs(Qa, k=k, sigma=0.0, which="LM", return_eigenvectors=False)
    lam = np.sort(-np.real(ev))
    return float(lam[0]), float(lam[1]), lam


def stage1_death(el, om):
    """Stage 1's own QSD decay rate on the SAME alive convention (n >= saddle)."""
    Q, m, _ = el._gen(om, 0.0, True)
    keep = np.where(np.arange(m) >= el.r2 * om)[0]
    A = Q[keep][:, keep].toarray()
    return float(np.sort(-np.real(np.linalg.eigvals(A)))[0])


def p9_exact(rep):
    """THE EXACT, WINDOW-FREE EFFECTIVE RATE. Predictions written before this was run.

    P8 found the position is a converging function of the window t, which no section swept, and
    rule 21 says expand the closed form before calling that a phenomenon: -ln(1 - q/pi_low)/t
    is window-free only if stage 2 is exactly two-state AND starts in the joint system's own
    quasi-stationary state. The quantity it is trying to measure has an exact definition that
    needs neither -- condition on stage 1 surviving (the Q-process) and stage 2's relaxation
    rate is lambda1 - lambda0 of the alive-restricted joint generator. It is the same kind of
    object as k(<x>) = escape_rate, which is itself a spectral gap.

      P9a  WIRING, exact: lambda0 must equal stage 1's own QSD decay rate (stage 2 does not feed
           back into stage 1), to eigensolver precision.
      P9b  the exact position approximates the t -> infinity limit of P8's survival-corrected
           sequence at Omega >= 30, where that sequence halves its increments cleanly: predict
           agreement with the Aitken limit to within 0.05. At Omega = 14 the two-state
           inversion's own window drift (+2.7% per doubling in the no-death control) makes
           that limit fuzzy, and no agreement is predicted there.
      P9c  THE DRIFT IS REAL. The survival term lowered the published curve's low end by ~0.30;
           the short window lowered its high end by ~0.3. Predict they roughly cancel, so the
           exact drift span sits within 0.15 of the published 1.176 -- and that the position
           still traverses the bracket. If instead the exact span is < 0.5, the arc's central
           "drift" was mostly instrument.
      P9d  §114.3's matched triple survives under the exact rate (outlier > 2x its pair).
      P9e  T-CASC-af's slope ratio under the exact rate: no prediction; reported.
    """
    print("\nP9 -- the exact window-free rate: lambda1 - lambda0 of the stage-1-alive joint generator")
    out = {"published": [], "rail": []}
    e = Elem("published", *PUBLISHED_RAILS, cap_mult=1.25)
    print(f"{'Om':>5}{'lam0':>12}{'stage1 QSD':>12}{'rel':>9}{'gap':>12}{'k_mean':>12}"
          f"{'pos exact':>11}{'pos t=2 pub':>13}")
    p8 = rep.get("p8", {}).get("positions", {})
    for om in OMEGAS_PUB:
        l0, l1, lam = conditioned_gap(e, om)
        d1 = stage1_death(e, om)
        km, kg, ka = e.candidate_averages(om)
        gap = l1 - l0
        pos = math.log(gap / km) / math.log(ka / km)
        seq = p8.get(str(om))
        aitken = float("nan")
        if seq and len(seq) >= 3:
            x0, x1, x2 = seq[-3:]
            den = (x2 - x1) - (x1 - x0)
            aitken = x2 - (x2 - x1) ** 2 / den if den != 0 else float("nan")
        out["published"].append({"omega": om, "lam0": l0, "lam1": l1, "stage1_death": d1,
                                 "gap": gap, "k_mean": km, "k_avg": ka, "pos_exact": pos,
                                 "aitken_p8": aitken, "spectrum": [float(x) for x in lam]})
        print(f"{om:>5}{l0:>12.4e}{d1:>12.4e}{abs(l0/d1-1):>9.1e}{gap:>12.4e}{km:>12.4e}"
              f"{pos:>11.4f}{'':>4}aitken {aitken:.4f}", flush=True)
    ps = [r["pos_exact"] for r in out["published"]]
    print(f"   exact drift span: {ps[-1]-ps[0]:.4f}   (published at t=2: 1.1756)")
    print(f"\n   rail cells, exact (cap 2.0):")
    rails = dict(ELEMENTS)
    for name, om in RAIL_CELLS:
        el = Elem(name, R1_FIX, rails[name], R3_FIX, cap_mult=CAP)
        l0, l1, lam = conditioned_gap(el, om)
        km, kg, ka = el.candidate_averages(om)
        pos = math.log((l1 - l0) / km) / math.log(ka / km)
        out["rail"].append({"cell": f"{name} Om={om}", "omega": om, "lam0": l0, "lam1": l1,
                            "gap": l1 - l0, "k_mean": km, "k_avg": ka, "pos_exact": pos,
                            "A_omega": el.action() * om,
                            "msq": el.margin_over_sigma(om) ** 2,
                            "pi_low": el.pi_low(om, el.operating_points(om, 2)[0][0])})
        print(f"   {name} Om={om:>3}  gap {l1-l0:.4e}  k_mean {km:.4e}  pos exact {pos:+.4f}",
              flush=True)
    rep["p9"] = out


def conditioned_qsd_low(el, om):
    """Stage 2's low-state occupancy in the Q-process's quasi-stationary law (left eigenvector
    of the alive-restricted generator at lambda0) -- the occupancy the conditioned system
    actually has, as opposed to `pi_low`, which pins the upstream at its mean."""
    import scipy.sparse.linalg as spla
    Qa = _joint_alive(el, om)
    _, v = spla.eigs(Qa.T.tocsc(), k=1, sigma=0.0, which="LM")
    q = np.abs(np.real(v[:, 0])); q /= q.sum()
    m = el.cap(om) + 1
    alive_n1 = np.arange(m)[np.arange(m) >= el.r2 * om]
    g = q.reshape(len(alive_n1), m)
    return float(g[:, np.arange(m) < el.r2 * om].sum())


def mode_check(el, om):
    """Is lambda1 stage 2's inter-well mode? Its right eigenvector should flip sign across
    stage 2's saddle. Returns |corr(eigvec, sign(n2 - saddle))| weighted by the QSD."""
    import scipy.sparse.linalg as spla
    Qa = _joint_alive(el, om)
    ev, vec = spla.eigs(Qa, k=3, sigma=0.0, which="LM")
    # smallest DECAY first. The first version had `[::-1]` here, which ordered decays
    # descending, so vec[:, 0] was a fast mode and the 'Doob-transformed' ratio divided the
    # lambda1 mode by it -- every correlation it printed (~0.05) was meaningless.
    order = np.argsort(-np.real(ev))
    lam = -np.real(ev[order]); vec = np.real(vec[:, order])
    m = el.cap(om) + 1
    alive_n1 = np.arange(m)[np.arange(m) >= el.r2 * om]
    n2 = np.tile(np.arange(m), len(alive_n1))
    sgn = np.where(n2 < el.r2 * om, -1.0, 1.0)
    _, lv = spla.eigs(Qa.T.tocsc(), k=1, sigma=0.0, which="LM")
    w = np.abs(np.real(lv[:, 0])); w /= w.sum()
    u = vec[:, 1] / vec[:, 0]                          # Doob-transformed mode
    mu_u, mu_s = (w * u).sum(), (w * sgn).sum()
    cov = (w * (u - mu_u) * (sgn - mu_s)).sum()
    corr = cov / np.sqrt((w * (u - mu_u) ** 2).sum() * (w * (sgn - mu_s) ** 2).sum())
    return float(abs(corr)), [float(x) for x in lam]


def p10_verify(rep):
    """VERIFYING A RETRACTION (rule 14). Predictions written before this was run.

    P9 put §114.3's 3.2x outlier (E-084, Omega = 14) at 0.493 under the exact gap, against 0.570
    and 0.551 for its bracketing pair -- which would make §114's refutation of A*Omega's
    sufficiency, its 'second variable' pi_low, and T-CASC-ac artifacts of the windowed inversion.
    A retraction of a refutation is a claim and gets checked three independent ways.

      V1  ONE BOX. P9 compared rails at cap 2.0 with the published element at cap 1.25. Redo the
          published element at 2.0: predict exact positions move < 0.03 and the triple's
          verdict -- E-084 within 0.1 of its pair -- is unchanged.
      V2  THE RIGHT MODE. lambda1 must be stage 2's inter-well relaxation, not a stage-1 mode
          or a stage-2 intra-well one: predict the Doob-transformed eigenvector flips sign across
          stage 2's saddle, |corr| > 0.9, in every cell tested, with lambda2 well separated.
      V3  AN INDEPENDENT ROUTE, and the mechanism's kill test. Suspect (rule 17): for k*t << 1
          the windowed inversion measures k_fwd,cond / pi_low,PINNED while the gap is
          k_fwd,cond / pi_low,COND, so the two differ by the occupancy ratio. Predict:
            (a) pi_low,cond differs from pinned pi_low, and the ratio is LARGEST for E-084
                (pinned 0.244 -- the cell that looked anomalous);
            (b) redoing the survival-corrected windowed inversion at t = 8 with pi_low,cond in
                place of the pinned value lands within 0.1 position units of the exact gap at
                Omega = 55 and 70, where k*t << 1 -- whereas with the pinned value it is 0.2-0.45
                off (P8 vs P9).
          If (b) fails, the exact gap and the windowed inversion disagree for a reason this does
          not name, and the retraction is NOT supported.
    """
    print("\nP10 -- verifying the retraction three ways (rule 14)")
    rails = dict(ELEMENTS)
    out = {}

    print("   V1 -- published element at cap 2.0, exact:")
    e2 = Elem("published", *PUBLISHED_RAILS, cap_mult=CAP)
    v1 = []
    old = {r["omega"]: r["pos_exact"] for r in rep["p9"]["published"]}
    for om in OMEGAS_PUB:
        l0, l1, _ = conditioned_gap(e2, om)
        km, kg, ka = e2.candidate_averages(om)
        pos = math.log((l1 - l0) / km) / math.log(ka / km)
        pos_fast = math.log(kg / km) / math.log(ka / km)
        v1.append({"omega": om, "pos_exact_cap2": pos, "pos_exact_cap125": old[om],
                   "A_omega": e2.action() * om, "pos_fast_limit": pos_fast})
        print(f"     Om={om:>3}  cap2.0 {pos:+.4f}   cap1.25 {old[om]:+.4f}"
              f"   (fast/geometric limit sits at {pos_fast:+.4f})", flush=True)
    out["v1"] = v1
    by = {r["cell"]: r["pos_exact"] for r in rep["p9"]["rail"]}
    pub30 = next(r for r in v1 if r["omega"] == 30)["pos_exact_cap2"]
    pair = (by["E-070 Om=30"], pub30)
    gap_ = abs(by["E-084 Om=14"] - sum(pair) / 2)
    print(f"     triple at one box: E-070 Om=30 {pair[0]:+.4f}  E-084 Om=14 "
          f"{by['E-084 Om=14']:+.4f}  PUB Om=30 {pair[1]:+.4f}   |outlier - pair mean| = {gap_:.4f}")
    out["v1_triple_gap"] = gap_

    print("   V2 -- is lambda1 stage 2's inter-well mode?")
    v2 = []
    for label, el, om in (("PUB", e2, 14), ("PUB", e2, 70),
                          ("E-084", Elem("E-084", R1_FIX, rails["E-084"], R3_FIX, cap_mult=CAP), 14),
                          ("E-070", Elem("E-070", R1_FIX, rails["E-070"], R3_FIX, cap_mult=CAP), 30)):
        corr, lam = mode_check(el, om)
        v2.append({"cell": f"{label} Om={om}", "corr": corr, "lam": lam})
        print(f"     {label} Om={om:>3}: |corr with sign(n2 - saddle)| = {corr:.4f}"
              f"   lam0..2 = {lam[0]:.3e}, {lam[1]:.3e}, {lam[2]:.3e}", flush=True)
    out["v2"] = v2

    print("   V3 -- conditioned occupancy, and the windowed inversion redone with it (t = 8):")
    v3 = []
    cells = [("PUB", e2, 30), ("PUB", e2, 55), ("PUB", e2, 70),
             ("E-084", Elem("E-084", R1_FIX, rails["E-084"], R3_FIX, cap_mult=CAP), 14),
             ("E-070", Elem("E-070", R1_FIX, rails["E-070"], R3_FIX, cap_mult=CAP), 30)]
    for label, el, om in cells:
        mus, _ = el.operating_points(om, 2)
        pin = el.pi_low(om, mus[0])
        cond = conditioned_qsd_low(el, om)
        km, kg, ka = el.candidate_averages(om)
        lb = math.log(ka / km)
        l0, l1, _ = conditioned_gap(el, om)
        exact = math.log((l1 - l0) / km) / lb
        v, s1 = el.joint_absorbing(om, 8.0, return_survival=True)
        w_pin = math.log(invert(v, s1, pin, 8.0) / km) / lb
        w_cond = math.log(invert(v, s1, cond, 8.0) / km) / lb
        v3.append({"cell": f"{label} Om={om}", "pi_pinned": pin, "pi_cond": cond,
                   "ratio": cond / pin, "exact": exact, "windowed_pinned": w_pin,
                   "windowed_cond": w_cond})
        print(f"     {label} Om={om:>3}: pi_low pinned {pin:.4f} cond {cond:.4f} (x{cond/pin:.3f})"
              f"   exact {exact:+.4f}   windowed: pinned {w_pin:+.4f}, cond {w_cond:+.4f}",
              flush=True)
    out["v3"] = v3
    rep["p10"] = out


def p11_mode_and_box(rep):
    """Two checks P10 left open. Predictions written before this was run.

      V2' THE RIGHT MODE, rerun with the sort fixed: the Doob-transformed lambda1 eigenvector
          must flip sign across stage 2's saddle, |corr| > 0.9, in every cell.
      V4  BOX CONVERGENCE of the exact position. P10's V1 found cap 1.25 -> 2.0 moves it 0.13 at
          Omega = 14, far more than predicted. §113 found escape rates converged by cap 2.0 at
          every Omega tested, so predict cap 2.0 -> 3.0 moves every exact position < 0.02 and
          leaves the matched triple's verdict (no 3x outlier) unchanged.
    """
    print("\nP11 -- the mode check (sort fixed) and box convergence of the exact position")
    rails = dict(ELEMENTS)
    mk = lambda name, cap: (Elem("published", *PUBLISHED_RAILS, cap_mult=cap) if name == "PUB"
                            else Elem(name, R1_FIX, rails[name], R3_FIX, cap_mult=cap))
    cells = (("PUB", 14), ("PUB", 30), ("E-084", 14), ("E-070", 30))
    out = {"mode": [], "box": []}
    for name, om in cells:
        corr, lam = mode_check(mk(name, 2.0), om)
        out["mode"].append({"cell": f"{name} Om={om}", "corr": corr, "lam": lam})
        print(f"   V2' {name} Om={om:>3}: |corr| = {corr:.4f}   decays {lam[0]:.3e}, "
              f"{lam[1]:.3e}, {lam[2]:.3e}", flush=True)
    for name, om in cells:
        pos = {}
        for cap in (2.0, 3.0):
            el = mk(name, cap)
            l0, l1, _ = conditioned_gap(el, om)
            km, kg, ka = el.candidate_averages(om)
            pos[cap] = math.log((l1 - l0) / km) / math.log(ka / km)
        out["box"].append({"cell": f"{name} Om={om}", "cap2": pos[2.0], "cap3": pos[3.0]})
        print(f"   V4  {name} Om={om:>3}: exact position cap2.0 {pos[2.0]:+.4f}   cap3.0 "
              f"{pos[3.0]:+.4f}   move {pos[3.0]-pos[2.0]:+.4f}", flush=True)
    rep["p11"] = out


def p12_redo_114(rep):
    """§114's analyses redone on the exact rate, ONE box (cap 2.0), stored gaps only.

    NOT A CLEAN PREDICTION, AND SAID SO. Reading P9/P10's numbers I estimated by hand that under
    the exact rate E-084's residual about the Omega-curve falls from 92% of the span to ~20%, and
    that T-CASC-af's rail/Omega slope ratio moves from 0.73 to ~1.1. Those estimates were made
    from the data this computes, so agreement here is arithmetic, not confirmation. What this
    adds is the full rule-15 table for every candidate, computed rather than eyeballed.
    """
    print("\nP12 -- §114's collapse and slope tests, redone on the exact rate (cap 2.0)")
    e2 = Elem("published", *PUBLISHED_RAILS, cap_mult=CAP)
    omg = []
    for r in rep["p10"]["v1"]:
        mus, _ = e2.operating_points(r["omega"], 2)
        omg.append({"cell": f"PUB Om={r['omega']}", "position": r["pos_exact_cap2"],
                    "A_omega": r["A_omega"], "msq": e2.margin_over_sigma(r["omega"]) ** 2,
                    "pi_low": e2.pi_low(r["omega"], mus[0])})
    rail = [{"cell": r["cell"], "position": r["pos_exact"], "A_omega": r["A_omega"],
             "msq": r["msq"], "pi_low": r["pi_low"]} for r in rep["p9"]["rail"]]
    hi = [r for r in rail if r["pi_low"] > 0.8]
    out = {"candidates": {}, "residuals": []}
    ox = np.array([r["A_omega"] for r in omg]); oy = np.array([r["position"] for r in omg])
    span = float(oy.max() - oy.min())
    print(f"   exact Omega-curve spans {span:.4f} (the windowed published curve spanned 1.1756)")
    print(f"   {'cell':>13}{'A*Omega':>9}{'exact':>9}{'curve':>9}{'resid':>9}{'% span':>8}")
    for r in rail:
        c = float(np.interp(r["A_omega"], ox, oy))
        out["residuals"].append({**r, "curve": c, "resid": r["position"] - c,
                                 "frac_span": (r["position"] - c) / span})
        print(f"   {r['cell']:>13}{r['A_omega']:>9.3f}{r['position']:>9.4f}{c:>9.4f}"
              f"{r['position']-c:>+9.4f}{(r['position']-c)/span*100:>7.1f}%")
    print(f"\n   {'candidate':>18}{'RMS/span (hi-occ)':>19}{'slope ratio':>13}")
    for key, lab in (("A_omega", "A*Omega"), ("msq", "(margin/sigma)^2")):
        sr, rs = slope_ratio([{**r, "p": r["position"]} for r in hi],
                             [{**r, "p": r["position"]} for r in omg], key, "p")
        out["candidates"][lab] = {"slope_ratio": sr, "rms_over_span": rs}
        print(f"   {lab:>18}{rs:>19.3f}{sr:>13.3f}")
    rep["p12"] = out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/the_survival_term.json"))
    ap.add_argument("--p4b-only", action="store_true",
                    help="run only the post-hoc control and merge it into the stored JSON")
    ap.add_argument("--p12-only", action="store_true",
                    help="run only the §114 redo and merge it into the stored JSON")
    ap.add_argument("--p11-only", action="store_true",
                    help="run only the mode/box checks and merge them into the stored JSON")
    ap.add_argument("--p10-only", action="store_true",
                    help="run only the verification of P9 and merge it into the stored JSON")
    ap.add_argument("--p9-only", action="store_true",
                    help="run only the exact-gap section and merge it into the stored JSON")
    ap.add_argument("--p8-only", action="store_true",
                    help="run only the window sweep and merge it into the stored JSON")
    args = ap.parse_args()
    if args.p12_only:
        rep = json.loads(args.out.read_text())
        p12_redo_114(rep)
        args.out.write_text(json.dumps(rep, indent=2, default=float))
        return
    if args.p11_only:
        rep = json.loads(args.out.read_text())
        p11_mode_and_box(rep)
        args.out.write_text(json.dumps(rep, indent=2, default=float))
        return
    if args.p10_only:
        rep = json.loads(args.out.read_text())
        p10_verify(rep)
        args.out.write_text(json.dumps(rep, indent=2, default=float))
        return
    if args.p9_only:
        rep = json.loads(args.out.read_text()) if args.out.exists() else {}
        p9_exact(rep)
        args.out.write_text(json.dumps(rep, indent=2, default=float))
        return
    if args.p8_only:
        rep = json.loads(args.out.read_text()) if args.out.exists() else {}
        p8_window(rep)
        args.out.write_text(json.dumps(rep, indent=2, default=float))
        return
    if args.p4b_only:
        rep = json.loads(args.out.read_text()) if args.out.exists() else {}
        p4b(rep)
        args.out.write_text(json.dumps(rep, indent=2, default=float))
        return
    rep = {}

    print("P1 -- wiring: survival ignored must reproduce §109 (published rails, cap 1.25)")
    pub109 = {r["omega"]: r for r in
              json.loads(pathlib.Path("results/it_was_the_seed.json").read_text())["p3"]}
    e125 = Elem("published", *PUBLISHED_RAILS, cap_mult=1.25)
    p1 = []
    for om in OMEGAS_PUB:
        c = cell(e125, om)
        rel = abs((c["k_uncorr"] / c["k_mean"]) / pub109[om]["qsd"] - 1)
        p1.append({**c, "rel_vs_109": rel})
        print(f"   Om={om:>3}  k_uncorr/k_mean = {c['k_uncorr']/c['k_mean']:.8f}"
              f"   §109 {pub109[om]['qsd']:.8f}   rel {rel:.1e}   S1 = {c['S1']:.4f}", flush=True)
    rep["p1_published_cap125"] = p1
    ok = max(r["rel_vs_109"] for r in p1) < 1e-5
    print(f"   -> P1 {'HOLDS' if ok else 'FAILS'}")
    if not ok:
        args.out.write_text(json.dumps(rep, indent=2, default=float))
        return

    print("\nP2/P3 -- S1 and the shift, published element (cap 1.25, as published)")
    print(f"{'Om':>5}{'S1':>9}{'exp(-l1 t)':>12}{'pos (pub)':>11}{'pos (corr)':>12}{'shift':>9}")
    for r in p1:
        l1 = e125.escape_rate(r["omega"], 0.0, first=True)
        r["S1_model"] = math.exp(-l1 * T_WINDOW)
        r["shift"] = r["pos_corr"] - r["pos_uncorr"]
        print(f"{r['omega']:>5}{r['S1']:>9.4f}{r['S1_model']:>12.4f}{r['pos_uncorr']:>11.4f}"
              f"{r['pos_corr']:>12.4f}{r['shift']:>9.4f}")
    span0 = p1[-1]["pos_uncorr"] - p1[0]["pos_uncorr"]
    span1 = p1[-1]["pos_corr"] - p1[0]["pos_corr"]
    print(f"   drift span: published {span0:.4f} -> corrected {span1:.4f}"
          f"   ({(1 - span1/span0)*100:.1f}% of the drift was the survival term)")
    rep["p3"] = {"span_published": span0, "span_corrected": span1}

    print("\nP4 -- kill test for the correction: is the corrected rate window-stable?")
    p4 = []
    for om in (14, 20, 30):
        a, b = cell(e125, om, 2.0), cell(e125, om, 4.0)
        ru, rc = b["k_uncorr"] / a["k_uncorr"], b["k_corr"] / a["k_corr"]
        p4.append({"omega": om, "ratio_uncorr": ru, "ratio_corr": rc,
                   "S1_t2": a["S1"], "S1_t4": b["S1"]})
        print(f"   Om={om:>3}  k(t=4)/k(t=2): uncorrected {ru:.4f}   corrected {rc:.4f}"
              f"   (S1: {a['S1']:.4f} -> {b['S1']:.4f})", flush=True)
    better = all(abs(r["ratio_corr"] - 1) < abs(r["ratio_uncorr"] - 1) for r in p4)
    print(f"   -> P4 {'HOLDS: the corrected extraction is more window-stable in every cell' if better else 'FAILS: dividing by S1 does not stabilise the rate -- P3 and P5 are unreadable'}")
    rep["p4"] = {"rows": p4, "corrected_more_stable": bool(better)}
    p4b(rep)

    print(f"\nP5/P6 -- the collapse, at ONE box (cap_mult = {CAP}) for both axes")
    ecap = Elem("published", *PUBLISHED_RAILS, cap_mult=CAP)
    omg = []
    for om in OMEGAS_PUB:
        omg.append({**cell(ecap, om), "cell": f"PUB Om={om}"})
        print(f"   PUB   Om={om:>3}  A*Om={omg[-1]['A_omega']:7.3f}  S1={omg[-1]['S1']:.4f}"
              f"  pos {omg[-1]['pos_uncorr']:+.4f} -> {omg[-1]['pos_corr']:+.4f}", flush=True)
    rails = dict(ELEMENTS)
    rail = []
    for name, om in RAIL_CELLS:
        el = Elem(name, R1_FIX, rails[name], R3_FIX, cap_mult=CAP)
        rail.append({**cell(el, om), "cell": f"{name} Om={om}"})
        print(f"   {name} Om={om:>3}  A*Om={rail[-1]['A_omega']:7.3f}  S1={rail[-1]['S1']:.4f}"
              f"  pos {rail[-1]['pos_uncorr']:+.4f} -> {rail[-1]['pos_corr']:+.4f}", flush=True)
    rep["omega_sweep"], rep["rail"] = omg, rail
    hi = [r for r in rail if r["pi_low"] > 0.8]

    print("\n   T-CASC-af -- slope ratio (rail / Omega) and RMS/span, high-occupancy cells:")
    p5 = {}
    for xkey, lab in (("A_omega", "A*Omega"), ("msq", "(margin/sigma)^2")):
        s0, r0 = slope_ratio(hi, omg, xkey, "pos_uncorr")
        s1, r1 = slope_ratio(hi, omg, xkey, "pos_corr")
        p5[lab] = {"slope_uncorr": s0, "slope_corr": s1, "rms_uncorr": r0, "rms_corr": r1}
        print(f"   {lab:>18}: slope ratio {s0:.3f} -> {s1:.3f}   RMS/span {r0:.3f} -> {r1:.3f}")
    toward = all(abs(v["slope_corr"] - 1) < abs(v["slope_uncorr"] - 1) for v in p5.values())
    print(f"   -> P5 {'HOLDS: the correction moves both slope ratios toward 1' if toward else 'FAILS: the survival bias is not what makes the rail axis shallower'}")
    rep["p5"] = p5

    by = {r["cell"]: r for r in rail + omg}
    lo, mid, hi_ = by["E-070 Om=30"], by["E-084 Om=14"], by["PUB Om=30"]
    print("\n   §114.3's matched triple under the correction:")
    for r in (lo, mid, hi_):
        print(f"     {r['cell']:>12}  A*Om={r['A_omega']:.3f}  pos {r['pos_uncorr']:+.4f}"
              f" -> {r['pos_corr']:+.4f}")
    ratio = mid["pos_corr"] / max(lo["pos_corr"], hi_["pos_corr"])
    print(f"   outlier / bracketing pair, corrected: {ratio:.2f}x"
          f"   -> P6 {'HOLDS: A*Omega is still refuted as sufficient' if ratio > 2 else 'FAILS: the refutation was partly the survival term'}")
    rep["p6"] = {"ratio_corrected": ratio}
    p8_window(rep)

    args.out.write_text(json.dumps(rep, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
