"""§133 -- where the committed falls cross Sigma, not stage 2's saddle line: a reading check on E24-E27.

E24-E27 (§124.12, §124.14) measured the committed never-dipped fall flux where stage 2 crosses its SADDLE LINE
x2 = s2 = 0.606, binned by stage 1's position x1 there. THEORIES then read that profile as "the region along
Sigma that carries the crossing flux". Sigma, the joint separatrix, coincides with the saddle line only near the
joint saddle (3, 0.606). Elsewhere it sits higher (x2 = 0.63 at x1 = 2.0, 1.14 at x1 = 1.18). A trajectory just
below Sigma does not fall straight down. It shadows Sigma toward the joint saddle and leaves along the saddle's
unstable manifold, by a distance set by how far below Sigma it is, which is O(Omega^-1/2). So the saddle-line
profile could drift toward x1 = 3 with Omega from geometry alone. If so, E24's "spread, drifting toward the
saddle" and §124.14's "slow localisation" suspect are partly readings of the wrong surface.

INDEPENDENT INPUT (computed before this script, Omega-free). The action along Sigma, eps(x1) = W_Sigma(x1) - W*,
from the ray family (`joint_action/sigma_profile.py`, lower envelope across a caustic near the saddle). It falls
monotonically toward the saddle: 7.6e-4 (x1 = 1.22), 1.5e-4 (1.6), 9e-5 (2.0), 5e-5 (2.5), 2.7e-5 (2.88),
1.0e-5 (2.99). W along Sigma is not flat; its slope is ~1e-4 per unit x1 in the interior.

METHOD. The same construction as `episode_committor.fall_profile`, for a general cut. On the Q-process stationary
law (stage 1 conditioned never to dip, stage 2 leaking at its low rail), the committed crossing of a cut is its
LAST downward crossing: the rate of each cut edge i -> j times the committor q(j) = P(reach the low rail before
re-crossing the cut upward), in the same measure as fall_profile. Two cuts: the saddle line (n2 = sad2, as
E24-E27), and Sigma's lattice image (n2 >= ceil(Sigma(x1) Omega) is "above"; both stage-2 and stage-1 moves can
cross it).

PREDICTIONS, WRITTEN BEFORE RUNNING (XK).
  S0  INSTRUMENT. With the cut at the saddle line, the generic code reproduces E26's stored flux_density at
      Omega = 200 and 400 to 1e-8 relative in every bin with flux above 1e-6 of the maximum.
  S1  CONSERVATION. The totals through the two cuts agree to 1e-6 relative at every Omega (each equals the
      committed fall rate).
  S2  THE READING CHECK. The Sigma-cut distribution's flux-weighted mean x1 is LOWER than the saddle-line
      distribution's by at least 0.2 at Omega = 400 and 800. Agreement within 0.05 would mean E24-E27's reading as
      "along Sigma" stands.
  S3  THE ACTION PROFILE, ABSOLUTELY (rule 16; needs Omega = 400, 800, 1600). Both windows have the same ln ratio
      (ln 2), so the prefactor power cancels in the second difference
          D2(x1) = ln F1600 - 2 ln F800 + ln F400   (Sigma-cut density per bin),
      leaving D2(x1) - D2(x_ref) = -400 (eps(x1) - eps(x_ref)) with x_ref = 2.5, and no free parameter. Predicted
      within 0.05 + 0.25 |prediction| in every bin of x1 in [1.2, 2.9] with flux above 1e-4 of the maximum. The
      bins where eps is large (x1 < 1.4: predicted -0.12 to -0.30) are the ones that can fail it.

OUTCOMES (appended after running; predictions above left as written).
  S0 HELD: the saddle-line cut reproduces E26's stored density to 3e-13 (Omega = 200, 400).
  S1 HELD: the two totals agree to 1e-13 at Omega = 200, 400, 800, 1600. The laptop (PARDISO) reproduces Omega = 400
     to 3e-13.
  S2 FAILED AS WRITTEN, and not the other way either. The mean-x1 gap (saddle line minus Sigma) is 0.060 / 0.102 /
     0.158 / 0.197 at Omega = 200 / 400 / 800 / 1600: under the 0.2 bar at 400 and 800, and far above the 0.05 at
     which E24's reading would stand. The criterion had no verdict for this middle. The reading: about a third
     of the saddle line's mean drift (+0.41 over 200 -> 1600) is geometry. On Sigma the mean moves +0.27, and the
     spread GROWS (sd 0.670 / 0.660 / 0.687 / 0.706), where on the saddle line it shrank (0.649 -> 0.577).
  S3 FAILED AS WRITTEN: 12 of 34 fine bins fall outside the tolerance. The failures are bin-to-bin alternation of
     about +-0.1 in the interior, where the predicted signal is < 0.05: lattice aliasing of the staircase cut,
     with ~20 columns per bin at Omega = 400. Where eps is large the fine bins track the prediction (1.225-1.375:
     -0.42 / -0.21 / -0.09 / -0.07 against -0.36 / -0.20 / -0.13 / -0.10). POST HOC, labelled as such: on 0.10-wide
     bins the interior (1.35-2.55) agrees within 0.03, and on 0.20-wide bins (1.5-2.3) within 0.012. Near the
     saddle (2.65-2.85) it is off by ~0.1, in the caustic region where the ray envelope is least certain.
     FINDINGS §133.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp

import experiments.episode_committor as ec
from experiments.weak_upstream import ARMS, element


_QP_CACHE = {}


def sigma_table(A, B, om_grid=None):
    xs = np.concatenate([np.linspace(A.r2 + 1e-3, 1.4, 80), np.linspace(1.4, 3.6, 120)[1:]])
    sg = ec.joint_separatrix(A, B, xs)
    return xs, sg


def cut_flux(A, B, om, cut="sad2", sig=None):
    m = A.cap(om) + 1
    sad1, sad2 = int(np.ceil(A.r2 * om)), int(np.ceil(B.r2 * om))
    key = (A.name, B.name, om, m, ec.SOLVER, ec.PIN_SCALE)
    if _QP_CACHE.get("key") != key:                 # both cuts share the stationary law: solve it once
        ec.qprocess_J(A, B, om, m, sad1, sad2, lower="rail")
        _QP_CACHE.clear()
        _QP_CACHE.update(key=key, last=ec.qprocess_J.last, ratio=ec.qprocess_J.last_ratio)
    L = _QP_CACHE["last"]
    P, n1s, base2, psi_n1, M1, M2 = L["P"], L["n1s"], L["base2"], L["psi_n1"], L["M1"], L["M2"]
    lt1, mt1 = L["lt1"].reshape(M1, -1)[:, 0], L["mt1"].reshape(M1, -1)[:, 0]
    n2s = base2 + np.arange(M2)
    x1 = n1s / om
    if cut == "sad2":
        n2cut = np.full(M1, sad2)
    else:
        xs, sg = sig
        sv = np.interp(x1, xs, np.where(np.isfinite(sg), sg, 99.0))
        sv[x1 < xs[np.isfinite(sg)][0]] = 99.0                   # everything falls there: the cut is above the box
        n2cut = np.clip(np.ceil(sv * om).astype(int), base2 + 1, base2 + M2)
    below = n2s[None, :] < n2cut[:, None]                          # region D (q unknown); q = 0 above
    I1, I2 = np.nonzero(below)
    idx = -np.ones((M1, M2), dtype=np.int64)
    idx[I1, I2] = np.arange(len(I1))
    n = len(I1)
    lam2, mu2 = B.rates(n2s[I2].astype(float), n1s[I1].astype(float), om, False)
    rows, cols, vals = [], [], []
    diag = np.zeros(n)
    rhs = np.zeros(n)
    for step, rate in ((+1, lt1[I1]), (-1, mt1[I1])):            # stage-1 moves
        J1 = I1 + step
        ok = (J1 >= 0) & (J1 < M1)
        tgt = np.where(ok, idx[np.clip(J1, 0, M1 - 1), I2], -1)
        inD = ok & (tgt >= 0)
        rows.append(np.nonzero(inD)[0]); cols.append(tgt[inD]); vals.append(rate[inD])
        diag -= np.where(ok, rate, 0.0)                            # to above-cut (q = 0) or within D
    up = I2 + 1
    tgt = np.where(up < M2, idx[I1, np.clip(up, 0, M2 - 1)], -1)
    inD = tgt >= 0
    rows.append(np.nonzero(inD)[0]); cols.append(tgt[inD]); vals.append(lam2[inD])
    diag -= np.where(up < M2, lam2, 0.0)
    dn = I2 - 1
    tgt = np.where(dn >= 0, idx[I1, np.clip(dn, 0, M2 - 1)], -1)
    inD = tgt >= 0
    rows.append(np.nonzero(inD)[0]); cols.append(tgt[inD]); vals.append(mu2[inD])
    diag -= mu2
    rhs -= np.where(dn < 0, mu2 / psi_n1[I1], 0.0)                 # the low rail: u = 1/psi (original measure)
    rows.append(np.arange(n)); cols.append(np.arange(n)); vals.append(diag)
    Mx = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsc()
    qv = ec.lu_solve(Mx, rhs, coords=(I1, I2))[0]
    q = np.zeros((M1, M2))
    q[I1, I2] = qv
    Pc = np.clip(P, 0.0, None)
    flux = np.zeros(M1)                                            # by stage 1's position after the crossing
    # stage-2 down-moves across the cut: (n1, n2cut) -> (n1, n2cut - 1)
    for i in range(M1):
        k = n2cut[i] - base2
        if 1 <= k < M2:
            _, mu_s = B.rates(np.array([float(n2s[k])]), np.array([float(n1s[i])]), om, False)
            flux[i] += Pc[i, k] * mu_s[0] * q[i, k - 1]
    # stage-1 moves across the cut: from above at n1 to below at n1 +- 1
    for i in range(M1):
        for step, rate in ((+1, lt1[i]), (-1, mt1[i])):
            j = i + step
            if 0 <= j < M1 and n2cut[j] > n2cut[i]:
                ks = np.arange(n2cut[i] - base2, min(n2cut[j] - base2, M2))
                ks = ks[ks >= 0]
                flux[j] += float((Pc[i, ks] * rate * q[j, ks]).sum())
    dens = [float(flux[(x1 >= 1.0 + 0.05 * k) & (x1 < 1.05 + 0.05 * k)].sum() / L["dip"] / 0.05) for k in range(44)]
    tot = float(flux.sum() / L["dip"])
    return {"omega": om, "cut": cut, "total": tot, "mean_x1": float((flux * x1).sum() / flux.sum()),
            "sd_x1": float(np.sqrt(max((flux * x1 ** 2).sum() / flux.sum() - ((flux * x1).sum() / flux.sum()) ** 2, 0))),
            "flux_bins_x1": [round(1.025 + 0.05 * k, 4) for k in range(44)], "flux_density": dens}


def main():
    A, B = element(ARMS["XK"][0]), element(ARMS["XK"][1])
    out = Path("results/sigma_flux.json")
    rep = json.loads(out.read_text()) if out.exists() else {}
    if len(sys.argv) > 2:
        ec.SOLVER = sys.argv[2]
    sig = None
    for om in [int(x) for x in sys.argv[1].split(",")]:
        for cut in ("sad2", "sigma"):
            from experiments.weak_upstream import CAP_MULT
            key = f"XK:{om}:{cut}" + ("" if CAP_MULT == 1.25 else f":cap{CAP_MULT:g}")
            if key in rep:
                continue
            if cut == "sigma" and sig is None:
                sig = sigma_table(A, B)
                rep["sigma_table"] = {"x1": list(map(float, sig[0])), "x2": [float(v) for v in sig[1]]}
            t = time.time()
            r = cut_flux(A, B, om, cut, sig)
            r["seconds"] = time.time() - t
            rep[key] = r
            print(key, {k: r[k] for k in ("total", "mean_x1", "sd_x1", "seconds")}, flush=True)
            out.write_text(json.dumps(rep, indent=2, default=float))


if __name__ == "__main__":
    main()
