"""§123 -- the exact p_d as a joint-chain committor: a route past the window to Omega > Omega_c.

§122's universal-floor suspect (T-CASC-ap) needs an EXACT p_d in a narrow-zone pair beyond its
crossover Omega_c = lam(s1)/(ls w^2) -- 227 for Z -- where evolving §119's whole window would need
~10^7 tagged states. A dip is an episode: it starts when stage 1 steps below its saddle, and ends
when stage 1 commits -- to its high rail or its low one. So, exactly on the joint chain,

    p_d_ep = sum_n2 J(n2) h(sad1-1, n2) / sum_n2 J(n2) c(sad1-1, n2),

  J(n2)    the stage-2 law AT the dip: the crossing flux sad1 -> sad1-1 out of the quasi-stationary
           law of the joint chain with stage 1 alive (n1 >= sad1) -- so the pre-dip sag is in it;
  c        P(stage 1 reaches its high rail before its low rail), from that state, on the joint chain;
  h        P(the same, AND stage 2 below its saddle at that moment).

c and h are each one sparse linear solve. Nothing is evolved, and no window enters.

PREDICTIONS, WRITTEN BEFORE RUNNING.

  E1  VALIDATION. At Omega = 40 and 55, p_d_ep reproduces §119's stored windowed exact p_d for X
      and Z within 10%. (They differ in definition: the window counts trajectories high at t = 16,
      which may dip again; the episode ends at commitment. If they disagree by more, the episode is
      not the window's quantity and E2 cannot be read against §122.)
  E2  THE KILL TEST (T-CASC-ap), conditional on E1. Z at Omega = 230 and 400 (past Omega_c = 227):
      the exact p_d_ep RISES from its Omega = 55 value, turning the decline that §119 measured
      (0.066 -> 0.054). If it keeps falling past Omega_c, the universal-floor suspect is refuted in
      the exact chain, whatever the reductions said.
  E3  phi_ep = sum J c / sum J -- the episode's upstream factor -- equals the stage-1 committor to
      within 2%: stage 2 does not feed back on stage 1, so it must.

  E4  Written after Z:800 returned p_d_ep = 0.44 with lam0 = -5e-15 and J_low_frac = 0.30 (against
      3e-8 at 400) -- the alive chain has two near-absorbing pockets (stage 2 high / low), and once
      both decay rates are below double precision the sigma = 0 eigenvector is an arbitrary mixture.
      Replacement with no eigenproblem: J from the OCCUPATION before the first dip, starting at the
      healthy product seed, x = (-Q_alive^T)^{-1} p_seed -- exactly stage 2's law at the first dip.
      Predict it reproduces the eigenvector-J p_d_ep within 1% at Omega <= 230, and returns
      J_low_frac at 800 to the ~1e-8 scale. Z:800 is read for E2 only from this J.

  E5  Written after E4's occupation solve was itself found broken at Omega >= 200 (relative
      residual ~1e2, total occupation FALLING from 200 to 400 where 1/k1 grows, J's shape jumping):
      stage 1 is autonomous, so the alive chain has TWO eigenvalues within a whisker of zero (stage 1
      dying, and that plus stage 2's own slow switching), and any solve near zero mixes them. All
      cells with lam0 <~ 1e-12 -- X at Omega >= 200, Z at Omega >= 400 -- are withdrawn, including
      E2's "still falling at 800" and X's non-monotone 0.596 / 0.538 / 0.674.
      Fix: make stage-2-low absorbing too. The chain alive AND stage 2 high has ONE tiny eigenvalue,
      O(1)-separated from the rest, so shift-invert is well-conditioned. It omits the fraction of
      stage 2 already low at the dip (<= 1.5% wherever it was reliably measured), carried as an
      error bar. Predict: matches the validated values within 1.5% at X <= 100, Z <= 230; and its
      own residual is <= 1e-8 at every Omega it is run at.

OUTCOMES (appended after running; every prediction above left as written). E1 held (window vs
episode within 6%). E3 held (phi_ep = stage-1 committor to 4 digits). E4's first half held at
Omega >= 55 (1.3% at 40); its second half was moot -- the occupation route itself broke at
lam0 <~ 1e-12. E5 held (within 1.2% of every validated cell; residuals <= 1.4e-10 at every Omega).
E2 FAILED as written: exact Z keeps falling through Omega = 400 = 1.76 Omega_c (0.065 -> 0.032);
Z:800 was not computed (stopped for system memory). Exact X rises monotonically to 0.807 at
Omega = 800 and phi * p_d to 0.380. FINDINGS §123.

  E6  Written before Z:800 ran (on the user's 64 GB laptop, after it reproduced Z:230/hi to 1.1e-10
      relative). T-CASC-ap needs exact Z to turn eventually; the reduction turned near Omega ~ 800
      and exact/reduction is falling. Predict Z:800 p_d >= 0.0306 (no more than 5% below Z:400's
      0.0322): flattening or rising. A fall of more than 10% (below 0.0290) weakens T-CASC-ap for
      narrow zones at accessible Omega; it would not refute a limit, but it would say Omega_c is
      not where the exact chain turns.
      -> HELD, and then some: Z:800 = 0.186 (from 0.032 at 400), eig_res 1.1e-10, on the laptop
         after it reproduced Z:230/hi to 1.1e-10. A 5.8x step is itself suspect (exact/reduction
         went 2.35 -> ~14), hence E7.

  E7  Written before running. If the turn is real it is smooth in Omega: Z:600 lies strictly
      between Z:400 and Z:800 (0.032 < p_d < 0.186), and Z:1000 lies above Z:800. Z:600 outside
      that band, or Z:1000 below Z:800, marks the 800 value as an isolated spike -- withdrawn
      until an independent route reproduces it.
      -> HELD: Z:600 = 0.0354, Z:1000 = 0.562 (laptop; eig_res <= 1.5e-10, |lambda_1| ~ 2.3).

  E8  Written before running: what makes the turn sharp? p_d = sum J(n2) f(n2) with f = h/c at
      (sad1 - 1, n2) -- stage 2's state at the dip (J) times what the episode does from there (f).
      Suspect (FINDINGS §123.4): pre-dip SAG -- J slides toward stage 2's separatrix as Omega grows.
      Cross-substitute on the common grid x2 = n2/Omega between Omega = 600 and 1000, both orders
      (rule 15). The suspect predicts: J_1000 with f_600 recovers >= 60% of the rise
      (p_d(1000) - p_d(600)) and J_600 with f_1000 recovers <= 40%. The reverse -- f carrying the
      rise -- means the damage is made after the dip, and the sag suspect is refuted.

  E9  Written after E8's run exposed that E5's route is NOT reproducible at Omega >= 600: the same
      cell, same code, same machine gave Z:1000 = 0.562 and then 0.188, Z:800 = 0.186 then 0.171
      (seed spread at Z:600 ~1%; Z:230 and Z:400 stable to 6 digits). The eigensolver starts from a
      random vector, and for this strongly non-normal generator a residual of 1e-10 does not bound
      the error in the eigenvector, even with |lambda_1| = 2.3. §123.4's turn is therefore unverified.
      Fix with no eigenproblem near zero: stage 1 is autonomous, so its conditioning function psi(n1)
      is one-dimensional -- psi = tau / tau(cap), tau the mean time to the dip, built in log domain
      (psi = lambda0 tau up to relative O(t_relax lambda0) ~ 1e-13). The Doob-transformed joint chain
      (stage 1's rates times psi ratios, stage 2 reflected at its saddle) is conservative with an O(1)
      gap; J's shape is its stationary law at n1 = sad1 -- one direct solve, no random start.
      Predict: within 1% of E5 at the seed-stable cells X:100, Z:230, Z:400; two different column
      orderings of the LU agree to 1e-6 at every Omega run. Z:800/1000 are read only from this route.

  E10 Written before running: the floor's other term. C = phi p_d + p0/(k1 t), and lim C = 1/2 only if
      p0/(k1 t) -> 0. From the same Q-process stationary law: r = (flux of stage 2 across its saddle
      while stage 1 is alive) / (flux of stage-1 dips), both read off phi = pi~/psi. Predict:
      (a) r reproduces §119's windowed p0/(k1 t) within 15% at X:30/40/55 and Z:40/55;
      (b) r(X) falls by >= 2x from Omega = 55 to 800 -- an approach that lingers near the saddle
          without crossing becomes rare relative to one that crosses. If r(X) instead stays within
          30% of 0.008, the floor's limit is above 1/2.
      -> (a) FAILED: r ran 21-30x the windowed value. The flux of stage 2 stepping below its saddle
         counts every recrossing, not committed falls -- the mistake §117.6 already recorded ("a
         one-step saddle current is not a committed rate"), repeated.

  E11 Written before running. Committed falls instead: in the Q-process, stage 2 now runs down to
      its low rail and is ABSORBED there (no reflection at the saddle); the fall rate is the flux
      into that rail, and J includes stage 2 already below its saddle at the dip (removing E5/E9's
      error bar). Still one near-zero mode, so the pinned solve stays well-conditioned. Predict:
      (a) the committed ratio reproduces windowed p0/(k1 t) within 15% at X:30/40/55, Z:40/55;
      (b) E10(b) unchanged; (c) p_d with the full J is >= the E9 value, by about the already-low
      fraction, at every cell.
      -> (a) HELD: 1.10-1.12 at all five cells. (c) HELD (difference <= 0.03%: the already-low
         fraction is negligible). (b) FAILED: p0/(k1 t) in X GROWS -- 0.009 / 0.018 / 0.073 / 0.304 at
         Omega = 55 / 100 / 200 / 400 -- so C passes 1/2 (0.66 at 400). The guess that it grows as
         exp(Omega dS), dS = V(s1) - V(edge2) = 0.0185 on stage 1's quasi-potential, does not fit: the
         local rate runs 0.0158 / 0.0137 / 0.0072, falling.

  E12 Written before running X:800 (laptop): the local growth rate of p0/(k1 t) from 400 to 800 keeps
      falling, below 0.005 -- sub-exponential, consistent with saturation. At or above 0.0185, the
      edge-approach yardstick is back.
      -> HELD by its letter (0.0037), but "consistent with saturation" was the wrong reading: local
         rate x Omega holds at 2.0-2.2, a power law ~Omega^2. The criterion could not tell (rule 19).

  E13 Written before running (T-CASC-ar's first measurement): WHERE is stage 1 when a never-dipped
      downstream fall commits? Reactive flux: stage 2 crossing its saddle (in the Q-process, stage 1
      alive), times q = P(stage 2 then reaches its low rail before re-crossing upward), by n1.
      Predict (a) its total equals the stored committed ratio (p0_rate_ratio) to 1e-6; (b) in X, >= 70%
      of it has stage 1 inside the zone [s1, edge2) at Omega = 100, 200, 400, not falling with Omega.
      The alternative -- falls with stage 1 ABOVE the edge -- would be stage 2's own noise over a
      lowered barrier, a different mechanism.
      -> (a) FAILED as built: 9-13% off. The fall rate is in the ORIGINAL process's QSD measure
         (phi = pi~/psi, stage 1 may still die); the crossings were weighted with the Q-PROCESS
         committor (stage 1 assumed to survive). psi varies sharply near the saddle, so this
         over-weights crossings made with stage 1 near its saddle -- exactly where the zone is.
         Fixed: u = E^Q[ 1/psi(n1 at the fall) ; low rail first ], boundary value 1/psi at the rail;
         then sum pi~ mu u = the fall flux exactly. (b) is judged on the corrected profile only.
"""


from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from experiments.weak_upstream import ARMS, element, pair_generator


def qsd_alive(Q, m, sad1):
    """Left principal eigenvector of the joint chain restricted to stage 1 alive (n1 >= sad1)."""
    N = Q.shape[0]
    alive = np.nonzero(np.arange(N) // m >= sad1)[0]
    Qa = Q[alive][:, alive].T.tocsc()
    vals, vecs = spla.eigs(Qa, k=1, sigma=0.0, which="LM")
    p = np.abs(np.real(vecs[:, 0]))
    full = np.zeros(N)
    full[alive] = p / p.sum()
    return full, float(-np.real(vals[0]))


def occupation_alive(Q, m, sad1, seed_vec):
    """Expected occupation of each alive state before stage 1 first steps below its saddle."""
    N = Q.shape[0]
    alive = np.nonzero(np.arange(N) // m >= sad1)[0]
    Qa = Q[alive][:, alive].T.tocsc()
    x = spla.spsolve(-Qa, seed_vec[alive])
    full = np.zeros(N)
    full[alive] = np.abs(x)
    return full / full.sum(), float("nan")


def qsd_alive_high(Q, m, sad1, sad2):
    """QSD of the chain with stage 1 alive AND stage 2 high: one tiny eigenvalue, well separated."""
    N = Q.shape[0]
    idx = np.arange(N)
    keep = np.nonzero((idx // m >= sad1) & (idx % m >= sad2))[0]
    Qk = Q[keep][:, keep].T.tocsc()
    vals, vecs = spla.eigs(Qk, k=2, sigma=0.0, which="LM")
    order = np.argsort(np.abs(vals))
    v = np.real(vecs[:, order[0]])
    v = v * np.sign(v.sum())
    res = float(np.abs(Qk @ v - np.real(vals[order[0]]) * v).max() / np.abs(v).max())
    full = np.zeros(N)
    full[keep] = np.clip(v, 0.0, None)
    gap = float(np.abs(vals[order[1]]) / max(np.abs(vals[order[0]]), 1e-300))
    return full / full.sum(), float(-np.real(vals[order[0]])), res, gap


def psi_stage1(A, om, sad1, cap):
    """Stage 1's conditioning function on n1 = sad1-1 .. cap: the mean time to first dip, normalised,
    in log domain. tau(n+1) - tau(n) = sum_{k>n} w_k / (lam_n w_n), w the product-formula weights."""
    n = np.arange(sad1 - 1, cap + 1, dtype=float)
    lam, mu = A.rates(n, np.zeros(len(n)), om, True)
    lw = np.zeros(len(n))
    for i in range(1, len(n)):
        lw[i] = lw[i - 1] + np.log(lam[i - 1]) - np.log(mu[i])
    # log of tail sums S_i = sum_{k>i} w_k
    ltail = np.full(len(n), -np.inf)
    acc = -np.inf
    for i in range(len(n) - 2, -1, -1):
        acc = np.logaddexp(acc, lw[i + 1])
        ltail[i] = acc
    ldelta = ltail[:-1] - (np.log(lam[:-1]) + lw[:-1])          # log increments, i = 0 .. len-2
    m0 = ldelta.max()
    tau = np.concatenate([[0.0], np.cumsum(np.exp(ldelta - m0))])
    return tau / tau[-1]                                         # psi on n1 = sad1-1 .. cap


def qprocess_J(A, B, om, m, sad1, sad2, permc="COLAMD", lower="sad"):
    """J's shape from the stationary law of the Doob-transformed (Q-process) joint chain, stage 2
    reflected at its saddle. Returns J over n2 (zero below sad2)."""
    cap = m - 1
    psi = psi_stage1(A, om, sad1, cap)                           # index j <-> n1 = sad1-1+j
    n1s = np.arange(sad1, cap + 1)
    lo2 = int(round(B.r1 * om))
    base2 = sad2 if lower == "sad" else lo2 + 1
    n2s = np.arange(base2, cap + 1)
    M1, M2 = len(n1s), len(n2s)
    I1, I2 = np.meshgrid(np.arange(M1), np.arange(M2), indexing="ij")
    I1, I2 = I1.ravel(), I2.ravel()
    idx = I1 * M2 + I2
    n1 = n1s[I1].astype(float); n2 = n2s[I2].astype(float)
    lam1, mu1 = A.rates(n1, np.zeros(len(n1)), om, True)
    lam2, mu2 = B.rates(n2, n1, om, False)
    ps = psi[I1 + 1]
    up1 = I1 < M1 - 1
    dn1 = I1 > 0
    lt1 = np.where(up1, lam1 * psi[np.minimum(I1 + 2, len(psi) - 1)] / ps, 0.0)
    mt1 = np.where(dn1, mu1 * psi[I1] / ps, 0.0)
    up2 = (I2 < M2 - 1) & (lam2 > 0)
    dn2 = (I2 > 0) & (mu2 > 0)
    leak = np.where((I2 == 0) & (mu2 > 0), mu2, 0.0) if lower == "rail" else np.zeros(len(I2))
    rows = [idx[up1], idx[dn1], idx[up2], idx[dn2]]
    cols = [idx[up1] + M2, idx[dn1] - M2, idx[up2] + 1, idx[dn2] - 1]
    vals = [lt1[up1], mt1[dn1], lam2[up2], mu2[dn2]]
    out = np.zeros(M1 * M2)
    for r_, v_ in zip(rows, vals):
        np.add.at(out, r_, v_)
    rows.append(idx); cols.append(idx); vals.append(-out - leak)
    Qt = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                       shape=(M1 * M2, M1 * M2)).tocsc()
    # pin one high-probability state (both stages at their rails) to 1 -- a single unit row keeps
    # the LU sparse; a dense normalisation row (the first draft) destroyed it at Omega >= 400
    k1 = int(np.clip(round(A.r3 * om) - sad1, 0, M1 - 1))
    k2 = int(np.clip(round(B.r3 * om) - base2, 0, M2 - 1))
    k = k1 * M2 + k2
    At = Qt.T.tocsr().tolil()
    At[k, :] = 0.0
    At[k, k] = 1.0
    b = np.zeros(M1 * M2); b[k] = 1.0
    pi = spla.splu(At.tocsc(), permc_spec=permc).solve(b)
    r_ = Qt.T @ pi
    r_[k] = 0.0                                                  # the pinned row is not an equation
    res = float(np.abs(r_).max() / np.abs(pi).max())
    P = pi.reshape(M1, M2)
    J = np.zeros(m)
    J[base2:] = np.clip(P[0, :], 0.0, None)                     # n1 = sad1 row
    # E10: phi = pi~/psi; dip flux = sum_n2 phi(sad1, n2) mu1(sad1); fall flux = sum_n1 phi(n1, sad2)
    # mu2(n1, sad2) -- stage 2 stepping below its saddle while stage 1 is alive
    psi_n1 = psi[np.arange(M1) + 1]
    _, mu1s = A.rates(np.array([float(sad1)]), np.zeros(1), om, True)
    dip = float(np.clip(P[0, :], 0, None).sum() / psi_n1[0] * mu1s[0])
    _, mu2s = B.rates(np.full(M1, float(base2)), n1s.astype(float), om, False)
    fall = float((np.clip(P[:, 0], 0, None) / psi_n1 * mu2s).sum())   # into the boundary below
    qprocess_J.last_ratio = fall / dip
    qprocess_J.last = dict(P=P, n1s=n1s, n2s=n2s, base2=base2, psi_n1=psi_n1, dip=dip, lt1=lt1, mt1=mt1,
                           M1=M1, M2=M2)
    return J / J.sum(), res


def fall_profile(A, B, om):
    """E13: reactive flux of committed never-dipped falls, by stage 1's position at the crossing."""
    from experiments.weak_upstream import edge as edge_of
    m = A.cap(om) + 1
    sad1, sad2 = int(np.ceil(A.r2 * om)), int(np.ceil(B.r2 * om))
    qprocess_J(A, B, om, m, sad1, sad2, lower="rail")
    L = qprocess_J.last
    P, n1s, base2, psi_n1, M1 = L["P"], L["n1s"], L["base2"], L["psi_n1"], L["M1"]
    lt1, mt1 = L["lt1"].reshape(M1, -1)[:, 0], L["mt1"].reshape(M1, -1)[:, 0]
    # committor q on n2 in (lo2, sad2): Q-process stage-1 moves + stage-2 moves; q = 1 at lo2, 0 at sad2
    lo2 = base2 - 1
    n2c = np.arange(lo2 + 1, sad2)
    K = len(n2c)
    I1, I2 = np.meshgrid(np.arange(M1), np.arange(K), indexing="ij")
    I1, I2 = I1.ravel(), I2.ravel()
    idx = I1 * K + I2
    lam2, mu2 = B.rates(n2c[I2].astype(float), n1s[I1].astype(float), om, False)
    l1, m1 = lt1[I1], mt1[I1]
    rows, cols, vals, diag = [], [], [], np.zeros(M1 * K)
    rhs = np.zeros(M1 * K)
    for ok, step, rate in ((I1 < M1 - 1, K, l1), (I1 > 0, -K, m1)):
        rows.append(idx[ok]); cols.append(idx[ok] + step); vals.append(rate[ok]); diag -= np.where(ok, rate, 0)
    up_in = I2 < K - 1
    rows.append(idx[up_in]); cols.append(idx[up_in] + 1); vals.append(lam2[up_in])
    diag -= lam2                                                  # up out of the last row hits sad2: q = 0
    dn_in = I2 > 0
    rows.append(idx[dn_in]); cols.append(idx[dn_in] - 1); vals.append(mu2[dn_in])
    diag -= mu2
    rhs -= np.where(~dn_in, mu2 / psi_n1[I1], 0.0)                # hits lo2: u = 1/psi(n1) (original measure)
    rows.append(idx); cols.append(idx); vals.append(diag)
    M = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(M1 * K, M1 * K)).tocsc()
    q = spla.splu(M).solve(rhs).reshape(M1, K)
    row = sad2 - base2                                            # P's column for n2 = sad2
    _, mu_s = B.rates(np.full(M1, float(sad2)), n1s.astype(float), om, False)
    flux = np.clip(P[:, row], 0, None) * mu_s * q[:, -1]   # (n1, sad2) -> (n1, sad2-1), committed, original measure
    x1 = n1s / om
    e2 = edge_of(B)
    tot = float(flux.sum() / L["dip"])
    return {"omega": om, "reactive_ratio": tot, "committed_ratio": qprocess_J.last_ratio,
            "frac_in_zone": float(flux[x1 < e2].sum() / flux.sum()),
            "frac_above_edge": float(flux[x1 >= e2].sum() / flux.sum()),
            "mean_x1": float((flux * x1).sum() / flux.sum()), "edge2": e2, "saddle1": A.r2}


def episode(A, B, om, method="qsd", dump=None, permc="COLAMD"):
    from experiments.weak_upstream import seed
    Q, m = pair_generator(A, B, om)
    N = m * m
    sad1, sad2 = int(np.ceil(A.r2 * om)), int(np.ceil(B.r2 * om))
    lo1, hi1 = int(round(A.r1 * om)), int(round(A.r3 * om))
    t0 = time.time()
    if method == "qsd":
        pi, lam0 = qsd_alive(Q, m, sad1)
    elif method == "occ":
        pi, lam0 = occupation_alive(Q, m, sad1, seed(A, B, om, m)[0])
    elif method == "hi":
        pi, lam0, eig_res, eig_gap = qsd_alive_high(Q, m, sad1, sad2)
    else:                                                        # "qp" (E9) / "qpr" (E11)
        Jq, qp_res = qprocess_J(A, B, om, m, sad1, sad2, permc, "rail" if method == "qpr" else "sad")
        pi, lam0 = None, float("nan")
    n2 = np.arange(m)
    _, mu1 = A.rates(np.full(m, float(sad1)), np.zeros(m), om, True)
    if method in ("qp", "qpr"):
        J = Jq
    else:
        J = pi.reshape(m, m)[sad1, :] * mu1                     # flux sad1 -> sad1-1, by n2
        J /= J.sum()
    # backward problems on interior n1 in (lo1, hi1); boundaries n1 = lo1 (0) and n1 = hi1 (target)
    idx = np.arange(N)
    n1 = idx // m
    inner = np.nonzero((n1 > lo1) & (n1 < hi1))[0]
    Qi = Q[inner][:, inner].tocsc()
    top = np.nonzero(n1 == hi1)[0]
    Qt = Q[inner][:, top]
    rhs_c = -np.asarray(Qt.sum(axis=1)).ravel()
    low2 = (top % m) < sad2
    rhs_h = -np.asarray(Qt[:, low2].sum(axis=1)).ravel()
    lu = spla.splu(Qi)
    c, h = lu.solve(rhs_c), lu.solve(rhs_h)
    pos = {s: i for i, s in enumerate(inner)}
    start = [pos[(sad1 - 1) * m + j] for j in range(m)]
    cs, hs = c[start], h[start]
    if dump is not None:
        f = np.where(cs > 0, hs / np.where(cs > 0, cs, 1.0), 0.0)
        np.savez_compressed(dump, x2=np.arange(m) / om, J=J, f=f, c=cs, sad2=sad2, om=om)
    return {"omega": om, "states": int(len(inner)), "secs": time.time() - t0,
            "phi_ep": float((J * cs).sum()), "p_d_ep": float((J * hs).sum() / (J * cs).sum()),
            "J_low_frac": float(J[:sad2].sum()), "lam0": lam0,
            **({"eig_res": eig_res, "eig_gap": eig_gap} if method == "hi" else {}),
            **({"qp_res": qp_res, "permc": permc, "p0_rate_ratio": qprocess_J.last_ratio}
               if method in ("qp", "qpr") else {})}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("results/episode_committor.json"))
    ap.add_argument("--cells", default="X:40,X:55,Z:40,Z:55")
    ap.add_argument("--dump", default=None, help="directory: save J and f per cell as .npz (E8)")
    ap.add_argument("--method", default="hi", choices=["qsd", "occ", "hi", "qp", "qpr"],
                    help="hi (default, E5): the only route valid at every Omega; qsd/occ kept as the record of E4's failure")
    args = ap.parse_args()
    rep = json.loads(args.out.read_text()) if args.out.exists() else {}
    for cell in args.cells.split(","):
        key = cell if args.method == "qsd" else f"{cell}/{args.method}"
        if key in rep:
            continue
        arm, om = cell.split(":")
        A, B = element(ARMS[arm][0]), element(ARMS[arm][1])
        dump = None if args.dump is None else str(pathlib.Path(args.dump) / f"{arm}_{om}_{args.method}.npz")
        if dump is not None:
            pathlib.Path(args.dump).mkdir(parents=True, exist_ok=True)
            key = f"{key}/dump"
        r = episode(A, B, int(om), args.method, dump)
        rep[key] = r
        print(key, json.dumps(r), flush=True)
        args.out.write_text(json.dumps(rep, indent=2, default=float))


if __name__ == "__main__":
    main()
