"""§118 (T-CASC-ak) -- do stage-1 near-misses CAUSE stage 2's excess falls?

§117 found the free two-stage chain's stage 2 falling ~1.2x as often at steady state as the
Q-process (stage 1 conditioned never to fail) predicts, while the ABSORBING chain -- which is
what the Q-process describes -- converges onto it exactly. The free chain's "stage 1 high at t"
is the union of trajectories where stage 1 never dipped below its saddle and trajectories where
it dipped and recovered, so the excess is carried by the second set. That much is accounting.
Whether the dips CAUSE the excess, rather than co-occurring with it, is what this tests.

TWO KILL TESTS, both exact CME, published element, cap_mult 1.25, matched seed, t up to 16.

  (1) TAGGING. Augment the state with a bit d = "stage 1 has been below its saddle". It starts
      at 0, flips to 1 on stage 1's first sub-saddle step, never resets. The state space doubles.
      At time t split the stage-1-high mass into d = 0 (never dipped) and d = 1 (dipped and
      recovered) and read P(stage 2 low) on each.
  (2) INTERVENTION. The same free chain, except stage 2's input is CLIPPED at the saddle:
      stage 2 sees max(n1, saddle), so stage 1's dips are invisible to it and nothing else
      changes (above the saddle the coupling is untouched). This is do(stage 2's input) rather
      than a conditioning, so it separates cause from co-occurrence.

PREDICTIONS, WRITTEN BEFORE RUNNING.

  P1  WIRING. The tagged chain's marginal (summed over d) reproduces §117's free-chain
      P(stage 2 low | stage 1 high) at t = 2 to the solver floor (~1e-6).
  P2  THE NEVER-DIPPED SET IS THE Q-PROCESS. Its model/measured ratio (pi_cond (1 - e^{-lambda t})
      over P(stage 2 low | d = 0, stage 1 high)) tracks the ABSORBING chain's §117 sequence --
      1.139, 1.066, 1.031, 1.015 at t = 2, 4, 8, 16 -- up to the seed difference (matched free
      QSD here vs the reflected law there), and in particular does NOT plateau near 0.84.
  P3  THE EXCESS SITS IN THE DIPPED SET. P(stage 2 low | d = 1, stage 1 high) is several times
      P(stage 2 low | d = 0, stage 1 high), and the dipped set's weight times that difference
      reproduces the free chain's excess over the never-dipped set.
  P4  CAUSATION, and it can print either verdict. With stage 2's input clipped at the saddle,
      the free chain's plateau ratio moves from ~0.836 to within 0.03 of the never-dipped value.
      If instead it stays near 0.84, the dips co-occur with stage 2's falls without causing them
      -- the suspect is then a common cause (e.g. stage 1 spending time NEAR its saddle, which
      clipping at the saddle does not remove), and T-CASC-ak's mechanism is refuted.
  P5  THE Omega-INDEPENDENCE, reported not predicted: dipped weight and per-trajectory effect at
      Omega = 14, 30, 55, separately. §117 found the product near-constant; if the two factors
      move in opposite directions, the constancy is a cancellation.

WHAT THIS CANNOT SETTLE. Clipping at the saddle removes the dips' effect on stage 2 but leaves
stage 1's excursions to just above the saddle visible. A partial reduction under P4 would mean
the channel is "stage 1 near its saddle", of which sub-saddle dips are only the deepest part.

OUTCOMES, appended after running (predictions above and in each function left as written).
P1-P3 held. P4 failed and was mis-designed (the saddle is below stage 2's edge). P6 failed, P7
failed, P8 failed its criterion while closing 70% of the gap, P9b failed (saturates at ~0.75),
P10b failed, P11 met its criterion by over-correction, P11b held (redraw alone ~0.24). Read in
FINDINGS §118: near-misses CAUSE the excess, three-quarters of it after stage 1 recovers.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from experiments.free_upstream_depth import build_gen, seed_matched
from experiments.two_axes_for_the_position import PUBLISHED_RAILS, Elem

MARKS = (2.0, 4.0, 8.0, 16.0)
EDGE = 1.5795                      # stage 2's bistability edge in x_up (§105)
DOSE = (1.3, 1.5795, 2.0, 2.5)


def tagged_mixed(Q0, Q1, m, sad):
    """Tagged chain whose d = 0 block evolves under Q0 and d = 1 block under Q1: intervene on
    trajectories only AFTER stage 1's first dip, leaving the never-dipped ones exactly as they
    are (so their law stays the Q-process's)."""
    N = Q0.shape[0]
    C0, C1 = Q0.tocoo(), Q1.tocoo()
    r, c, v = C0.row, C0.col, C0.data
    dip = (r != c) & (r // m >= sad) & (c // m < sad)
    rows = [r[~dip], r[dip], C1.row + N]
    cols = [c[~dip], c[dip] + N, C1.col + N]
    vals = [v[~dip], v[dip], C1.data]
    return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(2 * N, 2 * N)).tocsr()


def free_generator(e, om, clip=False, sub_saddle_as=None, below=None):
    """Two free stages. `clip`: stage 2 sees max(n1, saddle) (P4, and see P7 for why that was a
    mis-designed intervention). `sub_saddle_as`: while stage 1 is BELOW its saddle, stage 2 sees
    this concentration instead; above the saddle the coupling is untouched."""
    cap = e.cap(om)
    m = cap + 1
    N = m * m
    idx = np.arange(N)
    n1, n2 = (idx // m).astype(float), (idx % m).astype(float)
    nsad = float(np.ceil(e.r2 * om))
    rows, cols, vals = [], [], []
    diag = np.zeros(N)
    l1, u1 = e.rates(n1, np.zeros(N), om, True)
    up1 = (n1 < cap) & (l1 > 0)
    dn1 = (n1 > 0) & (u1 > 0)
    rows.append(idx[up1]); cols.append(idx[up1] + m); vals.append(l1[up1])
    rows.append(idx[dn1]); cols.append(idx[dn1] - m); vals.append(u1[dn1])
    diag -= np.where(up1, l1, 0.0) + np.where(dn1, u1, 0.0)
    seen = np.maximum(n1, nsad) if clip else n1
    if sub_saddle_as is not None:
        thresh = (below if below is not None else e.r2) * om
        seen = np.where(n1 < thresh, sub_saddle_as * om, n1)
    l2, u2 = e.rates(n2, seen, om, False)
    up2 = (n2 < cap) & (l2 > 0)
    dn2 = (n2 > 0) & (u2 > 0)
    rows.append(idx[up2]); cols.append(idx[up2] + 1); vals.append(l2[up2])
    rows.append(idx[dn2]); cols.append(idx[dn2] - 1); vals.append(u2[dn2])
    diag -= np.where(up2, l2, 0.0) + np.where(dn2, u2, 0.0)
    rows.append(idx); cols.append(idx); vals.append(diag)
    Q = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(N, N)).tocsr()
    return Q, m


def tagged(Q, m, sad):
    """Double the state space with a 'stage 1 has dipped' bit. A stage-1 step from n1 >= sad to
    n1 < sad moves d = 0 -> 1; d never resets. The d = 1 block is a copy of Q."""
    N = Q.shape[0]
    C = Q.tocoo()
    r, c, v = C.row, C.col, C.data
    off = r != c
    n1r, n1c = r // m, c // m
    dip = off & (n1r >= sad) & (n1c < sad)
    keep0 = ~dip                                   # d = 0 transitions that stay in d = 0
    rows = [r[keep0], r[dip], r + N]
    cols = [c[keep0], c[dip] + N, c + N]
    vals = [v[keep0], v[dip], v]
    return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(2 * N, 2 * N)).tocsr()


def conditional_series(Q, p0, m, sad, tagged_chain=False):
    """Evolve in 2-unit steps; at each mark return P(stage 2 low | stage 1 high) split by d."""
    QT = (Q.T * 2.0).tocsr()
    N = m * m
    s1hi = (np.arange(N) // m) >= sad
    s2lo = (np.arange(N) % m) < sad
    p, t, out = p0.copy(), 0.0, []
    while t < MARKS[-1] - 1e-9:
        p = spla.expm_multiply(QT, p)
        t += 2.0
        if any(abs(t - mk) < 1e-9 for mk in MARKS):
            blocks = [p[:N], p[N:]] if tagged_chain else [p]
            rec = {"t": t}
            hi_all = sum(b[s1hi].sum() for b in blocks)
            rec["all"] = float(sum(b[s1hi & s2lo].sum() for b in blocks) / hi_all)
            if tagged_chain:
                for d, b in zip((0, 1), blocks):
                    w = float(b[s1hi].sum())
                    rec[f"w{d}"] = w / hi_all
                    rec[f"d{d}"] = float(b[s1hi & s2lo].sum() / w) if w > 0 else float("nan")
            out.append(rec)
    return out


def run(om, e, pi_cond, lam):
    Q, m = free_generator(e, om)
    sad = int(np.ceil(e.r2 * om))
    _, ref, dims, strides, _, walled = build_gen(om, 2, 0, cap_mult=e.cap_mult)
    p0 = seed_matched(om, ref, dims, strides, walled)
    assert len(p0) == m * m
    T2 = tagged(Q, m, sad)
    tag = conditional_series(T2, np.concatenate([p0, np.zeros_like(p0)]), m, sad, True)
    Qc, _ = free_generator(e, om, clip=True)
    clip = conditional_series(Qc, p0, m, sad)
    model = lambda t: pi_cond * (1 - math.exp(-lam * t))
    rows = []
    for a, b in zip(tag, clip):
        t = a["t"]
        rows.append({"t": t, "model": model(t), "free_all": a["all"], "never": a["d0"],
                     "dipped": a["d1"], "w_dipped": a["w1"], "clipped": b["all"],
                     "r_free": model(t) / a["all"], "r_never": model(t) / a["d0"],
                     "r_clipped": model(t) / b["all"]})
    return rows


def approach(om, e, t_end=16.0):
    """Written before running. P4 refuted 'the sub-saddle input damages stage 2': clipping it
    changes nothing. Suspect instead: SELECTION. A stage-1 trajectory that crosses its saddle has
    first spent a long time approaching it -- mostly ABOVE it, where hill(x) already roughly
    halves stage 2's drive -- so stage 2 is damaged BEFORE the crossing, the crossing is only a
    marker, and the Q-process under-predicts because its conditioning deletes exactly those
    deepest excursions.

      P6  At the instant stage 1 first crosses (weighted by the crossing flux, taken from the
          never-dipped block near steady state), stage 2 is ALREADY low with probability >= 0.6
          at every Omega. Below 0.3 would mean the damage comes after the crossing -- which the
          clipping result would then leave unexplained.
      P6b Within the never-dipped block, P(stage 2 low | n1) rises steeply as n1 approaches the
          saddle from above: the gradient is where the damage is done.
    """
    Q, m = free_generator(e, om)
    sad = int(np.ceil(e.r2 * om))
    _, ref, dims, strides, _, walled = build_gen(om, 2, 0, cap_mult=e.cap_mult)
    p0 = seed_matched(om, ref, dims, strides, walled)
    T2 = tagged(Q, m, sad)
    QT = (T2.T * 2.0).tocsr()
    p = np.concatenate([p0, np.zeros_like(p0)])
    t = 0.0
    while t < t_end - 1e-9:
        p = spla.expm_multiply(QT, p)
        t += 2.0
    g0 = p[: m * m].reshape(m, m)                      # never-dipped block
    _, mu_sad = e.rates(float(sad), 0.0, om, True)
    J = g0[sad, :] * mu_sad                            # crossing flux sad -> sad-1, by n2
    at_cross = float(J[:sad].sum() / J.sum())
    rail = int(round(e.r3 * om))
    prof = {}
    for n1 in sorted(set([sad, sad + 1, sad + 2, sad + 4, (sad + rail) // 2, rail])):
        row = g0[n1, :]
        prof[int(n1)] = float(row[:sad].sum() / row.sum()) if row.sum() > 0 else float("nan")
    return {"omega": om, "p_low_at_crossing": at_cross, "profile": prof, "sad": sad, "rail": rail}


def intervene(om, e, pi_cond, lam, level):
    """P7, written before running. P4's clip held stage 2's input at stage 1's SADDLE (x = 1.0),
    which is below stage 2's own bistability edge (§105: x_up = 1.5795) -- so every dip still
    left stage 2 monostable-low, and the intervention could not remove the effect it tested.
    Its null was uninformative (rule 19, for an intervention rather than a criterion).

      P7  With stage 2 seeing stage 1 AS HEALTHY (at its operating point mu1) whenever stage 1 is
          below its saddle, and untouched otherwise, the free chain's ratio at t = 16 moves to
          within 0.03 of the never-dipped value: the sub-saddle phase CAUSES the excess.
      P7b Clipping at the edge itself (1.58) instead: predict most, not all, of the excess goes,
          since stage 2 is only marginally bistable there.
      If P7 fails -- the ratio stays near 0.84 -- the dips do not cause the excess even when the
      intervention is strong enough to matter, and T-CASC-ak's mechanism is refuted outright.
    """
    Q, m = free_generator(e, om, sub_saddle_as=level)
    sad = int(np.ceil(e.r2 * om))
    _, ref, dims, strides, _, walled = build_gen(om, 2, 0, cap_mult=e.cap_mult)
    p0 = seed_matched(om, ref, dims, strides, walled)
    ser = conditional_series(Q, p0, m, sad)
    return [{"t": r["t"], "ratio": pi_cond * (1 - math.exp(-lam * r["t"])) / r["all"]}
            for r in ser]


def after_first_dip(om, e, pi_cond, lam, edge, mu1):
    """P8, written before running. §118's zone check: in the dipped-and-recovered set stage 1 is
    back above stage 2's bistability edge 94% of the time, yet stage 2 is still low 83-89% --
    hysteresis, stage 2 trapped in its own low basin. P6 says it is mostly still high at the first
    crossing; P7 says hiding the sub-saddle time removes only 11-17%. Suspect: the damage is done
    while stage 1 lingers near its saddle AFTER first reaching it, mostly ABOVE it but below
    stage 2's edge (x_up < 1.58), where stage 2 has no high state.

      P8  Only AFTER stage 1's first dip (d = 1 block), stage 2 sees mu1 whenever stage 1 is below
          the EDGE; never-dipped trajectories are untouched. Predict the free chain's ratio comes
          within 0.03 of the never-dipped value at t = 16: the post-dip lingering below the edge
          CAUSES the excess. If it does not, something upstream of the first crossing -- which
          this cannot touch -- carries the rest.
    """
    Q0, m = free_generator(e, om)
    Q1, _ = free_generator(e, om, sub_saddle_as=mu1, below=edge)
    sad = int(np.ceil(e.r2 * om))
    _, ref, dims, strides, _, walled = build_gen(om, 2, 0, cap_mult=e.cap_mult)
    p0 = seed_matched(om, ref, dims, strides, walled)
    ser = conditional_series(tagged_mixed(Q0, Q1, m, sad),
                             np.concatenate([p0, np.zeros_like(p0)]), m, sad, True)
    return [{"t": r["t"], "ratio": pi_cond * (1 - math.exp(-lam * r["t"])) / r["all"]}
            for r in ser]


def dose(om, e, pi_cond, lam, levels, mu1):
    """P9, written before running, after P8 read 0.9640 / 0.9647 against never-dipped 1.0175 /
    1.0243 (Om = 30 / 55, t = 16): hiding post-dip stage 1 below the EDGE closes ~70% of the gap,
    at both volumes -- P8's 0.03 criterion FAILED, the effect is large and Om-flat.

      P9a Self-check: tagged_mixed(Q0, Q0) must reproduce the untagged free chain (r_free) to the
          solver floor (1e-9 relative) -- otherwise P8 measured the harness.
      P9b Dose-response in the post-dip hiding threshold, below = 1.3, 1.58, 2.0, 2.5 (x Omega).
          Predict the closed fraction (r - r_free) / (r_never - r_free) rises monotonically, through
          P7's ~0.17 at the saddle and P8's ~0.71 at the edge, and reaches >= 0.9 by 2.0: the
          remaining 30% is stage 1 still sagged just above the edge, where stage 2 is bistable but
          its high well is shallow. If it saturates well below 0.9 by 2.5, the remainder is not the
          input level at all and something else carries it (e.g. an earlier-fallen stage 2 that
          hiding cannot lift). Above 1 means over-correction: hiding pulls the dipped set BETTER
          than never-dipped, and the threshold has left the dip's own zone.
    """
    Q0, m = free_generator(e, om)
    sad = int(np.ceil(e.r2 * om))
    _, ref, dims, strides, _, walled = build_gen(om, 2, 0, cap_mult=e.cap_mult)
    p0 = seed_matched(om, ref, dims, strides, walled)
    p0t = np.concatenate([p0, np.zeros_like(p0)])
    model = lambda t: pi_cond * (1 - math.exp(-lam * t))
    out = {"self": conditional_series(tagged_mixed(Q0, Q0, m, sad), p0t, m, sad, True)[-1]["all"]}
    for lv in levels:
        Q1, _ = free_generator(e, om, sub_saddle_as=mu1, below=lv)
        last = conditional_series(tagged_mixed(Q0, Q1, m, sad), p0t, m, sad, True)[-1]
        out[str(lv)] = model(last["t"]) / last["all"]
    out["self_ratio"] = model(MARKS[-1]) / out["self"]
    return out


def tagged3(Q0, Q1, m, sad):
    """Three blocks: 0 never dipped (Q0); 1 dipped with stage 2 HIGH at the first dip (Q1);
    2 dipped with stage 2 already LOW at the first dip (Q1). A dip is a stage-1 step, so n2 is
    the same on both sides of it and the split is read off the transition's own row."""
    N = Q0.shape[0]
    C0, C1 = Q0.tocoo(), Q1.tocoo()
    r, c, v = C0.row, C0.col, C0.data
    dip = (r != c) & (r // m >= sad) & (c // m < sad)
    lo = (r % m) < sad
    rows = [r[~dip], r[dip & ~lo], r[dip & lo], C1.row + N, C1.row + 2 * N]
    cols = [c[~dip], c[dip & ~lo] + N, c[dip & lo] + 2 * N, C1.col + N, C1.col + 2 * N]
    vals = [v[~dip], v[dip & ~lo], v[dip & lo], C1.data, C1.data]
    return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(3 * N, 3 * N)).tocsr()


def already_low(om, e, mu1, below=2.5, t_end=16.0):
    """P10, written before running. P9b FAILED: the closed fraction saturates at ~0.75 (Om = 30:
    0.594, 0.706, 0.735, 0.746 hiding below 1.3, 1.58, 2.0, 2.5), not >= 0.9. P9's own alternative:
    part of the excess is a stage 2 that had ALREADY fallen before stage 1's first dip, which no
    post-dip intervention can lift. P6 measured P(stage 2 low | at the first crossing) = 0.20 at
    Om = 30 against 0.093 with stage 1 at its rail -- a selection effect, far smaller than P6's
    0.6 guess but not zero.

      P10a Wiring: blocks 1 + 2 together reproduce P9's dipped block (the 2.5 run's ratio) -- the
           marginal of tagged3 equals tagged_mixed's to the solver floor.
      P10b Under hiding below 2.5, split the dipped set by stage 2's state AT the first dip. The
           high-at-dip subgroup's P(stage 2 low | stage 1 high) is within 0.03 of the never-dipped
           set's (it saw a healthy input throughout, apart from stage 1 above 2.5 Omega), and the
           low-at-dip subgroup carries >= 80% of the dipped set's remaining excess. If instead the
           high-at-dip subgroup still sits well above never-dipped, the residual 25% is made AFTER
           the dip by something the input-hiding does not reach, and "already fallen" is refuted.
    """
    Q0, m = free_generator(e, om)
    Q1, _ = free_generator(e, om, sub_saddle_as=mu1, below=below)
    sad = int(np.ceil(e.r2 * om))
    _, ref, dims, strides, _, walled = build_gen(om, 2, 0, cap_mult=e.cap_mult)
    p0 = seed_matched(om, ref, dims, strides, walled)
    N = m * m
    QT = (tagged3(Q0, Q1, m, sad).T * 2.0).tocsr()
    p = np.concatenate([p0, np.zeros(2 * N)])
    t = 0.0
    while t < t_end - 1e-9:
        p = spla.expm_multiply(QT, p)
        t += 2.0
    s1hi = (np.arange(N) // m) >= sad
    s2lo = (np.arange(N) % m) < sad
    out = {"t": t}
    hi_all = sum(p[k * N:(k + 1) * N][s1hi].sum() for k in range(3))
    out["all"] = float(sum(p[k * N:(k + 1) * N][s1hi & s2lo].sum() for k in range(3)) / hi_all)
    for k, name in enumerate(("never", "hi_at_dip", "lo_at_dip")):
        b = p[k * N:(k + 1) * N]
        w = float(b[s1hi].sum())
        out[f"w_{name}"] = w / hi_all
        out[f"p_{name}"] = float(b[s1hi & s2lo].sum() / w)
    return out


def redraw(om, e, mu1, below=2.5, t_end=16.0):
    """P11, written before running. P10 at t = 16, hiding below 2.5: the LOW-at-dip subgroup
    carries 79% (Om = 30) and 64% (Om = 55) of the dipped set's remaining excess, and the
    HIGH-at-dip subgroup still falls far more often than never-dipped (0.133 vs 0.094; 0.089 vs
    0.0011). P10b FAILED at both volumes. P10's line -- stage 2 above or below its saddle at the
    dip -- is the same kind of line P4 drew for stage 1: a stage 2 that is "high" at the dip can be
    sagged just above its saddle after stage 1's long approach, and hiding the input afterwards
    does not undo that. Suspect: the whole residual is stage 2's STATE at the first dip, graded,
    not a binary fallen/not-fallen.

      P11 Counterfactual do(stage 2's state): at stage 1's first dip, redraw stage 2 from its
          healthy law (the matched seed's stage-2 marginal restricted to n2 >= saddle), then hide
          post-dip input below 2.5 as in P9/P10. Predict the closed fraction reaches >= 0.95 at
          both Omega: stage 2's state at the dip plus its input after it account for everything.
          If it stays near P9's 0.75, the residual is made after the dip by stage 1 above 2.5
          Omega, and the pre-dip account is refuted.

    P11 read 1.113 (Om = 30) and 1.001 (Om = 55): criterion met, but at Om = 30 by OVER-correction
    (the redrawn dipped set falls less often than never-dipped, 0.028 vs 0.094). The criterion
    could be met by any sledgehammer, so P11 shows the pair is SUFFICIENT and nothing about the
    split (rule 19). Hence, written before running:

      P11b Redraw only (below=None: no input hiding after the dip). Predict it closes ~0.25 at
           both Omega, the complement of P9's 0.75 -- the two shares are additive. If redraw-only
           plus hide-only (P9, 0.746 at 2.5) sum well above 1, the pre-dip sag and the post-dip
           input interact, the split is order-dependent, and both orders get reported (rule 15).
    """
    Q0, m = free_generator(e, om)
    Q1 = Q0 if below is None else free_generator(e, om, sub_saddle_as=mu1, below=below)[0]
    sad = int(np.ceil(e.r2 * om))
    _, ref, dims, strides, _, walled = build_gen(om, 2, 0, cap_mult=e.cap_mult)
    p0 = seed_matched(om, ref, dims, strides, walled)
    N = m * m
    q = p0.reshape(m, m).sum(axis=0)
    q[:sad] = 0.0
    q /= q.sum()
    sup = np.nonzero(q > 1e-300)[0]
    C0, C1 = Q0.tocoo(), Q1.tocoo()
    r, c, v = C0.row, C0.col, C0.data
    dip = (r != c) & (r // m >= sad) & (c // m < sad)
    dr, dc, dv = r[dip], c[dip], v[dip]
    rr = np.repeat(dr, len(sup))
    cc = (np.repeat(dc // m, len(sup)) * m + np.tile(sup, len(dr))) + N
    vv = np.repeat(dv, len(sup)) * np.tile(q[sup], len(dr))
    T = sp.coo_matrix((np.concatenate([v[~dip], vv, C1.data]),
                       (np.concatenate([r[~dip], rr, C1.row + N]),
                        np.concatenate([c[~dip], cc, C1.col + N]))), shape=(2 * N, 2 * N)).tocsr()
    ser = conditional_series(T, np.concatenate([p0, np.zeros(N)]), m, sad, True)
    return ser[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("results/near_misses.json"))
    args = ap.parse_args()
    from experiments.the_survival_term import conditioned_gap, conditioned_qsd_low

    e = Elem("published", *PUBLISHED_RAILS, cap_mult=1.25)
    rep = {"cells": {}}
    for om in (14, 30, 55):
        l0, l1, _ = conditioned_gap(e, om)
        pc = conditioned_qsd_low(e, om)
        rows = run(om, e, pc, l1 - l0)
        rep["cells"][str(om)] = rows
        print(f"\nOmega = {om}   (pi_cond {pc:.4f}, lambda {l1 - l0:.4e})")
        print(f"{'t':>5}{'free':>9}{'never-dipped':>14}{'CLIPPED':>10}"
              f"{'P(lo|dipped)':>14}{'P(lo|never)':>13}{'dipped wt':>11}")
        for r in rows:
            print(f"{r['t']:>5g}{r['r_free']:>9.4f}{r['r_never']:>14.4f}{r['r_clipped']:>10.4f}"
                  f"{r['dipped']:>14.4e}{r['never']:>13.4e}{r['w_dipped']:>11.3e}", flush=True)

    ref117 = json.loads(pathlib.Path("results/where_the_remainder_lives.json").read_text())
    r117 = {r["omega"]: r["los1"] for r in ref117["rows"]}
    f30 = rep["cells"]["30"][0]["free_all"]
    w = abs(f30 / r117[30] - 1)
    print(f"\nP1 -- tagged marginal vs §117's free chain at Om=30, t=2: rel {w:.1e}"
          f" -> {'HOLDS' if w < 1e-5 else 'FAILS'}")
    rep["p1_rel"] = w

    # P6-P11: the causal chain. Each was written after reading the one before; see docstrings.
    rep["approach"] = [approach(om, e) for om in (14, 30, 55)]
    rep["p7"], rep["p8"], rep["p9"], rep["p10"], rep["p11"], rep["p11b"] = {}, {}, {}, {}, {}, {}
    for om in (30, 55):
        l0, l1, _ = conditioned_gap(e, om)
        pc, lam = conditioned_qsd_low(e, om), l1 - l0
        mu1 = e.operating_points(om, 2)[0][0]
        model16 = pc * (1 - math.exp(-lam * MARKS[-1]))
        b = rep["cells"][str(om)][-1]
        closed = lambda r: (r - b["r_free"]) / (b["r_never"] - b["r_free"])
        rep["p7"][str(om)] = {"as_healthy": intervene(om, e, pc, lam, mu1),
                              "at_edge": intervene(om, e, pc, lam, EDGE),
                              "mu1": mu1, "edge": EDGE}
        rep["p8"][str(om)] = after_first_dip(om, e, pc, lam, EDGE, mu1)
        o = dose(om, e, pc, lam, DOSE, mu1)
        o["r_free"], o["r_never"] = b["r_free"], b["r_never"]
        rep["p9"][str(om)] = o
        rep["p10"][str(om)] = already_low(om, e, mu1)
        for key, below in (("p11", 2.5), ("p11b", None)):
            x = redraw(om, e, mu1, below=below)
            x["ratio"] = model16 / x["all"]
            x["closed"] = closed(x["ratio"])
            rep[key][str(om)] = x
        print(f"Om={om}: P8 {closed(rep['p8'][str(om)][-1]['ratio']):.3f}  "
              + "  ".join(f"P9@{lv} {closed(o[str(lv)]):.3f}" for lv in DOSE)
              + f"  P11 {rep['p11'][str(om)]['closed']:.3f}  P11b {rep['p11b'][str(om)]['closed']:.3f}",
              flush=True)
    args.out.write_text(json.dumps(rep, indent=2, default=float))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
