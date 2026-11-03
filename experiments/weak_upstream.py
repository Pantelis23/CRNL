"""§119 -- a weaker upstream: is near-miss transmission exponential in Omega, and does it need a zone?

§118 showed, by intervention, that stage-1 near-misses CAUSE stage 2's excess falls, that three-
quarters of the damage is done while a just-recovered stage 1 climbs between ITS saddle and stage
2's bistability EDGE (where stage 2 has no high state), and that the excess is Omega-flat. §118 used
two copies of ONE element, so stage 1's dip rate k1 and stage 2's own fall rate k2 share an action
and their ratio is flat by construction (T-CASC-am). This section breaks that.

THE ARGUMENT, stated before any run. Given stage 1 high at t, stage 2's low probability is
    P_free = (1 - w_d) p0 + w_d p_d,
with w_d the weight of dipped-and-recovered stage-1 trajectories (~ k1 t, a stage-1 quantity), p0
stage 2's low probability on never-dipped trajectories (~ k2 t: the Q-process, which every
conditioned model in this arc approximates) and p_d the same on dipped ones. §118 measured
p_d ~ 0.8-0.9 and said why: the fall happens where stage 2 has NO barrier, so p_d is O(1) and not
exponentially small. Then the excess over every conditioned model is
    E = P_free / p0 - 1  ~  (k1 / k2) p_d,
exponential in Omega at the difference of the two dip actions -- IF the zone [saddle_1, edge_2)
exists. If it does not, a near-miss that recrosses quickly never leaves stage 2 without a high
state, and p_d should itself be exponentially small.

ELEMENTS. The r1 = 0.15, r3 = 3.0 family (§114), so every pair shares its rails and the Hill
coupling (normalised by stage 2's own r3 = 3.0) is the published one. Dip rates k = decay rate of
the chain absorbed at the element's own saddle (stage 1 driven by nothing, stage 2 pinned at stage
1's operating point) -- measured below, not assumed. Edges from the deterministic cubic:
    E-062 saddle 1.233 (edge 1.770)   E-070 saddle 1.005 (edge 1.616)
    E-084 saddle 0.606, edge 1.362    F-025 (r2 = 0.25, new) edge 1.127

ARMS (stage 1 -> stage 2):
    X  E-070 -> E-084   zone [1.005, 1.362).   k1/k2 ~ 280, 2100, 45000 at Omega = 30, 40, 55.
    S  E-070 -> E-070   §118's shape (identical elements), same Omega: the flat reference.
    Z  E-062 -> E-084   zone [1.233, 1.362), narrow.
    N  E-062 -> F-025   NO zone: stage 1's saddle 1.233 sits above stage 2's edge 1.127.
Exact CME on the joint chain, cap_mult 1.25, stage-1 dip tagged as in §118, seed = product of the two
elements' quasi-stationary laws above their saddles (stage 2's pinned at stage 1's operating point).

PREDICTIONS, WRITTEN BEFORE RUNNING.

  P1  WIRING. Arm S at Omega = 30 through this module's generator reproduces §118's machinery on
      the same pair built by `near_misses.free_generator` (E-070 as a one-element `Elem`) to the
      solver floor, and the tagged marginal equals the untagged chain.
  P2  THE LAW, arm X. E(t = 16) grows with Omega at the rate of k1/k2: the product E * k2/k1
      changes by less than 3x over Omega = 30 -> 55 while k1/k2 itself changes ~150x. E exceeds 10
      at Omega = 55: every conditioned model under-predicts stage 2's failure-while-its-input-is-
      right by more than an order of magnitude. If E instead stays O(0.2) as in §118, the Omega-
      flatness was not the shared action, and the whole argument above is wrong.
      (Kill test for the criterion, rule 19: flat E gives E k2/k1 falling ~150x -> prints FAILS;
      E tracking k1/k2 gives a bounded product -> HOLDS. Either verdict is reachable.)
  P3  THE REFERENCE, arm S. E stays within a factor 1.5 across Omega = 30 -> 55 (identical elements).
  P4  THE ZONE IS NECESSARY. Same stage 1 (E-062), Omega = 40 and 55: p_d in arm N (no zone) is
      below 0.1 x p_d in arm Z (zone), at both Omega, and the gap widens with Omega (exponentially
      small p_d without a zone, O(1)-ish with one). If p_d(N) is comparable to p_d(Z), the zone
      reading of §118 is wrong: the damage would not need stage 2 to lose its high state.
  P5  ORDER OF MAGNITUDE, rule 16. Predict E from single-stage quantities alone, E_pred =
      w_d * p_d_ref / p0, taking w_d and p0 from the chain but p_d_ref = §118's 0.85, and report
      E / E_pred per cell. Not a pass/fail -- p_d is expected to depend on zone width -- but it is
      recorded so the ratio is never inferred from a fit.

OUTCOMES, appended after running (predictions above, in OUT_OF_SAMPLE and in CONVERGENCE left as
written). P1, P3, P4 held. P2 FAILED: E saturates in arm X (p0 is on stage 1's clock, not k2's).
P5 recorded: 0.34-0.36 in X, <= 0.08 in narrow zones. Q1 FAILED (equal-width pairs a, Z differ ~5x);
Q2, Q3, Q4 held -- Q3 is the result: P/(k1 t) flat in Omega out of sample. C1, C2 held. FINDINGS §119.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from experiments.two_axes_for_the_position import ELEMENTS, R1_FIX, R3_FIX, Elem

CAP_MULT = 1.25
MARKS = (4.0, 8.0, 16.0)


class ElemK(Elem):
    """An element whose input coupling has Hill threshold K instead of the published 1 (§124.9).
    The Hill factor stays normalised at r3, so a healthy input still gives hill = 1: the element's
    own rates at its operating point are untouched, and only its bistability edge moves."""

    def __init__(self, name, r1, r2, r3, hill_k, cap_mult=CAP_MULT):
        super().__init__(name, r1, r2, r3, cap_mult=cap_mult)
        self.hill_k = hill_k

    def hill(self, x_up):
        from experiments.two_axes_for_the_position import HILL_N
        z = (x_up / self.hill_k) ** HILL_N
        zr = (self.r3 / self.hill_k) ** HILL_N
        return (z / (1.0 + z)) / (zr / (1.0 + zr))


def element(name):
    """'E-084' etc.; 'E-084@K1.2' is the same element with its coupling threshold set to K = 1.2."""
    base, _, k = name.partition("@K")
    table = dict(ELEMENTS)
    table["F-025"] = 0.25
    table["F-040"] = 0.40
    table["F-050"] = 0.50
    if k:
        return ElemK(name, R1_FIX, table[base], R3_FIX, float(k), cap_mult=CAP_MULT)
    return Elem(name, R1_FIX, table[base], R3_FIX, cap_mult=CAP_MULT)


def dip_rate(el, om, x_up, first):
    """Decay rate of one element's chain absorbed at its own saddle."""
    Q, m, _ = el._gen(om, x_up, first)
    sad = int(np.ceil(el.r2 * om))
    keep = np.arange(sad, m)
    A = Q[keep][:, keep].toarray()
    return float(np.min(-np.real(np.linalg.eigvals(A))))


def pair_generator(A, B, om):
    """Stage 1 = A (undriven), stage 2 = B driven by stage 1's count. Common cap (shared r3)."""
    cap = A.cap(om)
    assert cap == B.cap(om)
    m = cap + 1
    N = m * m
    idx = np.arange(N)
    n1, n2 = (idx // m).astype(float), (idx % m).astype(float)
    rows, cols, vals, diag = [], [], [], np.zeros(N)
    l1, u1 = A.rates(n1, np.zeros(N), om, True)
    l2, u2 = B.rates(n2, n1, om, False)
    for lam, mu, step, n in ((l1, u1, m, n1), (l2, u2, 1, n2)):
        up = (n < cap) & (lam > 0)
        dn = (n > 0) & (mu > 0)
        rows += [idx[up], idx[dn]]
        cols += [idx[up] + step, idx[dn] - step]
        vals += [lam[up], mu[dn]]
        diag -= np.where(up, lam, 0.0) + np.where(dn, mu, 0.0)
    rows.append(idx); cols.append(idx); vals.append(diag)
    Q = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(N, N)).tocsr()
    return Q, m


def tag(Q, m, sad1):
    """§118's dip bit: a stage-1 step from n1 >= sad1 to n1 < sad1 moves d = 0 -> 1, never back."""
    N = Q.shape[0]
    C = Q.tocoo()
    r, c, v = C.row, C.col, C.data
    dip = (r != c) & (r // m >= sad1) & (c // m < sad1)
    rows = [r[~dip], r[dip], r + N]
    cols = [c[~dip], c[dip] + N, c + N]
    vals = [v[~dip], v[dip], v]
    return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(2 * N, 2 * N)).tocsr()


def seed(A, B, om, m):
    k1, p1 = A.qsd(om, 0.0, first=True)
    mu1 = A.operating_points(om, 1)[0][0]
    k2, p2 = B.qsd(om, mu1, first=False)
    a = np.zeros(m); a[k1] = p1
    b = np.zeros(m); b[k2] = p2
    return np.outer(a, b).ravel(), mu1


def run(A, B, om, tagged=True):
    Q, m = pair_generator(A, B, om)
    sad1, sad2 = int(np.ceil(A.r2 * om)), int(np.ceil(B.r2 * om))
    p0, mu1 = seed(A, B, om, m)
    N = m * m
    T = tag(Q, m, sad1) if tagged else Q
    QT = (T.T * 2.0).tocsr()
    p = np.concatenate([p0, np.zeros(N)]) if tagged else p0.copy()
    s1hi = (np.arange(N) // m) >= sad1
    s2lo = (np.arange(N) % m) < sad2
    out, t = [], 0.0
    while t < MARKS[-1] - 1e-9:
        p = spla.expm_multiply(QT, p)
        t += 2.0
        if not any(abs(t - mk) < 1e-9 for mk in MARKS):
            continue
        blocks = [p[:N], p[N:]] if tagged else [p]
        hi = sum(float(b[s1hi].sum()) for b in blocks)
        rec = {"t": t, "p_free": sum(float(b[s1hi & s2lo].sum()) for b in blocks) / hi}
        if tagged:
            w0, w1 = float(p[:N][s1hi].sum()), float(p[N:][s1hi].sum())
            rec.update(w_d=w1 / hi, p0=float(p[:N][s1hi & s2lo].sum()) / w0,
                       p_d=float(p[N:][s1hi & s2lo].sum()) / w1 if w1 > 0 else float("nan"))
            rec["E"] = rec["p_free"] / rec["p0"] - 1.0
            rec["E_pred"] = rec["w_d"] * 0.85 / rec["p0"]
        out.append(rec)
    return {"series": out, "mu1": mu1,
            "k1": dip_rate(A, om, 0.0, True), "k2": dip_rate(B, om, mu1, False)}


ARMS = {"X": ("E-070", "E-084"), "S": ("E-070", "E-070"),
        "Z": ("E-062", "E-084"), "N": ("E-062", "F-025"),
        "a": ("E-070", "F-025"), "b": ("E-070", "F-040"),
        "c": ("E-062", "F-050"), "d": ("E-062", "F-040"),
        "XK": ("E-070", "E-084@K1.2")}                     # §124.9: X with its coupling threshold at 1.2
ARM_OMEGAS = {"X": (30, 40, 55), "S": (30, 40, 55), "Z": (40, 55), "N": (40, 55),
              "a": (40, 55), "b": (40, 55), "c": (40, 55), "d": (40, 55)}


def edge(el, lo=0.2, hi=3.2, n=30000):
    """Stage-2 bistability edge: the least upstream level at which the deterministic cubic has
    three positive roots (§105's x_up for the published element)."""
    k1a, k1r, k2b, k2r = el.c
    for x in np.linspace(lo, hi, n):
        r = np.roots([-k1r, el.hill(x) * k1a, -k2r, k2b])
        if sum(1 for z in r if abs(z.imag) < 1e-9 and z.real > 0) == 3:
            return float(x)
    return float("nan")


OUT_OF_SAMPLE = """
§119b -- OUT OF SAMPLE, written after reading arms X, S, Z, N and BEFORE running a, b, c, d.

What arms X/S/Z/N showed and did not predict: once the zone [saddle_1, edge_2) exists, stage 2's
low probability given stage 1 high, divided by k1 t, is flat in Omega (X 0.122/0.123/0.136, Z
0.025/0.021) while k2 moves 200-7000x; it is phi * p_d + p0/(k1 t) with phi = w_d/(k1 t) ~ 0.37-0.41,
and p_d orders with the zone width w = edge_2 - saddle_1: S (+0.611) 0.84-0.89, X (+0.357)
0.29-0.31, Z (+0.129) 0.054-0.066, N (-0.106) 0.0014 -> 0.0005. A post-hoc pattern is not a result.
These arms use pairs never run, and the predictions are stated here first.

  Q1  ZONE WIDTH ALONE. Arm a = E-070 -> F-025, w = +0.122, shares NEITHER element with Z
      (w = +0.129). Predict p_d(a) within 1.5x of p_d(Z) at each of Omega = 40, 55. Kill test for
      the criterion: if p_d depends on which elements and not only on w, a and Z differ in both
      elements, so a mismatch beyond 1.5x prints FAILS; a match can only come from w (or from
      something that happens to track it -- say so if it holds).
  Q2  MONOTONE IN w. Arm b = E-070 -> F-040, w = +0.225: p_d(Z-like ~0.06) < p_d(b) < p_d(X ~0.29).
      Arm c = E-062 -> F-050, w = +0.061: p_d(N) < p_d(c) < p_d(Z). Arm d = E-062 -> F-040,
      w = -0.004 (the edge sits on the saddle): p_d(N) < p_d(d) < p_d(c).
  Q3  THE FLOOR. For every arm with w > 0 (a, b, c), P(stage 2 low | stage 1 high)/(k1 t) at
      t = 16 changes by less than 1.5x from Omega = 40 to 55. Against it: p0/(k1 t) and k1/k2
      move by orders of magnitude over the same step in X and Z, so a floor that is not slaved to
      k1 would fail this by a wide margin.
  Q4  THE SIGN OF w. For w < 0 or ~0 (arm d, and N) p_d FALLS from Omega = 40 to 55; for w > 0
      clearly (a, b) it does not fall by more than 1.5x.
"""


CONVERGENCE = """
§119c -- the instrument, written BEFORE running. §113 found cap_mult 1.25 unconverged for some
quantities; §118/§119 inherited 1.25 from §106. And E divides by p0, which reaches 1e-13 in arm N.

  C1  THE BOX. Rerun X:40, a:40 and c:40 at cap_mult 2.0. Predict P(stage 2 low | stage 1 high),
      p_d and w_d each move by < 5% -- well under the 1.37x drifts Q3 reads. Stage 1's dips happen
      at its saddle, far below the cap; the cap only trims the high-state tail. If any moves by
      more than 5%, every §119 cell is re-run at 2.0 before anything is written.
  C2  THE FLOOR. Rerun the same cells with 1-unit steps instead of 2 (expm_multiply's own
      tolerance is relative to the vector's norm, ~1): p0 agrees to < 1% wherever p0 > 1e-9, and
      cells with p0 below that are reported as floor-limited, their E not quoted.
"""


def check(A, B, om, cap_mult, step):
    A2, B2 = Elem(A.name, A.r1, A.r2, A.r3, cap_mult), Elem(B.name, B.r1, B.r2, B.r3, cap_mult)
    Q, m = pair_generator(A2, B2, om)
    sad1, sad2 = int(np.ceil(A2.r2 * om)), int(np.ceil(B2.r2 * om))
    p0, _ = seed(A2, B2, om, m)
    N = m * m
    QT = (tag(Q, m, sad1).T * step).tocsr()
    p = np.concatenate([p0, np.zeros(N)])
    for _ in range(int(round(MARKS[-1] / step))):
        p = spla.expm_multiply(QT, p)
    s1hi = (np.arange(N) // m) >= sad1
    s2lo = (np.arange(N) % m) < sad2
    hi = float(p[:N][s1hi].sum() + p[N:][s1hi].sum())
    w1 = float(p[N:][s1hi].sum())
    return {"p_free": float(p[:N][s1hi & s2lo].sum() + p[N:][s1hi & s2lo].sum()) / hi,
            "w_d": w1 / hi, "p_d": float(p[N:][s1hi & s2lo].sum()) / w1,
            "p0": float(p[:N][s1hi & s2lo].sum()) / float(p[:N][s1hi].sum())}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path, default=pathlib.Path("results/weak_upstream.json"))
    ap.add_argument("--arms", default="XSZNabcd")
    args = ap.parse_args()
    rep = json.loads(args.out.read_text()) if args.out.exists() else {}
    for arm in args.arms:
        a, b = ARMS[arm]
        A, B = element(a), element(b)
        for om in ARM_OMEGAS[arm]:
            key = f"{arm}:{om}"
            if key in rep:
                continue
            r = run(A, B, om)
            rep[key] = r
            last = r["series"][-1]
            print(f"{arm} {a}->{b} Om={om}: k1 {r['k1']:.3e} k2 {r['k2']:.3e} k1/k2 "
                  f"{r['k1'] / r['k2']:.3e} | t=16 w_d {last['w_d']:.3e} p0 {last['p0']:.3e} "
                  f"p_d {last['p_d']:.4f} E {last['E']:.4e} E/E_pred {last['E'] / last['E_pred']:.3f}",
                  flush=True)
            args.out.write_text(json.dumps(rep, indent=2, default=float))
    if "convergence" not in rep:                       # §119c, C1 and C2
        conv = {}
        for key in ("X:40", "a:40", "c:40"):
            arm, om = key.split(":")
            A, B = element(ARMS[arm][0]), element(ARMS[arm][1])
            ref = rep[key]["series"][-1]
            for name, cm, st in (("box2.0", 2.0, 2.0), ("step1", CAP_MULT, 1.0)):
                r = check(A, B, int(om), cm, st)
                r["rel"] = {q: r[q] / ref[q] - 1 for q in ("p_free", "w_d", "p_d", "p0")}
                conv[f"{key}/{name}"] = r
                print(key, name, "  ".join(f"{q} {v:+.2e}" for q, v in r["rel"].items()), flush=True)
        rep["convergence"] = conv
        args.out.write_text(json.dumps(rep, indent=2, default=float))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
