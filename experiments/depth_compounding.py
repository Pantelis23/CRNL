"""T-CASC-g: does the two-number description survive depth, and does the margin ERODE?

§91-§93 established that the composition penalty is a function of the noise margin and the
timescale ratio, both single-element quantities -- at D = 2. The reason that might not survive
depth is specific and predictable:

**the margin is measured in units of the UPSTREAM RAIL WIDTH, and each stage's output is wider
than its input.** Stage 2 inherits stage 1's fluctuation on top of its own, so stage 3 sees a
broader upstream, a SMALLER effective margin, and should fail more. If so the margin ERODES with
depth, the penalty compounds, and D_max is far below what a per-stage constant would give.

**That turns the two-number description into a PREDICTION rather than a fit** (rule 16): measure
stage 2's output width, convert it to stage 3's effective margin, put it through §91's measured law
log(penalty) = -0.952 x margin/sigma, and check against an exact D = 3 solve.

DESIGN. To compare stages apples-to-apples, every stage except the last is REFLECTED at its saddle
-- it fluctuates with whatever width it has inherited but cannot itself flip, so the last stage's
error is purely its own and nothing is conditioned (§92.1(b) is why this matters). Then:

    D = 2:  stage 1 reflected, stage 2 free   ->  penalty_2 = P(s2 low)/eps_iso
    D = 3:  stages 1,2 reflected, stage 3 free ->  penalty_3 = P(s3 low)/eps_iso

with the same eps_iso (the isolated element, upstream pinned at r3) in both.

PREDICTIONS, written before running.

  P1  GATE. The D = 3 chain marginalised over stage 3 must reproduce the D = 2 chain's joint law
      for stages 1-2 exactly -- the coupling is one-way, so adding a downstream stage cannot
      change anything upstream. **A violation means the generator is wired with back-action** and
      nothing below counts.
  P2  **THE WIDTH ERODES.** sigma(stage 2 output) > sigma(stage 1 output), both measured on the
      high rail from the exact joint law. **Predicted: yes**, since stage 2's rail fluctuation is
      its own plus what it inherits. Report the ratio, which is the erosion per stage.
  P3  **THE TEST, ABSOLUTE.** Predict penalty_3 from §91's law using the effective margin
      (r3 - x_crit)/sigma_2 measured in P2, and compare against the exact D = 3 value.
      **Predicted: penalty_3 > penalty_2, and the law predicts it within the ~18% scatter §91
      reported for its own collapse.** If the law predicts it, composition is computable to any
      depth from single-element quantities. If penalty_3 comes out EQUAL to penalty_2 the margin
      does not erode and the per-stage penalty is constant, which is a better world and a
      different result.
  P4  **WHAT IT MEANS FOR DEPTH -- and the naive reading is WRONG, worked out before running.**
      A first pass took the width ratio sigma_2/sigma_1 = 1.0855 and extrapolated it geometrically:
      the margin halves every 8.5 stages, so there is a hard maximum depth independent of Omega.
      **That cannot be right, because it contradicts what a restoring element IS.** In the LNA the
      output variance is

          sigma_out^2 = sigma_intr^2 + g^2 sigma_in^2

      with g the transfer gain at the operating point, and a RESTORING element attenuates input
      fluctuation, i.e. g < 1. So the widths do not grow geometrically -- they converge to a fixed
      point sigma_intr^2/(1-g^2), and the per-stage penalty converges to a CONSTANT. Extrapolating
      a ratio measured between the first two stages, where stage 1 uniquely has a noiseless
      chemostat for an input, is precisely rule 15's error.
      **The test:** take g^2 = (sigma_2/sigma_1)^2 - 1 from the first two stages, PREDICT
      sigma_3/sigma_1 = sqrt(1 + g^2 + g^4) and the limit 1/sqrt(1-g^2), and check both against
      the measured third stage. **Predicted: the widths converge, so there is no hard depth limit
      from erosion** -- and D_max is set by accumulated error at a constant per-stage penalty.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import experiments.chemical_cascade as cc
from experiments.cascade_schlogl import schlogl_consts
from experiments.margin_law import R1, R2, R3, stage1_stationary
from experiments.timescale_ratio import pinned_reference

C = schlogl_consts(R1, R2, R3)


def _hill_vec(x_up):
    f = lambda z: np.power(np.maximum(z, 0.0), cc.HILL_N) / (
        cc.HILL_K ** cc.HILL_N + np.power(np.maximum(z, 0.0), cc.HILL_N))
    return f(x_up) / f(R3)


def _rates_vec(n, n_up, om, first):
    """Vectorised propensities -- the scalar loop was too slow at 1.16M joint states."""
    k1a, k1r, k2b, k2r = C
    auto = k1a * n * (n - 1.0) / om
    mu = k1r * n * (n - 1.0) * (n - 2.0) / om ** 2 + k2r * n
    lam = (auto + k2b * om) if first else (_hill_vec(n_up / om) * auto + k2b * om)
    return np.maximum(lam, 0.0), np.maximum(mu, 0.0)


def build_chain(om, D, cap_mult=1.25, all_reflected=False):
    """D stages, hill coupling, every stage but the LAST reflected at its saddle.

    all_reflected=True reflects the last one too -- no free stage, so every width is measured
    without the escape-conditioning that truncates a free stage's low tail."""
    cap = int(np.ceil(cap_mult * R3 * om))
    nsad = int(np.ceil(R2 * om))
    ref = np.arange(nsad, cap + 1)
    nr = len(ref)
    mf = nr if all_reflected else cap + 1
    dims = [nr] * (D - 1) + [mf]
    strides = [int(np.prod(dims[i + 1:])) for i in range(D)]
    N = int(np.prod(dims))
    idx = np.arange(N)
    counts = []
    for i in range(D):
        ni = (idx // strides[i]) % dims[i]
        counts.append(ref[ni].astype(float) if (i < D - 1 or all_reflected)
                      else ni.astype(float))
    rows, cols, vals = [], [], []
    diag = np.zeros(N)
    for i in range(D):
        up = counts[i - 1] if i else np.zeros(N)
        lam, mu = _rates_vec(counts[i], up, om, i == 0)
        hi = float(cap)
        lo = float(ref[0]) if (i < D - 1 or all_reflected) else 0.0
        up_ok = (counts[i] < hi) & (lam > 0)
        dn_ok = (counts[i] > lo) & (mu > 0)
        rows.append(idx[up_ok]); cols.append(idx[up_ok] + strides[i]); vals.append(lam[up_ok])
        rows.append(idx[dn_ok]); cols.append(idx[dn_ok] - strides[i]); vals.append(mu[dn_ok])
        diag -= np.where(up_ok, lam, 0.0) + np.where(dn_ok, mu, 0.0)
    rows.append(idx); cols.append(idx); vals.append(diag)
    Q = sp.csr_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
    return Q, ref, dims, strides, cap


def seed(om, ref, dims, strides, pi1):
    """Stage 1 from its exact stationary law; every later stage starts at its rail."""
    p = np.zeros(int(np.prod(dims)))
    n_hi = int(round(R3 * om))
    rest = sum(strides[i] * (list(ref).index(n_hi) if i < len(dims) - 1 else n_hi)
               for i in range(1, len(dims)))
    for a, w in enumerate(pi1):
        p[a * strides[0] + rest] = w
    return p


def last_low(p, om, dims, strides):
    idx = np.arange(len(p))
    n_last = (idx // strides[-1]) % dims[-1]
    return float(p[n_last < R2 * om].sum())


def stage_stats(p, om, ref, dims, strides, k, all_reflected=False):
    """Mean and sd of stage k's count on its high side, from the exact joint law."""
    idx = np.arange(len(p))
    n = (idx // strides[k]) % dims[k]
    counts = ref[n] if (k < len(dims) - 1 or all_reflected) else n
    m = counts > R2 * om
    w = p[m] / p[m].sum()
    c = counts[m].astype(float)
    mu = float((w * c).sum())
    return mu / om, float(np.sqrt((w * (c - mu) ** 2).sum())) / om


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omega", type=int, default=30)
    ap.add_argument("--t", type=float, default=2.0)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/depth_compounding.json"))
    args = ap.parse_args()
    om, t = args.omega, args.t
    cc.HILL_N, cc.HILL_K = 4.0, 1.0
    _, pi1 = stage1_stationary(om)
    out = {}

    eps_iso = pinned_reference(om, 1.0, t)
    xs = np.linspace(R1, R3, 4001)
    xc = next(float(x) for x in xs[::-1] if len(cc.downstream_roots(x, C, R3, "hill")) < 3)
    print(f"eps_iso = {eps_iso:.5e};  collapse at x_crit = {xc:.4f}")

    res = {}
    for D in (2, 3):
        Q, ref, dims, strides, cap = build_chain(om, D)
        p0 = seed(om, ref, dims, strides, pi1)
        pT = spla.expm_multiply(Q.T * t, p0)
        res[D] = {"p": pT, "ref": ref, "dims": dims, "strides": strides, "N": len(p0)}
        print(f"  D={D}: {len(p0)} states")

    print("\n=== P1 GATE: does adding a downstream stage change anything upstream?")
    a2 = stage_stats(res[2]["p"], om, res[2]["ref"], res[2]["dims"], res[2]["strides"], 0)
    a3 = stage_stats(res[3]["p"], om, res[3]["ref"], res[3]["dims"], res[3]["strides"], 0)
    print(f"  stage 1 (mean, sd) at D=2: {a2[0]:.6f}, {a2[1]:.6f}")
    print(f"  stage 1 (mean, sd) at D=3: {a3[0]:.6f}, {a3[1]:.6f}")
    ok = abs(a2[0] - a3[0]) < 1e-9 and abs(a2[1] - a3[1]) < 1e-9
    print(f"  -> P1 {'HOLDS: the coupling is one-way, as built' if ok else 'FAILS: back-action -- the generator is wired wrong'}")
    assert ok

    print("\n=== P2: does the rail WIDTH erode from stage to stage?")
    print(f"{'stage':>7}{'mean':>11}{'sd':>11}{'margin in sd':>15}")
    widths = []
    for k in range(2):
        mu, sd = stage_stats(res[3]["p"], om, res[3]["ref"], res[3]["dims"],
                             res[3]["strides"], k)
        widths.append(sd)
        print(f"{k+1:>7}{mu:>11.5f}{sd:>11.5f}{(R3-xc)/sd:>15.3f}")
    erosion = widths[1] / widths[0]
    print(f"  width ratio stage2/stage1 = {erosion:.4f}")
    print(f"  -> P2 {'HOLDS: the output is wider than the input, so the effective margin ERODES' if erosion > 1.02 else ('the output is NARROWER -- the element is squeezing its input, and the margin IMPROVES with depth' if erosion < 0.98 else 'the width is unchanged to 2%, so the margin does not erode and the per-stage penalty should be constant')}")

    print("\n=== P3: the absolute test -- predict stage 3 from stage 2's measured width")
    pen2 = last_low(res[2]["p"], om, res[2]["dims"], res[2]["strides"]) / eps_iso
    pen3 = last_low(res[3]["p"], om, res[3]["dims"], res[3]["strides"]) / eps_iso
    m2, m3 = (R3 - xc) / widths[0], (R3 - xc) / widths[1]
    pred3 = pen2 * np.exp(-0.952 * (m3 - m2))
    print(f"  penalty_2 (measured)              = {pen2:.4f}   at margin {m2:.3f} sd")
    print(f"  penalty_3 PREDICTED from §91's law = {pred3:.4f}   at margin {m3:.3f} sd")
    print(f"  penalty_3 (measured)              = {pen3:.4f}")
    print(f"  predicted/measured = {pred3/pen3:.4f}")
    out.update({"eps_iso": eps_iso, "x_crit": xc, "widths": widths, "erosion": erosion,
                "pen2": pen2, "pen3": pen3, "pred3": pred3, "m2": m2, "m3": m3})
    within = abs(pred3 / pen3 - 1) < 0.20
    grew = pen3 > pen2
    print(f"  -> P3 {'HOLDS: the two-number law predicts the next stage from the previous one, within §91s own 18% scatter. Composition is computable to any depth from single-element quantities' if within else 'FAILS: the law does not carry to the next stage'}")
    print(f"     (penalty {'GREW' if grew else 'did not grow'} with depth: {pen2:.3f} -> {pen3:.3f})")

    print("\n=== P4: does the width CONVERGE, as a restoring element requires?")
    mu3, sd3 = stage_stats(res[3]["p"], om, res[3]["ref"], res[3]["dims"],
                           res[3]["strides"], 2)
    g2 = (widths[1] / widths[0]) ** 2 - 1.0
    pred3 = widths[0] * np.sqrt(1.0 + g2 + g2 ** 2)
    limit = widths[0] / np.sqrt(1.0 - g2) if g2 < 1 else float("inf")
    print(f"  transfer gain from stages 1-2:  g^2 = {g2:.4f}  (g = {np.sqrt(g2):.4f})")
    print(f"  sigma_1 = {widths[0]:.5f}   sigma_2 = {widths[1]:.5f}"
          f"   sigma_3 = {sd3:.5f} (measured)")
    print(f"  sigma_3 PREDICTED from the LNA recursion = {pred3:.5f}"
          f"   ratio {pred3/sd3:.4f}")
    print(f"  fixed point sigma_inf = {limit:.5f}"
          f"  ({limit/widths[0]:.4f} x stage 1, vs {widths[1]/widths[0]:.4f} at stage 2)")
    out.update({"sd3": sd3, "g2": g2, "sd3_pred": pred3, "sd_limit": limit})
    conv = g2 < 1.0 and abs(pred3 / sd3 - 1) < 0.05
    print(f"  -> P4 PROVISIONAL {'(g < 1 would give convergence to a fixed point)' if conv else '(the recursion already misses the third width)'}"
          f" -- but sigma_3 here is measured on the FREE stage, under escape-conditioning that")
    print(f"     truncates its low tail. P4b re-measures it unconditioned and DELIVERS the verdict.")
    print("\n  P4b: sigma_3 again with the LAST stage REFLECTED too, so no width is measured")
    print("       under escape-conditioning (a free stage's low tail is truncated, which biases")
    print("       its width DOWN -- and the free-stage sigma_3 above came out BELOW sigma_2,")
    print("       which the recursion cannot produce).")
    Qa, refa, dimsa, stridesa, _ = build_chain(om, 3, all_reflected=True)
    pa = seed(om, refa, dimsa, stridesa, pi1)
    pa = spla.expm_multiply(Qa.T * t, pa)
    sds = [stage_stats(pa, om, refa, dimsa, stridesa, k, all_reflected=True)[1]
           for k in range(3)]
    g2b = (sds[1] / sds[0]) ** 2 - 1.0
    pred3b = sds[0] * np.sqrt(1.0 + g2b + g2b ** 2)
    print(f"       sigma = " + ", ".join(f"{v:.5f}" for v in sds)
          + f"   (monotone: {all(sds[i+1] > sds[i] for i in range(2))})")
    print(f"       g^2 = {g2b:.4f};  sigma_3 predicted {pred3b:.5f} vs measured {sds[2]:.5f}"
          f"   ratio {pred3b/sds[2]:.4f}")
    out["all_reflected"] = {"sd": sds, "g2": g2b, "sd3_pred": pred3b}
    okb = abs(pred3b / sds[2] - 1) < 0.05 and all(sds[i + 1] > sds[i] for i in range(2))
    print(f"       -> P4b {'HOLDS: the widths rise monotonically toward the fixed point and the LNA recursion predicts the third from the first two, so the per-stage penalty tends to a constant' if okb else 'FAILS: the recursion does NOT predict the third width. The erosion ACCELERATES, so the widths do not converge to the LNA fixed point and BOTH simple readings -- geometric runaway and constant-gain convergence -- are refuted. The width sequence is not predicted by anything here.'}")
    print(f"  ** the geometric reading -- 'the margin halves every 8.5 stages, so depth is capped'")
    print(f"     -- is what a single ratio extrapolates to, and it is wrong (rule 15). Stage 1 is")
    print(f"     the ONLY stage whose input is a noiseless chemostat, so the 1->2 step is the")
    print(f"     least representative one in the chain. **")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
