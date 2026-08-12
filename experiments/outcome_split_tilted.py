"""T-TUR-e: §60 closed the reliability-dissipation door. Was it only closed by SYMMETRY?

§60 split §41's fluctuation-theorem identity by outcome and found it FACTORISES:
<e^(-S_tot) | o> = 1 for each outcome separately, so Phi_o = p_o and the identity carries no
information about the error rate. It closed the founding question's sharpest form.

**But §60's own explanation of why is symmetric-AM-specific.** It reasoned that the two
absorbing boundaries are EXCHANGE IMAGES carrying equal stationary weight, so the system term
cancels outcome by outcome. That premise is a fact about beta = 0. **For any tilted network --
that is, for every real device, since an inverter drives toward one rail and not the other --
the two boundaries are not exchange images and the cancellation has no reason to occur.**
§60's conclusion was stated for restoration in general. It was tested on a network whose
symmetry guarantees the answer.

This costs two solves per cell, on machinery that already exists: §60's tilted generator with
an outcome-selective boundary, run on `am_asymmetric` instead of `am_reversible`.

PREDICTIONS, written before running.

  P1  GATE. Phi_c + Phi_e = 1 at every beta, reproducing §41's identity. The IFT is an
      identity for ANY network, so a failure here is the pairing or the tilted generator,
      not physics, and nothing below counts. (`reverse_pairing` is derived by matching
      reactants to products, not hardcoded, and returns [3 4 5 0 1 2] on the tilted network.)
  P2  **THE TEST. Does the factorisation survive the tilt?** Report r_o = Phi_o / p_o for
      both outcomes against beta. **PREDICTED: it fails, and |r_o - 1| grows with beta**,
      because the mechanism §60 named is unavailable once the boundaries stop being exchange
      images. If instead r_o stays 1, the outcome-wise identity is DEEPER than exchange
      symmetry -- a theorem to prove rather than a number to file -- and that is the more
      interesting outcome of the two.
  P3  **THE MECHANISM IS A SUSPECT, NOT A RESULT (rule 17).** §60's account says the
      cancellation is about the two boundaries' stationary weights. That is checkable
      directly: compute w = pi(correct boundary) / pi(error boundary) and ask whether
      |r_o - 1| tracks ln w. **Reported as a correlation with its scatter, and the mechanism
      is named a suspect either way** -- this project has withdrawn three mechanisms proposed
      in the same breath as the measurement that motivated them.
  P4  **SCOPE GATE, and it is a hard one.** `am_asymmetric` is monostable above
      beta_c(gamma): past it the network "restores nothing, it just reports a constant",
      which would look like perfect fidelity to any metric that only ever feeds it X. Every
      cell must have beta < beta_c(gamma) strictly, the ratio beta/beta_c is reported per
      cell, and cells at or above it are excluded and counted, never fitted.
  P5  **RULE 9, an axis I did not choose.** beta is the axis the claim is about. Sweep Omega
      and the start offset too. A deviation that only appears at one Omega is a finite-size
      artifact; a deviation that is a property of the tilt must survive both.
  P6  **VERDICT RULE, and per §64's lesson it is unit-tested on synthetic data ENGINEERED to
      trigger each branch before this runs** (see tests/test_outcome_split_tilted.py). The
      rule must be able to print:
        (a) NO DEVIATION -- |r-1| at beta > 0 stays inside the beta = 0 residual, which §60
            measured at ~5%. Then §60 generalises and the factorisation is deeper than
            symmetry.
        (b) DEVIATION -- |r-1| grows monotonically in beta and exceeds the beta = 0 residual
            by more than the Omega-scatter. Then §60 is rescoped to symmetric elements and
            the founding question reopens where real devices live.
        (c) INCONCLUSIVE -- non-monotone, or growth inside the Omega-scatter.
      **The beta = 0 residual is measured in this same run, not taken from §60**, so the
      comparison is paired (rule 18) rather than against another run's number.

SECOND PASS, before the measurement is read. The first run failed P1 on 1 cell of 18
(beta = 0.40, Omega = 90: Phi_c + Phi_e = 0.813 against a median |sum - 1| of 4.2e-8), and
that one cell produced the entire apparent deviation. Diagnosis: **Phi_all itself is wrong
there, not the split**, and Phi_c tracks p_c to four digits in every cell while Phi_e degrades.
The reason is the instrument's own arithmetic: §41's convention weights the boundary by
pi(n)/pi(n0), and at beta = 0.40, Omega = 90 the error boundary carries
ln w = -38.9, i.e. e^-38.9 ~ 1e-17 of the correct boundary's stationary weight -- **below
double precision.** So:

  P1' PER-CELL GATE. |Phi_all - 1| < 1e-6 is now required cell by cell; failures are excluded
      and counted with their (beta, Omega, ln w), never averaged in.
  P7  **THE CONFOUND, and it is the dangerous kind: the nuisance grows with the cause.**
      Tilt is what makes the two boundaries' stationary weights differ, so |ln w| -- which
      is exactly what breaks the solve -- rises with beta, the very variable under test. A
      deviation rising with beta is therefore predicted by BOTH hypotheses. **Break it the
      way this project has broken confounds before: sweep the nuisance against the cause.**
      |ln w| also rises with Omega at fixed beta, so cells of matched |ln w| exist at
      different beta (large beta with small Omega against small beta with large Omega).
      **Compare |r - 1| across beta AT MATCHED |ln w|.** If it is flat there, the rise is
      conditioning and §60's closure survives as far as this instrument reaches; if it still
      rises with beta at matched |ln w|, the effect is the tilt.
      This also bounds the claim: the accessible beta range is limited by precision, and it
      is limited exactly where the predicted effect should be largest. That limit is part of
      the result, not a footnote.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import brentq

from crnl.approximations import propensities_batch
from crnl.cme import enumerate_states, generator, stationary
from crnl.networks.am_asymmetric import am_asymmetric, basin_boundary, beta_critical
from crnl.networks.am_reversible import reverse_pairing
from crnl.vectorized import compile_network


def slaved_b(net, delta, lo=1e-9, hi=1.0 - 1e-9):
    """b on the slow manifold at fixed delta = x - y, with x + y + b = 1.

    The `slaved` helper is written for the symmetric family; this is the same idea done
    generically, so no am_reversible closed form leaks into a tilted measurement.
    """
    S = net.stoichiometry_matrix()

    def db(b):
        x = (1.0 - b + delta) / 2.0
        y = (1.0 - b - delta) / 2.0
        if x < 0 or y < 0:
            return np.nan
        return float((S @ net.fluxes(np.array([x, y, b])))[2])

    xs = np.linspace(lo, hi, 400)
    vals = [db(b) for b in xs]
    for i in range(len(xs) - 1):
        a, c = vals[i], vals[i + 1]
        if np.isfinite(a) and np.isfinite(c) and a * c < 0:
            return float(brentq(db, xs[i], xs[i + 1], xtol=1e-14))
    return None


def cell(gamma, beta, omega, off, theta):
    """§60's outcome split, on the tilted network. Start offset from the SEPARATRIX."""
    bc = beta_critical(gamma)
    if beta >= bc:
        return {"skipped": "monostable", "beta_over_bc": beta / bc}
    net = am_asymmetric(gamma, beta)
    pairing = reverse_pairing(net)

    sep = basin_boundary(gamma, beta) if beta > 0 else 0.0
    x0 = sep - off                      # marginal start, on the Y side of the separatrix
    b0 = slaved_b(net, x0)
    if b0 is None:
        return None
    nb = int(round(b0 * omega))
    rest = omega - nb
    d0 = int(round(x0 * omega))
    if (rest - d0) % 2:
        d0 -= 1
    thr = max(2, int(round(theta * omega)))
    nx, ny = (rest + d0) // 2, (rest - d0) // 2
    if min(nx, ny, nb) < 0 or abs(d0) >= thr:
        return None
    n0 = np.array([nx, ny, nb], dtype=np.int64)

    states, index = enumerate_states(3, int(omega))
    dvec = np.array([int(s[0]) - int(s[1]) for s in states])
    absorb = np.abs(dvec) >= thr
    correct = absorb & (dvec <= -thr)      # started with Y ahead: Y is the correct answer
    error = absorb & (dvec >= thr)

    comp = compile_network(net, float(omega))
    S = net.stoichiometry_matrix().astype(np.int64)
    A = propensities_batch(comp, states.astype(float))
    a_tot = A.sum(axis=1)

    Q = generator(net, int(omega), float(omega))
    tr = np.where(~absorb)[0]
    tmap = -np.ones(len(states), dtype=np.int64)
    tmap[tr] = np.arange(len(tr))
    Qtt = Q[tr][:, tr].tocsr()
    Qta = Q[tr][:, np.where(absorb)[0]].tocsr()
    si = int(tmap[index[tuple(n0)]])
    p_e = float(spla.spsolve(Qtt, -(Qta @ error[absorb].astype(float)))[si])
    p_c = float(spla.spsolve(Qtt, -(Qta @ correct[absorb].astype(float)))[si])

    rows, cols, vals = [], [], []
    for j in range(net.n_reactions):
        rev = int(pairing[j])
        if rev < 0:
            continue
        n2 = states + S[:, j]
        ok = (n2 >= 0).all(axis=1) & (A[:, j] > 0)
        if not ok.any():
            continue
        idx = np.where(ok)[0]
        arev = propensities_batch(comp, n2[idx].astype(float))[:, rev]
        good = arev > 0
        src = idx[good]
        dst = np.array([index[tuple(x)] for x in n2[src]])
        rows.append(src); cols.append(dst); vals.append(arev[good])
    rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
    keep = ~absorb[rows]
    r, c, v = rows[keep], cols[keep], vals[keep]
    inner = ~absorb[c]
    M = (sp.coo_matrix((-v[inner], (tmap[r[inner]], tmap[c[inner]])),
                       shape=(len(tr), len(tr))).tocsr()
         + sp.diags(a_tot[tr])).tocsc()

    pi = stationary(net, int(omega), float(omega))
    ratio = pi / max(pi[index[tuple(n0)]], 1e-300)

    out, edge = {}, ~inner
    for name, mask in (("all", absorb), ("c", correct), ("e", error)):
        bnd = ratio * mask.astype(float)
        b = np.zeros(len(tr))
        np.add.at(b, tmap[r[edge]], v[edge] * bnd[c[edge]])
        out[name] = float(spla.spsolve(M, b)[si])

    # P3: the boundaries' stationary weights, which §60's mechanism says do the work
    wc = float(pi[correct].sum())
    we = float(pi[error].sum())
    return {"gamma": gamma, "beta": beta, "omega": omega, "off": off, "theta": theta,
            "beta_over_bc": beta / bc, "sep": float(sep), "d0": int(d0), "thr": int(thr),
            "p_c": p_c, "p_e": p_e, "Phi_all": out["all"], "Phi_c": out["c"],
            "Phi_e": out["e"], "w_c": wc, "w_e": we,
            "lnw": float(np.log(max(wc, 1e-300) / max(we, 1e-300)))}


def decide(betas, dev, dev0, omega_scatter):
    """P6's verdict. Unit-tested on engineered data in tests/test_outcome_split_tilted.py.

    betas/dev are paired and sorted; dev0 is |r-1| measured at beta = 0 IN THE SAME RUN;
    omega_scatter is the spread of dev across Omega at the largest beta.
    """
    d = np.asarray(dev, dtype=float)
    floor = max(dev0, omega_scatter)
    grew = bool(np.all(np.diff(d) > -1e-12)) and d[-1] > d[0]
    # The two floors are not the same kind of thing. dev0 is a physical baseline, so a
    # signal inside it means NO EFFECT; omega_scatter is instrument noise, so a signal
    # inside it means NOT MEASURABLE. Conflating them let the first draft print "no
    # deviation" on data that was merely too noisy to read (caught by the unit tests).
    noisy = omega_scatter > 2.0 * max(dev0, 1e-12)
    if noisy and d.max() <= floor:
        return "c", (f"INCONCLUSIVE: the Omega scatter ({omega_scatter:.4f}) dominates the "
                     f"beta=0 residual ({dev0:.4f}) and |r-1| peaks at only {d.max():.4f}. "
                     f"Too noisy to say either way.")
    if d.max() <= 2.0 * floor:
        return "a", (f"NO DEVIATION: |r-1| peaks at {d.max():.4f}, within 2x the floor "
                     f"{floor:.4f} (beta=0 residual {dev0:.4f}, Omega scatter "
                     f"{omega_scatter:.4f}). §60 generalises beyond exchange symmetry.")
    if grew and d[-1] > 2.0 * floor:
        return "b", (f"DEVIATION: |r-1| rises monotonically to {d[-1]:.4f}, clear of 2x the "
                     f"floor {floor:.4f}. §60 must be rescoped to symmetric elements.")
    return "c", (f"INCONCLUSIVE: |r-1| reaches {d.max():.4f} against a floor {floor:.4f} "
                 f"but is {'non-monotone' if not grew else 'not clear of it'}.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.25)
    ap.add_argument("--betas", type=float, nargs="+",
                    default=[0.0, 0.05, 0.10, 0.20, 0.30, 0.40])
    ap.add_argument("--omegas", type=int, nargs="+", default=[40, 60, 90])
    ap.add_argument("--off", type=float, default=0.06)
    ap.add_argument("--offs", type=float, nargs="+", default=[0.04, 0.06, 0.09])
    ap.add_argument("--theta", type=float, default=0.55)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/outcome_split_tilted.json"))
    args = ap.parse_args()
    t0 = time.time()
    g = args.gamma
    print(f"gamma = {g}, beta_c = {beta_critical(g):.4f}  (P4: every cell must stay below)")

    rows, skipped, bad = [], 0, []
    print(f"\n{'beta':>6}{'b/bc':>7}{'Om':>5}{'p_e':>11}{'Phi_c':>11}{'Phi_e':>10}"
          f"{'sum':>9}{'Phi_c/p_c':>11}{'Phi_e/p_e':>11}")
    for beta in args.betas:
        for om in args.omegas:
            try:
                r = cell(g, beta, om, args.off, args.theta)
            except Exception as e:
                print(f"{beta:>6.2f}{'':>7}{om:>5}   FAILED ({type(e).__name__})")
                continue
            if r is None or "skipped" in r:
                skipped += 1
                continue
            s = r["Phi_c"] + r["Phi_e"]
            r["sum"] = s
            if abs(r["Phi_all"] - 1.0) > 1e-6:        # P1' per-cell gate
                bad.append(r)
                print(f"{beta:>6.2f}{r['beta_over_bc']:>7.2f}{om:>5}"
                      f"   EXCLUDED by P1': |Phi_all - 1| = {abs(r['Phi_all']-1):.3e},"
                      f" ln w = {r['lnw']:.1f} (e^{r['lnw']:.0f} is past double precision)")
                continue
            r["rc"] = r["Phi_c"] / r["p_c"] if r["p_c"] > 0 else float("nan")
            r["re"] = r["Phi_e"] / r["p_e"] if r["p_e"] > 0 else float("nan")
            rows.append(r)
            print(f"{beta:>6.2f}{r['beta_over_bc']:>7.2f}{om:>5}{r['p_e']:>11.4e}"
                  f"{r['Phi_c']:>11.4e}{r['Phi_e']:>10.6f}{s:>9.6f}"
                  f"{r['rc']:>11.6f}{r['re']:>11.6f}")
    if not rows:
        print("no cells")
        return
    print(f"  P4: {skipped} cells skipped (monostable or unbuildable);"
          f" P1': {len(bad)} excluded on precision")

    print(f"\n=== P1 GATE: does Phi_c + Phi_e reconstruct §41's 1, on surviving cells?")
    dd = np.array([abs(r["sum"] - 1.0) for r in rows])
    print(f"  |sum - 1| over {len(rows)} cells: max {dd.max():.3e}, median {np.median(dd):.3e}")
    print(f"  -> P1 {'HOLDS' if dd.max() < 1e-6 else 'FAILS'}")

    print(f"\n=== P2: does the outcome-wise factorisation survive the tilt?")
    print(f"{'beta':>6}{'|rc-1|':>12}{'|re-1|':>12}{'cells':>7}")
    per_beta, betas_used = [], []
    for beta in args.betas:
        sel = [r for r in rows if r["beta"] == beta]
        if not sel:
            continue
        dc = float(np.mean([abs(r["rc"] - 1) for r in sel]))
        de = float(np.mean([abs(r["re"] - 1) for r in sel]))
        per_beta.append(max(dc, de))
        betas_used.append(beta)
        print(f"{beta:>6.2f}{dc:>12.6f}{de:>12.6f}{len(sel):>7}")

    dev0 = per_beta[0] if betas_used and betas_used[0] == 0.0 else 0.0
    top = [r for r in rows if r["beta"] == betas_used[-1]]
    om_scatter = float(np.ptp([max(abs(r["rc"] - 1), abs(r["re"] - 1)) for r in top])) \
        if len(top) > 1 else 0.0
    print(f"  beta=0 residual (measured HERE, paired): {dev0:.6f};"
          f" Omega scatter at beta={betas_used[-1]}: {om_scatter:.6f}")
    code, msg = decide(betas_used, per_beta, dev0, om_scatter)
    print(f"  -> P2 ({code}) {msg}")

    print(f"\n=== P3 (rule 17): does the deviation track the boundary weight ln w?")
    lw = np.array([r["lnw"] for r in rows])
    dv = np.array([max(abs(r["rc"] - 1), abs(r["re"] - 1)) for r in rows])
    ok = np.isfinite(lw) & np.isfinite(dv)
    if ok.sum() >= 4 and np.ptp(lw[ok]) > 0:
        cc = float(np.corrcoef(np.abs(lw[ok]), dv[ok])[0, 1])
        print(f"  |ln w| spans {np.abs(lw[ok]).min():.3f}..{np.abs(lw[ok]).max():.3f};"
              f" correlation with |r-1| = {cc:+.3f} over {ok.sum()} cells")
        print(f"  -> the boundary-weight account is a SUSPECT"
              f" ({'consistent' if abs(cc) > 0.7 else 'not supported'}), not a result")

    print(f"\n=== P5 (rule 9): the start offset, an axis the claim is not about")
    print(f"{'off':>7}{'beta':>7}{'|rc-1|':>12}{'|re-1|':>12}")
    offrows = []
    for off in args.offs:
        for beta in (0.0, args.betas[-1]):
            r = cell(g, beta, args.omegas[-1], off, args.theta)
            if r is None or "skipped" in r or r["p_c"] <= 0 or r["p_e"] <= 0:
                continue
            if abs(r["Phi_all"] - 1.0) > 1e-6:      # P1' applies here too
                print(f"{off:>7.3f}{beta:>7.2f}   EXCLUDED by P1'"
                      f" (|Phi_all-1| = {abs(r['Phi_all']-1):.2e}, ln w = {r['lnw']:.1f})")
                continue
            dc, de = abs(r["Phi_c"] / r["p_c"] - 1), abs(r["Phi_e"] / r["p_e"] - 1)
            offrows.append({"off": off, "beta": beta, "dc": dc, "de": de})
            print(f"{off:>7.3f}{beta:>7.2f}{dc:>12.6f}{de:>12.6f}")
    if offrows:
        hi = [r for r in offrows if r["beta"] > 0]
        lo = [r for r in offrows if r["beta"] == 0]
        if hi and lo:
            print(f"  mean |r-1|: beta=0 {np.mean([max(r['dc'],r['de']) for r in lo]):.6f}"
                  f" vs beta={args.betas[-1]}"
                  f" {np.mean([max(r['dc'],r['de']) for r in hi]):.6f}"
                  f"  -> {'survives the offset sweep' if np.mean([max(r['dc'],r['de']) for r in hi]) > 2*np.mean([max(r['dc'],r['de']) for r in lo]) else 'does NOT survive: offset-dependent'}")

    print(f"\n=== P7: the confound -- |ln w| rises with beta AND with Omega.")
    print(f"    Compare |r-1| across beta AT MATCHED |ln w| (nuisance swept against cause).")
    lwv = np.array([abs(r["lnw"]) for r in rows])
    dvv = np.array([abs(r["re"] - 1) for r in rows])
    bv = np.array([r["beta"] for r in rows])
    ov = np.array([r["omega"] for r in rows])
    print(f"  correlation |r-1| with |ln w|: {np.corrcoef(lwv, dvv)[0,1]:+.3f};"
          f" with beta: {np.corrcoef(bv, dvv)[0,1]:+.3f}"
          f"  (both predicted positive by BOTH hypotheses)")
    edges = np.quantile(lwv, [0.0, 1/3, 2/3, 1.0])
    verdicts = []
    for k in range(3):
        m = (lwv >= edges[k]) & (lwv <= edges[k + 1]) & (bv > 0)
        if m.sum() < 3 or np.ptp(bv[m]) <= 0:
            continue
        cc = float(np.corrcoef(bv[m], dvv[m])[0, 1])
        verdicts.append(cc)
        print(f"  |ln w| in [{edges[k]:.1f}, {edges[k+1]:.1f}]: {m.sum()} cells,"
              f" beta {bv[m].min():.2f}-{bv[m].max():.2f},"
              f" Omega {ov[m].min():.0f}-{ov[m].max():.0f},"
              f" corr(beta, |r-1|) = {cc:+.3f}")
    if verdicts:
        mc = float(np.mean(verdicts))
        ncell = int(((bv > 0)).sum())
        print(f"  mean within-band corr(beta, |r-1|) = {mc:+.3f} over ~{ncell//3} cells/band")
        # A bare threshold on a correlation is the rule-19 trap this project keeps falling
        # into. Report the strength, and say plainly what it does and does not settle.
        if mc > 0.8:
            print("  -> the rise clearly SURVIVES matched conditioning: the tilt is doing it")
        elif mc > 0.3:
            print(f"  -> SUGGESTIVE ONLY: the rise partly survives matched conditioning "
                  f"(r = {mc:+.3f} on ~{ncell//3} cells per band), which is too weak to "
                  f"separate tilt from solve error. NOT settled either way.")
        else:
            print("  -> the rise does NOT survive matched conditioning: it is the solve")
        print("  NOTE, and it does not depend on the above: §60 claimed a GENERAL closure "
              "and tested only beta = 0. Its scope must be corrected regardless of which "
              "hypothesis wins here.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"rows": rows, "per_beta": per_beta,
                                    "betas": betas_used, "dev0": dev0,
                                    "omega_scatter": om_scatter, "verdict": code,
                                    "excluded_precision": bad,
                                    "offrows": offrows}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
