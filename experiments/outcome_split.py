"""T-TUR-d: split the fluctuation theorem by OUTCOME — is reliability exactly dissipation?

§41 verified `<e^(-S_tot)> = 1` at the absorption time to 5.5e-14, validating sigma_local,
the reverse pairing and the sign convention against an exact identity. **It never split that
1 by outcome**, and the split is the founding question in its sharpest possible form.

Write the absorption-time average as a sum over the two outcomes,

    <e^(-S_tot)> = Phi_c + Phi_e,     Phi_o = p_o * <e^(-S_tot) | o>

Error trajectories run AGAINST the drive, so their entropy production is negative and
e^(-S) is exponentially LARGE; correct ones produce a lot of entropy and contribute almost
nothing. So Phi_e carries essentially the whole of the 1, and p_e is set by how much entropy
the error paths fail to produce. **If the split is exact in the natural way,**

    **Phi_e = p_c  and  Phi_c = p_e**,  equivalently  **<e^(-S_tot) | error> = p_c / p_e**

which says the ODDS OF BEING RIGHT are the exponentiated entropy along the error paths --
an exact identity between reliability and dissipation, which is what §37 set out to find and
§38 concluded was not a quantity in the form it was then asked.

**EVERYTHING NEEDED ALREADY EXISTS.** §41's tilted generator is built from reverse
propensities (`a_j e^(-ds_j) = a_rev(n')`), and §35's direct solve names an outcome at the
boundary to avoid `1 - split` cancellation. Combining them -- the tilted generator with an
OUTCOME-SELECTIVE boundary -- gives Phi_c and Phi_e as two more solves of the same system.
No new machinery, no Monte Carlo.

PREDICTIONS, written before running:

  P1  GATE. Phi_c + Phi_e = 1, reproducing §41 to solver precision. Computed as two separate
      solves whose sum must reconstruct the identity §41 measured in one. If it does not,
      the outcome-selective boundary is wrong and nothing below counts.
  P2  **THE TEST. Phi_e = p_c and Phi_c = p_e, exactly.** One relation, two ways of checking
      it, since Phi_c + Phi_e = 1 and p_c + p_e = 1 make them equivalent -- so agreement on
      one and disagreement on the other would itself indicate a coding error rather than a
      physical result.
  P3  If P2 holds, `<e^(-S_tot)|error> = p_c/p_e` and the log-odds of a correct decision are
      exactly ln<e^(-S)|error>. **That is an exact cost-of-reliability relation** and it
      would supersede §37's R = Sigma/L as the answer to the founding question.
  P4  **IF P2 FAILS, the ratios are the result.** Report Phi_e/p_c and Phi_c/p_e across
      gamma, Omega and eps. **A constant ratio is still an exact relation** with a
      coefficient to explain; a drifting one means outcome and entropy do not factorise this
      way and the identity is only the aggregate §41 already had.
  P5  Rule 9 applies to the confirming case too. If P2 holds it must hold along EVERY axis
      swept -- gamma, Omega and eps -- not just the one it was noticed on. An identity that
      holds on one axis and drifts on another is not an identity.
  P6  The Jensen direction is fixed and is a sanity check: <e^(-S)|e> >= e^(-<S|e>), so if
      P2 holds then p_c/p_e >= e^(-<S|e>), i.e. **<S|error> >= -ln(p_c/p_e)**. Reported, so
      that a violated inequality flags a sign error in sigma_local rather than being read as
      physics.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from crnl.approximations import propensities_batch
from crnl.cme import enumerate_states, generator, stationary
from crnl.networks.am_reversible import am_reversible, delta_star, reverse_pairing
from crnl.vectorized import compile_network
from experiments.slaving_axis import slaved


def cell(gamma, omega, eps, theta):
    net = am_reversible(gamma)
    pairing = reverse_pairing(net)
    ds = delta_star(gamma)
    x0 = eps * ds
    st = slaved(net, x0)
    if st is None:
        return None
    nb = int(round(st[2] * omega))
    rest = omega - nb
    d0 = max(1, int(round(x0 * omega)))
    if (rest - d0) % 2:
        d0 -= 1
    thr = max(2, int(round(theta * ds * omega)))
    if thr <= d0 or rest - d0 < 0:
        return None
    n0 = np.array([(rest + d0) // 2, (rest - d0) // 2, nb], dtype=np.int64)

    states, index = enumerate_states(3, int(omega))
    dvec = np.array([int(s[0]) - int(s[1]) for s in states])
    absorb = np.abs(dvec) >= thr
    correct = absorb & (dvec >= thr)          # started with X ahead
    error = absorb & (dvec <= -thr)

    comp = compile_network(net, float(omega))
    S = net.stoichiometry_matrix().astype(np.int64)
    A = propensities_batch(comp, states.astype(float))
    a_tot = A.sum(axis=1)

    # ---- splitting probabilities, named at the boundary (§35, no 1-split) -------
    Q = generator(net, int(omega), float(omega))
    tr = np.where(~absorb)[0]
    tmap = -np.ones(len(states), dtype=np.int64)
    tmap[tr] = np.arange(len(tr))
    Qtt = Q[tr][:, tr].tocsr()
    Qta = Q[tr][:, np.where(absorb)[0]].tocsr()
    si = int(tmap[index[tuple(n0)]])
    p_e = float(spla.spsolve(Qtt, -(Qta @ error[absorb].astype(float)))[si])
    p_c = float(spla.spsolve(Qtt, -(Qta @ correct[absorb].astype(float)))[si])

    # ---- tilted generator from REVERSE propensities (§41) ------------------------
    rows, cols, vals = [], [], []
    for j in range(net.n_reactions):
        rev = int(pairing[j])
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
    M = sp.coo_matrix((-v[inner], (tmap[r[inner]], tmap[c[inner]])),
                      shape=(len(tr), len(tr))).tocsr() + sp.diags(a_tot[tr])
    M = M.tocsc()

    pi = stationary(net, int(omega), float(omega))
    ratio = pi / max(pi[index[tuple(n0)]], 1e-300)     # §41's convention: pi(n)/pi(n0)

    out = {}
    edge = ~inner
    for name, mask in (("all", absorb), ("c", correct), ("e", error)):
        bnd = ratio * mask.astype(float)
        b = np.zeros(len(tr))
        np.add.at(b, tmap[r[edge]], v[edge] * bnd[c[edge]])
        out[name] = float(spla.spsolve(M, b)[si])

    return {"gamma": gamma, "omega": omega, "eps": eps, "theta": theta,
            "p_c": p_c, "p_e": p_e, "Phi_all": out["all"],
            "Phi_c": out["c"], "Phi_e": out["e"], "thr": int(thr)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.10, 0.20, 0.30])
    ap.add_argument("--omegas", type=int, nargs="+", default=[40, 60, 90])
    ap.add_argument("--epss", type=float, nargs="+", default=[0.35, 0.50])
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/outcome_split.json"))
    args = ap.parse_args()

    t0 = time.time()
    print("Phi_o = p_o * <e^(-S_tot) | o>;  §41 says Phi_c + Phi_e = 1")
    print(f"{'gamma':>7}{'Om':>5}{'eps':>6}{'p_e':>12}{'Phi_c':>12}"
          f"{'Phi_e':>10}{'sum':>10}{'Phi_e/p_c':>11}{'Phi_c/p_e':>11}")
    rows = []
    for g in args.gammas:
        for om in args.omegas:
            for eps in args.epss:
                try:
                    r = cell(g, om, eps, args.theta)
                except Exception as e:
                    print(f"{g:>7.2f}{om:>5}{eps:>6.2f}   SKIPPED ({type(e).__name__})")
                    continue
                if r is None:
                    continue
                s = r["Phi_c"] + r["Phi_e"]
                r["sum"] = s
                r["ratio_e"] = r["Phi_e"] / r["p_c"] if r["p_c"] else float("nan")
                r["ratio_c"] = r["Phi_c"] / r["p_e"] if r["p_e"] else float("nan")
                rows.append(r)
                print(f"{g:>7.2f}{om:>5}{eps:>6.2f}{r['p_e']:>12.4e}"
                      f"{r['Phi_c']:>12.4e}{r['Phi_e']:>10.6f}{s:>10.6f}"
                      f"{r['ratio_e']:>11.6f}{r['ratio_c']:>11.4e}")

    if not rows:
        print("no cells")
        return

    print(f"\n=== P1 GATE: does Phi_c + Phi_e reconstruct §41's 1?")
    d = np.array([abs(r["sum"] - 1.0) for r in rows])
    print(f"  |sum - 1| over {len(rows)} cells: max {d.max():.3e}, median {np.median(d):.3e}")
    print(f"  -> P1 {'HOLDS' if d.max() < 1e-6 else 'FAILS'}")

    print(f"\n=== P2: is Phi_e = p_c and Phi_c = p_e?")
    re = np.array([r["ratio_e"] for r in rows])
    rc = np.array([r["ratio_c"] for r in rows])
    print(f"  Phi_e/p_c: {re.min():.6f} .. {re.max():.6f}   (1 would be exact)")
    print(f"  Phi_c/p_e: {rc.min():.4e} .. {rc.max():.4e}")
    p2 = abs(re - 1).max() < 1e-6
    print(f"  -> P2 {'HOLDS -- reliability IS dissipation, exactly' if p2 else 'FAILS'}")

    print(f"\n=== P4: if P2 fails, are the ratios at least CONSTANT?")
    print(f"  Phi_e/p_c spread: {100*(re.max()-re.min())/abs(re.mean()):.2f}%")
    print(f"  Phi_c/p_e spread: {100*(rc.max()-rc.min())/abs(rc.mean()):.2f}%")

    print(f"\n=== P5 (rule 9): does the relation behave the same along EVERY axis?")
    for key in ("gamma", "omega", "eps"):
        vals = sorted({r[key] for r in rows})
        meds = [np.median([r["ratio_e"] for r in rows if r[key] == v]) for v in vals]
        print(f"  Phi_e/p_c by {key:>6}: "
              + "  ".join(f"{v}:{m:.6f}" for v, m in zip(vals, meds)))

    print(f"\n=== P6: Jensen sanity -- <e^-S|e> = Phi_e/p_e must exceed 1 for error paths")
    ee = np.array([r["Phi_e"] / r["p_e"] for r in rows])
    print(f"  <e^(-S)|error> spans {ee.min():.3e} .. {ee.max():.3e}"
          f"  -> {'all > 1, error paths consume entropy as they must' if ee.min() > 1 else 'SOME <= 1, suspect a sign error in sigma_local'}")
    oc = np.array([r["p_c"] / r["p_e"] for r in rows])
    print(f"  odds p_c/p_e spans {oc.min():.3e} .. {oc.max():.3e}")
    print(f"  ratio <e^-S|e> / (p_c/p_e): {(ee/oc).min():.6f} .. {(ee/oc).max():.6f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
