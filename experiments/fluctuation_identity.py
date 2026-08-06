"""T-TUR-b: an EXACT identity for the entropy machinery — the integral fluctuation theorem

§40 tested an inequality. THEORIES §5 records Neri's result that first-passage bounds
become EQUALITIES when the thresholded current is the entropy production, via the
martingale property of exp(-S). That suggested testing p_- = exp(-l_-), which needs an
augmented chain. **There is a cleaner form of the same physics that needs no new state
space.**

The integral fluctuation theorem, `<exp(-S_tot)> = 1`, holds at ANY stopping time --
including the absorption time this project already uses. And it is a LINEAR SOLVE, because
tilting each transition by its own entropy weight collapses:

    a_j(n) * exp(-Delta s_j) = a_j(n) * a_rev(n') / a_j(n) = a_rev(n')

**so the tilted generator is built from the REVERSE propensities.** With
psi(n) = E[exp(-S_med) | start n],

    a_tot(n) psi(n) = sum_j a_rev(n_j') psi(n_j'),     psi = boundary value on absorbing

which is the same sparsity and cost as every other solve here.

WHY THIS IS WORTH MORE THAN ANOTHER MEASUREMENT. §37-§40 all rest on one object:
`sigma_local(n) = sum_j a_j ln[a_j/a_rev]`, the reverse pairing that defines it, and the
sign convention. Nothing has ever checked that object against anything external -- §40
checked an inequality, which a wrong-by-a-constant entropy would still satisfy. **The IFT
is an equality that a wrong convention cannot satisfy**, and it holds for every gamma,
Omega, eps and threshold, so it is a far broader check than any single cell.

TOTAL vs MEDIUM. The IFT holds for the TOTAL entropy, medium plus system. The system term
is `ln[pi(n_0)/pi(n_T)]` with pi the stationary distribution of the UNABSORBED chain, which
`cme.stationary` computes. So two solves:

    psi_med : boundary 1            -> E[exp(-S_med)]
    psi_tot : boundary pi(n_0)/pi(n) -> E[exp(-S_med) * pi(n_0)/pi(n_T)] = E[exp(-S_tot)]

PREDICTIONS, written before running:

  P1  `E[exp(-S_med)] != 1` in general. The medium term alone does not satisfy the IFT,
      and quoting it as though it did would be the error this experiment exists to rule
      out.
  P2  THE TEST. `E[exp(-S_tot)] = 1` to solver precision, at every gamma, Omega and eps.
      **This validates sigma_local, the reverse pairing and the sign convention together
      against an exact external identity** -- the first equality, as opposed to
      inequality, that §37-§40's machinery has faced.
  P3  If P2 FAILS, the discrepancy is diagnostic rather than merely bad: a constant offset
      means a sign or factor convention, a gamma-dependent one means the pairing is wrong
      for some reaction, and an Omega-dependent one means a discretisation error in
      sigma_local. **§37-§40 all rest on this object, so a failure here would propagate to
      every cost result in the project.**
  P4  The identity is exact at ANY stopping time, so it must hold at every threshold
      theta as well. A theta-dependence would mean the absorbing set is interacting with
      the entropy accounting, which nothing in the construction allows.
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
from crnl.cme import enumerate_states, stationary
from crnl.networks.am_reversible import am_reversible, delta_star, reverse_pairing
from crnl.vectorized import compile_network
from experiments.slaving_axis import slaved


def identity_cell(gamma, omega, eps, theta):
    net = am_reversible(gamma)
    pairing = reverse_pairing(net)
    ds = delta_star(gamma)
    x0 = eps * ds
    st = slaved(net, x0)
    nb = int(round(st[2] * omega))
    rest = omega - nb
    d0 = max(1, int(round(x0 * omega)))
    if (rest - d0) % 2:
        d0 -= 1
    n0 = np.array([(rest + d0) // 2, (rest - d0) // 2, nb], dtype=np.int64)
    thr = max(2, int(round(theta * ds * omega)))

    states, index = enumerate_states(3, int(omega))
    absorb = np.array([abs(int(s[0]) - int(s[1])) >= thr for s in states])
    comp = compile_network(net, float(omega))
    S = net.stoichiometry_matrix().astype(np.int64)
    A = propensities_batch(comp, states.astype(float))
    a_tot = A.sum(axis=1)

    # tilted generator: weight of jump j out of n is a_rev(n + S_j), not a_j(n)
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
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    vals = np.concatenate(vals)

    tr = np.where(~absorb)[0]
    tmap = -np.ones(len(states), dtype=np.int64)
    tmap[tr] = np.arange(len(tr))
    keep = ~absorb[rows]
    r, c, v = rows[keep], cols[keep], vals[keep]

    pi = stationary(net, int(omega), float(omega))
    si = index[tuple(n0)]

    out = {}
    # e^{-S_tot} = e^{-S_med} * pi(n_T)/pi(n_0), so the absorbing boundary carries
    # pi(n)/pi(n_0) -- NOT its reciprocal, which is the convention error this
    # experiment exists to catch and which the first pass duly contained.
    for name, bnd in (("med", np.ones(len(states))),
                      ("tot", pi / max(pi[si], 1e-300))):
        inner = ~absorb[c]
        M = sp.coo_matrix((-v[inner], (tmap[r[inner]], tmap[c[inner]])),
                          shape=(len(tr), len(tr))).tocsr()
        M = M + sp.diags(a_tot[tr])
        b = np.zeros(len(tr))
        edge = ~inner
        np.add.at(b, tmap[r[edge]], v[edge] * bnd[c[edge]])
        psi = spla.spsolve(M.tocsc(), b)
        res = float(np.linalg.norm(M @ psi - b) / max(np.linalg.norm(b), 1e-300))
        out[name] = {"value": float(psi[tmap[si]]), "residual": res}
    return {"gamma": gamma, "omega": omega, "eps": eps, "theta": theta,
            "thr": int(thr), "E_exp_neg_S_med": out["med"]["value"],
            "E_exp_neg_S_tot": out["tot"]["value"],
            "residual": max(out["med"]["residual"], out["tot"]["residual"])}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.10, 0.20, 0.30])
    ap.add_argument("--omegas", type=int, nargs="+", default=[60, 100, 150])
    ap.add_argument("--epss", type=float, nargs="+", default=[0.35, 0.50])
    ap.add_argument("--thetas", type=float, nargs="+", default=[0.70, 0.80])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/fluctuation_identity.json"))
    args = ap.parse_args()

    t0 = time.time()
    print("integral fluctuation theorem: E[exp(-S_tot)] must equal 1 at ANY stopping time")
    print(f"{'gamma':>7}{'Omega':>7}{'eps':>6}{'theta':>7}{'E[e^-S_med]':>14}"
          f"{'E[e^-S_tot]':>14}{'|dev from 1|':>14}{'resid':>10}")
    rows = []
    for g in args.gammas:
        for om in args.omegas:
            for eps in args.epss:
                for th in args.thetas:
                    try:
                        r = identity_cell(g, om, eps, th)
                    except Exception as e:
                        print(f"{g:>7.2f}{om:>7}{eps:>6.2f}{th:>7.2f}   "
                              f"SKIPPED ({type(e).__name__})")
                        continue
                    rows.append(r)
                    dev = abs(r["E_exp_neg_S_tot"] - 1.0)
                    print(f"{g:>7.2f}{om:>7}{eps:>6.2f}{th:>7.2f}"
                          f"{r['E_exp_neg_S_med']:>14.6e}"
                          f"{r['E_exp_neg_S_tot']:>14.8f}{dev:>14.2e}"
                          f"{r['residual']:>10.1e}")

    med = np.array([r["E_exp_neg_S_med"] for r in rows])
    tot = np.array([r["E_exp_neg_S_tot"] for r in rows])
    print(f"\n=== P1: does the MEDIUM term alone satisfy the IFT?")
    print(f"  E[exp(-S_med)] spans {med.min():.3e} .. {med.max():.3e}"
          f"   -> {'NO, as predicted' if np.abs(med-1).max() > 1e-3 else 'unexpectedly yes'}")

    print(f"\n=== P2: does the TOTAL entropy satisfy it?")
    dev = np.abs(tot - 1.0)
    print(f"  |E[exp(-S_tot)] - 1| over {len(tot)} cells: max {dev.max():.3e},"
          f" median {np.median(dev):.3e}")
    print(f"  -> P2 {'HOLDS -- sigma_local, the pairing and the sign convention are '
                     'validated against an exact identity' if dev.max() < 1e-6 else 'FAILS'}")

    if dev.max() >= 1e-6:
        print(f"\n=== P3 diagnosis: what does the deviation depend on?")
        for key in ("gamma", "omega", "eps", "theta"):
            vals = sorted({r[key] for r in rows})
            meds = [np.median([abs(r["E_exp_neg_S_tot"] - 1)
                               for r in rows if r[key] == v]) for v in vals]
            print(f"  by {key:>6}: " + "  ".join(f"{v}:{m:.2e}"
                                                 for v, m in zip(vals, meds)))

    print(f"\n=== P4: is the identity threshold-independent?")
    for th in args.thetas:
        d = [abs(r["E_exp_neg_S_tot"] - 1) for r in rows if r["theta"] == th]
        if d:
            print(f"  theta={th}: median |dev| = {np.median(d):.3e}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
