"""T-COST-k: is the whole absorption correction 1-D discreteness? — one tridiagonal solve

§49 derived the absorption correction as O(1/Omega) and confirmed the EXPONENT (constant to
6% and 17% where the expansion applies), but its coefficient K -- the bulk diffusion term --
supplies only ~31% of the measured value and has the WRONG SIGN at gamma = 0.07, because K
is a near-cancellation: g rises then falls across the traversal.

Two same-order terms were named as missing: the absorbing boundary layer, and the
Kramers-Moyal truncation error (the CME is a jump process; its diffusion approximation is
itself wrong at O(1/Omega)). **Computing them separately is unnecessary. The 1-D slaved
BIRTH-DEATH CHAIN contains all three at once and needs no CME.**

On the slaved manifold n_delta hops +-1 with rates Omega*up(delta) and Omega*dn(delta),
both already returned by `updown`. That chain's exact MFPT to |n_delta| >= thr is a
TRIDIAGONAL solve. Comparing it against the same T_det = int ddelta/mu gives

    bd_coeff = (T_det/T_bd - 1) * Omega

which is the complete 1-D stochastic correction -- bulk, boundary layer and discreteness
together, with no expansion and no fitted parameter. The question is whether it equals the
correction actually measured against the 2-D CME,

    cme_coeff = (gap - <eps>_time) * Omega          [§48's numbers]

**If they agree, the entire absorption correction is one-dimensional discreteness and the
cost budget closes with no free parameter anywhere.** If the 1-D chain undershoots, the
remainder is genuine 2-D noise -- pool fluctuation feeding into delta -- which no 1-D
account can supply and which nothing in this project has yet measured.

PREDICTIONS, written before running:

  P1  GATE. bd_coeff must be Omega-INDEPENDENT, since §49 derived the exponent as 1. Checked
      over Omega = 300..2000, which the tridiagonal solve makes free. If bd_coeff drifts,
      the exponent is not 1 on the cheap instrument either and §49's confirmation was luck.
  P2  THE TEST, absolute. **bd_coeff = cme_coeff** in the two cells §49 found asymptotic
      (gamma=0.20 rho=1 at 4.21, gamma=0.20 rho=32 at 2.08).
  P3  It must beat K. K gave 0.313 of the measured value; the birth-death chain includes
      the two terms K omits, so **bd_coeff/cme_coeff must be much nearer 1 than 0.313**.
      Anything still near 0.31 would mean the extra terms are negligible and the deficit is
      2-D after all.
  P4  THE SIGN REPAIR. At gamma = 0.07 the bulk term K = -0.263 against a measured +3.65.
      **bd_coeff must be POSITIVE there.** A near-cancelling bulk integral going the wrong
      way is exactly what a complete calculation should fix, and if bd_coeff is also
      negative the problem is not the expansion but the slaved reduction itself.
  P5  THE PRE-ASYMPTOTIC CELLS. gamma=0.35 and rho=0.5 had cme_coeff still rising at
      Omega = 1000 (to 10.26 and 9.28). Their bd_coeff, being Omega-independent by P1, gives
      the asymptote those measurements are climbing toward -- so **bd_coeff should EXCEED
      the Omega=1000 cme_coeff in exactly those two cells** and match in the others. That is
      a per-cell prediction with a direction, not an average.
  P6  REFUTING OUTCOME. If bd_coeff undershoots cme_coeff systematically in the asymptotic
      cells, a genuine 2-D contribution exists. That is a real and interesting outcome --
      it would mean pool noise feeds the signal in a way the slaved reduction cannot carry,
      which is the stochastic counterpart of §47's deterministic lag and is not in evidence
      anywhere yet.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np
from scipy.linalg import solve_banded

from crnl.networks.am_reversible import reverse_pairing
from experiments.arrhenius_optimum import am_rho, delta_star_rho
from experiments.lag_endpoints import traversal
from experiments.slaving_axis import updown


def rates(net, delta):
    """(up, dn) at |delta|, using the exchange symmetry up(-d) = dn(d) for d < 0."""
    ud = updown(net, abs(delta))
    if ud is None:
        return None
    return (ud[0], ud[1]) if delta >= 0 else (ud[1], ud[0])


def bd_mfpt(net, omega, m0, thr, cache=None):
    """Exact MFPT of the 1-D slaved birth-death chain to |m| >= thr, tridiagonally.

    States m = -(thr-1) .. (thr-1). Birth rate Omega*up(m/Omega), death Omega*dn(m/Omega).
    """
    ms = np.arange(-(thr - 1), thr)
    n = len(ms)
    up = np.empty(n)
    dn = np.empty(n)
    for i, m in enumerate(ms):
        r = rates(net, m / omega) if cache is None else cache(m / omega)
        if r is None:
            return None
        up[i], dn[i] = omega * r[0], omega * r[1]
    ab = np.zeros((3, n))
    ab[1, :] = -(up + dn)          # diagonal
    ab[0, 1:] = up[:-1]            # super: coupling m -> m+1
    ab[2, :-1] = dn[1:]            # sub:   coupling m -> m-1
    t = solve_banded((1, 1), ab, -np.ones(n))
    i0 = int(np.where(ms == m0)[0][0])
    return float(t[i0])


def cell(gamma, rho, omega, eps_real, theta_real):
    ds = delta_star_rho(gamma, rho)
    net = am_rho(gamma, rho)
    pairing = reverse_pairing(net)
    lo, hi = eps_real * ds, theta_real * ds
    m0 = int(round(lo * omega))
    thr = int(round(hi * omega))
    if thr <= m0 + 1:
        return None
    # §48's lesson, applied here: T_det must integrate between the SAME endpoints the
    # chain runs between, not the unrounded ones. The rounding is O(1/Omega) in delta,
    # which after multiplying by Omega is exactly the size of the effect being measured --
    # and it is what made the first pass of this experiment scatter by up to 3029%.
    t_det = traversal(net, m0 / omega, thr / omega, pairing, n=4001)
    if t_det is None:
        return None
    memo = {}

    def cached(d):
        k = round(d, 12)
        if k not in memo:
            memo[k] = rates(net, d)
        return memo[k]

    t_bd = bd_mfpt(net, omega, m0, thr, cache=cached)
    if t_bd is None or not np.isfinite(t_bd) or t_bd <= 0:
        return None
    return {"gamma": gamma, "rho": rho, "omega": omega, "t_det": t_det, "t_bd": t_bd,
            "bd_coeff": float((t_det / t_bd - 1.0) * omega)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--omegas", type=int, nargs="+", default=[300, 500, 700, 1000, 2000])
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/birthdeath_absorption.json"))
    args = ap.parse_args()

    t0 = time.time()
    meas = json.load(open("results/lag_endpoints.json"))
    kk = {(c["gamma"], c["rho"]): c["K"]
          for c in json.load(open("results/absorption_coefficient.json"))}
    cells = sorted({(r["gamma"], r["rho"]) for r in meas})

    print("=== P1 GATE: is bd_coeff Omega-independent (exponent 1 on the cheap instrument)?")
    print(f"{'gamma':>6}{'rho':>6}" + "".join(f"{f'Om={o}':>10}" for o in args.omegas)
          + f"{'spread':>9}")
    rows, bd = [], {}
    for g, r in cells:
        ref = [c for c in meas if c["gamma"] == g and c["rho"] == r][0]
        line, vals = f"{g:>6.2f}{r:>6.1f}", []
        for om in args.omegas:
            c = cell(g, r, om, ref["eps_real"], ref["theta_real"])
            if c is None:
                line += f"{'--':>10}"
                continue
            rows.append(c)
            vals.append(c["bd_coeff"])
            line += f"{c['bd_coeff']:>10.3f}"
        if vals:
            v = np.array(vals)
            bd[(g, r)] = float(v[-1])
            line += f"{100*(v.max()-v.min())/abs(v.mean()):>8.1f}%"
        print(line)
    sp = [100 * (np.ptp([c["bd_coeff"] for c in rows if c["gamma"] == g and c["rho"] == r])
                 / abs(np.mean([c["bd_coeff"] for c in rows
                                if c["gamma"] == g and c["rho"] == r]))) for g, r in bd]
    print(f"  worst spread {max(sp):.1f}%  -> P1 {'HOLDS' if max(sp) < 10 else 'FAILS'}")

    print(f"\n=== P2/P3/P4/P5: the 1-D chain against the 2-D CME measurement")
    print(f"{'gamma':>6}{'rho':>6}{'bd_coeff':>10}{'cme@1000':>10}{'bd/cme':>9}"
          f"{'K/cme':>8}{'status':>13}")
    out = []
    for g, r in cells:
        if (g, r) not in bd:
            continue
        by = {c["omega"]: c for c in meas if c["gamma"] == g and c["rho"] == r}
        oms = sorted(by)
        cme = [(by[o]["gap_real"] - by[o]["pred_real"]) * o for o in oms]
        rising = cme[-1] > 1.25 * cme[0]
        ratio = bd[(g, r)] / cme[-1]
        kr = kk.get((g, r), float("nan")) / cme[-1]
        out.append({"gamma": g, "rho": r, "bd": bd[(g, r)], "cme": cme[-1],
                    "bd_over_cme": ratio, "K_over_cme": kr, "rising": bool(rising)})
        print(f"{g:>6.2f}{r:>6.1f}{bd[(g,r)]:>10.3f}{cme[-1]:>10.3f}{ratio:>9.3f}"
              f"{kr:>8.3f}{'PRE-ASYMPT' if rising else 'asymptotic':>13}")

    asym = [c for c in out if not c["rising"]]
    pre = [c for c in out if c["rising"]]

    print(f"\n=== P2/P3: the asymptotic cells -- does the 1-D chain supply the whole thing?")
    for c in asym:
        print(f"  gamma={c['gamma']:.2f} rho={c['rho']:<5}"
              f" bd/cme = {c['bd_over_cme']:.3f}   (bulk term K alone: {c['K_over_cme']:.3f})")
    if asym:
        v = np.array([c["bd_over_cme"] for c in asym])
        kv = np.array([c["K_over_cme"] for c in asym])
        print(f"  bd/cme {v.min():.3f}..{v.max():.3f} mean {v.mean():.3f}"
              f"   vs K/cme mean {kv.mean():.3f}")
        print(f"  -> P3 {'HOLDS, the 1-D chain beats the bulk term' if abs(v.mean()-1) < abs(kv.mean()-1) else 'FAILS, no better than K'}")
        print(f"  -> P2 {'HOLDS: the absorption correction is entirely 1-D discreteness' if abs(v.mean()-1) < 0.15 else 'FAILS: a residual remains'}")

    print(f"\n=== P4: the sign repair at gamma = 0.07")
    for c in out:
        if c["gamma"] == 0.07:
            print(f"  bulk K = {kk.get((c['gamma'],c['rho'])):.3f} (wrong sign);"
                  f" bd_coeff = {c['bd']:.3f}; measured {c['cme']:.3f}"
                  f"  -> {'REPAIRED' if c['bd'] > 0 else 'still wrong sign'}")

    print(f"\n=== P5: pre-asymptotic cells -- bd_coeff should EXCEED the rising measurement")
    for c in pre:
        print(f"  gamma={c['gamma']:.2f} rho={c['rho']:<5} bd = {c['bd']:.3f}"
              f"  vs cme@1000 = {c['cme']:.3f} (still rising)"
              f"  -> {'exceeds, AS PREDICTED' if c['bd'] > c['cme'] else 'BELOW it, contrary to P5'}")

    print(f"\n=== P6: verdict")
    if asym:
        v = np.array([c["bd_over_cme"] for c in asym])
        if abs(v.mean() - 1) < 0.15:
            print("  The cost budget closes with NO free parameter:")
            print("    T_det/MFPT - 1 = <eps>_time + bd_coeff/Omega")
            print("  the first term §47's deterministic lag, the second the exact 1-D")
            print("  slaved chain's own discreteness. Nothing is fitted anywhere.")
        else:
            print(f"  bd/cme = {v.mean():.3f}: a {100*(1-v.mean()):+.0f}% residual survives a")
            print("  COMPLETE 1-D account, so it is genuinely 2-D -- pool noise feeding the")
            print("  signal, the stochastic counterpart of §47's deterministic lag.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"sweep": rows, "compare": out}, indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
