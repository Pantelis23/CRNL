"""§24's kill test, on §21's own ladder: does the cliff track the COORDINATE or the NOISE?

§21 measured that every level keeping ANY noise recovers the restoration error
exponent to 2-12% while the ODE is categorically wrong, and concluded "having noise
at all" is what cannot be discarded. §24 argued that conclusion is under-determined:
every level on that ladder keeps every species, so it retains the noise AND the
coordinates together. Separating them on a stage-level kernel in the FUEL network
gave the coordinate axis +0.1000 of exponent against the noise axis's +0.0146, and
the refined claim -- **noise matters in the coordinate that carries the signal and
costs nothing in a bookkeeping coordinate** -- was recorded as a suspect on one
system, one channel, one Omega.

This runs that suspect's kill test in §21's own system, with §21's own observable
(P(error) from a start biased by eps, decided at |n_X - n_Y| >= theta*delta**Omega)
and §21's exact CME reference.

THE CONSTRUCTION. `am_reversible` conserves n_X + n_Y + n_B, so the CLE's noise
increment in species space sums to zero: `nz = (a, b, c)` with `a + b + c = 0`. Two
free coordinates, and they split cleanly by role:

    delta = n_X - n_Y     THE SIGNAL. The bit is its sign; the observable is
                          entirely a question about it crossing zero.
    s     = n_X + n_Y     BOOKKEEPING. It is `Omega - n_B`, the blank pool -- it
                          sets how much material is available to convert but no
                          decision happens in it.

Project the noise, keep the drift full in every arm:

    full        nz unchanged                                    = §21's CLE
    delta-only  a' = (a-b)/2, b' = -(a-b)/2, c' = 0     signal keeps its noise,
                                                        BLANK POOL goes deterministic
    s-only      a' = b' = (a+b)/2, c' = c               blank pool keeps its noise,
                                                        SIGNAL goes deterministic

Both projections conserve exactly. Note `s`-noise IS `n_B`-noise: with a+b+c = 0,
zeroing `c` forces `a+b = 0`, so "the blank pool evolves deterministically" and
"there is no noise along `s`" are the same statement.

WHICH ARM IS INFORMATIVE, stated before running so it cannot be claimed afterwards.
**The `delta-only` arm is close to definitional and proves little on its own**: the
observable is a question about `delta` crossing zero, so removing `delta`'s noise
should obviously kill the error probability, and the drift from a biased start runs
straight to the correct attractor. It is included as the anchor at one end, not as
a finding. **The `s-only` arm is the real test.** Nothing guarantees that blank-pool
noise is discardable -- every propensity in AM depends on `n_B`, so noise in `s`
feeds straight into `delta`'s dynamics, and it could easily matter as much as
`delta`'s own.

THE MAGNITUDE CONFOUND, and the control for it. Both projections REMOVE noise, so an
arm that fails might be failing on total amplitude rather than placement. The
realised removed variance is therefore reported per arm. The argument only works if
the two removals are comparable in size: if `delta-only` (less total noise than full)
still recovers the exponent while `s-only` (a comparable reduction) does not, then
placement is doing the work and amplitude is not.

PREDICTIONS, written before running:

  P1  THE ONE THAT MATTERS. `s-only` -- signal deterministic, blank pool noisy --
      fails like the ODE: p_error collapses toward 0 at every cell where the CME
      gives 1.5e-3 to 1.6e-1. If instead it recovers the exponent to within a few
      percent, then noise ANYWHERE suffices, §24's coordinate reading is wrong, and
      §21's "having noise at all" stands unqualified.
  P2  `delta-only` recovers the exponent to roughly the full CLE's accuracy despite
      carrying strictly less noise. This is the placement-over-amplitude half.
  P3  The full CLE reproduces §21's own numbers, which is the harness control. If it
      does not, nothing else here is readable.
  P4  THE MESSY OUTCOME I EXPECT IS POSSIBLE and am naming so it is not spun later:
      `s-only` fails PARTIALLY -- better than the ODE's categorical zero, worse than
      the CLE. That would mean the signal/bookkeeping distinction is real but GRADED
      rather than a cliff, and §24's "7x" would have to be restated as a continuum
      rather than a kind-difference. I would rather write that down now than discover
      the temptation afterwards.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from crnl.approximations import propensities_batch
from crnl.networks.am_reversible import am_reversible, delta_star
from crnl.vectorized import compile_network
from experiments.approximation_hierarchy import _setup, p_cme

MODES = ("full", "delta-only", "s-only", "uniform-11pct")

# the smoke run measured delta-only at 0.111 of the full noise variance; this arm
# keeps the SAME fraction spread over BOTH coordinates, which separates amplitude
# from placement outright: if 11% everywhere works while 88% in the wrong place
# fails, the total is not what the observable is sensitive to.
UNIFORM_FRAC = 0.111


def project(nz: np.ndarray, mode: str) -> np.ndarray:
    """Project the species-space noise increment onto one coordinate's subspace."""
    if mode == "full":
        return nz
    if mode == "uniform-11pct":
        return nz * np.sqrt(UNIFORM_FRAC)
    a, b = nz[:, 0], nz[:, 1]
    out = np.zeros_like(nz)
    if mode == "delta-only":
        h = 0.5 * (a - b)
        out[:, 0], out[:, 1], out[:, 2] = h, -h, 0.0
    elif mode == "s-only":
        h = 0.5 * (a + b)
        out[:, 0], out[:, 1], out[:, 2] = h, h, nz[:, 2]
    else:
        raise ValueError(mode)
    return out


def run_projected(comp, n0, rng, *, dt, thr, trials, t_max, mode,
                  max_steps=200_000) -> dict:
    """Batch CLE with the noise projected. Drift is untouched in every mode."""
    S = comp.S.astype(float)
    n = np.tile(np.asarray(n0, dtype=float), (trials, 1))
    live = np.ones(trials, bool)
    t = np.zeros(trials)
    rejected = 0
    var_full = var_kept = 0.0
    for _ in range(max_steps):
        if not live.any():
            break
        idx = np.where(live)[0]
        a = propensities_batch(comp, n[idx])
        mean = a * dt
        xi = np.sqrt(np.clip(mean, 0.0, None)) * rng.standard_normal(mean.shape)
        nz_full = xi @ S.T
        nz = project(nz_full, mode)
        var_full += float((nz_full ** 2).sum())
        var_kept += float((nz ** 2).sum())
        cand = n[idx] + (mean @ S.T) + nz
        ok = (cand >= 0.0).all(axis=1)
        rejected += int((~ok).sum())
        upd = idx[ok]
        n[upd] = cand[ok]
        t[upd] += dt
        done = (np.abs(n[idx, 0] - n[idx, 1]) >= thr) | (t[idx] >= t_max)
        live[idx[done]] = False
    fin = ~live
    wrong = int((n[fin, 0] <= n[fin, 1]).sum())
    nok = int(fin.sum())
    return {"p": wrong / nok if nok else float("nan"), "n_ok": nok,
            "unfinished": int(live.sum()), "rejected": rejected,
            "variance_kept_frac": var_kept / var_full if var_full > 0 else float("nan")}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gamma", type=float, default=0.30)
    ap.add_argument("--omegas", type=int, nargs="+", default=[40, 60, 80, 100])
    ap.add_argument("--eps-fracs", type=float, nargs="+", default=[0.20, 0.35])
    ap.add_argument("--theta", type=float, default=0.80)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--trials", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=20260730)
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/noise_placement.json"))
    args = ap.parse_args()

    t0 = time.time()
    print(f"gamma={args.gamma} theta={args.theta} dt={args.dt} "
          f"trials={args.trials}   (§21's cells and observable)")
    print(f"\n{'eps':>5} {'Omega':>6} {'CME (exact)':>12} {'CLE full':>11} "
          f"{'delta-only':>11} {'s-only':>11} {'unif-11%':>11} {'var kept d/s':>13}")
    rows = []
    for eps in args.eps_fracs:
        for om in args.omegas:
            n0, thr, _ = _setup(args.gamma, om, eps, args.theta)
            comp = compile_network(am_reversible(args.gamma), float(om))
            exact = p_cme(args.gamma, om, eps, args.theta)
            got = {}
            for mode in MODES:
                rng = np.random.default_rng(args.seed + 17 * om + int(1000 * eps))
                got[mode] = run_projected(comp, n0, rng, dt=args.dt, thr=thr,
                                          trials=args.trials, t_max=4000.0,
                                          mode=mode)
            print(f"{eps:>5.2f} {om:>6} {exact:>12.6f} "
                  f"{got['full']['p']:>11.6f} {got['delta-only']['p']:>11.6f} "
                  f"{got['s-only']['p']:>11.6f} {got['uniform-11pct']['p']:>11.6f} "
                  f"{got['delta-only']['variance_kept_frac']:>6.3f}/"
                  f"{got['s-only']['variance_kept_frac']:<6.3f}")
            rows.append({"eps_frac": eps, "omega": om, "p_cme": exact,
                         **{f"p_{m}": got[m]["p"] for m in MODES},
                         **{f"var_{m}": got[m]["variance_kept_frac"] for m in MODES},
                         **{f"unfinished_{m}": got[m]["unfinished"] for m in MODES}})

    print(f"\n=== exponent in Omega of P(error), per arm and per eps "
          f"(§21's comparison)")
    summary = {}
    for eps in args.eps_fracs:
        rs = [r for r in rows if abs(r["eps_frac"] - eps) < 1e-12]
        om = np.array([r["omega"] for r in rs], float)
        line = {}
        for key in ["p_cme"] + [f"p_{m}" for m in MODES]:
            y = np.array([r[key] for r in rs])
            good = np.isfinite(y) & (y > 0)
            if good.sum() < 2:
                line[key] = None
                print(f"  eps={eps}: {key:>14} -> no positive values "
                      f"(categorical failure, as the ODE gives)")
                continue
            p = float(np.polyfit(np.log(om[good]), np.log(y[good]), 1)[0])
            line[key] = p
            print(f"  eps={eps}: {key:>14} -> exponent {p:>8.4f}"
                  f"   ({good.sum()}/{len(rs)} cells positive)")
        summary[str(eps)] = line

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"rows": rows, "exponents": summary},
                                   indent=2, default=float))
    print(f"\n({time.time()-t0:.0f}s) wrote {args.out}")


if __name__ == "__main__":
    main()
