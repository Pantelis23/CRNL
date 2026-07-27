"""What a tilted landscape buys and what it costs -- FINDINGS 16.

Every network before `am_asymmetric` is symmetric under relabelling the symbols,
so both attractors are mirror images and one coefficient describes both. This
tilts the two autocatalytic branches by `beta`, keeping every reverse at gamma x
its OWN forward so the cycle affinity stays -3 ln gamma and beta costs no
thermodynamic force (see the module docstring; the dissipation RATE is a separate
question and is not claimed here).

Three measurements, `--part` selects:

  fold   the saddle-node beta_c(gamma) past which the network answers X whatever
         it is shown, and where the bias actually lives (the saddle, not the
         attractors)
  wells  both escape barriers from the exact quasipotential. Tests the prediction
         written before running -- "the barrier sum FALLS with tilt, so symmetric
         AM maximises restoration capacity" -- which is WRONG: the sum rises.
  info   mutual information through the tilted restorer for a symmetric source,
         which is the figure of merit the barrier sum fails to be. Falls
         monotonically in |beta|, and the penalty GROWS with Omega.

A NOTE ON THE INITIAL STATES, because the first version of `info` produced a
clean, plausible, impossible number. Building the two inputs as (rest +- d0)//2
gives biases of +9 and -11 counts when rest - d0 is odd, so at beta = 0 the two
symbols came out at P(ok) = 0.638 and 0.667 where symmetry forces them equal.
The artifact was 20% of the effect under test. `_inputs` forces the parity and
asserts the exact bias, which is why it asserts rather than rounds.
"""
from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

from crnl.cme import first_passage
from crnl.networks.am_asymmetric import (
    am_asymmetric,
    basin_boundary,
    beta_critical,
    interior_fixed_points,
)
from crnl.networks.am_reversible import delta_star
from crnl.quasipotential import barriers

THETA = 0.80          # absorbing at |n_X - n_Y| >= THETA * delta* * Omega


def _shannon(p: float) -> float:
    p = min(max(p, 1e-15), 1.0 - 1e-15)
    return float(-p * np.log2(p) - (1 - p) * np.log2(1 - p))


def _inputs(gamma: float, omega: int, eps_frac: float):
    """The two initial states, biased by EXACTLY +-d0 counts. See the docstring."""
    ds = delta_star(gamma)
    d0 = max(1, int(round(eps_frac * ds * omega)))
    nb = int(round(omega * gamma / (1.0 + gamma)))
    rest = omega - nb
    if (rest - d0) % 2:
        d0 -= 1
    if d0 < 1:
        raise ValueError(f"input bias underflows the lattice at Omega={omega}")
    out = []
    for s in (+1, -1):
        nx, ny = (rest + s * d0) // 2, (rest - s * d0) // 2
        assert nx - ny == s * d0 and nx + ny + nb == omega
        out.append(np.array([nx, ny, nb], dtype=np.int64))
    return out, max(2, int(round(THETA * ds * omega)))


def mutual_information(gamma: float, beta: float, omega: int,
                       eps_frac: float) -> float:
    """I(input bit ; absorbed state) in bits, symmetric source, exact."""
    net = am_asymmetric(gamma, beta)
    (n_x, n_y), thr = _inputs(gamma, omega, eps_frac)

    def absorbing(s, thr=thr):
        return abs(int(s[0]) - int(s[1])) >= thr

    outs = []
    for n0 in (n_x, n_y):
        fp = first_passage(net, int(omega), float(omega), n0, absorbing, None)
        if not fp["valid"]:
            return float("nan")
        outs.append(fp["split"])            # P(X wins)
    a, b = outs
    return _shannon(0.5 * (a + b)) - 0.5 * (_shannon(a) + _shannon(b))


def part_fold(gammas) -> list[dict]:
    print(f"{'gamma':>6} {'beta_c':>8} {'delta*':>8}   "
          "boundary (x-y)/delta* at beta = 0.25/0.50/0.75 of beta_c")
    rows = []
    for g in gammas:
        bc = beta_critical(g)
        if not np.isfinite(bc):
            print(f"{g:>6} {'--':>8}   monostable at beta=0")
            continue
        ds = delta_star(g)
        bnds = [basin_boundary(g, f * bc) / ds for f in (0.25, 0.50, 0.75)]
        att = [p for p in interior_fixed_points(g, 0.5 * bc) if p["kind"] == "attractor"]
        print(f"{g:>6} {bc:>8.4f} {ds:>8.4f}   "
              + "  ".join(f"{b:+.4f}" for b in bnds)
              + "   attractors at " + ", ".join(f"({p['x']:.3f},{p['y']:.3f})" for p in att))
        rows.append({"gamma": g, "beta_c": bc, "delta_star": ds,
                     "boundary_over_ds": bnds})
    return rows


def part_wells(gammas, omegas, fracs) -> list[dict]:
    rows = []
    for g in gammas:
        bc = beta_critical(g)
        print(f"\ngamma={g}  beta_c={bc:.4f}"
              "   PREDICTION UNDER TEST: the sum falls with tilt")
        print(f"{'beta':>7} {'b/bc':>6} " +
              " ".join(f"{'Om=' + str(o):>27}" for o in omegas))
        print(f"{'':>14} " +
              " ".join(f"{'dW_x':>9}{'dW_y':>9}{'sum':>9}" for _ in omegas))
        for f in fracs:
            be = f * bc
            cells = []
            for om in omegas:
                try:
                    br = barriers(am_asymmetric(g, be), int(om), float(om),
                                  split=basin_boundary(g, be))
                    cells.append(f"{br['dW_x']:>9.5f}{br['dW_y']:>9.5f}"
                                 f"{br['dW_x'] + br['dW_y']:>9.5f}")
                    rows.append({"gamma": g, "beta": be, "omega": om,
                                 "dW_x": br["dW_x"], "dW_y": br["dW_y"]})
                except (ValueError, RuntimeError):
                    cells.append(f"{'out of range':>27}")
            print(f"{be:>7.4f} {f:>6.2f} " + " ".join(cells))
    return rows


def part_info(gammas, omegas, fracs, eps_fracs) -> list[dict]:
    rows = []
    for g in gammas:
        bc = beta_critical(g)
        for ef in eps_fracs:
            print(f"\ngamma={g} beta_c={bc:.4f} input eps={ef}*delta*   I/I(beta=0)")
            print(f"{'beta/bc':>8} " + " ".join(f"{'Om=' + str(o):>10}" for o in omegas))
            base = {}
            for f in fracs:
                vals = []
                for om in omegas:
                    I = mutual_information(g, f * bc, int(om), ef)
                    if f == 0.0:
                        base[om] = I
                    vals.append(I / base[om] if base.get(om) else float("nan"))
                    rows.append({"gamma": g, "beta_frac": f, "omega": om,
                                 "eps_frac": ef, "I_bits": I})
                print(f"{f:>8.2f} " + " ".join(f"{v:>10.4f}" for v in vals))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--part", choices=["fold", "wells", "info", "all"],
                    default="all")
    ap.add_argument("--out", type=pathlib.Path,
                    default=pathlib.Path("results/asymmetric_landscape.json"))
    args = ap.parse_args()

    out = {}
    if args.part in ("fold", "all"):
        out["fold"] = part_fold([0.05, 0.20, 0.35, 0.42, 0.45])
    if args.part in ("wells", "all"):
        out["wells"] = part_wells([0.42, 0.45], [150, 250], [0.0, 0.5, 0.97])
    if args.part in ("info", "all"):
        out["info"] = part_info([0.42], [120, 200, 300, 400],
                                [0.0, 0.25, 0.5, 0.75, 0.95], [0.10, 0.25])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
