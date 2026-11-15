"""§134.9 (T-CASC-aw, measured part): the shape of the Sigma exit's core at the joint saddle, from the UNBINNED flux.

Why. §134.8: the rise of the total's a sits in the saddle bins (x1 in [2.95, 3.05]), and W's Hessian along Sigma at the
saddle is exactly zero on the instanton manifold. Suspect: eps ~ |delta|^3 along Sigma near B, so the core is
Omega^{-1/3} wide rather than Omega^{-1/2}, giving a_inf = 1/2 - 1/3 = +1/6. The stored 0.05 bins cannot resolve the core.

How. cut_flux(..., "sigma", raw=True) gives the committed flux per n1 (resolution 1/Omega), XK, cap 1.1. The flux
density per unit x1 is flux * Omega. Near x1 = 3: the peak, its half width at half maximum w (each side, by linear
interpolation), and the core's mass inside w.

Predictions, written before running:
  C1  The core's width scales as w ~ Omega^-nu. Suspect: nu = 1/3; a Gaussian core: nu = 1/2. Read per window
      (rule 21: nu_eff(window) is reported with its window, never as an exponent):
        every window's nu_eff in [0.25, 0.42] and not drifting toward 1/2 -> supports the cubic core;
        nu_eff >= 0.45 in the largest window                             -> the suspect is dead;
        otherwise UNRESOLVED.
  C2  (instrument) re-binning the raw flux into §133's 0.05 bins reproduces the stored flux_density at Omega = 1600
      to <= 1e-9 relative, on the bins that carry flux.
  My expectation, stated so it can be wrong: nu_eff ~ 0.3-0.4, i.e. the cubic suspect survives.
  -> C2 HELD: the re-binned raw flux equals the stored 1600 bins exactly (0.0).
  -> C1 AS WRITTEN MEASURED THE LATTICE, NOT THE CORE: FWHM = 1.29 / 1.25 / 1.24 / 1.25 x (1/Omega) at Omega = 400 /
     800 / 1200 / 1600, peak jumping 2.81 / 2.96 / 3.025 / 2.96. The Sigma cut is a staircase on the lattice and the
     stage-1 moves across each step land in one n1, so the raw per-n1 flux is spiky; its half width is one lattice
     spacing (rule 19: the criterion could only report 1/Omega).

REVISED BEFORE RE-RUNNING (C1'): the core's width from the density smoothed by a Gaussian kernel of width h = 0.01 and,
separately, 0.02 (>> the lattice spacing 1/Omega <= 0.0025, and << the 0.05 bins). Same nu bands as C1, read on each h;
if the two h disagree on the verdict, the width is UNRESOLVED at this resolution. Also reported: the flux within
|x1 - 3| <= 0.05 as a share of the total, against Omega (§134.8 read 28% at Omega = 5000 from bins).

Usage: CRNL_CAP_MULT=1.1 PYTHONPATH=. python experiments/sigma_core.py 400,800,1200,1600 [kr]
"""
import json
import sys
from pathlib import Path

import numpy as np

import experiments.episode_committor as ec
from experiments.sigma_flux import ARMS, cut_flux, element, sigma_table


def smoothed(x1, f, h):
    """Density per unit x1, Gaussian-smoothed at width h (flux conserved)."""
    g = np.linspace(2.4, 3.4, 2001)
    K = np.exp(-0.5 * ((g[:, None] - x1[None, :]) / h) ** 2) / (h * np.sqrt(2 * np.pi))
    return g, K @ f


def core(x1, dens):
    """Peak near x1 = 3 and its half widths at half maximum (left, right), by linear interpolation."""
    m = (x1 > 2.5) & (x1 < 3.3)
    xs, ds = x1[m], dens[m]
    k = int(np.argmax(ds)); pk = ds[k]; half = 0.5 * pk
    i = k
    while i > 0 and ds[i] > half:
        i -= 1
    xl = np.interp(half, [ds[i], ds[i + 1]], [xs[i], xs[i + 1]])
    j = k
    while j < len(ds) - 1 and ds[j] > half:
        j += 1
    xr = np.interp(half, [ds[j], ds[j - 1]], [xs[j], xs[j - 1]])
    inside = (xs >= xl) & (xs <= xr)
    return {"x_peak": float(xs[k]), "peak": float(pk), "w_left": float(xs[k] - xl), "w_right": float(xr - xs[k]),
            "fwhm": float(xr - xl), "core_mass": float(ds[inside].sum() * (x1[1] - x1[0]))}


if __name__ == "__main__":
    if len(sys.argv) > 2:
        ec.SOLVER = sys.argv[2]
    A, B = element(ARMS["XK"][0]), element(ARMS["XK"][1])
    sig = sigma_table(A, B)
    out_p = Path("results/sigma_core.json")
    out = json.loads(out_p.read_text()) if out_p.exists() else {}
    stored = json.loads(Path("results/sigma_flux.json").read_text())
    for om in [int(v) for v in sys.argv[1].split(",")]:
        x1, f = cut_flux(A, B, om, "sigma", sig, raw=True)
        dens = f * om                                              # per unit x1 (each n1 spans 1/Omega)
        rec = {"raw": core(x1, dens), "total": float(f.sum()),
               "share_within_0.05": float(f[np.abs(x1 - 3.0) <= 0.05].sum() / f.sum())}
        for h in (0.01, 0.02):
            g, d = smoothed(x1, f, h)
            rec[f"h{h:g}"] = core(g, d)
        key = f"XK:{om}:sigma:cap1.1"
        if key in stored:                                          # C2: the instrument against §133's stored bins
            st = np.array(stored[key]["flux_density"])
            rb = np.array([f[(x1 >= 1.0 + 0.05 * k) & (x1 < 1.05 + 0.05 * k)].sum() / 0.05 for k in range(44)])
            mk = st > 1e-6 * st.max()
            rec["rebin_vs_stored"] = float(np.max(np.abs(rb[mk] / st[mk] - 1)))
        out[str(om)] = rec
        print(om, {k: (round(v, 6) if isinstance(v, float) else v) for k, v in rec.items()}, flush=True)
        out_p.write_text(json.dumps(out, indent=1))

# OUTCOME (§134.9), Omega = 400-1600 (local; 2000/2800 were stopped by the host for low memory before any cell finished):
#   FWHM (h = 0.01): 0.3191 / 0.2365 / 0.2124 / 0.1797; (h = 0.02): 0.3369 / 0.2567 / 0.2186 / 0.1856.
#   nu_eff per window, h = 0.01 / 0.02: 400-800 0.43 / 0.39; 800-1200 0.27 / 0.40; 1200-1600 0.58 / 0.57;
#   400-1600 overall 0.41 / 0.43. Share within |x1 - 3| <= 0.05: 0.044 / 0.073 / 0.096 / 0.120.
#   C1' BY ITS REGISTERED RULE PRINTS "SUSPECT DEAD" (largest window 0.58 / 0.57 >= 0.45). But the window-to-window
#   scatter (0.27 -> 0.58, +-0.15) exceeds the 1/3-vs-1/2 gap the rule is meant to resolve, so the verdict rests on
#   the noisiest window (ratio 4/3). Recorded as printed and NOT acted on: UNRESOLVED (rule 19; rule 14 -- a
#   withdrawal is a claim too). Needs windows with ratio >= 2 at larger Omega.
