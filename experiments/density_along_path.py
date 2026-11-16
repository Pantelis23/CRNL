"""§134.12 D3 (registered in density_saddle.py before running): L(x) = ln[P(x)/psi(x1)] + Omega W(x) along the BVP orbit.
The orbit is saved by re-running bvp_gamma.py's continuation (HILL_K=1.2, tol 1e-6, drop 0) to /tmp/claude-1000/bvp_xk_orbit.npy."""
import json, math, sys
import numpy as np
import experiments.episode_committor as ec
from experiments.sigma_flux import ARMS, element
ec.SOLVER = "kr"
Yo = np.load("/tmp/claude-1000/bvp_xk_orbit.npy")          # rows x1, x2, p1, p2, S
A, B = element(ARMS["XK"][0]), element(ARMS["XK"][1])
# points along the orbit, by stage-2 coordinate (the path descends in x2)
targets = [2.5, 2.0, 1.6, 1.3, 1.0, 0.8, 0.7, 0.65, 0.62]
pts = []
for x2t in targets:
    k = int(np.argmax(Yo[1] < x2t)); pts.append((x2t, Yo[0, k], Yo[1, k], Yo[4, k]))
out = {}
for om in (200, 400, 800, 1200, 1600):
    m = A.cap(om) + 1
    sad1, sad2 = int(np.ceil(A.r2 * om)), int(np.ceil(B.r2 * om))
    ec.qprocess_J(A, B, om, m, sad1, sad2, lower="rail")
    L = ec.qprocess_J.last
    P = L["P"] / ec.PIN_SCALE / L["psi_n1"][:, None]          # the original measure
    n1s, base2 = L["n1s"], L["base2"]
    row = {}
    for x2t, x1, x2, W in pts:
        a, b = x1 * om - n1s[0], x2 * om - base2
        i, j = int(math.floor(a)), int(math.floor(b)); fa, fb = a - i, b - j
        lp = ((1 - fa) * (1 - fb) * math.log(P[i, j]) + fa * (1 - fb) * math.log(P[i + 1, j])
              + (1 - fa) * fb * math.log(P[i, j + 1]) + fa * fb * math.log(P[i + 1, j + 1]))
        row[f"{x2t}"] = {"x1": x1, "x2": x2, "W": W, "L": lp + om * W}
    out[om] = row
    print(om, " ".join(f"{k}:{v['L']:+.3f}" for k, v in row.items()), flush=True)
oms = sorted(out)
print("slope dL/dlnOmega per point, windows 200-400 400-800 800-1200 1200-1600:")
for x2t, x1, x2, W in pts:
    s = [(out[b][f"{x2t}"]["L"] - out[a][f"{x2t}"]["L"]) / math.log(b / a) for a, b in zip(oms, oms[1:])]
    print(f"  x=({x1:.3f},{x2:.3f}) W={W:.5f}: " + " ".join(f"{v:+.3f}" for v in s))
json.dump({str(k): v for k, v in out.items()}, open("results/density_along_path.json", "w"), indent=1)
