"""§134.25 (T-CASC-ba): where along the incoming optimal path does the saddle excess build up?

§134.12 read L = ln(P/psi) + Omega W along the BVP orbit at FIXED PHYSICAL points: regular (slope -> 0) everywhere it
looked. §134.21 read the excess across H along x2 at fixed stretched distance: it rises up the incoming side. Here: points
ON THE ORBIT at fixed STRETCHED distance from H, r = sqrt(Omega) |D^{-1/2} (x - H)| (D = diag(lam + mu) at H), so each
point moves toward H as Omega grows, keeping its place relative to the inner scale. W(x) is the orbit's own action S.

Predictions, written before running (XK, cap 1.1, Omega = 400 / 800 / 1200 / 1600; r = 1, 2, 5, 10, 20, 40):
  P1  the last-window slope dL/dlnOmega at fixed r is monotone in r: -> H's +0.37 as r -> 0 and -> 0 for large r.
  P2  the crossover r_c (the r where the slope is half of H's) is read per window. If r_c stays put across windows the
      excess is an inner-scale structure (built within r_c stretched units of H, so it scales as Omega^{-1/2}); if r_c
      grows with Omega it lives on an intermediate scale between Omega^{-1/2} and O(1).
  The verdict is P2's reading; P1 failing (non-monotone, or no approach to 0 by r = 40) is recorded as such.
  My expectation: r_c of order a few, roughly constant -- an inner-scale structure fed from the incoming side.

Usage: CRNL_CAP_MULT=1.1 PYTHONPATH=. python experiments/excess_along_path.py 400,800,1200,1600 kr
(regenerates the BVP orbit first: experiments/joint_action/bvp_gamma.py's continuation, HILL_K = 1.2)
"""
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PAIR = os.environ.get("PAIR", "XK")
HK = {"XK": "1.2", "X": "1.0", "XN": "1.2"}[PAIR]
HN = {"XN": "2"}.get(PAIR, "4")
ORBIT = Path(f"/tmp/claude-1000/bvp_{PAIR.lower()}_orbit.npy")


def orbit():
    if ORBIT.exists():
        return np.load(ORBIT)
    code = (
        "import sys; sys.argv=['x']; __file__='bvp_gamma.py'\n"
        "exec(open('bvp_gamma.py').read().split('if __name__')[0])\n"
        "tg, Yg = mam_guess(); t, Y = trim(tg, Yg, 3e-2); epsA = float(np.linalg.norm(Y[0, :2] - A[:2]))\n"
        "for step in range(40):\n"
        "    r, W, oa, ob, dr = solve(t, Y, 1e-6, epsA)\n"
        "    if max(oa, ob) < 3e-7 or r.status: break\n"
        "    t, Y, epsA = extend(r, 0.2 if oa > 3e-7 else 0.0, 1.0 if ob > 3e-7 else 0.0)\n"
        f"np.save('{ORBIT}', r.sol(np.linspace(0, r.x[-1], 200001)))\n"
        "print('orbit W', W)\n")
    env = dict(os.environ, HILL_K=HK, HILL_N=HN, BVP_DROP="0")
    subprocess.run([sys.executable, "-c", code], cwd=ROOT / "experiments" / "joint_action", env=env, check=True)
    return np.load(ORBIT)


if __name__ == "__main__":
    import experiments.episode_committor as ec
    from experiments.sigma_flux import ARMS, element
    if len(sys.argv) > 2:
        ec.SOLVER = sys.argv[2]
    Yo = orbit()
    A, B = element(ARMS[PAIR][0]), element(ARMS[PAIR][1])
    H = np.array([3.0, B.r2])
    sys.path.insert(0, str(ROOT / "experiments" / "joint_action"))
    os.environ["HILL_K"] = HK
    os.environ["HILL_N"] = HN
    from model import lam1, mu1, lam2, mu2
    Dh = np.sqrt(np.array([lam1(3.0) + mu1(3.0), lam2(3.0, B.r2) + mu2(B.r2)]))
    dist = np.sqrt(((Yo[0] - H[0]) / Dh[0]) ** 2 + ((Yo[1] - H[1]) / Dh[1]) ** 2)   # |D^{-1/2}(x - H)|
    rs = (1, 2, 5, 10, 20, 40)
    out = {}
    for om in [int(v) for v in sys.argv[1].split(",")]:
        m = A.cap(om) + 1
        sad1, sad2 = int(np.ceil(A.r2 * om)), int(np.ceil(B.r2 * om))
        ec.qprocess_J(A, B, om, m, sad1, sad2, lower="rail")
        L = ec.qprocess_J.last
        P = L["P"] / ec.PIN_SCALE / L["psi_n1"][:, None]
        n1s, base2 = L["n1s"], L["base2"]
        row = {}
        for r in rs:
            target = r / math.sqrt(om)
            k = int(np.argmin(np.abs(dist - target)))                 # the orbit point at that stretched distance
            x1, x2, W = Yo[0, k], Yo[1, k], Yo[4, k]
            a, b = x1 * om - n1s[0], x2 * om - base2
            i, j = int(math.floor(a)), int(math.floor(b)); fa, fb = a - i, b - j
            lp = ((1 - fa) * (1 - fb) * math.log(P[i, j]) + fa * (1 - fb) * math.log(P[i + 1, j])
                  + (1 - fa) * fb * math.log(P[i, j + 1]) + fa * fb * math.log(P[i + 1, j + 1]))
            row[str(r)] = {"x1": float(x1), "x2": float(x2), "W": float(W), "L": lp + om * W,
                           "dist_err": float(abs(dist[k] - target) / target)}
        out[str(om)] = row
        print(om, " ".join(f"r{k}:{v['L']:+.3f}" for k, v in row.items()), flush=True)
    name = "excess_along_path.json" if PAIR == "XK" else f"excess_along_path_{PAIR}.json"
    (ROOT / "results" / name).write_text(json.dumps(out, indent=1))
    ks = sorted(out, key=int)
    for r in rs:
        s = [(out[b][str(r)]["L"] - out[a][str(r)]["L"]) / math.log(int(b) / int(a)) for a, b in zip(ks, ks[1:])]
        print(f"r={r:>3}: slope per window " + " ".join(f"{v:+.3f}" for v in s))

# OUTCOME (§134.25): slopes dL/dlnOmega per window (400-800 / 800-1200 / 1200-1600) at fixed stretched r ON the orbit:
#   r=1: +0.176 / +0.249 / +0.288 (toward H's +0.37) | r=2: -0.147 / +0.029 / +0.107 | r=5: -0.916 / -0.958 / -0.682
#   r=10: -0.351 / -0.480 / -0.619 | r=20: -0.100 / -0.248 / -0.288 | r=40: only 1200-1600 valid (+0.404); at
#   Omega = 400/800 the r=40 point lies beyond the orbit's reach and snapped to the pin (L = 0) -- invalid, excluded.
# P1 FAILED: not monotone. A DEPLETION zone at r ~ 5-20 (density falls further below the WKB continuation as Omega
# grows, slopes to -0.9) separates the regular outer path from the excess inside r ~ 2. P2's half-crossover is
# ill-defined on a non-monotone profile -- recorded, not read. My expectation (monotone, r_c ~ few) was WRONG.

# REGISTERED BEFORE THE PAIR-X RUN (§134.26, suspect 6's kill test). Depletion depth = the most negative last-window
# (1200 -> 1600) slope over r = 5, 10, 20. XK: -0.682. X's orbit lies 5x closer to Sigma in tube units (§134.18), so
#   X's depth < -0.682 - 0.1 (deeper)            -> suspect 6 (effective grazing) SURVIVES;
#   X's depth >= -0.682 (not deeper)              -> suspect 6 DIES;
#   in between                                    -> UNRESOLVED.
# OUTCOME (§134.26), pair X, last window 1200 -> 1600: r=5 -0.856, r=10 -1.094, r=20 -0.565 (XK: -0.682 / -0.619 / -0.288).
# X's depth -1.094 < -0.782 -> suspect 6 SURVIVES as registered; deeper at every r in the zone. Caveat (rule 9): X and
# XK also differ in coupling, W* and Sigma's shape, so this supports suspect 6 but does not isolate it.

# REGISTERED BEFORE THE PAIR-XN RUN (§134.27, the isolating test). XN = XK with Hill exponent 2: same mu = 5.21 at H,
# different c and Sigma. Measured first (setup, graze_gap.py): the path-to-Sigma gap in tube units at Omega = 1600 is
# X 0.05-0.12 < XK 0.22-0.59 < XN 0.32-0.89. Suspect 6 needs depletion depth to order with the gap at FIXED mu:
#   XN's depth > -0.682 + 0.05 (shallower than XK)  -> suspect 6 SURVIVES, with mu held fixed;
#   XN's depth < -0.682 - 0.05 (deeper than XK)     -> suspect 6 DIES;
#   otherwise UNRESOLVED.
# OUTCOME (§134.27), XN, last window: r=5 -0.517, r=10 -0.717, r=20 -0.251 (XK -0.682 / -0.619 / -0.288).
# XN's depth -0.717 is within +-0.05 of XK's -0.682 -> UNRESOLVED as registered; mixed per r. Design flaw (rule 19):
# XN's gap is only ~1.5x XK's (vs X's 5x), too little lever to show an ordering -- the premise was not checked first.
