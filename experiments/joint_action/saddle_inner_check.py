"""§134.19 revision (registered in saddle_inner_grid.py): W1 normalisation and W2 Monte Carlo of saddle_inner's Gaussian."""
import json, math, os, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import saddle_inner as si

mu, c, _ = si.normal_form()
B = np.array([[1.0, 0.0], [c, -mu]])
e = si.incoming_direction(B)
out = {"mu": mu, "c": c}
# W1: normalisation of G_t over y, by quadrature in B's eigencoordinates (Y2 = X2 - k X1, Jacobian 1), where the
# Gaussian's long thin ridge is axis-aligned (an X1-X2 grid cannot resolve it at large t)
k = c / (1 + mu)
x0 = 2.0 * e
w1 = {}
for t in (0.3, 1.0, 3.0, 10.0):
    s11 = math.expm1(2 * t) / 2
    s12 = -k * (-math.expm1((1 - mu) * t)) / (mu - 1)
    s22 = (1 + k * k) * (-math.expm1(-2 * mu * t)) / (2 * mu)
    m1 = math.exp(t) * x0[0]
    m2Y = math.exp(-mu * t) * (x0[1] - k * x0[0])
    sd1 = math.sqrt(s11); sdc = math.sqrt(s22 - s12 * s12 / s11)
    g1 = np.linspace(m1 - 9 * sd1, m1 + 9 * sd1, 1201)
    tot = 0.0
    for a1 in g1:
        cm = m2Y + s12 / s11 * (a1 - m1)                           # conditional mean of Y2 given X1
        gY = np.linspace(cm - 9 * sdc, cm + 9 * sdc, 1201)
        ys = np.column_stack([np.full_like(gY, a1), gY + k * a1])  # back to X2 = Y2 + k X1
        lg = np.array([si.green_terms(B, y, x0, [t])[0][0] for y in ys])
        tot += np.trapezoid(np.exp(lg), gY)
    tot *= (g1[1] - g1[0])
    w1[str(t)] = tot
    print(f"W1 t={t}: int G dy = {tot:.8f}", flush=True)
out["W1"] = w1
# W2: Monte Carlo occupation density from the source, against the time-integrated Gaussian
rng = np.random.default_rng(12345)
ts = np.concatenate([np.geomspace(1e-5, 1.0, 600), np.linspace(1.0, 40.0, 40000)[1:]])
targets = {"0": np.zeros(2), "e1": 1.0 * e}
gauss = {k: math.exp(si.ln_v_and_S(B, y, x0, ts)[0]) for k, y in targets.items()}
w2 = {}
for dt in (1e-3, 5e-4):
    N = 1000000
    X = np.tile(x0, (N, 1)); alive = np.ones(N, bool)
    batch = np.arange(N) % 10
    occ = {k: np.zeros(10) for k in targets}
    half = 0.1
    steps = 0
    while alive.any() and steps < 200000:
        idx = np.where(alive)[0]
        Xa = X[idx]
        for k, y in targets.items():
            inb = np.all(np.abs(Xa - y) <= half, axis=1)
            occ[k] += np.bincount(batch[idx][inb], minlength=10) * dt
        X[idx] = Xa + (Xa @ B.T) * dt + rng.standard_normal((len(idx), 2)) * math.sqrt(dt)
        alive[idx] = np.linalg.norm(X[idx], axis=1) <= 12
        steps += 1
    rec = {}
    for k in targets:
        per = occ[k] / (N / 10) / (2 * half) ** 2                 # density estimate per batch
        dens, se = per.mean(), per.std(ddof=1) / math.sqrt(10)
        rec[k] = {"mc": dens, "mc_se": se, "gauss": gauss[k], "ratio": dens / gauss[k],
                  "z": (dens - gauss[k]) / se}
    w2[f"{dt:g}"] = rec
    print(f"W2 dt={dt}: " + " ".join(f"{k}: mc/gauss={v['ratio']:.4f} (z={v['z']:+.2f})" for k, v in rec.items()), flush=True)
out["W2"] = w2
(HERE.parents[1] / "results" / "saddle_inner_check_xk.json").write_text(json.dumps(out, indent=1))

# OUTCOME (§134.19): W1 HELD -- int G dy = 1.00000000 at t = 0.3 / 1 / 3 / 10 (in eigencoordinates; an X1-X2 grid had
# read 9.195 at t = 10, a quadrature artefact of the thin tilted ridge). W2 HELD by its registered dt trend: MC/Gauss at
# H = 0.932 (z -3.2) at dt = 1e-3 -> 0.992 (z -0.3) at 5e-4; on the ray 0.988 -> 0.994 (1e6 paths, 10-batch sigma).
# The Gaussian computation behind §134.17 is verified. My expectation (§134.17 stands) held.
