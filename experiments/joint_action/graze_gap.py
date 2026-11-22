"""Distance between each pair's BVP orbit and Sigma near the saddle (x1 in [1.9, 2.95]), in x2 and in units of the
orbit's transverse width at Omega = 1600 and 4000 (width = (Omega Z22)^{-1/2} with Z22 from the orbit's own p2/x2 slope
is not available; use the saddle's noise-stable Z22 from the linearisation as the scale)."""
import os, sys, json
sys.argv = ["x"]; __file__ = "bvp_gamma.py"
src = open("bvp_gamma.py").read()
exec(src.split("if __name__")[0])
import warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, "../..")
from experiments.sigma_flux import ARMS, element, sigma_table
pair = os.environ["PAIR"]
tg, Yg = mam_guess(); t, Y = trim(tg, Yg, 3e-2); epsA = float(np.linalg.norm(Y[0, :2] - A[:2]))
for step in range(40):
    r, W, oa, ob, dr = solve(t, Y, 1e-6, epsA)
    if max(oa, ob) < 3e-7 or r.status: break
    t, Y, epsA = extend(r, 0.2 if oa > 3e-7 else 0.0, 1.0 if ob > 3e-7 else 0.0)
Yo = r.sol(np.linspace(0, r.x[-1], 40001))
EA, EB = element(ARMS[pair][0]), element(ARMS[pair][1])
xs, sg = sigma_table(EA, EB); xs = np.array(xs); sg = np.array(sg, float)
w, V = np.linalg.eig(jac(B)); k = int(np.argmin(np.abs(w + w[w.real > 0].real.min())))
Z22 = abs(V[3, k].real / V[1, k].real)
gaps = []
for x1t in (1.9, 2.2, 2.5, 2.7, 2.85, 2.95):
    j = np.where((Yo[0] > x1t - 0.01) & (Yo[0] < x1t + 0.01) & (Yo[1] < R2_2 + 0.3))[0]
    if len(j) == 0: continue
    x2o = Yo[1, j].min(); s = float(np.interp(x1t, xs, sg))
    gaps.append((x1t, x2o, s, x2o - s))
print(pair, "W*", W, "Z22", Z22)
for x1t, x2o, s, g in gaps:
    print(f"  x1={x1t}: orbit x2={x2o:.4f} Sigma={s:.4f} gap={g:+.4f}  gap*sqrt(Om Z22): 1600 {g*np.sqrt(1600*Z22):.2f}  4000 {g*np.sqrt(4000*Z22):.2f}")
