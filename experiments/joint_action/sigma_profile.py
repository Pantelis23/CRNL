"""§133 input: the action along Sigma, eps(x1) = W_Sigma(x1) - W*, from the ray family.

Each ray from the linear-noise ellipse around (3, 3) is followed to its first hit of Sigma (x1 > s1). The hit
point (x1c, x2c) and its action S give W restricted to Sigma at x1c, so long as that ray is the cheapest one
reaching that point (no caustic crossing first). Rays near the limiting ray (graze.py) hit Sigma ever closer to
the joint saddle with S -> W*. The rest fan out along Sigma toward low x1.

Usage (from this directory): HILL_K=1.2 python sigma_profile.py EPS RTOL NSIG NRAYS
Writes sigma_profile_<tag>.npy with rows (d, x1c, x2c, S).
"""
import sys
import numpy as np
EPSa, RTOLa, NSIGa, NR = sys.argv[1:5]
sys.argv = ["x", EPSa, RTOLa, NSIGa]
exec(open('main_sigma.py').read().split("if __name__")[0])

nr = int(NR)
ds = -np.geomspace(1e-12, 1.0, nr)
rows = []
for d in ds:
    y, h = ray(d)
    if y is not None:
        rows.append((d, y[0], y[1], y[4]))
rows = np.array(rows)
tag = f"K{KH:g}_N{NH:g}_eps{EPSa}"
np.save(f"sigma_profile_{tag}.npy", rows)
o = np.argsort(rows[:, 1])
print(f"{len(rows)} hits of {nr} rays; W at hits from {rows[:,3].min():.9f} to {rows[:,3].max():.9f}")
for r in rows[o][:: max(1, len(rows) // 40)]:
    print(f"  d={r[0]:+.3e}  x1c={r[1]:.5f}  x2c={r[2]:.5f}  S={r[3]:.9f}")

# refine: bisect for the limiting ray d* (hit / no hit), then sample geometrically toward it
a, b = None, None
dd = ds[np.isin(ds, rows[:, 0])]
a = dd.min()                                    # the most negative d that still hits
nxt = ds[ds < a]
b = nxt.max() if len(nxt) else a * 1.5          # the first that misses
for _ in range(80):
    m = 0.5 * (a + b)
    if ray(m)[0] is not None:
        a = m
    else:
        b = m
dstar = a
rows2 = []
for f in np.geomspace(1e-9, 0.6, int(NR)):
    d = dstar * (1 - f)                         # |d| just inside the limit, moving away from it
    y, h = ray(d)
    if y is not None:
        rows2.append((d, y[0], y[1], y[4]))
rows2 = np.array(rows2)
allr = np.vstack([rows, rows2])
np.save(f"sigma_profile_{tag}.npy", allr)
o = np.argsort(rows2[:, 1])
print(f"d* = {dstar:.15e}; refined hits {len(rows2)}")
for r in rows2[o][:: max(1, len(rows2) // 40)]:
    print(f"  d={r[0]:+.6e}  x1c={r[1]:.5f}  x2c={r[2]:.5f}  S={r[3]:.9f}")
