"""W* = min over the joint deterministic separatrix Sigma of the constrained quasi-potential
W_c(x1,x2) from (3,3), x1 > s1, by characteristic (ray) shooting.

Predictions (written before running, after computing Sigma and V1 only):
  * Sigma rises above x2=3.4 for x1 < ~1.03, so the free-x2 (p2=0) ray, which drives x1
    down the reversed stage-1 path, crosses Sigma at x1c ~ 1.02-1.04 with cost V1(x1c).
  * Near s1, V1(s1)-V1(x) ~ 0.183 (x-s1)^2, so the p2=0 ray alone gives Gamma ~ 1e-4.
  * 2-D optimisation (p2 < 0 pushes x2 down) can only lower W*; predict Gamma in [1e-4, 5e-3],
    argmin x1 in [1.02, 1.15], i.e. far below the CME's mean crossing x1 ~ 1.5.
  * The minimum does NOT hug x1 -> s1 exactly (Gamma > 0 strictly), but is close to it.
Usage: main_sigma.py EPS RTOL NSIG
"""
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from model import *
import shoot as sh
from sigma import sigma

EPS = float(sys.argv[1]); RTOL = float(sys.argv[2]); NSIG = int(sys.argv[3])

# Sigma table: dense near s1
xa = np.concatenate([S1 + np.geomspace(1e-4, 0.4, NSIG), np.linspace(S1 + 0.41, 3.4, NSIG // 4)])
sg = np.array([sigma(x, tol=1e-10) for x in xa])
sg[np.isinf(sg)] = 50.0      # everything fallen there
sg[np.isnan(sg)] = -50.0


def sig(x1):
    if x1 <= xa[0]:
        return 50.0
    return np.interp(x1, xa, sg)


def ev_sigma(t, y):
    if not sh.inside(y):
        return 1.0
    return y[1] - sig(y[0])
ev_sigma.terminal = True
ev_sigma.direction = -1


def ray(eta, tmax=3000.0):
    # eta: momentum direction p ~ (-1, eta); eta = 0 is the p2 = 0 (free-x2) ray
    p = np.array([-1.0, eta])
    dx = sh.Cm @ p
    dx *= np.sqrt(2 * EPS / (dx @ sh.Z @ dx))
    x0 = sh.xs + dx; p0 = sh.Z @ dx
    s = solve_ivp(sh.rhs, (0, tmax), [x0[0], x0[1], p0[0], p0[1], EPS], method='DOP853',
                  rtol=RTOL, atol=RTOL * 1e-3, events=[ev_sigma, sh.ev_s1])
    if len(s.t_events[0]):
        y = s.y_events[0][0]
        return y, H(y[:2], y[2:4])
    return None, None


def cost(eta):
    y, _ = ray(eta)
    return 10.0 if y is None else y[4]


if __name__ == '__main__':
    V1s = V1(S1)
    print(f'EPS={EPS:g} RTOL={RTOL:g} NSIG={NSIG}  V1(s1)={V1s:.10f}')
    # scan around phi = pi (p2 = 0 ray) on log offsets, both signs
    offs = np.concatenate([-np.geomspace(1e-9, 1.0, 46), [0.0], np.geomspace(1e-9, 1.0, 46)])
    rows = []
    for d in np.sort(offs):
        y, h = ray(-d)
        if y is not None:
            rows.append((d, *y, h))
    rows = np.array(rows)
    i = np.argmin(rows[:, 5])
    for r in rows[max(0, i - 4): i + 5]:
        print(f'  d={r[0]:+.3e} x1={r[1]:.6f} x2={r[2]:.6f} p=({r[3]:+.5f},{r[4]:+.5f}) S={r[5]:.9f} H={r[6]:.1e}')
    # p2 = 0 ray exactly
    y0, h0 = ray(0.0)
    print(f'p2=0 ray: hits Sigma at x1={y0[0]:.6f} x2={y0[1]:.6f} S={y0[4]:.10f}  V1(x1c)={V1(y0[0]):.10f}')
    # refine
    lo, hi = rows[max(0, i - 1), 0], rows[min(len(rows) - 1, i + 1), 0]
    res = minimize_scalar(lambda d: cost(-d), bounds=(lo, hi), method='bounded',
                          options={'xatol': 1e-13})
    y, h = ray(-res.x)
    print(f'OPT d={res.x:+.4e} x1*={y[0]:.6f} x2*={y[1]:.6f} W*={y[4]:.10f} H={h:.1e} '
          f'Gamma={V1s - y[4]:.4e}  V1(x1*)={V1(y[0]):.10f}')
    # global check: coarse full-circle scan
    full = [(e, cost(e)) for e in np.concatenate([-np.geomspace(1e-9, 1e3, 120), np.geomspace(1e-9, 1e3, 120)])]
    full = sorted(full, key=lambda t: t[1])[:3]
    print('full-circle best 3:', [(round(a, 4), round(b, 8)) for a, b in full])
