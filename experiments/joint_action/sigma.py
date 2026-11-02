"""Joint deterministic separatrix Sigma: x2 = sigma(x1), by bisection in x2 at each x1.
Fallen <=> deterministic pair ends with x2 < 0.606 (stage-2 middle root at x1=3).
"""
import numpy as np
from scipy.integrate import solve_ivp
from model import *


def det_rhs(t, y):
    x1, x2 = y
    return [lam1(x1) - mu1(x1), lam2(x1, x2) - mu2(x2)]


def fallen(x1, x2, T=400.0):
    s = solve_ivp(det_rhs, (0, T), [x1, x2], method='LSODA', rtol=1e-11, atol=1e-13)
    assert abs(s.y[0, -1] - 3.0) < 1e-3, (x1, x2, s.y[:, -1])
    xe = s.y[1, -1]
    assert xe < 0.3 or xe > 2.5, xe  # unambiguous end state
    return xe < R2_2


def sigma(x1, tol=1e-12):
    lo, hi = 0.02, 3.4   # lo fallen, hi not
    if not fallen(x1, lo):
        return np.nan
    if fallen(x1, hi):
        return np.inf
    while hi - lo > tol:
        m = 0.5 * (lo + hi)
        if fallen(x1, m):
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


if __name__ == '__main__':
    x1g = np.concatenate([S1 + np.logspace(-10, -1, 19), np.linspace(1.11, 3.3, 45)])
    out = []
    for x1 in x1g:
        s = sigma(x1)
        out.append((x1, s))
        print(f'{x1:.10f} {s:.10f}  mid/ghost={mid_or_ghost(x1):.6f}', flush=True)
    np.save('sigma.npy', np.array(out))
