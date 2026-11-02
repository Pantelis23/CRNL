"""Robust ray shooting. Rays launched from ellipse dx^T Z dx = 2*EPS around (3,3),
p = Z dx, parametrised by momentum angle phi (phi = pi is the pure stage-1 ray, p2 == 0).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_continuous_lyapunov
from model import *

xs = np.array([3.0, 3.0])
f1p = dlam1(3.0) - dmu1(3.0)
Amat = np.array([[f1p, 0.0], [dhill(3.0) * A2 * 9.0, hill(3.0) * A2 * 6.0 - (27 + C2)]])
D = np.diag([lam1(3.0) + mu1(3.0), lam2(3.0, 3.0) + mu2(3.0)])
Cm = solve_continuous_lyapunov(Amat, -D)
Z = np.linalg.inv(Cm)


def inside(y):
    return (np.all(np.isfinite(y)) and 0.9 < y[0] < 3.5 and 0.05 < y[1] < 3.6
            and abs(y[2]) < 3 and abs(y[3]) < 3)


def rhs(t, y):
    if not inside(y):
        return [0.0] * 5
    return ham_rhs(t, y)


def ev_target(t, y):
    if not inside(y):
        return 1.0
    return y[1] - mid_or_ghost(y[0])
ev_target.terminal = True
ev_target.direction = -1


def ev_s1(t, y):
    if not np.all(np.isfinite(y)):
        return 1.0
    return y[0] - S1
ev_s1.terminal = True
ev_s1.direction = -1


def shoot(phi, eps, rtol=1e-11, tmax=3000.0, dense=False):
    p = np.array([np.cos(phi), np.sin(phi)])
    dx = Cm @ p
    dx *= np.sqrt(2 * eps / (dx @ Z @ dx))
    x0 = xs + dx; p0 = Z @ dx
    sol = solve_ivp(rhs, (0, tmax), [x0[0], x0[1], p0[0], p0[1], eps], method='DOP853',
                    rtol=rtol, atol=rtol * 1e-3, events=[ev_target, ev_s1], dense_output=dense)
    if len(sol.t_events[0]):
        return 'hit', sol.y_events[0][0], sol
    if len(sol.t_events[1]):
        return 's1', sol.y_events[1][0], sol
    return ('out' if not inside(sol.y[:, -1]) else 'tmax'), sol.y[:, -1], sol
