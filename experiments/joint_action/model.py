"""Shared model definitions for the E-070 -> E-084 cascade WKB problem."""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


def kconst(r):
    r1, r2, r3 = r
    return r1 + r2 + r3, r1 * r2 * r3, r1 * r2 + r1 * r3 + r2 * r3


import os
# pair selection by env (defaults = original E-070 -> E-084 pair)
R2_1 = float(os.environ.get('R2_STAGE1', '1.005'))
R2_2 = float(os.environ.get('R2_STAGE2', '0.606'))
TAG = os.environ.get('PAIR_TAG', '')
A1, B1, C1 = kconst((0.15, R2_1, 3.0))   # stage 1: k1a, k2b, k2r
A2, B2, C2 = kconst((0.15, R2_2, 3.0))   # stage 2
S1 = R2_1
KH = float(os.environ.get('HILL_K', '1.0'))      # coupling threshold (§124.9); 1 = published
NH = float(os.environ.get('HILL_N', '4.0'))      # Hill exponent (§124.13); 4 = published
H3 = (3.0 / KH)**NH / (1 + (3.0 / KH)**NH)


def hill(x1):
    z = (x1 / KH)**NH
    return (z / (1 + z)) / H3


def dhill(x1):
    z = (x1 / KH)**NH
    return (NH * z / x1 / (1 + z) ** 2) / H3


def lam1(x1): return A1 * x1 * x1 + B1
def mu1(x1): return x1**3 + C1 * x1
def dlam1(x1): return 2 * A1 * x1
def dmu1(x1): return 3 * x1 * x1 + C1
def lam2(x1, x2): return hill(x1) * A2 * x2 * x2 + B2
def mu2(x2): return x2**3 + C2 * x2


def V1(x):
    return quad(lambda y: np.log(lam1(y) / mu1(y)), x, 3.0, epsabs=1e-14, epsrel=1e-13)[0]


def roots2(x1):
    r = np.roots([-1.0, hill(x1) * A2, -C2, B2])
    return np.sort(r.real[np.abs(r.imag) < 1e-9])


def _disc(x1):
    # number of real roots changes; use discriminant of cubic x^3 - h a x^2 + c x - b
    a = -hill(x1) * A2; b = C2; c = -B2
    return 18 * a * b * c - 4 * a**3 * c + a * a * b * b - 4 * b**3 - 27 * c * c


_g = np.linspace(0.3, 3.0, 2701); _d = np.array([_disc(x) for x in _g])
_k = np.where(np.sign(_d[:-1]) != np.sign(_d[1:]))[0][-1]
EDGE = brentq(_disc, _g[_k], _g[_k + 1], xtol=1e-15)
X2SN = hill(EDGE) * A2 / 3.0  # approx; refine below
# saddle-node point: double root of f(x)=x^3 - hA x^2 + C x - B, f'=0 -> larger root of 3x^2-2hAx+C
_h = hill(EDGE) * A2
X2SN = (2 * _h + np.sqrt(4 * _h * _h - 12 * C2)) / 6.0


def mid_or_ghost(x1):
    """Target curve: middle root of stage-2 cubic for x1 >= EDGE, saddle-node x2 for x1 < EDGE."""
    if x1 < EDGE:
        return X2SN
    r = roots2(x1)
    if len(r) == 3:
        return r[1]
    return X2SN


def H(x, p):
    x1, x2 = x; p1, p2 = p
    return (lam1(x1) * np.expm1(p1) + mu1(x1) * np.expm1(-p1)
            + lam2(x1, x2) * np.expm1(p2) + mu2(x2) * np.expm1(-p2))


def ham_rhs(t, y):
    x1, x2, p1, p2, S = y
    l1, m1 = lam1(x1), mu1(x1)
    l2, m2 = lam2(x1, x2), mu2(x2)
    e1, f1 = np.exp(p1), np.exp(-p1)
    e2, f2 = np.exp(p2), np.exp(-p2)
    dx1 = l1 * e1 - m1 * f1
    dx2 = l2 * e2 - m2 * f2
    dl2dx1 = dhill(x1) * A2 * x2 * x2
    dl2dx2 = hill(x1) * A2 * 2 * x2
    dm2dx2 = 3 * x2 * x2 + C2
    dp1 = -(dlam1(x1) * (e1 - 1) + dmu1(x1) * (f1 - 1) + dl2dx1 * (e2 - 1))
    dp2 = -(dl2dx2 * (e2 - 1) + dm2dx2 * (f2 - 1))
    return [dx1, dx2, dp1, dp2, p1 * dx1 + p2 * dx2]


def Lag1d(v, l, m):
    """Legendre transform of l(e^p-1)+m(e^-p-1)."""
    s = np.sqrt(v * v + 4 * l * m)
    return v * np.log((v + s) / (2 * l)) - s + l + m
