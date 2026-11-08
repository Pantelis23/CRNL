"""§124.11 -- the Eyring-Kramers prefactor of the never-dipped fall, from the rate equations alone.

WKB for the master equation: pi ~ K(x) exp(-Omega W(x)). Along a characteristic (x' = H_p, p' = -H_x):
  Z = Hess W:   dZ/dt   = -(H_xx + H_xp Z + (H_xp Z)^T + Z H_pp Z)            (Riccati)
  ln K:         dlnK/dt = -Tr H_xp - (1/2) Tr(H_pp Z) + H1(x, p)               (transport)
with (H_xp)_ik = d2H/dx_i dp_k and H1 = sum_r w1_r (e^{r.p} - 1), where w1 are the O(1/Omega) parts of
the CME propensities (n(n-1)/Omega = Omega(x^2 - x/Omega), n(n-1)(n-2)/Omega^2 = Omega(x^3 - 3x^2/Omega
+ ...)). Started at the stable point m = (3, 3) with Z = inverse linear-noise covariance, K = 1; read
at the joint saddle s = (3, r2_stage2) as the limiting ray grazes it. Then, as for a diffusion near s
(where p = 0 and the dynamics are locally Gaussian),

  k_fall = (lambda_+ / 2 pi) sqrt(det Z_m / |det Z_s|) (K_s / K_m) exp(-Omega W*),

lambda_+ the unstable eigenvalue of the rate equations at s. The dip rate needs no theory: it is
the exact 1-D stage-1 rate (chain killed at its first step below s1), so R_pred = k_fall / k_dip.

PREDICTIONS, WRITTEN BEFORE RUNNING.
  K1  VALIDATION (no coupling: hill == 1). The 2-D formula must reduce to stage 2's own escape:
      k_EK / k_exact (stage 2's 1-D chain, absorbed at its low rail) -> 1 with a deviation shrinking
      roughly like 1/Omega over Omega = 100 .. 1600 -- tests the Riccati, transport, H1 and det code.
  K2  THE TEST, absolute, rule 16. Pair XK (E-070 -> E-084@K1.2): R_meas / R_pred over Omega =
      800 .. 1600 approaches 1 -- its deviation from 1 shrinks with Omega and, extrapolated in 1/Omega,
      lies within 0.05 of 1. Equivalently C_pred = lim R_pred e^{-Gamma Omega} matches the measured
      plateau 0.63. A ratio settling at another constant means the formula misses an O(1) factor
      (e.g. the factor 2 between reaching the saddle and committing), and that is reported as such.
Usage: HILL_K=1.2 prefactor.py EPS RTOL NSIG [--uncoupled]

OUTCOMES (appended after running; predictions above left as written).
  An instrument failure first: the Riccati equation is UNSTABLE in forward time near the saddle (its
  physical root there repels; an error grows like 1/distance^2). Z is therefore taken at the saddle
  from the algebraic Riccati (a Lyapunov equation, A C + C A^T = -D, Z_s = C^{-1}), and ln K is read
  in a window |p| in [1e-3, 6e-3] on the descent into the saddle and extrapolated linearly.
  K1 HELD: uncoupled, k_EK / k_exact = 1.032 / 1.018 / 1.007 / 1.0024 / 1.0000 / 0.9989 at Omega =
     50 .. 1600 (W* = 0.35395995 vs exact 0.35395993; ln K_s 2.39699 vs exact 2.39923: ~0.2%).
  K2 NOT EVALUABLE as designed: on the coupled instanton (XK) the Hessian component along the
     separatrix collapses (Z11: 0.15 -> ~3e-4) while the path slides along Sigma into the saddle, and
     the Riccati solution diverges just before it (a caustic). The non-degenerate formula's Z_s and K_s
     do not exist along this path; the exit is degenerate. FINDINGS §124.11.
"""
import sys
import numpy as np
from scipy.integrate import solve_ivp

UNCOUPLED = "--uncoupled" in sys.argv
args = [a for a in sys.argv[1:] if not a.startswith("--")]
import model
if UNCOUPLED:                                   # patch before anything copies the bindings
    model.hill = lambda x1: 1.0
    model.dhill = lambda x1: 0.0
sys.argv = ["x"] + args[:3]
exec(open("main_sigma.py").read().split("if __name__")[0])
from model import A1, B1, C1, A2, B2, C2, KH, H3, R2_2, S1, V1


def h2(x1):
    if UNCOUPLED:
        return 0.0
    z = (x1 / KH) ** 4
    return 4 * z * (3 - 5 * z) / (x1 * x1 * (1 + z) ** 3) / H3


def aug_rhs(t, y):
    x1, x2, p1, p2, S, z11, z12, z22, lK = y
    if not sh.inside(y[:5]):
        return [0.0] * 9
    hh, dh, ddh = model.hill(x1), model.dhill(x1), h2(x1)
    e1, f1, e2, f2 = np.exp(p1), np.exp(-p1), np.exp(p2), np.exp(-p2)
    l1, m1 = A1 * x1 * x1 + B1, x1 ** 3 + C1 * x1
    l2, m2 = hh * A2 * x2 * x2 + B2, x2 ** 3 + C2 * x2
    l1_1, m1_1 = 2 * A1 * x1, 3 * x1 * x1 + C1
    l2_1, l2_2, m2_2 = dh * A2 * x2 * x2, 2 * hh * A2 * x2, 3 * x2 * x2 + C2
    Hp = np.array([l1 * e1 - m1 * f1, l2 * e2 - m2 * f2])
    Hx = np.array([l1_1 * (e1 - 1) + m1_1 * (f1 - 1) + l2_1 * (e2 - 1), l2_2 * (e2 - 1) + m2_2 * (f2 - 1)])
    Hpp = np.diag([l1 * e1 + m1 * f1, l2 * e2 + m2 * f2])
    Hxp = np.array([[l1_1 * e1 - m1_1 * f1, l2_1 * e2], [0.0, l2_2 * e2 - m2_2 * f2]])
    Hxx = np.array([[2 * A1 * (e1 - 1) + 6 * x1 * (f1 - 1) + ddh * A2 * x2 * x2 * (e2 - 1), 2 * dh * A2 * x2 * (e2 - 1)],
                    [2 * dh * A2 * x2 * (e2 - 1), 2 * hh * A2 * (e2 - 1) + 6 * x2 * (f2 - 1)]])
    H1 = -A1 * x1 * (e1 - 1) - 3 * x1 * x1 * (f1 - 1) - hh * A2 * x2 * (e2 - 1) - 3 * x2 * x2 * (f2 - 1)
    Z = np.array([[z11, z12], [z12, z22]])
    M = Hxp @ Z
    dZ = -(Hxx + M + M.T + Z @ Hpp @ Z)
    dlK = -np.trace(Hxp) - 0.5 * np.trace(Hpp @ Z) + H1
    return [Hp[0], Hp[1], -Hx[0], -Hx[1], p1 * Hp[0] + p2 * Hp[1], dZ[0, 0], dZ[0, 1], dZ[1, 1], dlK]


def ray_aug(eta):
    p = np.array([-1.0, eta])
    dx = sh.Cm @ p
    dx *= np.sqrt(2 * EPS / (dx @ sh.Z @ dx))
    x0 = sh.xs + dx; p0 = sh.Z @ dx
    y0 = [x0[0], x0[1], p0[0], p0[1], EPS, sh.Z[0, 0], sh.Z[0, 1], sh.Z[1, 1], 0.0]
    ev = lambda t, y: ev_sigma(t, y[:5]); ev.terminal = True; ev.direction = -1
    e1 = lambda t, y: sh.ev_s1(t, y[:5]); e1.terminal = True; e1.direction = -1
    s = solve_ivp(aug_rhs, (0, 3000.0), y0, method="DOP853", rtol=RTOL, atol=RTOL * 1e-3, events=[ev, e1])
    return s.y_events[0][0] if len(s.t_events[0]) else None


def approach_samples(pdir, tmax=120.0):
    """Integrate one ray (initial momentum direction pdir) and return rows (|p|, W, Z11, Z12, Z22, lnK)
    along its approach to the saddle (x1 > 2.9, before the closest point)."""
    p = np.asarray(pdir, float)
    dx = sh.Cm @ p
    dx *= np.sqrt(2 * EPS / (dx @ sh.Z @ dx))
    x0 = sh.xs + dx; p0 = sh.Z @ dx
    y0 = [x0[0], x0[1], p0[0], p0[1], EPS, sh.Z[0, 0], sh.Z[0, 1], sh.Z[1, 1], 0.0]
    s_ = solve_ivp(aug_rhs, (0, tmax), y0, method="DOP853", rtol=RTOL, atol=RTOL * 1e-3, dense_output=True)
    Y = s_.sol(np.linspace(0, s_.t[-1], 200001))
    pn = np.hypot(Y[2], Y[3])
    near = (Y[0] > 2.9) & (np.abs(Y[1] - R2_2) < 0.3)
    idx = np.nonzero(near)[0]
    j = idx[np.argmin(pn[idx])]                          # closest approach to the saddle
    seg = idx[idx <= j]
    i0 = seg[np.argmax(pn[seg])]                         # descent into the saddle starts at the |p| peak
    print(f"  closest approach: |p| = {pn[j]:.2e} at x = ({Y[0, j]:.6f}, {Y[1, j]:.6f})")
    return np.array([[pn[i], Y[4, i], Y[5, i], Y[6, i], Y[7, i], Y[8, i]] for i in seg if i >= i0])


WIN = (1e-3, 6e-3)        # |p| window: close enough to be linear, far enough that the Riccati is accurate


def saddle_lyapunov():
    """Z_s = C_s^{-1}, A C + C A^T = -D at the saddle (p = 0): the algebraic Riccati's full-rank root."""
    from scipy.linalg import solve_continuous_lyapunov
    xs1, xs2 = 3.0, R2_2
    f1p = 2 * A1 * xs1 - (3 * xs1 * xs1 + C1)
    hh, dh = model.hill(xs1), model.dhill(xs1)
    f2p = 2 * hh * A2 * xs2 - (3 * xs2 * xs2 + C2)
    A = np.array([[f1p, 0.0], [dh * A2 * xs2 * xs2, f2p]])
    D = np.diag([2 * (A1 * xs1 * xs1 + B1), 2 * (hh * A2 * xs2 * xs2 + B2)])
    C = solve_continuous_lyapunov(A, -D)
    return np.linalg.inv(C), f2p


def report(H, label):
    w = H[(H[:, 0] >= WIN[0]) & (H[:, 0] <= WIN[1])]
    Zs, lam_plus = saddle_lyapunov()
    W = np.polyfit(w[:, 0] ** 2, w[:, 1], 1)[1]
    zf = [np.polyfit(w[:, 0], w[:, j], 1)[1] for j in (2, 3, 4)]
    lK = np.polyfit(w[:, 0], w[:, 5], 1)[1]
    detZm, detZs = np.linalg.det(sh.Z), np.linalg.det(Zs)
    pref = lam_plus / (2 * np.pi) * np.sqrt(detZm / abs(detZs)) * np.exp(lK)
    print(f"{label}: EPS={EPS:g} RTOL={RTOL:g} window samples={len(w)}")
    print(f"  Z_s Lyapunov=({Zs[0,0]:.6f},{Zs[0,1]:.6f},{Zs[1,1]:.6f})  Riccati-window fit=({zf[0]:.6f},{zf[1]:.6f},{zf[2]:.6f})")
    print(f"  W*={W:.10f} Gamma={V1(S1) - W:.6e} lnK_s={lK:.6f} detZm={detZm:.5e} detZs={detZs:.5e} lambda+={lam_plus:.5f}")
    print(f"PREFACTOR {pref:.8e} W {W:.12f} LNK {lK:.8f}")


if __name__ == "__main__" and UNCOUPLED:
    report(approach_samples([0.0, -1.0]), "uncoupled")
    sys.exit()

if __name__ == "__main__":
    ds = -np.geomspace(1e-10, 1e-1, 181)
    hits = [ray_aug(d) is not None for d in ds]
    k = max(i for i in range(180) if hits[i] and not hits[i + 1])
    a, b = ds[k], ds[k + 1]
    for _ in range(70):
        m_ = 0.5 * (a + b)
        if ray_aug(m_) is not None:
            a = m_
        else:
            b = m_
    report(approach_samples([-1.0, a]), "coupled")
