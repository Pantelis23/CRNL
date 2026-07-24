"""AM phase portrait: the landscape, made visible (design.md §2.3, §7).

Draw the flow of the deterministic mass-action ODE on the (X, Y) simplex, mark
the four fixed points with their stability, trace the separatrix through the
saddle, and overlay a handful of exact Gillespie trajectories so the noise sits
next to the smooth flow it fluctuates around. The point: the saddle's role as a
restoring *threshold* shows up as basin structure -- a wall the flow is pushed
away from -- not as a slogan.

    python -m experiments.phase_portrait            # writes experiments/phase_portrait.png
    python -m experiments.phase_portrait --omega 60 --trajectories 8
"""

from __future__ import annotations

import argparse
import os

import numpy as np

from crnl import approximate_majority, integrate, gillespie, seed_for


FIXED_POINTS = {
    "all-X (stable)": (1.0, 0.0),
    "all-Y (stable)": (0.0, 1.0),
    "all-B (repeller)": (0.0, 0.0),
    "saddle": (1 / 3, 1 / 3),
}


def make_figure(omega, n_traj, out_path, seed_base):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    net = approximate_majority()

    fig, ax = plt.subplots(figsize=(7.5, 7.0))

    # -- simplex boundary x + y <= 1, x,y >= 0 --
    ax.plot([0, 1, 0, 0], [0, 0, 1, 0], color="#333333", lw=1.2, zorder=1)

    # -- flow field on a grid inside the simplex --
    g = 22
    xs = np.linspace(0.001, 0.999, g)
    ys = np.linspace(0.001, 0.999, g)
    X, Y = np.meshgrid(xs, ys)
    U = np.zeros_like(X)
    V = np.zeros_like(Y)
    mask = (X + Y) < 1.0
    for i in range(g):
        for j in range(g):
            if not mask[i, j]:
                continue
            x, y = X[i, j], Y[i, j]
            r = net.rhs([x, y, 1 - x - y])
            U[i, j], V[i, j] = r[0], r[1]
    speed = np.hypot(U, V)
    ax.streamplot(X, Y, np.where(mask, U, np.nan), np.where(mask, V, np.nan),
                  color=np.where(mask, speed, np.nan), cmap="Blues",
                  density=1.3, linewidth=0.7, arrowsize=0.8, zorder=2)

    # -- separatrix: integrate BACKWARD from just off the saddle along the
    #    stable eigenvector to trace the basin boundary (X=Y line for AM) --
    # AM's separatrix is exactly the diagonal x = y; draw it as the threshold.
    ax.plot([0, 0.5], [0, 0.5], "--", color="#d62728", lw=2,
            label="separatrix (decision threshold)", zorder=3)

    # -- overlay exact stochastic trajectories from a slightly biased start --
    x0 = np.array([0.52, 0.48, 0.0])
    n0 = np.round(x0 * omega).astype(int)
    n0[2] = int(omega) - n0[0] - n0[1]
    for t in range(n_traj):
        res, ts, ns = gillespie(net, n0, omega, seed_for(omega, t, base=seed_base),
                                record=True)
        conc = ns / omega
        ax.plot(conc[0], conc[1], color="#7f7f7f", lw=0.6, alpha=0.55,
                zorder=4)
    ax.plot([], [], color="#7f7f7f", lw=0.8, alpha=0.8,
            label=f"Gillespie trajectories (Ω={omega}, 52/48 start)")

    # -- deterministic trajectory from the same start (the smooth 'lie') --
    traj = integrate(net, x0, t_span=(0, 60))
    ax.plot(traj.x[0], traj.x[1], color="#2ca02c", lw=2.4,
            label="deterministic ODE (same start)", zorder=6)

    # -- fixed points --
    styles = {
        "all-X (stable)": ("o", "#000000", 90),
        "all-Y (stable)": ("o", "#000000", 90),
        "all-B (repeller)": ("^", "#d62728", 90),
        "saddle": ("X", "#d62728", 120),
    }
    for label, (x, y) in FIXED_POINTS.items():
        m, c, s = styles[label]
        ax.scatter([x], [y], marker=m, c=c, s=s, zorder=7,
                   edgecolors="white", linewidths=1.0)
        ax.annotate(label, (x, y), textcoords="offset points",
                    xytext=(8, 8), fontsize=9, zorder=8)

    ax.set_xlabel("[X]")
    ax.set_ylabel("[Y]")
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.set_aspect("equal")
    ax.set_title("Approximate Majority: logic in, landscape out\n"
                 "flow rolls downhill to a rail; the saddle is the restoring "
                 "threshold")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    print(f"wrote figure -> {out_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--omega", type=int, default=80)
    p.add_argument("--trajectories", type=int, default=10)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default=os.path.join(
        os.path.dirname(__file__), "phase_portrait.png"))
    args = p.parse_args()
    make_figure(args.omega, args.trajectories, args.out, args.seed)


if __name__ == "__main__":
    main()
