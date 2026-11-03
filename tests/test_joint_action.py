"""§124.7 -- the least joint action of a never-dipped fall, and its agreement with the exact CME."""

import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
J = json.loads((ROOT / "results/joint_action.json").read_text())


def test_free_x2_bound_reproduces():
    """The free-x2 ray for X: x_dag ~ 1.085, Gamma bound ~ 1.1e-3 (coarse Sigma grid, 2%)."""
    out = subprocess.run([sys.executable, "-W", "ignore", "p20.py", "1e-6", "1e-12", "40"],
                         cwd=ROOT / "experiments/joint_action", capture_output=True, text=True,
                         timeout=600).stdout
    line = next(l for l in out.splitlines() if l.startswith("p2=0 ray"))
    g = float(line.split("Gamma_bound=")[1].split()[0])
    assert abs(g / 1.105e-3 - 1) < 0.02, line


def test_gamma_ordering_and_bounds():
    """Converged Gamma exceeds its free-x2 bound in every pair, and falls as the edge sinks."""
    for k in ("X", "b", "Z"):
        assert J[k]["Gamma"] > J[k]["free_x2_bound"]["Gamma"]
    assert J["X"]["Gamma"] > J["b"]["Gamma"] > J["Z"]["Gamma"] > J["a"]["Gamma_estimate"]


def test_stationary_law_extrapolates_to_the_action():
    """X's exact stationary-law Gamma_hat on Sigma, extrapolated with Gamma + c ln(Om)/Om from
    800/1600, lies within 15% of the minimum-action Gamma (the c/Om ansatz is reported too)."""
    import math
    e = json.loads((ROOT / "results/episode_committor.json").read_text())["E13"]
    g = {om: e[str(om)]["gamma_hat_sigma"] for om in (800, 1600)}
    L = {om: math.log(om) / om for om in g}
    c = (g[800] - g[1600]) / (L[800] - L[1600])
    gam = g[1600] - c * L[1600]
    assert abs(gam / J["X"]["Gamma"] - 1) < 0.15


def test_b_versus_X_discriminates():
    """E19(a) / E20: b's never-dipped local exponent falls at large Omega while X's rises."""
    import math
    e = json.loads((ROOT / "results/episode_committor.json").read_text())

    def slope(arm):
        r = {om: e[f"{arm}:{om}/qpr"]["p0_rate_ratio"] for om in (800, 1200, 1600)}
        e1 = math.log(r[1200] / r[800]) / math.log(1.5)
        e2 = math.log(r[1600] / r[1200]) / math.log(4 / 3)
        return (e2 - e1) / (math.sqrt(1200 * 1600) - math.sqrt(800 * 1200))

    assert slope("X") > 0 > slope("b")
    assert slope("b") <= 0.3 * slope("X")


def test_no_joint_path_means_decay():
    """E21: pairs d and N (no x_dag) -- the never-dipped term decays exponentially; verified across
    three solvers at 1e-21 and 1e-45."""
    e = json.loads((ROOT / "results/episode_committor.json").read_text())
    for arm in "dN":
        r = [e[f"{arm}:{om}/qpr"]["p0_rate_ratio"] for om in (200, 400, 800, 1200, 1600)]
        assert all(b < a for a, b in zip(r, r[1:])), arm
    for v in e["E21_crosscheck"].values():
        assert v["rel"] < 1e-10 and v["rel_vs_laptop"] < 1e-10
