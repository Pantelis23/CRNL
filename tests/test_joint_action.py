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


def test_large_gamma_pair_converges_onto_the_action():
    """E22: XK (Hill threshold 1.2), Gamma = 7.34e-3 fixed beforehand; the local-exponent slope over
    the two largest windows lies in [0.6, 1.3] Gamma and increases toward it."""
    import math
    e = json.loads((ROOT / "results/episode_committor.json").read_text())
    oms = (100, 200, 400, 600, 800)
    r = {om: e[f"XK:{om}/qpr"]["p0_rate_ratio"] for om in oms}
    ex = [(math.sqrt(a * b), math.log(r[b] / r[a]) / math.log(b / a)) for a, b in zip(oms, oms[1:])]
    sl = [(e2 - e1) / (m2 - m1) for (m1, e1), (m2, e2) in zip(ex, ex[1:])]
    G = J["XK"]["Gamma"]
    assert all(0.6 * G <= v <= 1.3 * G for v in sl[-2:])
    assert sl == sorted(sl)


def test_hill_k_one_is_the_published_coupling():
    """ElemK at K = 1 is the published element: same Hill factor everywhere."""
    from experiments.weak_upstream import element
    a, b = element("E-084"), element("E-084@K1.0")
    for x in (0.5, 1.0, 1.362, 2.0, 3.0):
        assert abs(a.hill(x) - b.hill(x)) < 1e-14


def test_prefactor_is_a_constant_past_the_crossover():
    """E23: XK to Omega = 1600 -- |a| <= 0.05 in every window past u ~ 6 (H0), not 0.4-0.5 (H1/2), and
    ln R - Gamma Omega flat to < 0.15 over 1000 -> 1600."""
    import math
    e = json.loads((ROOT / "results/episode_committor.json").read_text())
    G = J["XK"]["Gamma"]
    oms = (800, 1000, 1200, 1400, 1600)
    r = {om: e[f"XK:{om}/qpr"]["p0_rate_ratio"] for om in oms}
    for a_, b_ in zip(oms, oms[1:]):
        L = (b_ - a_) / math.log(b_ / a_)
        assert abs(math.log(r[b_] / r[a_]) / math.log(b_ / a_) - G * L) <= 0.05
    off = [math.log(r[om]) - G * om for om in (1000, 1600)]
    assert abs(off[1] - off[0]) < 0.15


def test_prefactor_machinery_reproduces_an_exact_rate():
    """K1 (§124.11): with the coupling removed the 2-D Eyring-Kramers formula reproduces stage 2's
    exact 1-D escape rate, approaching 1 like 1/Omega to the method's ~0.2%; rerun the prefactor."""
    P = json.loads((ROOT / "results/prefactor.json").read_text())["K1_uncoupled"]
    r = [P["ratio_EK_over_exact"][k] for k in ("50", "100", "200", "400", "800")]
    assert all(abs(b - 1) < abs(a - 1) for a, b in zip(r, r[1:]))
    assert abs(P["ratio_EK_over_exact"]["1600"] - 1) < 0.003
    out = subprocess.run([sys.executable, "-W", "ignore", "prefactor.py", "1e-6", "1e-12", "40", "--uncoupled"],
                         cwd=ROOT / "experiments/joint_action", capture_output=True, text=True, timeout=900).stdout
    pref = float(next(l for l in out.splitlines() if l.startswith("PREFACTOR")).split()[1])
    assert abs(pref / P["prefactor"] - 1) < 1e-3


def test_exit_is_spread_at_accessible_omega():
    """E24: the committed-fall crossing distribution keeps an O(1) width over Omega = 200-1600
    (ratio 400 -> 1600 within 25%), while its mean drifts toward the saddle."""
    e = json.loads((ROOT / "results/episode_committor.json").read_text())["E24"]
    sd = {om: e[str(om)]["sd_x1"] for om in (200, 400, 800, 1600)}
    assert 0.75 <= sd[1600] / sd[400] <= 1.25
    mu = [e[str(om)]["mean_x1"] for om in (400, 800, 1600)]
    assert mu == sorted(mu)
