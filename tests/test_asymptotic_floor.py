"""§122 -- phi is the committor; p_d tends to a single-element erf law; the finite-Omega goal failed."""

import json
import pathlib

STORE = pathlib.Path("results/asymptotic_floor.json")


def _d():
    return json.loads(STORE.read_text())


def test_phi_is_the_committor():
    """F2: within 1.5% wherever the window is not saturating (all but Z:40, k1 t = 0.43)."""
    phi = _d()["phi"]
    for k in ("X:30", "X:40", "X:55", "Z:55"):
        assert abs(phi[k]["ratio"] - 1) < 0.015, k
    assert abs(phi["Z:40"]["ratio"] - 1) < 0.05


def test_phi_tends_to_one_half():
    for name in ("E-070", "E-062"):
        seq = [v for _, v in sorted(_d()["phi"][f"{name}/limit"].items(), key=lambda kv: int(kv[0]))]
        assert seq == sorted(seq) and 0.48 < seq[-1] < 0.5


def test_u_star_is_horizon_stable_and_positive_for_every_zone():
    for p, v in _d()["u_star"].items():
        assert v["u_star"] > 0 and v["w"] > 0
        assert abs(v["u_star_h200"] / v["u_star"] - 1) < 1e-6


def test_reduction_converges_to_the_erf_law():
    """F4 held: within 15% at Omega = 400 for X and b, and every pair within O(sigma/w) of 1 at
    Omega = 3200 -- not monotone to the end (b undershoots to 0.94)."""
    conv = _d()["conv"]
    at = lambda p, om: next(r for r in conv[p] if r["omega"] == om)["ratio"]
    assert abs(at("X", 400) - 1) < 0.15 and abs(at("b", 400) - 1) < 0.15
    for p in conv:
        r = [row["ratio"] for row in conv[p] if row["omega"] <= 1600]
        assert r == sorted(r, reverse=True), p
        assert abs(at(p, 3200) - 1) < 0.08, p


def test_the_crossover_orders_the_convergence():
    """Omega_c = lam(s1)/(ls w^2) orders how late each pair converges: X < b < Z."""
    u = _d()["u_star"]
    assert u["X"]["omega_c"] < u["b"]["omega_c"] < u["Z"]["omega_c"]
    conv = _d()["conv"]
    at = lambda p, om: next(r for r in conv[p] if r["omega"] == om)["ratio"]
    assert at("X", 100) < at("b", 100) < at("Z", 100)


def test_the_goal_failed_at_finite_omega():
    """G3: the single-element formula is 3x (X) to ~1000x (c) below the exact floor at 30-55, and
    zero where the edge sits on the saddle (d) -- recorded, not hidden."""
    g = _d()["goal"]
    assert all(v["ratio"] > 2.5 for v in g.values())
    assert g["d:40"]["C_formula"] == 0.0


def test_exact_chain_at_omega_100_moves_the_predicted_way():
    """G4 held: exact p_d rises to 0.385 (predicted [0.33, 0.45]), C rises past 0.136, phi equals
    the committor in a cell never used, and exact/formula keeps closing (3.6 -> 2.8 -> 2.6)."""
    e = _d()["exact100"]
    assert 0.33 <= e["p_d"] <= 0.45 and e["C"] > 0.136
    assert abs(e["phi"] / e["committor"] - 1) < 0.01
    assert e["p_d"] / e["formula_p_d"] < 0.3082 / 0.1089
