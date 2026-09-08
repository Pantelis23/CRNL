"""§114 — the rail axis for the position: reproduction, validity, and the two refutations."""

import json
import math
import pathlib

import numpy as np
import pytest

from experiments.two_axes_for_the_position import (
    BRACKET_MIN,
    ELEMENTS,
    PUBLISHED_RAILS,
    R1_FIX,
    R3_FIX,
    SURV_MIN,
    Elem,
    validity,
)

STORE = pathlib.Path("results/two_axes_for_the_position.json")


def _grid():
    return [r for r in json.loads(STORE.read_text())["grid"] if "position" in r]


# ---------------------------------------------------------------- the parameterised chain

def test_every_element_is_neutral_at_its_own_rail():
    """§103.1's fixed-point structure: F(r3) = r3 deterministically, for every element.

    The whole section rests on the constructed elements being real cascade elements, and this
    is the property §91 built the Hill coupling to satisfy.
    """
    for name, r2 in ELEMENTS:
        el = Elem(name, R1_FIX, r2, R3_FIX)
        k1a, k1r, k2b, k2r = el.c
        drift = el.hill(el.r3) * k1a * el.r3 ** 2 - k1r * el.r3 ** 3 + k2b - k2r * el.r3
        assert abs(drift) < 1e-12, f"{name}: drift at its own rail is {drift}"


@pytest.mark.parametrize("om", (14, 30, 55))
def test_parameterised_chain_reproduces_section_109(om):
    """P1. Below the escape-rate solver's own reproducibility, which is what it must be."""
    pub = {r["omega"]: r for r in
           json.loads(pathlib.Path("results/it_was_the_seed.json").read_text())["p3"]}
    e = Elem("published", *PUBLISHED_RAILS, cap_mult=1.25)
    pos, k_eff, km, kg, ka = e.position(om)
    assert k_eff / km == pytest.approx(pub[om]["qsd"], rel=3e-6)
    assert ka / km == pytest.approx(pub[om]["bracket_top"], rel=3e-6)


def test_the_reproduction_check_cannot_be_tighter_than_the_published_number():
    """P1's gate was first written as 1e-6, which §109's own numbers do not reach at Om=70.

    The escape rate there is ~3e-6 out of a generator of norm ~1e4; a dense QR eigensolve
    resolves it to about 1e-6. The gate has to be the instrument's measured noise.
    """
    e = Elem("published", *PUBLISHED_RAILS, cap_mult=1.25)
    mus, _ = e.operating_points(70, 2)
    vals = [e.escape_rate(70, mus[0] * (1.0 + 1e-15 * i)) for i in range(6)]
    floor = (max(vals) - min(vals)) / min(vals)
    assert floor > 1e-7, f"if the solver were this good the 1e-6 gate would have been fine: {floor}"


# ---------------------------------------------------------------- validity

def test_validity_rejects_a_degenerate_bracket():
    """E-021's two limits collapse together and the position reads -304."""
    rows = _grid()
    e21 = [r for r in rows if r["element"] == "E-021"]
    assert e21, "E-021 missing from the stored grid"
    for r in e21:
        assert math.exp(r["bracket"]) < BRACKET_MIN
        assert not validity(r)[0]
        assert abs(r["position"]) > 10.0, "this is what a degenerate denominator prints"


def test_validity_rejects_an_inverted_bracket():
    """k_avg < k_mean for the deepest element: there is no bracket to hold a position in."""
    rows = _grid()
    inv = [r for r in rows if math.exp(r["bracket"]) < 1.0]
    assert len(inv) == 2, [r["element"] for r in inv]
    assert all(r["element"] == "E-084" for r in inv)
    for r in inv:
        assert not validity(r)[0]
        assert "INVERTED" in validity(r)[1]


def test_validity_rejects_cells_where_stage_one_dies():
    rows = _grid()
    dead = [r for r in rows if r["surv"] < SURV_MIN]
    assert dead, "expected some cells where stage 1 does not survive the window"
    for r in dead:
        assert not validity(r)[0]


def test_the_filter_makes_the_conclusion_harder_not_easier():
    """The raw grid would have 'refuted' everything. This pins that the filter removes the
    wild cells rather than the inconvenient ones."""
    rows = _grid()
    valid = [r for r in rows if validity(r)[0]]
    invalid = [r for r in rows if not validity(r)[0]]
    assert len(valid) == 5, [r["element"] + str(r["omega"]) for r in valid]
    # every rejected cell is wilder than every accepted one
    assert max(abs(r["position"]) for r in valid) < min(abs(r["position"]) for r in invalid
                                                        if abs(r["position"]) > 2.0)


# ---------------------------------------------------------------- the two refutations

def test_a_omega_collapses_the_two_axes_in_the_high_occupancy_regime():
    """§114.2: the out-of-sample confirmation, measured as a collapse over ALL valid cells.

    A first version of this asserted every matched high-occupancy pair agreed to 0.01. Two of
    them do not (gaps 0.176 and 0.209, both where the curve is steepest), and the section text
    had quoted only the three that agree to 0.004 -- the flattering-subset error. The claim is
    now the RMS over every valid cell, with the tight pairs as a separate, weaker statement.
    """
    c = json.loads(STORE.read_text())["collapse"]
    assert c["n_high"] >= 4, c
    assert c["rms_high_occupancy"] / c["curve_span"] < 0.12, c
    assert c["rms_low_occupancy"] / c["curve_span"] > 0.5, "the outlier must stay an outlier"
    assert c["rms_low_occupancy"] > 8.0 * c["rms_high_occupancy"]


def test_the_tightest_matched_pairs_agree_to_four_decimals():
    """The weaker, separate statement: where the curve is flat the collapse is much tighter."""
    d = json.loads(STORE.read_text())["p4"]
    tight = [p for p in d["pairs"]
             if min(p["other_a"], p["other_b"]) > 0.8
             and min(p["val_a"], p["val_b"]) > 4.5]
    assert len(tight) == 3, [(p["a"], p["b"]) for p in tight]
    assert max(p["gap"] for p in tight) < 0.005, [(p["a"], p["b"], p["gap"]) for p in tight]


def test_a_omega_is_refuted_as_sufficient_by_a_matched_triple():
    """§114.3: E-084's A*Omega is bracketed by two cells that agree to four decimals."""
    d = json.loads(STORE.read_text())["p4"]
    pool = {p["cell"]: p for p in d["pool"]}
    lo, mid, hi = pool["E-070 Om=30"], pool["E-084 Om=14"], pool["PUB Om=30"]
    assert lo["A_omega"] < mid["A_omega"] < hi["A_omega"], "the outlier must be bracketed"
    assert abs(lo["position"] - hi["position"]) < 0.001, "the two brackets must agree"
    assert mid["position"] / lo["position"] > 3.0, "and the outlier must be far off"
    assert d["refuted"] is True


def test_pi_low_is_refuted_the_same_way():
    """§114.4, rule 15: the distinguishing candidate gets the identical test, not a kinder one."""
    d = json.loads(STORE.read_text())["p4"]
    assert d["pi_worst"] is not None
    w = d["pi_worst"]
    assert abs(w["val_a"] / w["val_b"] - 1) < 0.05, "pi_low really is matched in this pair"
    assert w["gap"] > 0.25 * d["span"], "and the positions really do disagree"
    assert d["pi_refuted"] is True


def test_the_position_moves_at_fixed_omega():
    """§114.1: kills every volume-only candidate, and it is the one thing a second axis
    is guaranteed to settle."""
    d = json.loads(STORE.read_text())["p2"]
    assert d, "no Omega had two valid cells"
    assert max(v["span"] for v in d.values()) > 0.05


# ---------------------------------------------------------------- the framings that did not transfer

def test_the_bracket_inverts_exactly_where_the_rate_stops_being_convex():
    """§114.5: Jensen is the whole content of the frozen/fast bracket."""
    curv = json.loads(STORE.read_text())["p6"]["curvature"]
    inverted = [c for c in curv if c["ratio"] < 1.0]
    upright = [c for c in curv if c["ratio"] >= 1.0]
    assert inverted and upright
    assert max(c["frac_convex"] for c in inverted) < 0.5
    assert min(c["frac_convex"] for c in upright) > 0.9


def test_section_103s_intrinsic_depression_changes_sign():
    """§114.6: d_intr > 0 for shallow wells -- the QSD mean sits ABOVE the rail."""
    d = json.loads(STORE.read_text())["p5"]["d_intr"]
    assert any(v > 0 for _, v in d["E-021"]), d["E-021"]
    assert all(v < 0 for _, v in d["E-070"]), d["E-070"]


def test_the_stored_grid_matches_a_fresh_recomputation_of_one_cell():
    """Guard against the cached grid drifting from the code that made it."""
    rows = _grid()
    r = next(x for x in rows if x["element"] == "E-070" and x["omega"] == 30)
    el = Elem("E-070", R1_FIX, 1.0050, R3_FIX, cap_mult=2.0)
    assert el.position(30)[0] == pytest.approx(r["position"], rel=1e-6)
    assert el.action() == pytest.approx(r["A"], rel=1e-9)


# ---------------------------------------------------------------- §114.9's correction

def test_proximity_to_the_curve_is_not_lying_on_it():
    """§114.9. The three 'tight' pairs are NOT confirmations of A*Omega: the curve predicts
    they should differ by 0.05-0.11 and they differ by 0.000-0.004."""
    d = json.loads(STORE.read_text())
    om = sorted((r["A_omega"], r["position"]) for r in d["p3"]["omega_sweep"])
    cx = np.array([a for a, _ in om]); cy = np.array([p for _, p in om])
    tight = [p for p in d["p4"]["pairs"]
             if min(p["other_a"], p["other_b"]) > 0.8 and min(p["val_a"], p["val_b"]) > 4.5]
    assert len(tight) == 3
    for p in tight:
        predicted = abs(float(np.interp(p["val_a"], cx, cy))
                        - float(np.interp(p["val_b"], cx, cy)))
        assert predicted > 5.0 * max(p["gap"], 1e-4), (p["a"], p["b"], predicted, p["gap"])


def test_no_candidate_makes_the_two_axes_one_curve():
    """§114.9's headline: every slope ratio is well below 1, so the rail axis is shallower."""
    p7 = json.loads(STORE.read_text())["p7"]
    ratios = [c["slope_ratio"] for c in p7]
    assert all(0.4 < r < 0.8 for r in ratios), {c["candidate"]: c["slope_ratio"] for c in p7}
    assert max(ratios) < 0.85, "if any reached ~1 the two axes would be one curve in it"


def test_tau_cross_is_refuted_and_the_other_two_are_not_separated():
    """§114.9: the one candidate the second axis eliminates, and the one it does not."""
    p7 = {c["candidate"]: c for c in json.loads(STORE.read_text())["p7"]}
    a, m, t_ = p7["A*Omega"], p7["(margin/sigma)^2"], p7["tau_cross"]
    assert t_["rms_over_span"] > 2.0 * m["rms_over_span"], "tau_cross must be distinctly worse"
    assert t_["slope_ratio"] < 0.6 < min(a["slope_ratio"], m["slope_ratio"])
    # and A*Omega vs (margin/sigma)^2 are NOT separated -- if this ever fails, T-CASC-ab moved
    assert abs(a["rms_over_span"] - m["rms_over_span"]) < 0.05
    assert m["rms_over_span"] <= a["rms_over_span"], "(margin/sigma)^2 is the marginally better one"


def test_pi_low_axes_do_not_overlap_so_the_slope_test_cannot_run():
    p7 = {c["candidate"]: c for c in json.loads(STORE.read_text())["p7"]}
    pl = p7["pi_low"]
    assert pl["n_overlap"] == 0
    assert math.isnan(pl["rms_over_span"])
    assert pl["rail_range"][0] > pl["omega_range"][1] - 0.01
