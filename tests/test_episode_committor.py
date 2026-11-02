"""§123 -- exact p_d as an episode committor; two instruments broke; the floor climbs toward 1/2."""

import json
import pathlib

STORE = pathlib.Path("results/episode_committor.json")
WINDOW = pathlib.Path("results/weak_upstream.json")


def _e():
    return json.loads(STORE.read_text())


def test_episode_reproduces_the_window():
    """E1: within 6-7% of §119's windowed exact p_d at Omega = 40, 55."""
    w = json.loads(WINDOW.read_text())
    for k in ("X:40", "X:55", "Z:40", "Z:55"):
        r = _e()[k]["p_d_ep"] / w[k]["series"][-1]["p_d"]
        assert 0.97 < r < 1.07, k


def test_episode_upstream_factor_is_the_committor():
    """E3, against §122's independently computed committor."""
    a = json.loads(pathlib.Path("results/asymptotic_floor.json").read_text())["phi"]
    for k in ("X:40", "X:55", "Z:40", "Z:55"):
        assert abs(_e()[k]["phi_ep"] / a[k]["committor"] - 1) < 1e-3, k


def test_well_conditioned_route_matches_where_others_were_valid():
    """E5: within 1.5% of the occupation route wherever that route was reliable."""
    e = _e()
    for k in ("X:40", "X:55", "X:100", "Z:40", "Z:55", "Z:100", "Z:230"):
        assert abs(e[f"{k}/hi"]["p_d_ep"] / e[f"{k}/occ"]["p_d_ep"] - 1) < 0.015, k
    assert all(v["eig_res"] < 1e-8 for k, v in e.items() if k.endswith("/hi"))


def test_the_broken_instruments_are_on_record():
    """The occupation route at X:400 and the QSD eigenvector at Z:800 printed confident wrong
    numbers; both are kept, labelled, and visibly far from the valid route."""
    e = _e()
    assert abs(e["X:400/occ"]["p_d_ep"] / e["X:400/hi"]["p_d_ep"] - 1) > 0.2
    assert e["Z:800"]["J_low_frac"] > 0.1 and e["Z:800"]["lam0"] < 0


def test_exact_wide_zone_floor_climbs_toward_one_half():
    """§123.5 (E9 route): X's p_d rises monotonically to 0.94 at Omega = 800; phi * p_d to 0.443."""
    e = _e()
    oms = (40, 100, 200, 400, 800)
    pd = [e[f"X:{om}/qp"]["p_d_ep"] for om in oms]
    lead = [e[f"X:{om}/qp"]["p_d_ep"] * e[f"X:{om}/qp"]["phi_ep"] for om in oms]
    assert pd == sorted(pd) and pd[-1] > 0.9
    assert lead == sorted(lead) and 0.42 < lead[-1] < 0.5


def test_narrow_zone_has_not_turned_E2_failed():
    """E2 as written required Z to rise above 0.057 by Omega = 230-400. It did not."""
    e = _e()
    z = [e[f"Z:{om}/hi"]["p_d_ep"] for om in (40, 55, 100, 230, 400)]
    assert z == sorted(z, reverse=True)
    assert z[-1] < 0.057
    assert (z[-2] - z[-1]) / z[-2] < (z[-3] - z[-2]) / z[-3], "but the decline is flattening"


def test_narrow_zone_turns_gently_on_the_deterministic_route():
    """E6 / E7 hold on the E9 numbers: monotone 400 -> 600 -> 800 -> 1000, but gently (0.032 ->
    0.0415), not §123.4's withdrawn 0.186 / 0.562."""
    e = _e()
    z = [e[f"Z:{om}/qp"]["p_d_ep"] for om in (400, 600, 800, 1000)]
    assert z == sorted(z) and z[2] >= 0.0306
    assert z[-1] < 0.05


def test_the_withdrawn_route_is_on_record_and_irreproducible():
    """§123.5: the same Z:1000 cell, same code, gave 0.562 and 0.188 on E5's random-start
    eigenvector; the deterministic route gives 0.0415. All three are kept, labelled."""
    e = _e()
    assert "withdrawn" in e["Z:1000/hi"]
    assert e["Z:1000/hi"]["p_d_ep"] / e["Z:1000/hi_rerun"]["p_d_ep"] > 2
    assert e["Z:1000/hi_rerun"]["p_d_ep"] / e["Z:1000/qp"]["p_d_ep"] > 4


def test_deterministic_route_validated_and_order_independent():
    """E9: within 0.6% of E5 at its seed-stable cells; LU orderings agree to 1e-12."""
    e = _e()
    for k in ("X:40", "X:100", "Z:55", "Z:230", "Z:400"):
        assert abs(e[f"{k}/qp"]["p_d_ep"] / e[f"{k}/hi"]["p_d_ep"] - 1) < 0.006, k
        assert e[f"{k}/qp"]["order_diff"] < 1e-12, k


def test_the_narrow_zone_rise_is_mostly_pre_dip_sag():
    """E8 on the E9 dumps: J_1000 with f_600 recovers >= 60% of the 600 -> 1000 rise, J_600 with
    f_1000 <= 40%."""
    import numpy as np
    D = pathlib.Path("results/episode_dumps")
    L = {om: np.load(D / f"Z_{om}_qp.npz") for om in (600, 1000)}

    def pd(a, b):
        A, B = L[a], L[b]
        w = A["J"] * A["c"]
        return float((w * np.interp(A["x2"], B["x2"], B["f"])).sum() / w.sum())

    lo, hi = pd(600, 600), pd(1000, 1000)
    assert (pd(1000, 600) - lo) / (hi - lo) >= 0.6
    assert (pd(600, 1000) - lo) / (hi - lo) <= 0.4


def test_committed_never_dipped_term_matches_the_window():
    """E11(a): within 15% of §119's windowed p0/(k1 t) at five cells, and of §122's at X:100."""
    e = _e()
    w = json.loads(WINDOW.read_text())
    for k in ("X:30", "X:40", "X:55", "Z:40", "Z:55"):
        s = w[k]["series"][-1]
        assert abs(e[f"{k}/qpr"]["p0_rate_ratio"] / (s["p0"] / (w[k]["k1"] * s["t"])) - 1) < 0.15, k
    assert abs(e["X:100/qpr"]["p0_rate_ratio"] / e["X:100/window_p0"]["p0_over_k1t"] - 1) < 0.15


def test_never_dipped_term_grows_like_omega_squared_in_X():
    """E11(b) failed, E12 held by letter, reading corrected: a steady local exponent ~2, so C
    passes 1/2 and reaches ~1.8 at Omega = 800. 'C -> 1/2' withdrawn."""
    import math
    e = _e()
    r = {om: e[f"X:{om}/qpr"]["p0_rate_ratio"] for om in (100, 200, 400, 800)}
    expo = [math.log(r[2 * om] / r[om]) / math.log(2) for om in (100, 200, 400)]
    assert all(1.8 < x < 2.4 for x in expo)
    v = e["X:800/qpr"]
    assert v["phi_ep"] * v["p_d_ep"] + v["p0_rate_ratio"] > 1.5
    assert max(e[f"Z:{om}/qpr"]["p0_rate_ratio"] for om in (230, 400, 800, 1000)) < 1e-5


def test_never_dipped_falls_two_channels():
    """E13: the original-measure reactive flux equals the committed fall rate to 1e-12; the in-zone
    share rises 0.12 -> 0.53 over Omega = 55 -> 400 (E13(b)'s >= 70% failed)."""
    e = _e()["E13"]
    for v in e.values():
        assert abs(v["reactive_ratio"] / v["committed_ratio"] - 1) < 1e-12
    z = [e[str(om)]["frac_in_zone"] for om in (55, 100, 200, 400)]
    assert z == sorted(z) and z[0] < 0.2 and 0.5 < z[-1] < 0.7


def test_to_omega_1600_on_pardiso():
    """E15 / E14 (laptop, MKL PARDISO): the dip channel reaches 0.478 -> 1/2, the never-dipped term
    grows x5.8 from 800 to 1600; the factorisation is verified by phi = the exact 1-D committor."""
    e = _e()
    for om in (1200, 1600):
        v = e[f"X:{om}/qpr"]
        assert v["solver"] == "pardiso"
        assert abs(v["phi_ep"] / v["phi_check_1d_committor"] - 1) < 1e-10
    assert e["X:1600/qpr"]["committor_nnz"] > 2 ** 31, "X:1600's factors exceed any 32-bit index"
    lead = [e[f"X:{om}/qpr"]["phi_ep"] * e[f"X:{om}/qpr"]["p_d_ep"] for om in (800, 1200, 1600)]
    assert lead == sorted(lead) and 0.443 < lead[-1] < 0.5
    r = e["X:1600/qpr"]["p0_rate_ratio"] / e["X:800/qpr"]["p0_rate_ratio"]
    assert 3.5 <= r <= 6.0


def test_channels_hand_off_and_the_exponent_drifts():
    """E16 failed: past Omega ~ 800 both channels grow alike; the in-zone share peaks at 800; the
    total local exponent drifts upward ('~Omega^2' was a window)."""
    import math
    e = _e()["E13"]
    share = {om: e[str(om)]["frac_in_zone"] for om in (400, 800, 1200, 1600)}
    assert share[800] > share[400] and share[800] > share[1600]
    tot = {om: e[str(om)]["committed_ratio"] for om in (200, 400, 800, 1200, 1600)}
    ex = [math.log(tot[b] / tot[a]) / math.log(b / a) for a, b in ((200, 400), (400, 800), (800, 1200), (1200, 1600))]
    assert ex == sorted(ex) and ex[-1] > 2.6
    zone = lambda om: tot[om] * share[om]
    above = lambda om: tot[om] * (1 - share[om])
    ez = math.log(zone(1600) / zone(800)) / math.log(2)
    ea = math.log(above(1600) / above(800)) / math.log(2)
    assert abs(ez - ea) < 0.3 and ea > ez
