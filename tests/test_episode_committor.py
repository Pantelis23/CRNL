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
    e = _e()
    oms = (40, 55, 100, 200, 400, 800)
    pd = [e[f"X:{om}/hi"]["p_d_ep"] for om in oms]
    lead = [e[f"X:{om}/hi"]["p_d_ep"] * e[f"X:{om}/hi"]["phi_ep"] for om in oms]
    assert pd == sorted(pd) and pd[-1] > 0.8
    assert lead == sorted(lead) and 0.35 < lead[-1] < 0.5


def test_narrow_zone_has_not_turned_E2_failed():
    """E2 as written required Z to rise above 0.057 by Omega = 230-400. It did not."""
    e = _e()
    z = [e[f"Z:{om}/hi"]["p_d_ep"] for om in (40, 55, 100, 230, 400)]
    assert z == sorted(z, reverse=True)
    assert z[-1] < 0.057
    assert (z[-2] - z[-1]) / z[-2] < (z[-3] - z[-2]) / z[-3], "but the decline is flattening"
