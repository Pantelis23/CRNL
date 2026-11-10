"""§127: W-bit's joint minimum action, by grid, collocation and the closed form for N."""

import json
from pathlib import Path

import numpy as np
import pytest

from experiments import wbit_action as wa
from experiments.wbit_cascade import PAIRS

RES = Path(__file__).resolve().parents[1] / "results" / "wbit_action.json"


def test_escape_action_at_zero_drive_is_the_dip_action():
    assert wa.escape_action(0.0) == pytest.approx(wa.S_DIP, abs=2e-5)


def test_n_closed_form_matches_collocation():
    """W*_N = S + S2(b2): cell 1 parked on its saddle, cell 2 escaping against the drive it sees there."""
    c, b2 = PAIRS["N"]
    _, _, fld = wa.joint_action(81, c, b2, field=True)
    P = wa.grid_path(fld)
    P0 = np.vstack([np.repeat(P[:1], 40, axis=0), P])
    W, feas, _ = wa.path_action_states(c, b2, len(P0) - 1, P0)
    assert feas
    assert W - wa.S_DIP == pytest.approx(wa.escape_action(b2), abs=1e-3)


def test_z_path_is_interior_and_gamma_positive():
    c, b2 = PAIRS["Z"]
    _, _, fld = wa.joint_action(81, c, b2, field=True)
    P = wa.grid_path(fld)
    P0 = np.vstack([np.repeat(P[:1], 20, axis=0), P])
    W, feas, r = wa.path_action_states(c, b2, len(P0) - 1, P0)
    assert feas
    assert W == pytest.approx(0.31646, abs=5e-5)
    assert r.x[0::2].min() > 0.55                       # never near cell 1's saddle
    assert wa.S_DIP - W > 0


def test_grid_is_an_upper_bound_converging_down():
    rep = json.loads(RES.read_text())
    for name, best in (("Z", 0.31646), ("N", 0.43552)):
        ws = [rep[f"{name}:{n}"]["W_star"] for n in (81, 121, 161, 201, 241)]
        assert all(a > b for a, b in zip(ws, ws[1:]))
        assert ws[-1] > best
        assert rep[f"colloc:{name}"][-1]["W_star"] == pytest.approx(best, abs=2e-5)
