"""§106 — where §103's closure degrades, and why."""

import json
import pathlib

import numpy as np
import pytest

from experiments.depth_compounding import R3
from experiments.escape_accounts_for_it import escape_rate
from experiments.what_reflection_costs import spectral_gap

RESULTS = pathlib.Path("results/where_the_expansion_frays.json")


@pytest.fixture(scope="module")
def data():
    return json.loads(RESULTS.read_text())


def test_operating_points_converge_super_algebraically(data):
    """P2 refuted: not 1/Omega. The local exponent steepens without bound."""
    rows = data["cells"]
    oms = np.array([r["omega"] for r in rows], float)
    y = np.array([abs(r["rel"][0]) for r in rows])
    loc = [np.log(y[i + 1] / y[i]) / np.log(oms[i + 1] / oms[i]) for i in range(len(y) - 1)]
    assert all(l < -2.0 for l in loc), "every local exponent must beat 1/Omega"
    assert loc == sorted(loc, reverse=True), "and must steepen monotonically"
    assert loc[-1] < -8.0


def test_residuals_are_one_signed(data):
    """P3 holds: the model always predicts a less degraded chain."""
    assert all(v > 0 for r in data["cells"] for v in r["rel"])


def test_closure_ratio_diverges(data):
    """P4 refuted, and this is §106's finding."""
    rows = data["cells"]
    seq = [r["pred_over_meas"] for r in rows]
    assert seq == sorted(seq), "must move monotonically away from 1"
    assert abs(seq[-1] - 1) > 4 * abs(seq[0] - 1)


def test_the_convicted_term_crosses_one_where_it_was_validated(data):
    """The one-way escape model spans 1.71x and is accidentally right near Omega=30."""
    esc = data["escape_model"]
    r = {e["omega"]: e["model_over_meas"] for e in esc}
    assert max(r.values()) / min(r.values()) > 1.6
    assert r[30] == pytest.approx(0.9746, abs=5e-3), "accidentally accurate where §102 checked"
    assert r[14] < 1.0 < r[70], "it crosses 1 inside the swept range"


def test_pure_is_bounded_while_contam_diverges(data):
    """§106.1's localisation: the contaminated channel is the divergent one."""
    dec = data["decomposition"]
    contam = [d["contam_meas_over_model"] for d in dec]
    pure = [d["pure_meas_over_model"] for d in dec]
    assert contam == sorted(contam, reverse=True), "contam falls steadily"
    assert contam[-1] < 0.7
    assert max(pure) / min(pure) < 1.3, "pure stays bounded"
    assert pure != sorted(pure), "and is non-monotone"


@pytest.mark.parametrize("om", (14, 20, 24, 30, 40, 55, 70))
def test_stage_one_rate_is_the_gap_not_its_own_operating_point(om):
    """§106.3: the corrected indexing reproduces the free gap exactly; as-coded does not."""
    gap, _ = spectral_gap(om, False)
    assert escape_rate(om, R3) == pytest.approx(gap, rel=1e-12)


def test_the_as_coded_indexing_is_measurably_off(data):
    """It is a real error, not a notational one -- 12.7% at Omega=14."""
    rows = {r["omega"]: r for r in data["cells"]}
    gap14, _ = spectral_gap(14, False)
    as_coded = escape_rate(14, rows[14]["meas_mus"][0])
    assert as_coded / gap14 > 1.10
