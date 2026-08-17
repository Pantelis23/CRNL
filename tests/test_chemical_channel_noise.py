"""FINDINGS §75: the intrinsic rail noise, and the exact molecules-at-the-rail collapse."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.chemical_channel_noise import am_rail_width, rail_width

B = (0.5, 1.0, 1.5)


def _ratio(lam, om):
    r = rail_width(om, lam * B[0], lam * B[1], lam * B[2])
    assert r is not None
    return (lam * (B[2] - B[0]) / 2.0) / r["sd_exact"]


def test_delta_over_sigma_depends_only_on_the_molecule_count():
    """§75.1: the headline. Scaling rails and scaling Omega are the SAME act."""
    for prod, cells in ((800, [(0.5, 1600), (2.0, 400), (1.0, 800)]),
                        (3200, [(0.5, 6400), (2.0, 1600), (8.0, 400), (1.0, 3200)])):
        vals = [_ratio(lam, om) for lam, om in cells]
        assert np.ptp(vals) / np.mean(vals) < 1e-9, (prod, vals)
    # and it genuinely grows with the count, so the test is not vacuous
    assert _ratio(1.0, 3200) > 2 * _ratio(1.0, 200)


def test_rail_noise_is_poissonian_in_lambda_and_omega():
    """§75: sigma_x ~ sqrt(lambda/Omega), i.e. sigma_n ~ sqrt(n̄). NOT sigma ~ lambda."""
    lams = [0.5, 1.0, 2.0, 4.0, 8.0]
    sds = [rail_width(6400, l * B[0], l * B[1], l * B[2])["sd_exact"] for l in lams]
    e_lam = float(np.polyfit(np.log(lams), np.log(sds), 1)[0])
    assert 0.40 < e_lam < 0.60, e_lam            # 1/2, and emphatically not 1

    oms = [400, 1600, 6400]
    sdo = [rail_width(o, *B)["sd_exact"] for o in oms]
    e_om = float(np.polyfit(np.log(oms), np.log(sdo), 1)[0])
    assert -0.60 < e_om < -0.40, e_om


def test_am_rail_noise_is_also_poissonian_so_74s_ceiling_dissolves():
    """§75.2: delta*/sigma ~ sqrt(Omega), so the conservative ceiling was a convention."""
    oms = [60, 120, 240]
    sds = [am_rail_width(0.05, o)[0] for o in oms]
    e = float(np.polyfit(np.log(oms), np.log(sds), 1)[0])
    assert e == pytest.approx(-0.5, abs=0.05), e
    assert sds[0] > sds[-1]
