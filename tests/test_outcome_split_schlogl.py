"""§69's verdict rule, tested before the experiment ran (§66's convention)."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.outcome_split_schlogl import consts, decide


def test_verdict_reaches_all_three_branches():
    assert decide([1, 2, 3], [1e-14, 1.1e-14, 0.9e-14], 1e-14, 1e-15)[0] == "a"
    assert decide([1, 2, 3], [1e-3, 1e-2, 1e-1], 1e-3, 1e-12)[0] == "b"
    assert decide([1, 2, 3], [1e-3, 9e-2, 2e-2], 1e-3, 1e-12)[0] == "c"


def test_verdict_does_not_fire_on_solver_noise():
    """A rise entirely inside the noise floor must not read as a deviation."""
    assert decide([1, 2, 3], [1e-9, 2e-9, 3e-9], 1e-9, 1e-8)[0] == "a"


def test_consts_reproduce_the_placed_roots():
    r1, r2, r3 = 0.5, 1.0, 2.5
    k1a, k1r, k2b, k2r = consts(r1, r2, r3)
    roots = np.sort(np.roots([-k1r, k1a, -k2r, k2b]).real)
    assert roots == pytest.approx([r1, r2, r3], abs=1e-12)


def test_factorisation_holds_without_any_symmetry():
    """FINDINGS §69: Phi_o = p_o on a one-species element at 4:1 skew, ln w over tens.

    §60 attributed the outcome-wise identity to the boundaries being exchange images. There
    is one species here, so that account cannot apply -- and the identity holds anyway.
    """
    from experiments.outcome_split_schlogl import cell

    for r3, om in ((1.5, 400), (2.0, 400), (3.0, 200)):
        r = cell(om, 0.5, 1.0, r3, 0.35, 0.80)
        assert r is not None
        assert r["Phi_lo"] + r["Phi_hi"] == pytest.approx(1.0, abs=1e-8)
        assert r["Phi_lo"] == pytest.approx(r["p_lo"], rel=1e-9)
        assert r["Phi_hi"] == pytest.approx(r["p_hi"], rel=1e-9)
    # and the boundary weights really do differ, so the test is not vacuous
    assert abs(cell(200, 0.5, 1.0, 3.0, 0.35, 0.80)["lnw"]) > 3.0
