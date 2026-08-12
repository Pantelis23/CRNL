"""§66's verdict rule, tested on engineered data BEFORE the experiment ran.

§63 and §64 lost four verdict rules to defects that a five-line test would have caught: a
threshold no data could cross, a statistic whose argmax was pinned to the window edge, a rule
that read a 0.001 wobble as a trend, and a three-way rule with an unreachable branch. All four
were mechanically testable. This file makes that the convention rather than the exception:
every branch of `decide` must be reachable, and the rule must not fire on noise.
"""
from __future__ import annotations

import numpy as np
import pytest

from experiments.outcome_split_tilted import decide


def test_branch_a_no_deviation_is_reachable():
    """Flat |r-1| inside the beta=0 residual must read as 'no deviation'."""
    code, _ = decide([0.0, 0.1, 0.2, 0.3], [0.04, 0.041, 0.039, 0.042],
                     dev0=0.04, omega_scatter=0.01)
    assert code == "a"


def test_branch_b_deviation_is_reachable():
    """Monotone growth clear of the floor must read as a real deviation."""
    code, _ = decide([0.0, 0.1, 0.2, 0.3], [0.05, 0.20, 0.45, 0.80],
                     dev0=0.05, omega_scatter=0.02)
    assert code == "b"


def test_branch_c_inconclusive_is_reachable():
    """Non-monotone, or growth inside the floor, must NOT be called a deviation."""
    assert decide([0.0, 0.1, 0.2, 0.3], [0.05, 0.30, 0.09, 0.28],
                  dev0=0.05, omega_scatter=0.02)[0] == "c"
    assert decide([0.0, 0.1, 0.2, 0.3], [0.05, 0.06, 0.07, 0.08],
                  dev0=0.05, omega_scatter=0.30)[0] == "c"


def test_the_rule_does_not_fire_on_omega_scatter_alone():
    """A rise smaller than the Omega spread is not evidence -- §63's failure mode."""
    code, _ = decide([0.0, 0.2, 0.4], [0.05, 0.09, 0.12],
                     dev0=0.05, omega_scatter=0.20)
    assert code != "b"


def test_the_rule_does_not_fire_on_the_beta_zero_residual_alone():
    """If beta=0 already sits where beta>0 sits, there is no tilt effect."""
    code, _ = decide([0.0, 0.2, 0.4], [0.30, 0.31, 0.32],
                     dev0=0.30, omega_scatter=0.01)
    assert code == "a"


def test_slaved_b_solves_the_slow_manifold_generically():
    """The start point must come from the tilted network, not an am_reversible closed form."""
    from crnl.networks.am_asymmetric import am_asymmetric
    from experiments.outcome_split_tilted import slaved_b

    for beta in (0.0, 0.2):
        net = am_asymmetric(0.25, beta)
        b = slaved_b(net, -0.05)
        assert b is not None and 0 < b < 1
        S = net.stoichiometry_matrix()
        x = (1.0 - b - 0.05) / 2.0
        y = (1.0 - b + 0.05) / 2.0
        assert float((S @ net.fluxes(np.array([x, y, b])))[2]) == pytest.approx(0, abs=1e-10)
