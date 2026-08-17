"""FINDINGS §82: the escape action against the thermodynamic force.

The verdict rules are tested FIRST, on engineered data, against the branch each must be able
to print -- the convention that has caught a defect in four consecutive sections.
"""
from __future__ import annotations

import numpy as np
import pytest

from experiments.action_is_not_priced import (
    am_affinity, am_gammas, both_ways_verdict, r3_at_affinity, schlogl_A,
    schlogl_affinity, spread_verdict,
)


class TestVerdictRules:
    def test_spread_verdict_prints_the_other_branch_on_a_constant_A(self):
        """Rule 19: name the data that flips it. A constant A on the level set means
        reliability HAS a thermodynamic price, and the rule must say so."""
        assert spread_verdict([0.5, 0.5, 0.5])[1] is False
        assert spread_verdict([0.5, 0.5001, 0.4999])[1] is False
        assert spread_verdict([0.01, 0.5, 2.0])[1] is True

    def test_spread_verdict_refuses_a_family_too_small_to_be_evidence(self):
        """An empty or singleton family gives max/min = 1.0, which reads as 'constant' --
        §71 P2 printed HOLDS off a None. It must return None instead."""
        assert spread_verdict([])[1] is None
        assert spread_verdict([0.4])[1] is None
        assert spread_verdict([0.4, 0.9])[1] is None
        assert spread_verdict([None, np.nan, -1.0, 0.4])[1] is None

    def test_both_ways_verdict_refuses_a_single_point(self):
        assert both_ways_verdict([3.5])[1] is None
        assert both_ways_verdict([3.5, 3.5001])[1] is False
        assert both_ways_verdict([3.5, 4.9])[1] is True


class TestAffinityIsHeldFixed:
    def test_the_schlogl_level_set_really_is_a_level_set(self):
        F0 = schlogl_affinity(0.1, 1.9)
        for r1 in (0.08, 0.15, 0.3, 0.5):
            assert schlogl_affinity(r1, r3_at_affinity(r1, F0)) == pytest.approx(F0, abs=1e-10)

    def test_schlogl_affinity_is_invariant_under_rescaling_the_roots(self):
        """This is WHY the saddle is pinned at r2 = 1: without it, A moves at fixed affinity
        for free, by a change of units."""
        base = schlogl_affinity(0.1, 1.9, 1.0)
        for lam in (2.0, 5.0):
            assert schlogl_affinity(lam * 0.1, lam * 1.9, lam) == pytest.approx(base, abs=1e-12)
        # ... while A does NOT share that invariance, so the rescaling route is real and excluded
        assert schlogl_A(0.2, 3.8, 2.0)[0] / schlogl_A(0.1, 1.9, 1.0)[0] == pytest.approx(2.0, rel=1e-6)

    def test_am_level_set_holds_the_cycle_affinity_exactly(self):
        C = 0.35 ** 3
        vals = [am_affinity(C / g2 ** 2, g2) for g2 in (0.30, 0.35, 0.45, 0.55)]
        assert max(vals) - min(vals) < 1e-12
        assert vals[0] == pytest.approx(-np.log(C), abs=1e-12)

    def test_am_gammas_keeps_exchange_symmetry_when_g3_defaults(self):
        net = am_gammas(0.5, 0.3)
        x = np.array([0.4, 0.4, 0.2])
        S = net.stoichiometry_matrix()
        f = S @ net.fluxes(x)
        assert f[0] == pytest.approx(f[1], abs=1e-12)   # saddle stays on the diagonal


class TestTheMeasurement:
    def test_A_moves_by_orders_of_magnitude_on_the_affinity_level_set(self):
        F0 = schlogl_affinity(0.1, 1.9)
        As = [schlogl_A(r1, r3_at_affinity(r1, F0))[0] for r1 in (0.08, 0.2, 0.5)]
        assert max(As) / min(As) > 100

    def test_equal_delta_nearly_but_not_quite_fixes_A(self):
        """§82.1. The PREDICTION was that Delta would not determine A; it is refuted in its
        strong form and the refutation is the result. On the narrow sweep r1 in 0.05..0.35 --
        a 1.4x move in the saddle's relative position -- A spans 3% and the verdict printed
        'Delta alone reproduces A'. Over the FULL admissible range it spans 1.75x. Both
        numbers are kept: the narrow one is what a too-short sweep says (rule 9)."""
        narrow = [schlogl_A(r1, r1 + 2 * 1.5164)[0] for r1 in (0.05, 0.15, 0.35)]
        assert max(narrow) / min(narrow) < 1.05          # what the short sweep saw

        wide = [schlogl_A(r1, r1 + 2 * 0.9)[0] for r1 in (0.02, 0.5, 0.97)]
        assert max(wide) / min(wide) > 1.5               # what the full range says

    def test_delta_constrains_A_far_more_tightly_than_the_affinity_does(self):
        """The section's actual claim: BOTH fail to determine A, by margins three orders
        of magnitude apart -- 1.75x against 926x."""
        F0 = schlogl_affinity(0.1, 1.9)
        at_fixed_affinity = [schlogl_A(r1, r3_at_affinity(r1, F0))[0]
                             for r1 in (0.08, 0.2, 0.5)]
        at_fixed_delta = [schlogl_A(r1, r1 + 2 * 0.9)[0] for r1 in (0.02, 0.5, 0.97)]
        s_aff = max(at_fixed_affinity) / min(at_fixed_affinity)
        s_del = max(at_fixed_delta) / min(at_fixed_delta)
        assert s_aff > 100 * s_del
