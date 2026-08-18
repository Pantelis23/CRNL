"""FINDINGS §86: the quasipotential ridge sits off the deterministic slow manifold."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.off_the_manifold import FLOOR, A_on_curve, ridge_points
from experiments.where_the_deficit_lives import MEASURED, slow_manifold, u_star


class TestTheDisplacementIsResolved:
    def test_ridge_sits_above_the_manifold_and_survives_in_physical_units(self):
        """P1(b): a quantisation artifact is a FIXED number of lattice units. This is not --
        it holds in physical units while growing in lattice units."""
        phys, units = [], []
        for om in (200, 300, 450):
            rows, margin, _ = ridge_points(0.40, om, n=4)
            d = np.mean([r["b_ridge"] - r["b_det"] for r in rows])
            assert d > 0
            phys.append(d)
            units.append(d * om)
        assert units[-1] > units[0] * 1.5          # grows in lattice units
        assert max(phys) / min(phys) < 1.5         # holds in physical units

    def test_the_underflow_floor_is_far_from_every_point_used(self):
        """§81.1's trap. The floor must not reach the region being read."""
        for om in (200, 450):
            _, margin, _ = ridge_points(0.40, om, n=4)
            assert margin > 100.0

    def test_a_slice_at_the_floor_is_dropped_not_ridged(self):
        """§86.1(1): at M=8 a slice had every point at -690.8 but one, and that survivor
        became the 'argmax'. Such slices must be excluded and counted."""
        rows, _, dropped = ridge_points(0.40, 400, M=8.0, n=6)
        assert len(dropped) >= 1
        for r in rows:
            assert abs(r["b_ridge"] - r["b_det"]) < 0.05      # no survivors of that kind


class TestTheCurveNotTheProjection:
    @pytest.mark.parametrize("g", [0.40, 0.44])
    def test_integrating_along_the_ridge_closes_the_deficit(self, g):
        rows, _, _ = ridge_points(g, 600, n=10)
        uu = np.array([0.0] + [r["u"] for r in rows] + [u_star(g)])
        bb = np.array([slow_manifold(0.0, g)] + [r["b_ridge"] for r in rows] + [g / (1 + g)])
        o = np.argsort(uu)
        uu, bb = uu[o], bb[o]
        m = MEASURED[g]
        a_det = A_on_curve(g, lambda u: slow_manifold(u, g))[0]
        a_rid = A_on_curve(g, lambda u: float(np.interp(u, uu, bb)))[0]
        assert abs(1 - a_det / m) > 0.05                   # the manifold misses
        assert abs(1 - a_rid / m) < 0.02                   # the ridge does not

    def test_interp_must_be_anchored_at_both_ends(self):
        """§86.1(3): rule 19's np.interp trap. Unanchored, the whole interval below the first
        traced slice sits at a constant wrong b."""
        g = 0.40
        rows, _, _ = ridge_points(g, 600, n=10)
        uu = np.array([r["u"] for r in rows])
        bb = np.array([r["b_ridge"] for r in rows])
        bare = A_on_curve(g, lambda u: float(np.interp(u, uu, bb)))[0]
        anch_u = np.array([0.0] + list(uu) + [u_star(g)])
        anch_b = np.array([slow_manifold(0.0, g)] + list(bb) + [g / (1 + g)])
        anch = A_on_curve(g, lambda u: float(np.interp(u, anch_u, anch_b)))[0]
        m = MEASURED[g]
        assert abs(1 - anch / m) < abs(1 - bare / m)
