"""Tests for the cascade experiment (experiments/cascade.py)."""

from __future__ import annotations

import numpy as np

from crnl.networks import approximate_majority
from crnl.vectorized import compile_network
from crnl import seed_for
from experiments.cascade import _am_restore, run_cascade


def test_am_restore_returns_full_rail():
    net = approximate_majority()
    comp = compile_network(net, 50)
    names = list(net.species)
    # a strong bias must restore to the corresponding full-magnitude rail
    out = [_am_restore(0.6, 50, comp, names, seed_for(50, t)) for t in range(20)]
    assert set(out) <= {1.0, -1.0}          # always a clean +-1
    assert sum(o == 1.0 for o in out) >= 18  # strong +bias -> almost always +1


def test_restoration_beats_passthrough_over_depth():
    # the founding claim: restoration keeps the bit alive deeper than analog drift.
    depth, sigma, trials = 20, 0.35, 800
    non = run_cascade("nonrestoring", depth, sigma, 80, trials, 0)
    res = run_cascade("restoring", depth, sigma, 80, trials, 0)
    # both start comparable (stage 1), but by the end restoration is well ahead
    assert res[-1] > non[-1] + 0.1
    # analog passthrough decays toward the coin flip; restoration stays high
    assert non[-1] < 0.62
    assert res[-1] > 0.72


def test_restoration_is_monotone_in_population():
    # bigger Omega -> stronger restoration -> flatter survival curve
    depth, sigma, trials = 20, 0.35, 800
    small = run_cascade("restoring", depth, sigma, 20, trials, 0)[-1]
    big = run_cascade("restoring", depth, sigma, 80, trials, 0)[-1]
    assert big >= small - 0.02              # larger population never worse
