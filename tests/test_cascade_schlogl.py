"""§71's verdict rule and instrument, tested before the experiment ran (§66's convention)."""
from __future__ import annotations

import numpy as np
import pytest

from experiments.cascade_schlogl import d_max, stage_kernel, verdict


def test_verdict_reaches_both_branches():
    assert verdict(3.2, 3.0)[0] == "transfers"
    assert verdict(0.2, 3.0)[0] == "does-not"
    assert verdict(None, 3.0)[0] == "none"


def test_verdict_compares_ratios_not_raw_depths():
    """AM's own measured/predicted is ~3, so agreement with 1 would be the WRONG target."""
    assert verdict(1.0, 3.0)[0] == "does-not"
    assert verdict(3.0, 3.0)[0] == "transfers"


def test_d_max_interpolates_and_returns_none_when_never_crossed():
    assert d_max(np.array([1.0, 0.9, 0.6, 0.4])) == pytest.approx(3.5, abs=0.01)
    assert d_max(np.array([1.0, 0.99, 0.98])) is None


def test_stage_kernel_is_stochastic():
    K, cap = stage_kernel(60, 0.5, 1.0, 1.5, 6.0)
    assert K.shape == (cap + 1, cap + 1)
    assert np.abs(K.sum(axis=1) - 1.0).max() < 1e-12
    assert (K >= -1e-14).all()


def test_fast_and_exact_cascades_agree():
    """FINDINGS §71.2: run_fast is validated against the dense instrument.

    It is 50-90x SLOWER (||Qt|| ~ Omega, so Krylov cost scales with Omega), so it is kept
    for correctness, not speed. Two bugs were caught by exactly this comparison: a globally
    rather than per-source normalised channel, and channel/chemistry applied in the wrong
    order.
    """
    from experiments.cascade_schlogl import run, run_fast, d_max

    a, _ = run(120, 0.1, 1.0, 1.9, 2.0, 25, 0.45 * 0.9)
    b, _ = run_fast(120, 0.1, 1.0, 1.9, 2.0, 25, 0.45 * 0.9)
    assert np.abs(a - b).max() < 1e-4
    assert d_max(a) == pytest.approx(d_max(b), rel=1e-3)
