"""The Lean proofs in formal/ build, rest on no unproved step, and the code implements what they prove.

The proofs say nothing about the Python unless something ties the two together. These tests do it
in both directions: the Lean library must build with only Lean's standard axioms (no `sorryAx`),
and each Python helper that a theorem is about must agree with the proved statement on real inputs.
"""

import math
import pathlib
import shutil
import subprocess

import pytest

from experiments.cascade_schlogl import schlogl_consts
from experiments.where_pooling_stops_winning import ELEMENTS, Element

FORMAL = pathlib.Path(__file__).resolve().parent.parent / "formal"
LAKE = shutil.which("lake") or str(pathlib.Path.home() / ".elan" / "bin" / "lake")
THEOREMS = (
    "MergeLaw.below_saddle_iff_threshold_le", "MergeLaw.threshold_eq_k_iff",
    "MergeLaw.threshold_eq_one_iff", "Schlogl.stability", "Schlogl.drift_eq",
    "Schlogl.consts_pos", "EffectiveExponent.ν_eq", "EffectiveExponent.ν_close",
    "EffectiveExponent.ν_tenth", "EffectiveExponent.A_half_sub",
    "BirthDeath.detailed_balance", "BirthDeath.global_balance",
)


@pytest.mark.skipif(not pathlib.Path(LAKE).exists(), reason="Lean toolchain (elan) not installed")
def test_lean_library_builds_on_standard_axioms_only():
    run = dict(cwd=FORMAL, capture_output=True, text=True, timeout=1800)
    b = subprocess.run([LAKE, "build"], **run)
    assert b.returncode == 0, b.stdout[-3000:] + b.stderr[-3000:]
    a = subprocess.run([LAKE, "env", "lean", "scripts/axioms.lean"], **run)
    assert a.returncode == 0, a.stderr[-3000:]
    lines = [ln for ln in a.stdout.splitlines() if "depends on axioms" in ln]
    for name in THEOREMS:
        line = next((ln for ln in lines if f"'CRNL.{name}'" in ln), None)
        assert line is not None, f"{name} not audited"
        assert "sorryAx" not in line, line
        assert line.endswith("[propext, Classical.choice, Quot.sound]"), line


def test_code_merge_threshold_is_the_proved_one():
    """MergeLaw.below_saddle_iff_threshold_le: the least j with x_merged(k, j) < r2 is
    floor(f k) + 1. Brute-force the least j from the code's own x_merged and compare with m_of_k."""
    for name, r2 in ELEMENTS:
        el = Element(name, 0.15, r2, 3.0)
        for k in range(1, 40):
            least = next(j for j in range(k + 1) if el.x_merged(k, j) < el.r2)
            assert el.m_of_k(k) == least, (name, k)


def test_regime_boundaries_match_the_table_conditions():
    """threshold_eq_k_iff / threshold_eq_one_iff, read on the six elements over k = 1..39."""
    for name, r2 in ELEMENTS:
        el = Element(name, 0.15, r2, 3.0)
        for k in range(1, 40):
            assert (el.m_of_k(k) == k) == (el.f >= 1 - 1 / k), (name, k)
            assert (el.m_of_k(k) == 1) == (el.f < 1 / k), (name, k)


def test_code_constants_are_the_vieta_construction():
    """Schlogl.drift_eq / drift'_rails: the drift vanishes at the rails, and the code's tau is
    1 / (k (r3 - r1)(r3 - r2))."""
    for name, r2 in ELEMENTS:
        r1, r3 = 0.15, 3.0
        k1a, k1r, k2b, k2r = schlogl_consts(r1, r2, r3)
        for r in (r1, r2, r3):
            assert abs(k1a * r * r + k2b - k1r * r ** 3 - k2r * r) < 1e-12
        el = Element(name, r1, r2, r3)
        assert el.tau == pytest.approx(1.0 / ((r3 - r1) * (r3 - r2)), rel=1e-12)


def test_effective_exponent_printed_values_and_bound():
    """EffectiveExponent.ν_eq at the three ε §85.1 printed (1.8651 / 1.9254 / 1.9837), and the
    proved bound |ν − (2 − 5ε/3)| ≤ 7ε² on (0, 1/4], checked on a grid."""
    nu = lambda e: 2 - e / (1 - e) - e / (0.5 + e) + 2 * e / (1.5 - e)
    for e, printed in ((0.1, 1.8651), (0.05, 1.9254), (0.01, 1.9837)):
        assert abs(nu(e) - printed) < 6e-5
    assert nu(0.1) == pytest.approx(235 / 126, rel=1e-14)
    for i in range(1, 251):
        e = i / 1000
        assert abs(nu(e) - (2 - 5 * e / 3)) <= 7 * e * e
    # and the closed form for A itself: A(1/2 - e) = 4e^2(1-e)/((1/2+e)(3/2-e)^2)
    for e in (0.01, 0.1, 0.3):
        g = 0.5 - e
        A = (1 + 2 * g) * (1 - 2 * g) ** 2 / (2 * (1 - g) * (1 + g) ** 2)
        assert A == pytest.approx(4 * e * e * (1 - e) / ((0.5 + e) * (1.5 - e) ** 2), rel=1e-12)
        assert math.isfinite(A)


def test_stage1_stationary_is_the_proved_product_formula():
    """BirthDeath.global_balance: the product formula is stationary for the reflecting box. The
    code's `stage1_stationary` must therefore zero pi Q on the chain it says it describes -- stage 1
    restricted to [saddle, cap], no death out of the saddle, no birth out of the cap."""
    import numpy as np

    from experiments.two_axes_for_the_position import PUBLISHED_RAILS, Elem

    for el in (Elem("pub", *PUBLISHED_RAILS, cap_mult=1.25), Elem("E-062", 0.15, 1.233, 3.0)):
        for om in (14, 30):
            up, pi = el.stage1_stationary(om)
            n = up.astype(float)
            lam, mu = el.rates(n, np.zeros(len(n)), om, True)
            lam[-1], mu[0] = 0.0, 0.0
            Q = np.diag(lam[:-1], 1) + np.diag(mu[1:], -1)
            Q -= np.diag(Q.sum(axis=1))
            flux = np.abs(pi @ Q).max() / (pi * (lam + mu)).max()
            assert flux < 1e-12, (el.name, om, flux)
