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
    "Committor.h_zero", "Committor.h_top", "Committor.backward",
    "EdgeInterval.F_pos_of_nonpos", "EdgeInterval.disc_eq", "EdgeInterval.bistable_levels_convex",
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


def test_committor_code_is_the_proved_harmonic_function():
    """Committor.backward / h_zero / h_top: the code's closed form must be 0 at the low rail, 1 at
    the high rail, and solve lam (h+ - h) + mu (h- - h) = 0 in between -- checked on the code's own
    `asymptotic_floor.committor`, state by state."""
    import numpy as np

    from experiments.asymptotic_floor import committor
    from experiments.weak_upstream import element

    for name in ("E-070", "E-062"):
        A = element(name)
        for om in (30, 55):
            lo, hi = int(round(A.r1 * om)), int(round(A.r3 * om))
            h = np.array([committor(A, om, n) for n in range(lo, hi + 1)])
            assert h[0] == 0.0 and abs(h[-1] - 1.0) < 1e-12
            n = np.arange(lo + 1, hi, dtype=float)
            lam, mu = A.rates(n, np.zeros(len(n)), om, True)
            res = lam * (h[2:] - h[1:-1]) + mu * (h[:-2] - h[1:-1])
            assert np.abs(res).max() < 1e-10 * lam.max(), (name, om)


def test_edge_is_the_left_end_of_the_proved_interval():
    """EdgeInterval: three positive roots <=> discriminant > 0 on the code's own elements (the
    classical half the proof does not cover, checked numerically), the bistable inputs are one
    contiguous run, and weak_upstream.edge() is its left end."""
    import numpy as np

    from experiments.weak_upstream import edge, element

    for name in ("E-021", "E-045", "E-050", "E-062", "E-070", "E-084", "F-025", "F-040", "F-050"):
        B = element(name)
        k1a, k1r, k2b, k2r = B.c
        xs = np.linspace(0.2, 3.2, 3001)
        three, dpos = [], []
        for x in xs:
            h = B.hill(x)
            r = np.roots([-k1r, h * k1a, -k2r, k2b])
            real = [z.real for z in r if abs(z.imag) < 1e-9]
            three.append(len(real) == 3 and all(z > 0 for z in real))
            D = (-4 * k1a ** 3 * k2b * h ** 3 + k1a ** 2 * k2r ** 2 * h ** 2
                 + 18 * k1r * k1a * k2r * k2b * h - (4 * k1r * k2r ** 3 + 27 * k1r ** 2 * k2b ** 2))
            dpos.append(D > 0)
        assert three == dpos, name
        idx = np.nonzero(three)[0]
        assert len(idx) > 0 and np.all(np.diff(idx) == 1), f"{name}: bistable inputs not contiguous"
        assert abs(xs[idx[0]] - edge(B)) < 2e-3, name
