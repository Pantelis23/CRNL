"""§101 — the free-upstream channel with depth."""

import numpy as np
import pytest

from experiments import free_upstream_depth as fud
from experiments.depth_compounding import build_chain
from experiments.depth_compounding import seed as seed94
from experiments.margin_law import stage1_stationary
from experiments.what_reflection_costs import build_free


@pytest.mark.parametrize("om", (14, 30))
@pytest.mark.parametrize("D", (2, 3))
def test_generalised_builder_reproduces_section_94_exactly(om, D):
    """Same vectorised propensities, so this must be bit-exact -- generator and seed."""
    _, pi1 = stage1_stationary(om)
    for nr, all_ref in ((D - 1, False), (D, True)):
        Qg, ref, dims, strides, _, walled = fud.build_gen(om, D, nr)
        Q9, ref9, dims9, str9, _ = build_chain(om, D, all_reflected=all_ref)
        assert Qg.shape == Q9.shape
        assert abs(Qg - Q9).max() == 0.0
        sg = fud.seed_gen(om, ref, dims, strides, walled, pi1)
        s9 = seed94(om, ref9, dims9, str9, pi1, all_reflected=all_ref)
        assert np.abs(sg - s9).max() == 0.0


@pytest.mark.parametrize("om", (14, 30))
def test_generalised_builder_matches_section_100_to_float_tolerance(om):
    """§100 uses a scalar loop, this uses a vectorised expression -- agreement is not bitwise.

    §101's P1 asked for max|dQ| == 0 here and printed FAILS off a 1.8e-12 residual. That
    criterion cannot be met by two different formulations of the same arithmetic, and rule 19
    is about exactly that: a criterion must be satisfiable by the thing it claims to test.
    """
    Qg, *_ = fud.build_gen(om, 2, 0)
    Qf, _, _ = build_free(om)
    assert Qg.shape == Qf.shape
    d = abs(Qg - Qf).max()
    assert 0.0 <= d < 1e-9
    assert abs(Qg).max() > 1.0                      # the tolerance is not trivially satisfied


def test_n_reflected_is_range_checked():
    with pytest.raises(ValueError):
        fud.build_gen(14, 2, 3)
    with pytest.raises(ValueError):
        fud.build_gen(14, 2, -1)


def test_channel_split_is_exhaustive():
    """pure + contaminated must equal P(last low) -- no probability may go missing."""
    p, ref, dims, strides, walled = fud.solve(14, 2, 0, 1.0)
    tot, pure, contam = fud.channel_split(p, 14, dims, strides, ref, walled)
    assert pure + contam == pytest.approx(tot, rel=1e-12)
    assert tot == pytest.approx(fud.last_low(p, 14, dims, strides, walled, ref), rel=1e-12)
    assert 0.0 < pure < tot and 0.0 < contam < tot


def test_reflected_arm_reproduces_section_100_at_D2():
    """The D=2 free/reflected ratio must match §100's independently computed value."""
    pf, ref, dims, strides, walled = fud.solve(30, 2, 0, 2.0)
    pr, refr, dimr, strr, walr = fud.solve(30, 2, 1, 2.0)
    ratio = (fud.last_low(pf, 30, dims, strides, walled, ref)
             / fud.last_low(pr, 30, dimr, strr, walr, refr))
    assert ratio == pytest.approx(1.1219, abs=2e-3)


def test_the_wall_understates_total_and_overstates_the_fluctuation_channel():
    """The two errors have opposite signs; the modest net is a cancellation, not soundness."""
    om, D, t = 14, 3, 2.0
    pf, ref, dims, strides, walled = fud.solve(om, D, 0, t)
    tot, pure, contam = fud.channel_split(pf, om, dims, strides, ref, walled)
    pr, refr, dimr, strr, walr = fud.solve(om, D, D - 1, t)
    refl = fud.last_low(pr, om, dimr, strr, walr, refr)
    assert refl < tot, "the wall must understate the total (it deletes a channel)"
    assert refl > pure, "the wall must overstate the fluctuation channel (it reflects mass back)"
    assert refl / pure > 2.0                        # 2.64 measured
    assert tot / refl == pytest.approx(1.2121, abs=5e-3)


def test_equilibration_gate_fires_at_short_windows():
    """§101.1: below t~2 a DEEPER reflected chain reports LESS error, which is impossible.

    Run at Omega=14; the same inversion holds at Omega=30 (1.5353e-4 against 1.5149e-3) but a
    D=3 solve there costs ~20 minutes, so that cell is covered by the stored-run test below.
    """
    def refl_err(om, D, t):
        p, ref, dims, strides, walled = fud.solve(om, D, D - 1, t)
        return fud.last_low(p, om, dims, strides, walled, ref)

    assert refl_err(14, 3, 0.5) < refl_err(14, 2, 0.5), "the gate must exclude t=0.5"
    assert refl_err(14, 3, 2.0) > refl_err(14, 2, 2.0), "the gate must admit t=2.0"


def test_last_stage_is_cleaner_than_its_input_at_short_windows():
    """The mechanism behind the gate: later stages seed as a delta and have not spread."""
    from experiments.depth_compounding import stage_stats
    p, ref, dims, strides, _ = fud.solve(14, 3, 2, 0.5)
    m2, s2 = stage_stats(p, 14, ref, dims, strides, 1, all_reflected=False)
    m3, s3 = stage_stats(p, 14, ref, dims, strides, 2, all_reflected=False)
    assert m3 > m2 and s3 < s2, "at t=0.5 stage 3 must read cleaner than stage 2"
    p, ref, dims, strides, _ = fud.solve(14, 3, 2, 4.0)
    m2b, _ = stage_stats(p, 14, ref, dims, strides, 1, all_reflected=False)
    assert m2b < m2, "stage 2 must degrade as the seeded delta spreads"


def test_one_added_stage_moves_the_majority():
    """§101's P3, computed live at the cheap barrier: minority at D=2, majority at D=3."""
    shares = {}
    for D in (2, 3):
        p, ref, dims, strides, walled = fud.solve(14, D, 0, 2.0)
        tot, pure, contam = fud.channel_split(p, 14, dims, strides, ref, walled)
        assert tot <= 0.5, "saturated cells are excluded from this reading"
        shares[D] = contam / tot
    assert shares[2] < 0.5 < shares[3], f"no crossing: {shares}"
    assert shares[2] == pytest.approx(0.473, abs=5e-3)
    assert shares[3] == pytest.approx(0.688, abs=5e-3)


def test_stored_run_shows_the_crossing_in_every_admitted_cell():
    """The Omega=30 cells cost hours to recompute; guard the recorded values instead."""
    import json
    import pathlib
    rows = json.loads(pathlib.Path("results/free_upstream_depth.json").read_text())
    idx = {(r["omega"], r["D"], r["t"]): r for r in rows}
    admitted = ((30, 2.0), (30, 8.0), (14, 2.0))
    for om, t in admitted:
        d2, d3 = idx[(om, 2, t)], idx[(om, 3, t)]
        assert d3["refl"] > d2["refl"], f"equilibration gate: Om={om} t={t}"
        assert max(d2["free"], d3["free"]) <= 0.5, f"saturated: Om={om} t={t}"
        assert d2["contam_share"] < 0.5 < d3["contam_share"], f"no crossing at Om={om} t={t}"
        assert d3["ratio"] > d2["ratio"], "the ratio must grow with depth"
        # the wall overstates the fluctuation channel and understates the total
        assert d3["refl"] / d3["pure"] > 2.0
        assert d3["refl"] < d3["free"]
    # sub-additive: the D=3 excess is well under twice the D=2 excess
    for om, t in admitted:
        d2, d3 = idx[(om, 2, t)], idx[(om, 3, t)]
        frac = (d3["ratio"] - 1) / (2 * (d2["ratio"] - 1))
        assert 0.55 < frac < 0.75, f"{om} {t}: {frac}"
