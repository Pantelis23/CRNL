"""§100 — what the reflecting boundary costs, and the T-CASC-l discriminator."""

import numpy as np
import pytest

from experiments import what_reflection_costs as wrc


@pytest.fixture(scope="module")
def scope():
    return wrc.run_scope()


def test_relax_rate_has_no_omega_in_it():
    """The whole discriminator rests on this being a macroscopic quantity."""
    assert wrc.relax_rate() == pytest.approx(6.6195, abs=1e-3)


def test_p1_reproduces_section_93(scope):
    """Same code path as §93; a mismatch means the harness moved under it."""
    assert scope["pen_refl"] == pytest.approx(4.4419, abs=1e-3)


def test_free_upstream_has_more_error(scope):
    """Reflection can only delete trajectories, so it can only lower stage-2 error."""
    assert scope["p_free"] > scope["p_refl"]
    assert scope["ratio"] == pytest.approx(1.1219, abs=2e-3)


def test_the_wall_inflates_the_surviving_branch(scope):
    """The net -10.9% hides two opposite effects, each far larger than the net."""
    both = scope["product"]
    hi_free = scope["p_free"] - both
    inflation = scope["p_refl"] - hi_free
    assert both > 0 and inflation > 0, "the two effects must have opposite signs"
    assert inflation / hi_free == pytest.approx(0.49, abs=0.02)
    # each component dwarfs the net
    assert both > 3 * abs(scope["p_refl"] - scope["p_free"])
    assert inflation > 2 * abs(scope["p_refl"] - scope["p_free"])


def test_saturation_does_not_protect_against_failure(scope):
    """P3 refuted: a failed upstream transmits almost perfectly through the Hill map."""
    assert scope["p_s2lo_given_s1lo"] == pytest.approx(0.9376, abs=5e-3)
    assert scope["p_s2lo_given_s1lo"] > 0.9


@pytest.mark.parametrize("om", wrc.GAP_OMEGAS)
def test_free_gap_tracks_the_escape_action(om):
    """The free stage-1 gap is an escape rate: same order as exp(-A*Omega)."""
    gap, _ = wrc.spectral_gap(om, False)
    assert 0.8 < gap / np.exp(-wrc.A_UP * om) < 2.0


def test_the_two_timescales_are_separated_by_the_omega_axis():
    """§93's speed knob could not do this; Omega separates them by 3 orders."""
    refl = [wrc.spectral_gap(om, True)[0] for om in wrc.GAP_OMEGAS]
    free = [wrc.spectral_gap(om, False)[0] for om in wrc.GAP_OMEGAS]
    assert max(refl) / min(refl) < 1.3, "a relaxation time must be flat in Omega"
    assert max(free) / min(free) > 100, "an escape time must be exponential in Omega"


def test_reflected_gap_is_not_the_rail_relaxation_rate():
    """P4 refuted: the wall installs a box-scale timescale 4.6x slower than |f'(r3)|."""
    gap, _ = wrc.spectral_gap(30, True)
    assert gap == pytest.approx(1.4327, abs=5e-3)
    assert wrc.relax_rate() / gap == pytest.approx(4.62, abs=0.05)
    assert gap < 0.5 * wrc.relax_rate(), "it does not approach the rail rate"


def test_the_headline_numbers_are_window_dependent():
    """§100.2: both P2's and P3's numbers are one cell of a strong trend."""
    conds, ratios = [], []
    for t in (0.5, 2.0, 8.0):
        s = wrc.run_scope(t0=t)
        conds.append(s["p_s2lo_given_s1lo"])
        ratios.append(s["ratio"])
    # co-occurrence rises monotonically toward 1 -- it is not a transmission probability
    assert conds == sorted(conds)
    assert conds[0] == pytest.approx(0.7254, abs=5e-3)
    assert conds[-1] == pytest.approx(0.9830, abs=5e-3)
    # the wall hides MORE at short windows, not less
    assert ratios == sorted(ratios, reverse=True)
    assert ratios[0] == pytest.approx(1.9615, abs=1e-2)
    assert ratios[0] / ratios[-1] > 1.5, "the effect is not absorbable into a constant"


def test_saturation_gives_no_protection_at_any_window():
    """The claim that survives §100.2: co-occurrence is already high at the shortest window."""
    assert wrc.run_scope(t0=0.5)["p_s2lo_given_s1lo"] > 0.7
