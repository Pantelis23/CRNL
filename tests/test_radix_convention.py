"""Which convention isolates the effect of alphabet size (FINDINGS open Q2)."""

import numpy as np
import pytest

from experiments.radix_convention import (
    counts_margin, counts_share, p_champion_wins, pairwise_lead,
)


# -- the two conventions differ in what they hold fixed ---------------------

@pytest.mark.parametrize("n", [2, 3, 4, 8, 16, 24])
def test_margin_convention_holds_the_pairwise_lead_fixed(n):
    omega = 120
    lead = pairwise_lead(counts_margin(n, omega, 0.10), n, omega)
    assert lead == pytest.approx(0.10, abs=0.02), (n, lead)


def test_share_convention_lets_the_pairwise_lead_grow():
    """The reason the two conventions disagree. THEORIES T3 predicted the lead
    would SHRINK under fixed share; it is s - (1-s)/(n-1), which grows toward s
    (measured 0.100 -> 0.525 over n=2..24, a 5.2x increase)."""
    omega = 120
    leads = [pairwise_lead(counts_share(n, omega, 0.55), n, omega)
             for n in (2, 3, 4, 8, 16, 24)]
    assert leads == sorted(leads)
    assert leads[-1] / leads[0] > 4.0, leads


@pytest.mark.parametrize("n", [3, 8, 24])
def test_share_convention_holds_the_share_fixed(n):
    omega = 120
    counts = counts_share(n, omega, 0.55)
    assert counts[0] / omega == pytest.approx(0.55, abs=0.02)


def test_margin_convention_champion_share_falls_toward_delta():
    """Under fixed margin the champion's SHARE is what shrinks -- toward delta.
    That is the quantity §3's penalty is really about."""
    omega = 120
    shares = [counts_margin(n, omega, 0.10)[0] / omega
              for n in (2, 3, 4, 8, 16, 24)]
    assert shares == sorted(shares, reverse=True)
    assert shares[-1] == pytest.approx(0.10, abs=0.06)


# -- both conventions start from the same place -----------------------------

def test_the_two_conventions_coincide_at_n_equals_two():
    """Anchored to 0.55/0.45, so any divergence at n>2 is the convention and
    not a different starting point."""
    omega = 120
    assert np.array_equal(counts_margin(2, omega, 0.10),
                          counts_share(2, omega, 0.55))


# -- the guard that keeps 'champion wins' meaningful ------------------------

@pytest.mark.parametrize("n,share", [(2, 0.50), (3, 0.50), (4, 0.25)])
def test_champion_always_starts_strictly_ahead(n, share):
    """A tie at t=0 would make 'champion wins' a coin flip that reads as
    physics. Both builders must guarantee a strict lead."""
    counts = counts_share(n, 120, share)
    assert counts[0] > counts[1:n].max()


def test_counts_conserve_omega():
    for n in (2, 5, 13):
        for builder, arg in ((counts_margin, 0.10), (counts_share, 0.55)):
            counts = builder(n, 120, arg)
            assert counts.sum() == 120
            assert counts[-1] == 0            # B starts empty, per design.md 4


# -- the measured answer ----------------------------------------------------

def test_penalty_exists_under_fixed_margin_and_vanishes_under_fixed_share():
    """The result, at reduced trials. Measured at 3000 trials, Omega=120:
    margin 0.9710 -> 0.3703 across n=2..24; share 0.9710 -> 1.0000."""
    omega, trials = 120, 600
    m2 = p_champion_wins(2, omega, counts_margin(2, omega), trials)["p_win"]
    m16 = p_champion_wins(16, omega, counts_margin(16, omega), trials)["p_win"]
    s16 = p_champion_wins(16, omega, counts_share(16, omega), trials)["p_win"]
    assert m16 < m2 - 0.3, (m2, m16)          # a real penalty under fixed margin
    assert s16 > 0.99, s16                    # and none under fixed share


def test_share_alone_does_not_determine_difficulty():
    """The other half of why fixed share is not a fair comparator: at the SAME
    champion share of 0.50, P(win) is 0.606 at n=2 and 0.997 at n=3, because the
    same share is a pairwise lead of ~0 at n=2 and 0.25 at n=3."""
    omega, trials = 120, 800
    p2 = p_champion_wins(2, omega, counts_share(2, omega, 0.50), trials)["p_win"]
    p3 = p_champion_wins(3, omega, counts_share(3, omega, 0.50), trials)["p_win"]
    assert p3 - p2 > 0.25, (p2, p3)
