"""Tests for the cooling drive (FINDINGS 19).

`test_reduces_to_expanding_at_w_zero` is the anchor: with no cooling the new
integrator must reproduce `gillespie_expanding` exactly, the same bit-for-bit
check FINDINGS 5.1 used on the time change itself. Everything else in the module
is a generalisation of that case, so if this drifts, nothing else is meaningful.
"""
from __future__ import annotations

import numpy as np
import pytest

from crnl.cooling import CoolingSchedule, gillespie_cooling, reverse_mask
from crnl.expanding import gillespie_expanding
from crnl.networks.am_reversible import am_reversible, reverse_pairing
from crnl.vectorized import compile_network

N0 = np.array([30, 30, 30])
SPECIES = ["X", "Y", "B"]


def _setup(gamma, w):
    net = am_reversible(gamma)
    mask = reverse_mask(net, reverse_pairing(net))
    return compile_network(net, 90.0), CoolingSchedule.build(gamma, w), mask


def test_reverse_mask_selects_the_reverses():
    net = am_reversible(0.2)
    mask = reverse_mask(net, reverse_pairing(net))
    assert list(mask) == [False, False, False, True, True, True]


@pytest.mark.parametrize("hubble", [0.02, 0.05, 0.2])
def test_reduces_to_expanding_at_w_zero(hubble):
    """w = 0 is no cooling, so it must BE the FINDINGS 5.1 integrator."""
    comp, sched, mask = _setup(0.2, 0.0)
    for seed in range(40):
        a = gillespie_expanding(comp, N0, np.random.default_rng(seed), hubble,
                                species=SPECIES)
        b = gillespie_cooling(comp, N0, np.random.default_rng(seed), hubble,
                              sched, mask, species=SPECIES)
        assert np.array_equal(a.n_final, b.n_final), (seed, a.n_final, b.n_final)
        assert a.steps == b.steps, (seed, a.steps, b.steps)


def test_schedule_is_monotone_and_spans_the_full_drive():
    sched = CoolingSchedule.build(0.55, 1 / 3)
    assert sched.g[0] == pytest.approx(1.0)
    assert np.all(np.diff(sched.g) <= 1e-12)          # never increases
    assert sched.gamma_at(0.0) == pytest.approx(0.55)
    # pin the closed form rather than a guessed threshold: at s = 0.999,
    # (1-s)^(-1/3) = 10, so gamma = 0.55**10 = 0.00253
    assert sched.gamma_at(0.999) == pytest.approx(0.55 ** 10, rel=1e-3)
    assert sched.gamma_at(1.0 - 1e-9) < 1e-4          # cools essentially to zero
    assert np.all(np.diff(sched.cum) >= 0)


def test_schedule_has_no_h_dependence():
    """The whole point: gamma(s) is universal, H only sets how much fits in it."""
    a = CoolingSchedule.build(0.55, 1 / 3)
    b = CoolingSchedule.build(0.55, 1 / 3)
    assert np.array_equal(a.g, b.g)
    # and the pitchfork crossing is a property of (gamma0, w) alone
    assert a.s_of_gamma(0.5) == pytest.approx(0.358, abs=0.01)


def test_no_cooling_means_no_bifurcation_time():
    """At w = 0 the drive never moves, so asking where it crosses is undefined."""
    assert np.isnan(CoolingSchedule.build(0.55, 0.0).s_of_gamma(0.5))


@pytest.mark.parametrize("bad", [0.0, 1.0, 1.5, -0.1])
def test_schedule_rejects_gamma0_outside_the_open_unit_interval(bad):
    with pytest.raises(ValueError):
        CoolingSchedule.build(bad, 1 / 3)


def test_cooling_rejects_zero_hubble():
    comp, sched, mask = _setup(0.55, 1 / 3)
    with pytest.raises(ValueError, match="hubble must be > 0"):
        gillespie_cooling(comp, N0, np.random.default_rng(0), 0.0, sched, mask)


def test_cooling_freezes_and_reports_the_drive_it_froze_at():
    comp, sched, mask = _setup(0.55, 1 / 3)
    r = gillespie_cooling(comp, N0, np.random.default_rng(1), 0.05, sched, mask,
                          species=SPECIES)
    assert r.status == "frozen" and r.frozen
    assert 0.0 < r.gamma_final < 0.55           # cooled below where it started
    assert r.n_final.sum() == N0.sum()          # every reaction is 2 -> 2


def test_cooling_decides_where_the_fixed_drive_cannot():
    """gamma0 = 0.55 has NO landscape at s=0; only cooling creates one.

    With w = 0 the drive sits above gamma_c forever and the population stays
    mixed; with w > 0 it crosses the pitchfork and resolves. This is the
    experiment's premise, so it is pinned here.
    """
    comp, cold, mask = _setup(0.55, 1 / 3)
    _, flat, _ = _setup(0.55, 0.0)
    minority = {}
    for name, sched in (("cooled", cold), ("flat", flat)):
        vals = []
        rng = np.random.default_rng(7)
        for _ in range(60):
            r = gillespie_cooling(comp, N0, rng, 0.01, sched, mask, species=SPECIES)
            nx, ny = int(r.n_final[0]), int(r.n_final[1])
            if nx + ny:
                vals.append(min(nx, ny) / (nx + ny))
        minority[name] = float(np.median(vals))
    assert minority["cooled"] < 0.05, minority
    assert minority["flat"] > 0.30, minority
