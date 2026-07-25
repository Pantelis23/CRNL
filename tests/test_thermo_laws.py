"""The two laws, in the forms that survive measurement.

Each test carries the counterexample that killed the naive version:

  * "at detailed balance entropy production is zero" -- true for the stationary
    RATE, false for a trajectory's cumulative medium EP (+28.5 k_B from a rail
    at Omega=30, which is just the entropy of mixing).
  * "a biased start at gamma=1 never reaches a rail" -- FALSE. The exact MFPT to
    |delta| >= 0.5 at gamma=1 is 16.8 / 344 / 7.8e4 at Omega = 12 / 30 / 60. So
    the assertion is a TIME BUDGET at a stated Omega, never impossibility.
  * "mean medium EP is non-negative" -- false at a CONDITIONED stopping time, but
    ONLY for the start the spec rejects. Measured, gamma=1:
        (1/3,1/3,1/3) start, theta=0.5*Omega : -6.095 (Om=30), -12.014 (Om=60)
        B(0)=0 start,       theta=0.35*Omega : +16.3
    The negative series belongs to the (1/3,1/3,1/3) start. Under the mandated
    B(0)=0 convention the medium EP is POSITIVE, because an n_B=0 state sits far
    below the multinomial mode and mixing dominates ordering. An earlier draft
    asserted the negative sign for the B(0)=0 start and would have failed while
    instructing the implementer that Plan 1 was broken.

Note decompose()["total"] is MEDIUM (heat) EP over the trajectory, not the
quantity the second law bounds -- the system term is computed nowhere here.
"""

import numpy as np
import pytest

from crnl.cme import first_passage
from crnl.networks.am_reversible import (
    am_reversible, cycle_affinity, delta_star, initial_counts, reverse_pairing,
)
from crnl.stochastic import seed_for
from crnl.thermo import decompose, gillespie_instrumented, ln_multinomial
from crnl.vectorized import compile_network


def test_cumulative_ep_from_a_rail_is_positive_even_at_detailed_balance():
    """Detailed balance kills the stationary RATE, not the relaxation."""
    omega = 30
    net = am_reversible(1.0)
    pairing = reverse_pairing(net)
    comp = compile_network(net, float(omega))
    n0 = np.array([omega, 0, 0], dtype=np.int64)          # a rail
    A = cycle_affinity(net, pairing)
    assert A == pytest.approx(0.0, abs=1e-12)             # gamma=1: no drive

    totals = [decompose(n0, gillespie_instrumented(
                  comp, n0, seed_for(omega, t, base=41), pairing, t_max=20.0,
                  species=list(net.species)).n_final,
              0, A)["total"] for t in range(40)]
    mean = float(np.mean(totals))
    # with A = 0 the total IS the boundary term, bounded above by the mode of
    # lnW (29.345 at N=30, the uniform state) since lnW(rail) = 0.
    assert 20.0 < mean <= ln_multinomial(np.array([10, 10, 10])) + 1e-9, mean


@pytest.mark.parametrize("omega,budget", [(12, 40.0), (30, 800.0)])
def test_biased_start_at_detailed_balance_does_reach_a_rail(omega, budget):
    """A time budget, not an impossibility claim (see module docstring).
    Measured: 16.83 at Omega=12, 344.3 at Omega=30."""
    net = am_reversible(1.0)
    theta = int(round(0.5 * omega))
    n_x = (omega + 2) // 2
    n0 = np.array([n_x, omega - n_x, 0], dtype=np.int64)
    fp = first_passage(net, omega, float(omega), n0,
                       lambda n: abs(int(n[0]) - int(n[1])) >= theta,
                       reverse_pairing(net))
    assert fp["valid"]
    assert 0.0 < fp["mean_time"] < budget


@pytest.mark.parametrize("gamma", [0.1, 0.3, 0.45])
def test_cycle_term_is_positive_when_driven(gamma):
    """For gamma < 1 the drive does net forward cycles, so the CYCLE term is
    positive. Deliberately asserted on the cycle term rather than on `total`:
    total = boundary + cycle and the boundary term is start-dependent, so a
    "total > 0" assertion is the one that silently depends on protocol."""
    omega = 60
    net = am_reversible(gamma)
    pairing = reverse_pairing(net)
    theta = max(2, int(round(0.7 * delta_star(gamma) * omega)))
    n0 = initial_counts(omega, gamma,
                        count_diff=max(1, int(round(0.2 * delta_star(gamma) * omega))))
    fp = first_passage(net, omega, float(omega), n0,
                       lambda n: abs(int(n[0]) - int(n[1])) >= theta, pairing)
    assert fp["valid"]
    dec = decompose(n0, None, fp["net_reaction_firings"],
                    cycle_affinity(net, pairing), boundary=fp["boundary"])
    assert dec["cycle"] > 0.0
    assert dec["total"] > 0.0


@pytest.mark.parametrize("omega,expected", [(30, -6.095), (60, -12.014)])
def test_medium_ep_at_a_conditioned_stop_is_negative_from_the_symmetric_start(
        omega, expected):
    """The counterexample, with its initial condition STATED.

    Heat flows out of the medium because the stopping rule selects trajectories
    that have ordered themselves; the second law bounds medium + system, and the
    system entropy fell by more. This holds for the (1/3,1/3,1/3) start with
    theta = 0.5*Omega -- NOT for the B(0)=0 start, where the same quantity is
    +16.3. Restricted to Omega <= 60: first_passage rejects Omega = 90 and 120
    here on residual.
    """
    net = am_reversible(1.0)
    pairing = reverse_pairing(net)
    n = omega // 3
    n0 = np.array([n, n, omega - 2 * n], dtype=np.int64)
    theta = int(round(0.5 * omega))
    fp = first_passage(net, omega, float(omega), n0,
                       lambda s: abs(int(s[0]) - int(s[1])) >= theta, pairing)
    assert fp["valid"]
    dec = decompose(n0, None, fp["net_reaction_firings"],
                    cycle_affinity(net, pairing), boundary=fp["boundary"])
    assert dec["total"] == pytest.approx(expected, abs=0.05)


def test_the_same_protocol_from_the_mandated_start_is_positive():
    """The other half of the pair: the sign depends on the START, not on gamma.
    Asserted so nobody 'fixes' the test above by changing the initial condition."""
    omega = 60
    net = am_reversible(1.0)
    pairing = reverse_pairing(net)
    n0 = initial_counts(omega, 0.49, count_diff=2)        # B(0) = 0
    theta = max(2, int(round(0.35 * omega)))
    fp = first_passage(net, omega, float(omega), n0,
                       lambda s: abs(int(s[0]) - int(s[1])) >= theta, pairing)
    assert fp["valid"]
    dec = decompose(n0, None, fp["net_reaction_firings"],
                    cycle_affinity(net, pairing), boundary=fp["boundary"])
    assert dec["total"] > 10.0
