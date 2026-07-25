"""The exact cascade kernel.

Every test here guards an invariant that a previous version of this experiment
got wrong, or that only an independent instrument can check. Two full designs
were discarded before this one:

  * design 1 stopped a stage at 0.7*delta_star and emitted +-1, so the stop
    predicate fired on the INITIAL state and 83-96% of stages ran zero
    reactions -- the harness did the restoring, for free;
  * design 2 fixed that but compared against a control confined to +-1 while
    scaling the noise by delta_star(gamma), which manufactured an entire
    headline ("restoration needs a minimum Omega") out of the mismatch.

So the tests below check the CONTROL as carefully as the chemistry.
"""

import numpy as np
import pytest

from crnl.networks.am_reversible import GAMMA_C
from crnl.cascade_exact import (
    TIE_BAND, channel_matrix, input_alphabet, run_cascade,
    run_control, stage_kernel, verdict,
)
from crnl.networks.am_reversible import (
    am_reversible, delta_star, reverse_pairing,
)
from crnl.stochastic import seed_for
from crnl.thermo import decompose, gillespie_instrumented
from crnl.vectorized import compile_network


# -- alphabets and the parity trap -----------------------------------------

@pytest.mark.parametrize("omega,expected", [(30, 31), (60, 61), (120, 121)])
def test_input_alphabet_has_omega_plus_one_sites(omega, expected):
    """Omega+1, NOT Omega/2+1. The miscount put a memory table 2x low."""
    din = input_alphabet(omega)
    assert len(din) == expected
    assert (np.abs(np.diff(din)) == 2).all()
    assert ((omega + din) % 2 == 0).all()


def test_kernel_rows_are_probability_distributions():
    k = stage_kernel(0.3, 30, 2.0)
    assert np.allclose(k["K"].sum(axis=1), 1.0, rtol=0.0, atol=1e-9)
    assert (k["K"] >= -1e-12).all()
    assert len(k["dout"]) == 61          # outputs span BOTH parities


def test_about_half_the_output_mass_lands_on_the_opposite_parity():
    """The quantitative form of the parity trap. A loose '> 0.05' bound would
    pass while 90% of the mass was wrongly discarded."""
    omega = 30
    k = stage_kernel(0.3, omega, 2.0)
    opposite = [(omega + int(d)) % 2 == 1 for d in k["dout"]]
    mass = k["K"][:, opposite].sum(axis=1)
    assert 0.30 < mass.min(), mass.min()
    assert 0.35 < mass.mean() < 0.65, mass.mean()


def test_collapsing_the_alphabets_destroys_the_signal_entirely():
    """Documents the ACTUAL failure mode. The spec claimed collapsing drives
    fidelity to 0.5; it drives it to 0, because unrenormalised rows lose ~half
    their mass every stage (0.5^depth), which is a different symptom to look
    for."""
    omega, depth = 30, 12
    k = stage_kernel(0.3, omega, 2.0)
    din, dout, K = k["din"], k["dout"], k["K"]
    idx = {int(d): i for i, d in enumerate(din)}
    collapsed = np.zeros((len(din), len(din)))
    for o, d in enumerate(dout):
        if int(d) in idx:                      # silently drop odd parity
            collapsed[:, idx[int(d)]] += K[:, o]
    p = np.zeros(len(din)); p[np.argmin(np.abs(din / omega - 0.7))] = 1.0
    for _ in range(depth):
        p = p @ collapsed
    assert p.sum() < 0.01, p.sum()             # mass gone, not redistributed


# -- <M> against an independent instrument ---------------------------------

def test_exact_M_matches_an_unbiased_ssa_oracle():
    """The augmented generator vs the integer counter.

    halt_before_tmax=True is REQUIRED: the default loop applies the jump that
    crosses t_max (gillespie_fast's convention), which biases the counter +6.7%
    here -- z = +7.3 at 6000 trials against the exact value, versus z = +0.3
    with the exact stop. A test written without it passes at 800 trials and
    fails as trials rise, incriminating the innocent kernel.
    """
    gamma, omega, t_stage, trials = 0.3, 30, 2.0, 4000
    k = stage_kernel(gamma, omega, t_stage)
    din = list(k["din"])
    d = int(k["din"][np.argmin(np.abs(k["din"] / omega - 0.4))])
    exact = float(k["M"][din.index(d)])

    n_x = (omega + d) // 2
    n0 = np.array([n_x, omega - n_x, 0], dtype=np.int64)
    net = am_reversible(gamma)
    comp = compile_network(net, float(omega))
    pairing = reverse_pairing(net)
    got = [gillespie_instrumented(comp, n0, seed_for(omega, t, base=91), pairing,
                                  t_max=t_stage, halt_before_tmax=True,
                                  species=list(net.species)).net_firings
           for t in range(trials)]
    mean = float(np.mean(got))
    sem = float(np.std(got, ddof=1) / np.sqrt(trials))
    assert abs(mean - exact) < 4.0 * sem, f"{mean=} +-{sem=} vs {exact=}"


def test_stage_cost_goes_through_thermo_decompose():
    """The /3 lives in thermo.decompose and nowhere else (Global Constraint).
    Writing (A/3)*M inline in run_cascade would be invisible until the two
    silently diverged."""
    k = stage_kernel(0.3, 30, 2.0)
    i = len(k["din"]) // 2
    expected = decompose(None, None, float(k["M"][i]), k["affinity"],
                         boundary=float(k["boundary"][i]))["total"]
    assert expected == pytest.approx(
        k["boundary"][i] + (k["affinity"] / 3.0) * k["M"][i], rel=1e-12)


# -- the channel ------------------------------------------------------------

def test_channel_rows_are_normalised_on_the_two_alphabet_path():
    """Exercises the production path (output lattice -> input lattice), not the
    square special case, so a bin-edge bug on the odd rows is reachable."""
    omega = 30
    C = channel_matrix(np.arange(-omega, omega + 1), input_alphabet(omega),
                       omega, 0.25)
    assert C.shape == (61, 31)
    assert np.allclose(C.sum(axis=1), 1.0, atol=1e-10)


def test_odd_outputs_split_evenly_between_their_two_neighbours():
    omega = 30
    C = channel_matrix(np.array([1]), input_alphabet(omega), omega, 1e-7)
    assert C.max() == pytest.approx(0.5, abs=1e-6)
    assert sorted(C[0])[-2:] == pytest.approx([0.5, 0.5], abs=1e-6)


# -- guards -----------------------------------------------------------------

@pytest.mark.parametrize("gamma", [0.0, GAMMA_C, 0.7, 1.0])
def test_no_landscape_raises_rather_than_returning_nan(gamma):
    """delta_star = 0 above GAMMA_C makes sigma_ch = 0 and the channel NaN.
    theta_decide already refuses here; so must this."""
    with pytest.raises(ValueError, match="bistable|landscape"):
        stage_kernel(gamma, 30, 2.0)


def test_zero_or_negative_stage_time_raises():
    with pytest.raises(ValueError, match="positive"):
        stage_kernel(0.3, 30, 0.0)


# -- symmetry and conventions ----------------------------------------------

def test_cascade_is_antisymmetric_under_bit_flip():
    r = run_cascade(0.15, 30, 4.0, 8)
    flipped = run_cascade(0.15, 30, 4.0, 8, start_delta=-delta_star(0.15))
    for a, b in zip(r["p_correct"], flipped["p_correct"]):
        assert a + b == pytest.approx(1.0, abs=1e-9)


def test_start_delta_reproduces_the_findings_7_convention():
    """FINDINGS 7 starts from a weak s_init=0.3 against sigma=0.35 and so opens
    at Phi(0.857) = 0.80; starting from the rail opens near 1. Both are
    defensible and they are DIFFERENT experiments -- pinned so the comparison is
    never quietly claimed to be quantitative."""
    weak = run_cascade(0.05, 60, 8.0, 3, start_delta=0.3)
    rail = run_cascade(0.05, 60, 8.0, 3)
    assert weak["p_correct"][0] == pytest.approx(0.808, abs=0.01)
    assert rail["p_correct"][0] > 0.99


def test_strict_and_half_credit_differ_only_near_gamma_c():
    """The half-credit tie convention is cosmetically load-bearing near
    gamma_c: it turns a strict 0.4938 into a tidier 0.5005."""
    far = run_cascade(0.05, 60, 8.0, 20)
    near = run_cascade(0.45, 60, 8.0, 20)
    assert abs(far["p_correct"][-1] - far["p_correct_strict"][-1]) < 1e-3
    assert near["p_correct"][-1] - near["p_correct_strict"][-1] > 3e-3


# -- the control IS an axis: the regression guard on the withdrawn headline --

def test_control_rail_changes_the_verdict_and_so_must_be_reported():
    """THE guard on the artifact that produced a withdrawn finding.

    At gamma=0.30, Omega=30 the chemistry LOSES to a control given the full
    +-1 lattice and WINS against one railed to the chemistry's own +-delta_star.
    The chemistry arm is identical in both. If this test ever stops
    discriminating, the two conventions have converged and the experiment may
    quote one control; until then it must quote both.
    """
    omega, gamma, depth = 30, 0.30, 30
    r = run_cascade(gamma, omega, 8.0, depth)
    d = r["delta_star"]
    wide = run_control(omega, depth, r["sigma_ch"], d, rail=1.0)
    matched = run_control(omega, depth, r["sigma_ch"], d, rail=d)
    assert verdict(r["p_correct"][-1], wide["p_correct"][-1]) == "loses"
    assert verdict(r["p_correct"][-1], matched["p_correct"][-1]) == "wins"


def test_matched_rail_control_is_gamma_independent():
    """The property that makes it a fair comparator: with the rail matched to
    the landscape the control's difficulty no longer varies with gamma (spread
    < 0.01), whereas the +-1 control swings by ~0.15 across the same sweep --
    monotonically, in the direction of the claimed effect."""
    depth = 30
    matched, wide = [], []
    for gamma in (0.05, 0.15, 0.30, 0.45):
        d = delta_star(gamma)
        s = 0.35 * d
        matched.append(run_control(60, depth, s, d, rail=d)["p_correct"][-1])
        wide.append(run_control(60, depth, s, d, rail=1.0)["p_correct"][-1])
    assert max(matched) - min(matched) < 0.01, matched
    assert max(wide) - min(wide) > 0.10, wide


def test_verdict_tie_band_suppresses_meaningless_differences():
    assert verdict(0.5000001, 0.5000000) == "tie"
    assert verdict(0.60, 0.50) == "wins"
    assert verdict(0.50, 0.60) == "loses"
    assert TIE_BAND > 0.0


# -- physics that survives every convention --------------------------------

def test_restoration_degrades_monotonically_toward_gamma_c():
    """The robust finding: the chemistry arm alone, no control involved."""
    f = [run_cascade(g, 60, 8.0, 30)["p_correct"][-1]
         for g in (0.05, 0.15, 0.30, 0.45)]
    assert f == sorted(f, reverse=True), f
    assert f[0] > 0.85 and f[-1] < 0.55


def test_fidelity_improves_with_population_at_fixed_gamma():
    f = [run_cascade(0.30, om, 8.0, 30)["p_correct"][-1] for om in (30, 60, 120)]
    assert f == sorted(f), f


def test_cost_is_extensive_in_omega():
    """dS/stage scales ~linearly in Omega (the O(Omega) scaling of FINDINGS 9.2).
    Measured at depth 30, gamma=0.15, t=8: 19.7 / 40.3 / 81.8 for Om=30/60/120.
    Uses the SAME depth as the quoted numbers -- an earlier version ran depth 5
    while its docstring quoted depth-30 values."""
    depth = 30
    a = run_cascade(0.15, 30, 8.0, depth)["ds_per_stage"]
    b = run_cascade(0.15, 60, 8.0, depth)["ds_per_stage"]
    assert 1.6 < b / a < 2.6, b / a


def test_ds_per_stage_is_depth_dependent_and_labelled_as_such():
    """cum[-1]/depth averages a non-stationary sequence, so it is only
    meaningful with its depth quoted. Pinned so nobody compares two runs at
    different depths."""
    short = run_cascade(0.15, 30, 8.0, 1)["ds_per_stage"]
    long_ = run_cascade(0.15, 30, 8.0, 30)["ds_per_stage"]
    assert abs(long_ - short) > 1.0
