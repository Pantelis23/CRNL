"""Reversible n-winner AM and the affinity floor (FINDINGS 13)."""

import numpy as np
import pytest

from crnl.networks.n_winner_reversible import (
    affinity_critical, breaking_mode, gamma_critical, lambda_breaking,
    n_winner_reversible, symmetric_state,
)
from crnl.deterministic import jacobian


# -- the n=2 case must reproduce FINDINGS 9.1 exactly ----------------------

@pytest.mark.parametrize("gamma", [0.0, 0.1, 0.25, 0.4, 0.5, 0.6, 0.9])
def test_n2_reproduces_the_known_closed_form(gamma):
    """lambda_antisym = (1-2 gamma)/3 from am_reversible. If this drifts, the
    generalisation is not a generalisation."""
    assert lambda_breaking(2, gamma) == pytest.approx((1 - 2 * gamma) / 3, abs=1e-12)


def test_n2_recovers_gamma_c_one_half_and_three_ln_two():
    assert gamma_critical(2) == pytest.approx(0.5, abs=1e-12)
    assert affinity_critical(2) == pytest.approx(3 * np.log(2), abs=1e-11)


def test_n2_symmetric_point_is_gamma_independent_at_one_third():
    """The degeneracy that is SPECIAL to n=2 and does not survive n>=3."""
    for gamma in (0.02, 0.3, 0.8):
        x, b = symmetric_state(2, gamma)
        assert x == pytest.approx(1 / 3, abs=1e-12)
        assert b == pytest.approx(1 / 3, abs=1e-12)


# -- structure -------------------------------------------------------------

def test_network_shape_and_reversibility():
    n = 4
    net = n_winner_reversible(n, 0.3)
    assert net.n_species == n + 1
    assert net.n_reactions == 2 * (n * (n - 1) // 2 + n)   # forwards + reverses
    fwd = net.n_reactions // 2
    for j in range(fwd):
        f, r = net.reactions[j], net.reactions[j + fwd]
        assert f.reactants == r.products and f.products == r.reactants
        assert r.k == pytest.approx(0.3 * f.k)


@pytest.mark.parametrize("n", [2, 3, 4, 6])
def test_cycle_dimension_is_n_choose_2(n):
    """Counting each reversible PAIR as one edge -- the counting under which
    FINDINGS 9.1 called AM's cycle space one-dimensional. C(n,2)+n edges against
    rank(S)=n leaves exactly C(n,2): 1, 3, 6, 15.

    Mixing this with the count that treats forward and reverse as separate edges
    (4, 9, 16) is how an earlier draft of this module claimed n=3 had a
    9-dimensional cycle space while comparing it to 9.1's 1.
    """
    net = n_winner_reversible(n, 0.3)
    forward = net.n_reactions // 2
    S_forward = net.stoichiometry_matrix()[:, :forward]
    rank = np.linalg.matrix_rank(S_forward)
    assert rank == n
    assert forward - rank == n * (n - 1) // 2


def test_symmetric_state_is_a_fixed_point():
    for n in (2, 3, 5, 9):
        for gamma in (0.01, 0.2, 0.7):
            x, b = symmetric_state(n, gamma)
            net = n_winner_reversible(n, gamma)
            state = np.concatenate([np.full(n, x), [b]])
            assert np.abs(net.rhs(state)).max() < 1e-12
            assert x * n + b == pytest.approx(1.0, abs=1e-12)


def test_symmetric_state_moves_with_gamma_for_n_at_least_three():
    """Unlike n=2. Measured at n=3: x = 0.2050 / 0.2269 / 0.2431."""
    xs = [symmetric_state(3, g)[0] for g in (0.02, 0.2, 0.6)]
    assert xs == sorted(xs)
    assert xs[-1] - xs[0] > 0.03
    assert not any(x == pytest.approx(1 / 4, abs=1e-3) for x in xs)  # not 1/(n+1)


def test_breaking_mode_is_a_true_eigenvector():
    """Otherwise the Rayleigh quotient would be a bound, not the eigenvalue."""
    for n in (3, 8, 16):
        gamma = 0.05
        x, b = symmetric_state(n, gamma)
        J = jacobian(n_winner_reversible(n, gamma),
                     np.concatenate([np.full(n, x), [b]]))
        v = breaking_mode(n)
        Jv = J @ v
        lam = v @ Jv
        assert np.linalg.norm(Jv - lam * v) / np.linalg.norm(Jv) < 1e-10


# -- the measured law ------------------------------------------------------

def test_gamma_c_falls_much_faster_than_one_over_n():
    """Kills the pre-run prediction gamma_c = 1/n (which would give A_c=3 ln n).
    Measured at n=32: 3.58e-5 against 1/n = 0.031, ~870x below."""
    assert gamma_critical(32) < 0.031 / 100


def test_affinity_floor_approaches_nine_ln_n():
    """A_c/ln n climbs 3.0 -> ~9. Measured 8.9894 at n=256."""
    ratios = [affinity_critical(n) / np.log(n) for n in (2, 8, 32, 128)]
    assert ratios == sorted(ratios)
    assert ratios[0] == pytest.approx(3.0, abs=1e-6)     # exact at n=2
    assert 8.9 < ratios[-1] < 9.0


def test_the_landscape_really_dies_above_gamma_c():
    """The eigenvalue crossing is a bifurcation, not just a sign change:
    below gamma_c a broken state survives, above it symmetry is restored."""
    from scipy.integrate import solve_ivp

    n = 4
    gc = gamma_critical(n)
    for factor, expect_broken in ((0.7, True), (1.4, False)):
        gamma = gc * factor
        x, b = symmetric_state(n, gamma)
        net = n_winner_reversible(n, gamma)
        y0 = np.concatenate([np.full(n, x), [b]])
        y0[0] += 1e-6
        y0[1] -= 1e-6
        sol = solve_ivp(lambda t, y: net.rhs(y), (0, 3000.0), y0,
                        method="LSODA", rtol=1e-11, atol=1e-14)
        spread = sol.y[:n, -1].max() - sol.y[:n, -1].min()
        # bool(): numpy comparisons return np.bool_, and np.bool_(True) is not
        # the Python singleton, so `is expect_broken` fails on a passing case
        assert bool(spread > 0.5) is expect_broken, (gamma, spread)


def test_rejects_degenerate_inputs():
    with pytest.raises(ValueError, match="n >= 2"):
        n_winner_reversible(1, 0.3)
    with pytest.raises(ValueError, match="gamma"):
        n_winner_reversible(3, -0.1)


# -- the quasipotential ingredients (FINDINGS 14) --------------------------

@pytest.mark.parametrize("n", [2, 3, 4, 6, 8, 16, 32, 64])
def test_closed_forms_match_the_numeric_route(n):
    """The closed forms exist only for speed -- n=4096 would need 8.4M
    reactions -- so they must be pinned against the network computation."""
    from crnl.networks.n_winner_reversible import (
        breaking_diffusion, diffusion_closed, lambda_closed,
    )
    assert lambda_closed(n) == pytest.approx(lambda_breaking(n, 0.0), rel=1e-9)
    assert diffusion_closed(n) == pytest.approx(breaking_diffusion(n, 0.0), rel=1e-9)


def test_diffusion_at_n2_is_the_design_doc_value():
    """design.md section 9 derives D = 1/(9 Omega) for irreversible AM."""
    from crnl.networks.n_winner_reversible import breaking_diffusion
    assert breaking_diffusion(2, 0.0) == pytest.approx(1 / 9, rel=1e-12)


def test_predicted_barrier_reproduces_the_design_doc_at_n2():
    """c(eps) = (3/2) eps^2, the one first-principles prediction in the project."""
    from crnl.networks.n_winner_reversible import predicted_barrier
    for delta in (0.04, 0.10, 0.20):
        assert predicted_barrier(2, delta) == pytest.approx(1.5 * delta**2, rel=1e-12)


def test_lambda_and_diffusion_vanish_at_the_same_rate():
    """THE mechanism behind FINDINGS 3's saturation: both go like 1/(2n-1), so
    their ratio (2n-1)/(2n-3) tends to 1 instead of diverging."""
    from crnl.networks.n_winner_reversible import diffusion_closed, lambda_closed
    ratios = [lambda_closed(n) / diffusion_closed(n) for n in (2, 8, 64, 1024)]
    assert ratios == sorted(ratios, reverse=True)
    assert ratios[0] == pytest.approx(3.0)         # n=2
    assert ratios[-1] == pytest.approx(1.0, abs=0.002)


def test_predicted_barrier_saturates_like_the_measurement():
    """Predicted floor delta^2/2. The measured floor is 2.3x lower -- a constant
    offset, not a different shape (FINDINGS 14)."""
    from crnl.networks.n_winner_reversible import predicted_barrier
    d = 0.10
    assert predicted_barrier(1024, d) == pytest.approx(d * d / 2, rel=0.002)
    assert predicted_barrier(64, d) / predicted_barrier(1024, d) < 1.02


# -- the pairwise multiplicative identity (FINDINGS 30, T15-a) -------------

def _pair_bracket(n, gamma, omega, counts, i, j):
    others = counts[:n].sum() - counts[i] - counts[j]
    return (counts[n] - others - gamma * (counts[i] + counts[j] - 1.0)) / omega


@pytest.mark.parametrize("n", [2, 3, 4, 6])
@pytest.mark.parametrize("gamma_frac", [0.1, 0.6, 0.9])
def test_pairwise_drift_is_exactly_multiplicative(n, gamma_frac):
    """d(n_i - n_j)/dt = (n_i - n_j) * [n_B - sum_{l!=i,j} n_l - gamma(n_i+n_j-1)]/Omega.

    FINDINGS 30. This is why `bookkeeping-only` returns a categorical zero at every n
    (24.2) and why `decision-only` froze tied rivals (24.2's Omega-parity trap): the
    identity makes sign(n_i - n_j) a conserved quantity once that difference direction
    carries no noise. It holds only because the reverse disagreement reaction
    `2B -> X_i + X_j` reaches X_i and X_j through the same (n-1) pairs, so a change to
    that reaction's stoichiometry or rate convention breaks the theorem silently --
    the CLE arms would simply stop returning zero and nothing else would complain.
    """
    from crnl.approximations import propensities_batch
    from crnl.vectorized import compile_network

    gamma = gamma_frac * gamma_critical(n)
    omega = 97
    comp = compile_network(n_winner_reversible(n, gamma), float(omega))
    rng = np.random.default_rng(4242 + n)
    for _ in range(12):
        cuts = np.sort(rng.integers(0, omega + 1, size=n))
        counts = np.diff(np.concatenate([[0], cuts, [omega]])).astype(np.int64)
        assert counts.sum() == omega
        a = propensities_batch(comp, counts[None, :].astype(float))[0]
        b = comp.S @ a
        for i in range(n):
            for j in range(i + 1, n):
                d = float(counts[i] - counts[j])
                rhs = d * _pair_bracket(n, gamma, omega, counts, i, j)
                traffic = float(np.abs(comp.S[i] - comp.S[j]) @ a)
                assert abs(float(b[i] - b[j]) - rhs) <= 1e-12 * max(traffic, 1.0)


def test_tied_species_have_exactly_equal_drift():
    """The n_i = n_j corner of the identity, where it forces b_i - b_j == 0.

    This is the case that made 24.2's `decision-only` arm return 0: rivals that start
    tied cannot be separated by drift, only by noise in their own difference.
    """
    from crnl.approximations import propensities_batch
    from crnl.vectorized import compile_network

    n, omega = 3, 90
    gamma = 0.6 * gamma_critical(n)
    comp = compile_network(n_winner_reversible(n, gamma), float(omega))
    for tied in (5, 12, 20, 28):
        counts = np.array([omega - 2 * tied - 10, tied, tied, 10], dtype=np.int64)
        assert (counts > 0).all()
        a = propensities_batch(comp, counts[None, :].astype(float))[0]
        b = comp.S @ a
        assert b[1] == pytest.approx(b[2], abs=1e-12)


def test_champion_margin_is_eroded_quadratically_by_rival_spread():
    """du/dt = u*Gbar - (1-gamma)*delta_23^2/(4 Omega) at n = 3 (FINDINGS 30.1).

    The additive term is why `rivals-only` is NOT protected by a conservation law even
    though it never failed in 440,000 trajectories: the champion's mean margin has no
    noise but it does have a sink, quadratic in how far the rivals have spread.
    """
    from crnl.approximations import propensities_batch
    from crnl.vectorized import compile_network

    n, omega = 3, 97
    for gamma_frac in (0.1, 0.6, 0.9):
        gamma = gamma_frac * gamma_critical(n)
        comp = compile_network(n_winner_reversible(n, gamma), float(omega))
        rng = np.random.default_rng(909)
        for _ in range(12):
            cuts = np.sort(rng.integers(0, omega + 1, size=n))
            c = np.diff(np.concatenate([[0], cuts, [omega]])).astype(np.int64)
            a = propensities_batch(comp, c[None, :].astype(float))[0]
            b = comp.S @ a
            u = float(c[0]) - 0.5 * float(c[1] + c[2])
            d23 = float(c[1] - c[2])
            gbar = 0.5 * (_pair_bracket(n, gamma, omega, c, 0, 1)
                          + _pair_bracket(n, gamma, omega, c, 0, 2))
            rhs = u * gbar - (1.0 - gamma) * d23 ** 2 / (4.0 * omega)
            assert float(b[0] - 0.5 * (b[1] + b[2])) == pytest.approx(rhs, abs=1e-11)


# --- FINDINGS 42/43: the identity needs exchange symmetry, not conservation laws ---


def test_pairwise_identity_survives_extra_conservation_laws():
    """FINDINGS 42: cofactor pairs add conservation laws and the identity is unmoved.

    B (two laws) and C (four laws) keep the constant-ratio property; the conservation
    structure is not what the identity depends on.
    """
    from experiments.conservation_identity import am_cofactor, split_spread

    omega = 60
    for net in (am_cofactor(0.25), am_cofactor(0.25, double=True)):
        counts = [17, 23, 20] + [13, 47] * ((len(net.species) - 3) // 2)
        spread, _ = split_spread(net, counts, omega)
        assert spread < 1e-12


def test_pairwise_identity_dies_when_only_a_rate_constant_breaks_symmetry():
    """FINDINGS 42: D differs from B in ONE rate constant and the identity fails.

    Same species, same reactions, same orders, same two conservation laws -- so nothing
    structural is available as an alternative explanation for the failure.
    """
    from experiments.conservation_identity import am_cofactor, split_spread

    omega = 60
    counts = [17, 23, 20, 13, 47]
    assert split_spread(am_cofactor(0.25, beta=0.0), counts, omega)[0] < 1e-12
    assert split_spread(am_cofactor(0.25, beta=0.2), counts, omega)[0] > 1e-3


def test_exchange_symmetry_forces_divisibility_by_the_lead():
    """FINDINGS 43: an antisymmetric polynomial is divisible by the difference.

    Symmetrised random networks -- including ones with NO conservation law at all --
    have b_X - b_Y exactly zero at n_X = n_Y. Unsymmetrised ones generically do not.
    """
    from experiments.exchange_theorem import probe, random_network

    rng = np.random.default_rng(4242)
    omega, failures = 60, 0
    for sym in (True, False):
        rng = np.random.default_rng(4242)
        for _ in range(40):
            net = random_network(rng, n_extra=2, n_rx=5, max_order=4, symmetrise=sym)
            counts = [int(rng.integers(6, omega)) for _ in net.species]
            div = probe(net, counts, omega)[0]
            if sym:
                assert div < 1e-12
            else:
                failures += int(div > 1e-12)
    assert failures > 20, "probe has no power if unsymmetrised networks pass too"


def test_cubic_pair_term_separates_divisibility_from_the_constant_ratio():
    """FINDINGS 43: am_cubic is divisible but its ratio is affine in delta^2.

    The delta^2 coefficient is k/(4 Omega^2) EXACTLY, derived by hand -- an absolute
    check, not a fit (rule 16).
    """
    from experiments.exchange_theorem import am_cubic, probe

    omega = 60
    div, spread, deltas, ratios = probe(am_cubic(0.25), [21, 19, 20], omega)
    assert div < 1e-12
    assert spread > 1e-3
    A = np.vstack([np.ones_like(deltas), deltas ** 2]).T
    coef, *_ = np.linalg.lstsq(A, ratios, rcond=None)
    assert np.max(np.abs(A @ coef - ratios)) < 1e-12 * np.abs(ratios).max()
    assert coef[1] == pytest.approx(1.0 / (4.0 * omega ** 2), rel=1e-9)
