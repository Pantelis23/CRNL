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


# --- FINDINGS 44: Arrhenius kinetics and the rho lever ---


def test_delta_star_rho_reduces_to_delta_star_at_rho_one():
    """FINDINGS 44: the rho-generalised attractor separation must match the closed form."""
    from crnl.networks.am_reversible import delta_star
    from experiments.arrhenius_optimum import delta_star_rho

    for gamma in (0.05, 0.10, 0.25, 0.40, 0.49):
        assert delta_star_rho(gamma, 1.0) == pytest.approx(delta_star(gamma), abs=1e-14)


def test_delta_star_rho_refuses_the_spurious_landscape_below_rho_c():
    """FINDINGS 44: for rho < gamma the closed form returns a PLAUSIBLE positive value.

    0.806 at (0.25, 0.05), and 6.37 and 11.4 near rho = gamma where delta* <= 1 is a hard
    bound -- in a region where the ODE nullcline says there is no landscape at all. The
    guard rho > rho_c must return exactly 0 there.
    """
    from experiments.arrhenius_optimum import delta_star_rho, rho_critical

    for gamma in (0.25, 0.40):
        rc = rho_critical(gamma)
        assert rc > gamma, "rho_c > gamma is what lets one guard cover the double flip"
        for rho in (0.01, 0.05, 0.15, gamma, rc * 0.999):
            assert delta_star_rho(gamma, rho) == 0.0
        assert 0.0 < delta_star_rho(gamma, rc * 1.2) < 1.0 / (1.0 + gamma)


def test_delta_star_rho_matches_the_ode_nullcline():
    """FINDINGS 44: the closed form is checked against the network's own fluxes."""
    from scipy.optimize import brentq

    from experiments.arrhenius_optimum import am_rho, delta_star_rho, rho_critical

    for gamma in (0.10, 0.25, 0.40):
        for rho in (0.7, 1.0, 2.0, 5.0, 20.0):
            if rho <= rho_critical(gamma):
                continue
            net = am_rho(gamma, rho)
            S = net.stoichiometry_matrix()
            b = gamma / (1.0 + gamma)
            s = 1.0 - b

            def f(d, net=net, S=S, b=b, s=s):
                x, y = 0.5 * (s + d), 0.5 * (s - d)
                return float((S @ net.fluxes(np.array([x, y, b])))[:2].sum())

            assert f(1e-9) * f(s - 1e-9) < 0
            root = brentq(f, 1e-9, s - 1e-9, xtol=1e-15)
            assert delta_star_rho(gamma, rho) == pytest.approx(root, abs=1e-9)


def test_rho_leaves_the_cycle_affinity_at_three_ln_gamma():
    """FINDINGS 44: rho is a catalyst, so it must not touch the drive.

    A catalyst accelerates both directions equally. §16 pins the affinity of the cycle
    X+Y->2B, B+X->2X, B+Y->2Y at -3 ln gamma, and it must be rho-independent.
    """
    from experiments.arrhenius_optimum import am_rho

    for gamma in (0.10, 0.25, 0.40):
        for rho in (0.5, 1.0, 4.0):
            net = am_rho(gamma, rho)
            by = {r.name.split(":")[0]: r.k for r in net.reactions}
            aff = np.log((by["f1"] * by["f2"] * by["f3"])
                         / (by["r1"] * by["r2"] * by["r3"]))
            assert aff == pytest.approx(-3.0 * np.log(gamma), abs=1e-12)


def test_sep_of_matches_the_closed_form_at_rho_one():
    """FINDINGS 45: sep_of is anchored to 3(1+2g)/(1-2g) before being used off rho = 1.

    The whole of §45 reads a residual against sep on a network family where sep has no
    published closed form, so the instrument is pinned where one exists.
    """
    from experiments.arrhenius_optimum import am_rho
    from experiments.slaving_axis import sep_of

    for gamma in (0.07, 0.12, 0.20, 0.28, 0.35):
        sep, _ = sep_of(am_rho(gamma, 1.0))
        assert sep == pytest.approx(3.0 * (1.0 + 2 * gamma) / (1.0 - 2 * gamma), rel=1e-12)


def test_sep_is_non_monotone_in_rho():
    """FINDINGS 45: the control that made the test discriminating.

    sep dips near rho ~ 1.5 while §44.2's cost falls monotonically, so "residual tracks
    sep" and "residual tracks rho" predict opposite shapes. If sep ever became monotone in
    rho this test loses its power and §45's P2 would need rebuilding.
    """
    from experiments.arrhenius_optimum import am_rho
    from experiments.slaving_axis import sep_of

    rhos = [0.5, 1.0, 1.5, 2.0, 4.0, 8.0, 32.0]
    seps = [sep_of(am_rho(0.20, r))[0] for r in rhos]
    i = int(np.argmin(seps))
    assert 0 < i < len(seps) - 1, "sep must have an INTERIOR minimum in rho"
    assert seps[i] < seps[0] and seps[i] < seps[-1]
    assert seps[-1] > 5 * seps[i], "sep must still grow strongly at large rho"


def test_slow_eigenvalue_crosses_zero_mid_path():
    """FINDINGS 46: why only a HARMONIC path average of sep is defined.

    The slow eigenvalue changes sign near delta/delta* ~ 0.57, where the drift peaks, so
    sep(delta) diverges there. An arithmetic mean of sep is meaningless; <1/sep> is not.
    """
    from experiments.arrhenius_optimum import am_rho, delta_star_rho
    from experiments.path_separation import eigs_at

    net, ds = am_rho(0.20, 1.0), delta_star_rho(0.20, 1.0)
    slow = [eigs_at(net, f * ds)[1] for f in (0.35, 0.45, 0.55, 0.65, 0.75)]
    assert slow[0] > 0 and slow[-1] < 0, "the slow eigenvalue must change sign on the path"
    assert all(a > b for a, b in zip(slow, slow[1:])), "and do so monotonically"


def test_rho_scales_both_eigenvalues_together():
    """FINDINGS 46: the structural reason no ratio can explain the residual.

    rho = 0.5 -> 4 moves fast and slow by nearly the SAME factor, so the ratio barely
    changes -- which is why the matched pair stayed matched under every convention.
    """
    from experiments.arrhenius_optimum import am_rho, delta_star_rho
    from experiments.path_separation import eigs_at

    f = []
    for rho in (0.5, 4.0):
        ds = delta_star_rho(0.20, rho)
        f.append(eigs_at(am_rho(0.20, rho), 0.35 * ds))
    fast_ratio = f[1][0] / f[0][0]
    slow_ratio = f[1][1] / f[0][1]
    assert 3.0 < fast_ratio < 4.5 and 3.0 < slow_ratio < 4.5
    assert abs(fast_ratio - slow_ratio) / fast_ratio < 0.15, "they scale together"


def test_lag_epsilon_is_parameter_free_and_h_converged():
    """FINDINGS 47: the closed-form lag correction, with no fitted constant.

    eps = (dmu/ds)(ds*/ddelta)/(dnu/ds) is built from three finite differences, so its
    step h is a second axis (rule 13). It must be converged before any comparison.
    """
    from crnl.networks.am_reversible import reverse_pairing
    from experiments.arrhenius_optimum import am_rho, delta_star_rho
    from experiments.lag_absolute import predict

    net = am_rho(0.20, 1.0)
    ds = delta_star_rho(0.20, 1.0)
    vals = [predict(net, 0.35 * ds, 0.80 * ds, reverse_pairing(net), h)
            for h in (4e-4, 1e-4, 2.5e-5)]
    assert all(v is not None for v in vals)
    assert abs(vals[-1] - vals[0]) / abs(vals[-1]) < 1e-6
    assert vals[-1] == pytest.approx(0.090945, rel=1e-4)


def test_lag_prediction_explains_the_fitted_constant():
    """FINDINGS 47: <eps>*sep reproduces the C that §39.2 and §46 fitted separately.

    §39.2 fitted 0.6465 on the T axis and §46 fitted 0.5963 on the gamma axis, and the
    non-transfer was a published puzzle. The computed value lands between them.
    """
    from crnl.networks.am_reversible import reverse_pairing
    from experiments.arrhenius_optimum import am_rho, delta_star_rho
    from experiments.lag_absolute import predict
    from experiments.slaving_axis import sep_of

    vals = []
    for gamma in (0.07, 0.12, 0.20, 0.28):
        net = am_rho(gamma, 1.0)
        ds = delta_star_rho(gamma, 1.0)
        p = predict(net, 0.35 * ds, 0.80 * ds, reverse_pairing(net), 1e-4)
        vals.append(p * sep_of(net)[0])
    assert 0.50 < min(vals) and max(vals) < 0.75
    assert 0.55 < float(np.mean(vals)) < 0.70


def test_traversal_endpoints_match_the_lattice_the_cme_uses():
    """FINDINGS 48: rule 11 -- the control must share its endpoints with its arm.

    T_det integrated over the NOMINAL range differs from the realised lattice range the
    first-passage solve actually runs between. The offset is small (<1% in T_det) but it
    is what made the nominal ratios bounce with Omega instead of converging.
    """
    from crnl.networks.am_reversible import reverse_pairing
    from experiments.arrhenius_optimum import am_rho, delta_star_rho
    from experiments.lag_endpoints import traversal
    from experiments.slaving_axis import slaved

    gamma, omega = 0.35, 300
    ds = delta_star_rho(gamma, 1.0)
    net = am_rho(gamma, 1.0)
    st = slaved(net, 0.35 * ds)
    nb = int(round(st[2] * omega))
    d0 = max(1, int(round(0.35 * ds * omega)))
    if (omega - nb - d0) % 2:
        d0 -= 1
    thr = max(2, int(round(0.80 * ds * omega)))

    assert abs(d0 / omega / ds - 0.35) > 1e-3, "this cell must actually be off-lattice"
    t_nom = traversal(net, 0.35 * ds, 0.80 * ds, reverse_pairing(net))
    t_real = traversal(net, d0 / omega, thr / omega, reverse_pairing(net))
    assert t_nom is not None and t_real is not None
    assert 0 < abs(t_real - t_nom) / t_nom < 0.02, "small, but not zero"


def test_birthdeath_chain_needs_matched_lattice_endpoints():
    """FINDINGS 50: rule 11 again -- T_det must integrate between the chain's own endpoints.

    Integrating the unrounded limits while the chain runs between the rounded ones is an
    O(1/Omega) error in delta, which after multiplying by Omega is the size of the whole
    effect. It produced sign-flipping coefficients before it was fixed.
    """
    from crnl.networks.am_reversible import reverse_pairing
    from experiments.arrhenius_optimum import am_rho, delta_star_rho
    from experiments.birthdeath_absorption import cell
    from experiments.lag_endpoints import traversal

    gamma, rho = 0.07, 1.0
    ds = delta_star_rho(gamma, rho)
    net = am_rho(gamma, rho)
    vals = [cell(gamma, rho, om, 0.35, 0.80)["bd_coeff"] for om in (1000, 2000)]
    assert all(v > 0 for v in vals), "the exact chain must give a positive correction"
    # 2.1% between Omega = 1000 and 2000 at these nominal endpoints; the run's own
    # realised endpoints converge to 0.2%. Either way it is settling, not drifting.
    assert abs(vals[1] - vals[0]) / vals[1] < 0.03, "and be Omega-independent once converged"

    om = 1000
    m0, thr = int(round(0.35 * ds * om)), int(round(0.80 * ds * om))
    t_match = traversal(net, m0 / om, thr / om, reverse_pairing(net), n=4001)
    t_un = traversal(net, 0.35 * ds, 0.80 * ds, reverse_pairing(net), n=4001)
    # the mismatch, times Omega, is comparable to the coefficient being measured
    assert abs(t_match - t_un) * om > 0.1 * vals[0]


def test_signal_drift_is_exactly_linear_in_the_pool():
    """FINDINGS 51: mu = k*delta*(1 - (1+gamma)s), which is §30's bracket in concentrations.

    This is why the pool-fluctuation Jensen term is identically zero: the drift has no
    curvature in s. It is also rho-independent, because the disagreement channel moves
    delta by exactly zero -- §30's first cancellation.
    """
    from experiments.arrhenius_optimum import am_rho
    from experiments.lag_absolute import field

    for gamma in (0.07, 0.20, 0.35):
        for rho in (0.5, 1.0, 32.0):
            net = am_rho(gamma, rho)
            for d in (0.1, 0.3, 0.5):
                for s in (0.55, 0.70, 0.85):
                    if s <= d:
                        continue
                    mu = field(net, d, s)[0]
                    assert mu == pytest.approx(d * (1.0 - (1.0 + gamma) * s), abs=1e-14)


def test_pool_jensen_term_vanishes_identically():
    """FINDINGS 51: J = 0 exactly, so §39.1's candidate (iii) is dead for the time too."""
    from crnl.networks.am_reversible import reverse_pairing
    from experiments.arrhenius_optimum import am_rho, delta_star_rho
    from experiments.pool_jensen import J_of

    for gamma, rho in ((0.07, 1.0), (0.20, 1.0), (0.20, 32.0), (0.35, 1.0)):
        ds = delta_star_rho(gamma, rho)
        net = am_rho(gamma, rho)
        J, lo, hi = J_of(net, 0.35 * ds, 0.80 * ds, reverse_pairing(net))
        # The floor here is the second difference's own roundoff, eps/h^2 ~ 2e-6 in
        # d2mu/ds2 at h = 1e-5, which propagates to ~1e-8 in J. J would have to be
        # ~2.0 to explain §50's missing 79%, so this is seven orders short of mattering.
        assert abs(J) < 1e-6, f"Jensen term must vanish, got {J}"
        assert abs(lo) < 1e-4 and abs(hi) < 1e-4, "d2mu/ds2 must be identically zero"


def test_P_equals_the_symmetry_breaking_eigenvalue():
    """FINDINGS 53: §43's P and T7/§14's lambda are the same object.

    At a symmetric state P is the antisymmetric directional derivative of b_i - b_j; at a
    symmetric FIXED point that is the Jacobian eigenvalue along (1,-1,0).
    """
    from crnl.networks.am_reversible import am_reversible
    from experiments.amplification_sign import P_at, antisym_eig
    from experiments.slaving_axis import slaved

    for gamma in (0.05, 0.20, 0.35, 0.49):
        net = am_reversible(gamma)
        st = slaved(net, 0.0)
        assert P_at(net, st) == pytest.approx(antisym_eig(net, st), rel=1e-8)


def test_P_at_the_symmetric_point_is_one_minus_two_gamma_over_three():
    """FINDINGS 53: P = (1-2*gamma)/3 exactly, and 1/3 = 1/(2n-1) at n = 2 (T7/§14)."""
    from crnl.networks.am_reversible import am_reversible
    from experiments.amplification_sign import P_at
    from experiments.slaving_axis import slaved

    for gamma in (0.0, 0.10, 0.20, 0.35, 0.49, 0.55, 0.90):
        net = am_reversible(gamma)
        st = slaved(net, 0.0)
        assert P_at(net, st) == pytest.approx((1.0 - 2.0 * gamma) / 3.0, abs=1e-10)
    assert P_at(am_reversible(1e-12), slaved(am_reversible(1e-12), 0.0)) == \
        pytest.approx(1.0 / 3.0, abs=1e-9)


def test_attractors_lie_on_P_zero_set():
    """FINDINGS 53: d(delta)/dt = delta*P, so an off-symmetric fixed point forces P = 0."""
    from crnl.networks.am_reversible import am_reversible, delta_star
    from experiments.amplification_sign import P_at
    from experiments.slaving_axis import slaved

    for gamma in (0.05, 0.20, 0.35, 0.45):
        net = am_reversible(gamma)
        st = slaved(net, delta_star(gamma))
        assert abs(P_at(net, st)) < 1e-13


def test_P_decomposition_is_exact():
    """FINDINGS 54: P = sum_pairs d_r * c_r * (non-negative bracket), exactly."""
    from crnl.networks.am_reversible import am_reversible
    from experiments.amplification_sign import P_at
    from experiments.amplification_signature import P_decomposed
    from experiments.exchange_theorem import am_cubic

    for gamma in (0.05, 0.25, 0.55):
        for net in (am_reversible(gamma), am_cubic(gamma)):
            for s, d in ((0.4, 0.1), (0.7, 0.25), (0.9, 0.05)):
                x = [(s + d) / 2, (s - d) / 2, 1 - s]
                dec, brackets = P_decomposed(net, x)
                assert dec == pytest.approx(P_at(net, x), rel=1e-11)
                assert all(b >= 0 for b in brackets), "brackets must be non-negative"


def test_am_decomposes_into_the_section_30_bracket():
    """FINDINGS 54: AM is MIXED -- recruitment d=+1, its reverse d=-1 weighted by gamma.

    The disagreement channel is self-mirror and contributes exactly zero, which is §30's
    first cancellation and §51's rho-independence of mu, one structural cause.
    """
    from crnl.networks.am_reversible import am_reversible
    from experiments.amplification_signature import classify, mirror_pairs

    net = am_reversible(0.25)
    pairs = mirror_pairs(net)
    ds = sorted(d for d, _ in pairs)
    assert ds == [-1, 1], f"AM must have exactly d = -1 and +1, got {ds}"
    assert classify(net) == "mixed"
    ks = {d: net.reactions[i].k for d, i in pairs}
    assert ks[1] == pytest.approx(1.0) and ks[-1] == pytest.approx(0.25)


def test_unanimous_topologies_cannot_flip_sign():
    """FINDINGS 54: all d_r <= 0 forbids amplification whatever the rate constants."""
    import numpy as _np

    from experiments.amplification_sign import P_at
    from experiments.amplification_signature import classify
    from experiments.exchange_theorem import random_network

    rng = _np.random.default_rng(31337)
    checked = 0
    for _ in range(120):
        net = random_network(rng, n_extra=2, n_rx=5, max_order=3, symmetrise=True)
        if net is None or classify(net) != "all<=0":
            continue
        checked += 1
        for _ in range(8):
            x = rng.uniform(0.05, 1.5, len(net.species))
            assert P_at(net, x) <= 1e-9, "an all-d<=0 network amplified"
    assert checked >= 5, "need some all<=0 networks for this test to have power"


def test_P_at_symmetric_point_is_one_over_2n_minus_1():
    """FINDINGS 55: P(symmetric, gamma=0) = 1/(2n-1) = T7/§14's lambda(n), no fit."""
    from experiments.amplification_n import P_symmetric

    for n in (2, 3, 4, 5, 6):
        p, _ = P_symmetric(n, 0.0)
        assert p == pytest.approx(1.0 / (2 * n - 1), abs=1e-10)


def test_P_vanishes_at_gamma_critical():
    """FINDINGS 55: the decomposition predicts gamma_c(n) with no fitted parameter."""
    from experiments.amplification_n import P_symmetric

    for n in (2, 3, 4, 5, 6):
        p, _ = P_symmetric(n, gamma_critical(n))
        assert abs(p) < 1e-9, f"P should vanish at gamma_c(n={n}), got {p}"


def test_n_winner_term_counts_are_one_against_n_minus_one():
    """FINDINGS 55: exactly one amplifying term against n-1 contracting ones.

    p == q pairs have an EMPTY bracket sum and contribute zero; mirror_pairs records
    their d as an arbitrary +-1, which is what made §55's first pass miscount.
    """
    from experiments.amplification_signature import mirror_pairs

    for n in (2, 3, 4, 5, 6):
        net = n_winner_reversible(n, 0.5 * gamma_critical(n))
        i, j = net.species[0], net.species[1]
        pos = neg = zero = 0
        for d, idx in mirror_pairs(net, i, j):
            r = net.reactions[idx]
            if r.reactants.get(i, 0) == r.reactants.get(j, 0):
                zero += 1
            elif d > 0:
                pos += 1
            else:
                neg += 1
        assert (pos, neg, zero) == (1, n - 1, n - 2), f"n={n} gave {(pos, neg, zero)}"


def test_nonrestoring_rates_form_a_convex_cone():
    """FINDINGS 56: P is linear in c, so {c : P <= 0 everywhere} is closed under addition.

    c1, c2 both non-restoring must imply c1 + c2 non-restoring. One counterexample would
    refute the theorem, not merely this test.
    """
    import numpy as _np

    from experiments.amplification_signature import classify
    from experiments.exchange_theorem import random_network
    from experiments.restoration_cone import orthant_states, restores, v_of

    rng = _np.random.default_rng(777)
    checked = 0
    for _ in range(150):
        net = random_network(rng, n_extra=2, n_rx=5, max_order=3, symmetrise=True)
        if net is None or classify(net) != "mixed":
            continue
        states = orthant_states(rng, len(net.species), 30)
        v0, _ = v_of(net, states[0])
        if v0 is None or v0.size == 0:
            continue
        bad = []
        for _ in range(30):
            c = _np.exp(rng.uniform(-4, 4, v0.size))
            if restores(net, c, states) is False:
                bad.append(c)
            if len(bad) == 2:
                break
        if len(bad) < 2:
            continue
        checked += 1
        assert not restores(net, bad[0] + bad[1], states), "cone theorem violated"
    assert checked >= 3, "need some non-restoring pairs for this test to have power"


def test_capability_is_combinatorial():
    """FINDINGS 56: some d_r > 0 => restores for SOME c; all d_r <= 0 => never."""
    import numpy as _np

    from experiments.amplification_signature import classify, mirror_pairs
    from experiments.exchange_theorem import random_network
    from experiments.restoration_cone import orthant_states, restores, v_of

    rng = _np.random.default_rng(4242)
    seen = {"mixed": 0, "all<=0": 0}
    for _ in range(150):
        net = random_network(rng, n_extra=int(rng.integers(1, 4)),
                             n_rx=int(rng.integers(3, 8)), max_order=4, symmetrise=True)
        if net is None:
            continue
        cls = classify(net)
        if cls not in seen:
            continue
        states = orthant_states(rng, len(net.species), 60)
        v0, _ = v_of(net, states[0])
        if v0 is None or v0.size == 0:
            continue
        ds = _np.array([d for d, _ in mirror_pairs(net)], dtype=float)
        c = _np.where(ds > 0, 1e4, _np.where(ds < 0, 1e-4, 1.0))
        got = bool(restores(net, c, states))
        seen[cls] += 1
        if cls == "mixed":
            assert got, "a mixed network must be capable for some rates"
        else:
            assert not got, "an all-d<=0 network must never restore"
    assert seen["mixed"] >= 3 and seen["all<=0"] >= 3


def test_classify_ignores_zero_bracket_pairs():
    """FINDINGS 56: p == q pairs contribute nothing and must not enter the classification.

    Counting them misclassified 148 of 232 networks as 'mixed' when they are all-d<=0.
    """
    import numpy as _np

    from experiments.amplification_signature import classify, mirror_pairs
    from experiments.exchange_theorem import random_network

    rng = _np.random.default_rng(20260812)
    found = 0
    for _ in range(200):
        net = random_network(rng, n_extra=2, n_rx=5, max_order=4, symmetrise=True)
        if net is None:
            continue
        pairs = mirror_pairs(net)
        if pairs is None:
            continue
        raw = [d for d, _ in pairs]
        contrib = [d for d, i in pairs
                   if net.reactions[i].reactants.get("X", 0)
                   != net.reactions[i].reactants.get("Y", 0)]
        if len(raw) == len(contrib):
            continue
        found += 1
        # the classification must depend only on the contributing pairs
        expect = ("trivial" if not contrib or all(d == 0 for d in contrib)
                  else "all<=0" if all(d <= 0 for d in contrib)
                  else "all>=0" if all(d >= 0 for d in contrib) else "mixed")
        assert classify(net) == expect
    assert found >= 5, "need networks with p==q pairs for this test to have power"


def test_landscape_bracket_survives_large_delta_star():
    """FINDINGS 57: the screening bug that discarded the best cells.

    slaving_axis.delta_star_of returns None where delta* > ~0.95 because `slaved` dies as
    delta -> 1. AM's Q minimum sits at gamma ~ 0.05 where delta* = 0.952, so the first
    pass of §57's search threw away exactly the best cells.
    """
    from crnl.networks.am_reversible import am_reversible, delta_star
    from experiments.optimal_element import landscape
    from experiments.slaving_axis import delta_star_of

    assert delta_star_of(am_reversible(0.05)) is None, "the bug this test pins"
    for gamma in (0.05, 0.10, 0.20, 0.35, 0.45):
        assert landscape(am_reversible(gamma)) == pytest.approx(delta_star(gamma), abs=1e-9)


def test_optimal_element_family_is_the_expected_size():
    """FINDINGS 57: 16 exchange-symmetric classes, and AM is built from two of them."""
    from experiments.optimal_element import build, symmetric_classes

    cls = symmetric_classes()
    assert len(cls) == 16
    assert sum(1 for c in cls if len(c) == 1) == 2
    net = build([c for c in cls
                 if c[0] == (("X", "Y"), ("B", "B")) or c[0] == (("B", "X"), ("X", "X"))],
                0.2)
    assert len(net.reactions) == 6, "AM has 3 forward + 3 reverse"
