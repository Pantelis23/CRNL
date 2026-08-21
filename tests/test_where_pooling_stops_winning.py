"""§112 — the merge threshold as rail geometry, and the two instrument replacements."""

import math

import numpy as np
import pytest

from experiments.where_pooling_stops_winning import ELEMENTS, R1, R3, Element


def _el(name):
    return Element(name, R1, dict(ELEMENTS)[name], R3)


ALL = [_el(n) for n, _ in ELEMENTS]


# ---------------------------------------------------------------- the constructed elements

@pytest.mark.parametrize("el", ALL, ids=[e.name for e in ALL])
def test_rails_round_trip_exactly(el):
    """The whole section rests on being able to ask for three roots and get them."""
    rts = el.roots()
    assert len(rts) == 3, f"{el.name} is not bistable: {rts}"
    for got, want in zip(rts, (el.r1, el.r2, el.r3)):
        assert got == pytest.approx(want, abs=1e-9)
    assert all(v > 0 for v in el.c), "a rate constant went non-positive"


def test_the_elements_are_not_one_element_reparametrised():
    """c differs by an order of magnitude across the set, so they are distinct chemistry."""
    cs = [el.fit(el.grid(2.0, 8.0, 4))[0] for el in ALL]
    assert max(cs) / min(cs) > 9.0, f"escape actions too similar to be independent: {cs}"


# ---------------------------------------------------------------- instrument (i): the MFPT

@pytest.mark.parametrize("N", (10, 14, 20, 30, 42, 60))
def test_log_domain_mfpt_reproduces_the_dense_solve_where_it_is_sound(N):
    """P1's gate. §111's dense solve is correct below its ceiling and must be matched."""
    from experiments.depth_compounding import R2 as S_R2, R3 as S_R3
    from experiments.is_the_vote_a_majority import hold_lifetime as dense

    # box held at §111's 1.25: this test compares SOLVERS, not box widths
    sch = Element("S", 0.15, S_R2, S_R3, cap_mult=1.25)
    assert abs(sch.log_hold_lifetime(N) - math.log(dense(N))) < 1e-9


def test_log_domain_mfpt_survives_past_the_dense_solve_ceiling():
    """§111 capped at k*N <= 150 because the dense solve returned NEGATIVE lifetimes."""
    from experiments.depth_compounding import R2 as S_R2, R3 as S_R3
    from experiments.is_the_vote_a_majority import MFPT_CEILING, hold_lifetime as dense

    sch = Element("S", 0.15, S_R2, S_R3, cap_mult=1.25)
    vals = [sch.log_hold_lifetime(N) for N in (200, 400, 700, 1000)]
    assert all(np.isfinite(v) and v > 0 for v in vals)
    assert vals == sorted(vals)
    # the local slope must converge, not wander
    sl = [(vals[i + 1] - vals[i]) for i in range(len(vals) - 1)]
    sl = [s / d for s, d in zip(sl, (200, 300, 300))]
    assert max(sl) / min(sl) < 1.01, f"escape action not converging: {sl}"
    # and the dense solve really does fail up there -- this is why the replacement exists
    with pytest.raises(FloatingPointError):
        dense(3 * MFPT_CEILING)


def test_it_reproduces_section_111s_published_headline():
    """Rule 7: §111's number stands, so the new instrument must land on it, not near it."""
    from experiments.depth_compounding import R2 as S_R2, R3 as S_R3
    from experiments.is_the_vote_a_majority import VOLUMES, m_of_k

    sch = Element("S", 0.15, S_R2, S_R3, cap_mult=1.25)
    c, a, _ = sch.fit(VOLUMES)
    assert 2.0 * (math.log(sch.tau) - a) == pytest.approx(-3.4131, abs=5e-4)
    first = sch.log_hold_lifetime(30) - sch.ln_remerge(sch.log_hold_lifetime(10), 3, m_of_k(3))
    assert first == pytest.approx(-3.4236, abs=5e-4)


def test_remerge_lifetime_matches_section_111s_closed_form():
    """The tau term is in the DENOMINATOR. Getting this sign wrong flipped P1 by 226%."""
    from experiments.is_the_vote_a_majority import remerge_lifetime

    el = _el("E-070")
    T, k, m = 12.5, 3, el.m_of_k(3)
    assert m == 3
    ref = T ** m / (math.comb(k, m) * el.tau ** (m - 1))
    assert math.exp(el.ln_remerge(math.log(T), k, m)) == pytest.approx(ref, rel=1e-12)
    # and against §111's own function, on §111's own rails
    from experiments.depth_compounding import R2 as S_R2, R3 as S_R3
    from experiments.is_the_vote_a_majority import m_of_k

    sch = Element("S", 0.15, S_R2, S_R3, cap_mult=1.25)
    mm = m_of_k(3)
    assert math.exp(sch.ln_remerge(math.log(T), 3, mm)) == pytest.approx(
        remerge_lifetime(T, 3, sch.tau), rel=1e-12)


# ------------------------------------------------------- instrument (ii): the commit probability

@pytest.mark.parametrize("el", ALL, ids=[e.name for e in ALL])
@pytest.mark.parametrize("k", (3, 5))
def test_commit_probability_is_exact_at_the_two_rails(el, k):
    """j=0 puts the pool ON the high rail and j=k on the low one: no dynamics needed."""
    N = 20
    assert el.commit_high(k, 0, N) == pytest.approx(1.0, abs=1e-12)
    assert el.commit_high(k, k, N) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("el", ALL, ids=[e.name for e in ALL])
def test_commit_probability_is_monotone_in_the_number_of_failures(el):
    ps = [el.commit_high(5, j, 24) for j in range(6)]
    # monotone up to float rounding: the partial sums can overshoot 1 by ~1e-16
    assert all(ps[i] >= ps[i + 1] - 1e-12 for i in range(len(ps) - 1)), ps
    assert ps[0] > ps[-1]


def test_commit_probability_has_no_settle_time_knob():
    """§111 evolved for 20 tau; a splitting probability has no time in it at all."""
    import inspect

    src = inspect.signature(Element.commit_high).parameters
    assert set(src) == {"self", "k", "j", "N"}, "a free time parameter crept back in"


# ---------------------------------------------------------------- P2b: the limit, not a tolerance

@pytest.mark.parametrize("name,k,j", (("E-070", 7, 5), ("E-084", 7, 6), ("E-021", 5, 1)))
def test_soft_cells_converge_toward_the_sign_of_their_own_margin(name, k, j):
    """T16-c. §111 could not resolve these; the criterion is the limit (rule 20)."""
    el = _el(name)
    c, _, _ = el.fit(el.grid(2.0, 8.0, 4))
    N0 = max(10, min(60, int(round(6.0 / c))))
    series = [el.commit_high(k, j, om) for om in (N0, 2 * N0, 4 * N0, 8 * N0)]
    margin = el.x_merged(k, j) - el.r2
    assert abs(margin) < 0.11, "this cell is supposed to be a near-saddle one"
    if margin > 0:
        assert series == sorted(series) and series[-1] > series[0]
    else:
        assert series == sorted(series, reverse=True) and series[-1] < series[0]


# ---------------------------------------------------------------- the law itself

@pytest.mark.parametrize("el", ALL, ids=[e.name for e in ALL])
@pytest.mark.parametrize("k", (3, 5, 7))
def test_m_is_where_the_merged_mean_crosses_the_saddle(el, k):
    """m = floor(f k) + 1 must be the first j whose pooled mean sits below r2 -- definitional,
    and worth pinning because the whole section hangs on it."""
    m = el.m_of_k(k)
    assert el.x_merged(k, m) < el.r2
    assert el.x_merged(k, m - 1) > el.r2


def test_unanimity_and_pooling_boundaries_are_pure_rail_geometry():
    for el in ALL:
        for k in (3, 5, 7, 9, 11):
            m = el.m_of_k(k)
            assert (m == k) == (el.f >= 1.0 - 1.0 / k), f"{el.name} k={k}"
            assert (m == 1) == (el.f < 1.0 / k), f"{el.name} k={k}"


def test_e084_straddles_the_unanimity_boundary_within_one_element():
    """P3: same rails, same c and a -- unanimity at k=3,5 and a finite crossover at k=7."""
    el = _el("E-084")
    assert el.m_of_k(3) == 3 and el.m_of_k(5) == 5, "must be unanimous below the boundary"
    assert el.m_of_k(7) == 6, "and must NOT be unanimous above it"
    assert 1 - 1 / 5 <= el.f < 1 - 1 / 7


def test_k_independence_holds_exactly_on_the_symmetric_band():
    """P6, arithmetic: (m-1)/(k-m) = 1 iff |f - 1/2| < 1/(2k). AM sits at f = 1/2."""
    for k in (3, 5, 7, 9, 11):
        for f in np.linspace(0.05, 0.95, 91):
            m = int(math.floor(f * k)) + 1
            flat = (m != k) and (m - 1) == (k - m)
            assert flat == (0.5 - 0.5 / k <= f < 0.5 + 0.5 / k), f"k={k} f={f}"


@pytest.mark.parametrize("name,k", (("E-070", 3), ("E-084", 3), ("E-084", 5)))
def test_unanimity_gives_a_ratio_flat_in_volume(name, k):
    """The measured ln(L_hold/L_remerge) must be near-flat where m = k -- slope small
    against the c(k-m) that a NON-unanimous k of the same element shows."""
    el = _el(name)
    gr = el.grid(0.6, 10.0, 14)
    lr = [el.ln_ratio(N, k) for N in gr]
    slope = float(np.polyfit(np.array(gr, float), np.array(lr), 1)[0])
    lr7 = [el.ln_ratio(N, 7) for N in gr]
    slope7 = float(np.polyfit(np.array(gr, float), np.array(lr7), 1)[0])
    assert el.m_of_k(k) == k
    # The claim is that re-merging wins at EVERY volume; that is the assertion below and it is
    # exact. The residual drift is ln T's curvature and is reported in §112.3 rather than gated
    # -- the factor here only has to separate "flat" from "rising", within one element (rule 18).
    assert abs(slope) < 0.25 * abs(slope7), f"{name} k={k}: {slope} vs {slope7}"
    assert all(v < 0 for v in lr), "re-merging must win at every volume, not just on average"


def test_one_failure_regime_means_pooling_wins_at_every_volume():
    """P4. E-021 has f < 1/3, so m = 1 at k = 3 and there is no crossover on that side."""
    el = _el("E-021")
    assert el.m_of_k(3) == 1
    gr = el.grid(0.6, 10.0, 14)
    lr = [el.ln_ratio(N, 3) for N in gr]
    assert all(v > 0 for v in lr), lr
    assert lr[-1] > lr[0]


def test_the_straddling_pair_behaves_qualitatively_differently():
    """P5: E-062 and E-070 differ in r2 by 0.228 and land on opposite sides of 2/3."""
    lo, hi = _el("E-062"), _el("E-070")
    assert lo.m_of_k(3) == 2 and hi.m_of_k(3) == 3
    assert abs(lo.r2 - hi.r2) == pytest.approx(0.228, abs=1e-9)
    lr_lo = [lo.ln_ratio(N, 3) for N in lo.grid(0.6, 10.0, 14)]
    lr_hi = [hi.ln_ratio(N, 3) for N in hi.grid(0.6, 10.0, 14)]
    assert all(v > 0 for v in lr_lo), "below the boundary, pooling takes over"
    assert all(v < 0 for v in lr_hi), "above it, re-merging never loses"


def test_predicted_crossover_is_where_the_ln_ratio_actually_vanishes():
    """Guard the algebra of §34's formula against the ln-ratio it is supposed to zero."""
    el = _el("E-070")
    gr = el.grid(0.6, 10.0, 14)
    c, a, _ = el.fit(gr[len(gr) // 2:])
    for k in (5, 7):
        xo = el.predicted_crossover(c, a, k)
        m = el.m_of_k(k)
        model = c * xo * (k - m) + a * (1 - m) + math.log(math.comb(k, m)) \
            + (m - 1) * math.log(el.tau)
        assert model == pytest.approx(0.0, abs=1e-9)


def test_the_absolute_crossover_is_window_limited_not_exact():
    """P7. Refitting over a window covering k*Omega moves E-084 k=7 by 18 points, so no
    cell is resolved better than the window ambiguity. Guards against over-claiming."""
    el = _el("E-084")
    gr = el.grid(0.6, 10.0, 14)
    c_a, a_a, _ = el.fit(gr[len(gr) // 2:])
    cov = sorted({int(round(17 * (117 / 17) ** (i / 7))) for i in range(8)})
    c_c, a_c, _ = el.fit(cov)
    p_a, p_c = el.predicted_crossover(c_a, a_a, 7), el.predicted_crossover(c_c, a_c, 7)
    assert abs(p_a / p_c - 1) > 0.15, "the window really does move this prediction"


@pytest.mark.parametrize("el", ALL, ids=[e.name for e in ALL])
def test_the_inherited_box_width_was_not_converged(el):
    """Rule 13, and the criterion is CONVERGENCE, not a tolerance (rule 20).

    A first version of this test asserted |ln T(2.0) - ln T(3.0)| < 1e-9, which is a fixed gate
    on a converging quantity and failed on the shallowest element at 3.8e-5. What the box has to
    do is converge: successive differences must shrink, and the last one must vanish.
    """
    vs = [Element(el.name, el.r1, el.r2, el.r3, cap_mult=c).log_hold_lifetime(20)
          for c in (1.25, 1.5, 2.0, 3.0, 4.0)]
    d = [abs(vs[i + 1] - vs[i]) for i in range(len(vs) - 1)]
    assert all(d[i] > d[i + 1] for i in range(len(d) - 1)), f"not converging: {d}"
    assert d[-1] == 0.0, "3.0 and 4.0 must agree exactly -- that is the box being irrelevant"
    assert d[0] > 0.04, "and 1.25 must be visibly off, or there was nothing to fix"
    assert vs == sorted(vs), "a wider box can only lengthen the escape"


def test_the_commit_probability_is_exactly_box_independent():
    """Why §112.1 and §112.2 survive the box correction untouched."""
    worst = 0.0
    for el in ALL:
        for k in (3, 5, 7):
            for j in range(k + 1):
                a = Element(el.name, el.r1, el.r2, el.r3, cap_mult=1.25).commit_high(k, j, 30)
                b = Element(el.name, el.r1, el.r2, el.r3, cap_mult=3.0).commit_high(k, j, 30)
                worst = max(worst, abs(a - b))
    assert worst == 0.0, worst


@pytest.mark.parametrize("N", (11, 30, 77))
def test_vectorised_rates_match_the_shared_helper(N):
    """_rates restates two expressions for speed; it must not restate the chemistry."""
    import experiments.chemical_cascade as cc

    for el in ALL:
        lam, mu, cap = el._rates(N)
        for s in range(0, cap + 1, max(1, cap // 25)):
            l, m = cc.rates_stage(float(s), 0.0, N, el.c, el.r3, True, "hill")
            assert mu[s] == pytest.approx(m, rel=1e-12, abs=1e-12)
            if s != cap:
                assert lam[s] == pytest.approx(l, rel=1e-12, abs=1e-12)
        assert lam[cap] == 0.0, "the top of the box must be reflecting"
