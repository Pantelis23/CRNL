"""§104 — predicting p_transmit, and §103 without a free parameter."""

import numpy as np
import pytest

from experiments.chain_without_a_joint_solve import (
    MEASURED as CHAIN_MEASURED,
    chain_operating_points,
    split_from,
)
from experiments.depth_compounding import R1, R2, R3
from experiments.predicting_transmission import (
    MEASURED,
    descent_rate,
    p_transmit,
    pinned_roots,
)


def test_downstream_is_monostable_once_the_upstream_reaches_its_saddle():
    """P1, and its second half: the high rail is gone BEFORE the upstream formally fails."""
    assert len(pinned_roots(R3)) == 3
    assert len(pinned_roots(R2)) == 1
    assert len(pinned_roots(R1)) == 1
    assert pinned_roots(R1)[0] < 0.2


def test_k_low_is_macroscopic():
    """P4: a descent rate must be nearly Omega-independent; an escape rate moves 3 orders."""
    ks = [descent_rate(om)[0] for om in (14, 30, 55)]
    assert max(ks) / min(ks) < 1.15
    assert all(4.0 < k < 7.0 for k in ks)


def test_absorbing_boundary_matches_the_measured_definition():
    """§104.1(a): last_low counts n < R2*om, so absorption is one site BELOW ceil(R2*om)."""
    import experiments.predicting_transmission as pt
    import inspect
    src = inspect.getsource(pt.descent_rate)
    assert "int(np.ceil(R2 * om)) - 1" in src, "the off-by-one fix must stay fixed"


def test_predicted_curve_is_low_everywhere_and_worst_at_the_shortest_window():
    """P3: the residual is one-signed and monotone in the window."""
    k, _ = descent_rate(30)
    res = {t: (p_transmit(k, t) - m) / m for t, m in MEASURED.items()}
    assert all(r < 0 for r in res.values()), "every residual must be negative"
    ts = sorted(res)
    assert [abs(res[t]) for t in ts] == sorted((abs(res[t]) for t in ts), reverse=True)
    assert res[0.5] == pytest.approx(-0.0895, abs=2e-3)
    assert res[8.0] == pytest.approx(-0.0058, abs=2e-3)


def test_p2_was_refuted():
    """Recorded as refuted, not quietly widened: the computed rate misses the interval."""
    k, _ = descent_rate(30)
    assert k == pytest.approx(5.5169, abs=5e-3)
    assert not (7.071 <= k <= 8.771)


def test_operating_point_start_moves_it_only_three_percent():
    """§104.1(b): the candidate cause was tested and is largely refuted."""
    import experiments.predicting_transmission as pt
    k_rail, _ = pt.descent_rate(30)
    # start the descent from the degraded operating point instead of the pristine rail
    import experiments.chemical_cascade as cc
    from experiments.depth_compounding import C

    def from_point(om, start_x, x_up=R1, cap_mult=1.25):
        cap = int(np.ceil(cap_mult * R3 * om))
        a = int(np.ceil(R2 * om)) - 1
        idx = np.arange(a, cap + 1)
        A = np.zeros((len(idx), len(idx)))
        rhs = -np.ones(len(idx))
        for i, s in enumerate(idx):
            if s == a:
                A[i, i] = 1.0
                rhs[i] = 0.0
                continue
            lam, mu = cc.rates_stage(float(s), x_up * om, om, C, R3, False, "hill")
            if s == cap:
                lam = 0.0
            A[i, i] = -(lam + mu)
            if s + 1 <= cap:
                A[i, i + 1] = lam
            if s - 1 >= a:
                A[i, i - 1] = mu
        T = np.linalg.solve(A, rhs)
        return 1.0 / float(T[list(idx).index(int(round(start_x * om)))])

    k_op = from_point(30, 2.9759)
    assert k_op > k_rail, "starting lower must speed the descent"
    assert k_op / k_rail < 1.06, "but only marginally -- the gap is 22-37%"


def test_chain_stays_inside_the_gate_with_no_free_parameter():
    """§104.1: the whole point. Predicted p_transmit, four cells, §102's gate."""
    k, _ = descent_rate(30)
    pt_pred = p_transmit(k, 2.0)
    assert pt_pred == pytest.approx(0.9094, abs=2e-3)
    for (om, D), ref in CHAIN_MEASURED.items():
        mus, _ = chain_operating_points(om, D)
        _, c, p = split_from(om, mus, 2.0, p_transmit=pt_pred, legacy=True)
        ratio = (c / p) / (ref["contam"] / ref["pure"])
        assert 0.5 < ratio < 2.0, f"Om={om} D={D}: {ratio}"


def test_the_improvement_is_a_cancellation_not_an_improvement():
    """Predicted p_t is ~3% low and the model runs high; the errors partly cancel."""
    k, _ = descent_rate(30)
    pt_pred = p_transmit(k, 2.0)
    assert pt_pred < MEASURED[2.0], "the predicted p_transmit is LOW"
    for (om, D), ref in CHAIN_MEASURED.items():
        mus, _ = chain_operating_points(om, D)
        meas = ref["contam"] / ref["pure"]
        _, c1, p1 = split_from(om, mus, 2.0, legacy=True)
        _, c2, p2 = split_from(om, mus, 2.0, p_transmit=pt_pred, legacy=True)
        assert (c1 / p1) / meas > 1.0, "the model runs high with the measured p_t"
        assert (c2 / p2) / meas < (c1 / p1) / meas, "and lowering p_t moves it toward 1"
