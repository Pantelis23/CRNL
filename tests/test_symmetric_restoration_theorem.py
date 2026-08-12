"""The Symmetric Restoration Theorem (FINDINGS §65), as one machine-checked object.

Each clause is tested on the SAME randomly drawn exchange-symmetric networks, so this file
asserts the theorem rather than five separate facts that happen to appear in one document.
Scope is asserted too: clause 1 is checked to FAIL when the exchange symmetry is broken, and
clause 4 is applied only where the symmetric steady state is unique.

Prior art and what is and is not claimed: see FINDINGS §65.1. Clause 1's invariance half is
standard (Golubitsky-Stewart-Schaeffer); the novelty claim rests on clause 2 and is recorded
as UNVERIFIED pending T-THM-a.
"""
from __future__ import annotations

import numpy as np
import pytest

from crnl.cme import generator
from crnl.networks.am_reversible import am_reversible
from experiments.amplification_sign import P_at
from experiments.amplification_signature import mirror_pairs
from experiments.exchange_theorem import probe, random_network
from experiments.free_rate_optimum import build_free
from experiments.optimal_element import symmetric_classes
from experiments.restoration_boundary import P_closed, symmetric_fixed_points
from experiments.restoration_cone import restores, orthant_states, v_of
from experiments.threshold_sharpness import antisym_block, lambda_A


def _nets(n, seed, symmetrise=True):
    """Random exchange-symmetric networks with both a float state and integer counts."""
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < n:
        net = random_network(rng, n_extra=2, n_rx=5, max_order=3, symmetrise=symmetrise)
        if net is None:
            continue
        x = np.exp(rng.uniform(-2, 2, len(net.species)))
        counts = rng.integers(6, 40, size=len(net.species))
        out.append((net, x, counts))
    return out


def test_clause1_invariance_holds_and_needs_the_symmetry():
    """b_X - b_Y = (n_X - n_Y) P(n): divisible, and NOT divisible without the symmetry."""
    worst = 0.0
    for net, _, counts in _nets(40, 20260812):
        div = probe(net, counts, 1.0, 0, 1)[0]          # [0] is the divisibility residual
        worst = max(worst, abs(div))
    assert worst < 1e-10, worst

    bad = 0
    for net, _, counts in _nets(40, 20260812, symmetrise=False):
        div = probe(net, counts, 1.0, 0, 1)[0]
        if abs(div) > 1e-8:
            bad += 1
    assert bad > 0, "scope: clause 1 must fail somewhere without exchange symmetry"


def test_clause2_decomposition_is_exact_and_every_bracket_is_nonnegative():
    """P = sum_r c_r d_r B_r with B_r >= 0, and d_r = 0 contributes identically zero."""
    checked = 0
    for net, x, _ in _nets(40, 7):
        v, c = v_of(net, x, "X", "Y")
        if v is None or v.size == 0:
            continue
        assert float(np.dot(c, v)) == pytest.approx(
            P_at(net, x, list(net.species).index("X"), list(net.species).index("Y")),
            rel=1e-7, abs=1e-12)
        pairs = mirror_pairs(net, "X", "Y")
        for slot, (d, idx) in enumerate(pairs):
            r = net.reactions[idx]
            p, q = r.reactants.get("X", 0), r.reactants.get("Y", 0)
            # p == q records an arbitrary +-1 for d (§55's counting bug); d == 0 is a
            # self-mirror term. Both contribute identically zero and carry no bracket.
            if p == q or d == 0:
                assert abs(float(v[slot])) < 1e-12 or p == q
                continue
            B = float(v[slot]) / d
            assert B >= -1e-12, (B, d)      # every bracket is nonnegative
        checked += 1
    assert checked >= 20


def test_clause3_capability_is_combinatorial():
    """Some d_r > 0 (with p != q) => restores for SOME c on the open orthant."""
    rng = np.random.default_rng(99)
    capable_tested = 0
    for net, _, _c in _nets(60, 11):
        pairs = mirror_pairs(net, "X", "Y")
        if not pairs:
            continue
        ds = [(d, i) for d, i in pairs
              if net.reactions[i].reactants.get("X", 0)
              != net.reactions[i].reactants.get("Y", 0)]
        if not ds:
            continue
        states = orthant_states(rng, len(net.species), 40)
        # load the rates onto the positive-d reactions, starve the rest
        c = np.array([50.0 if d > 0 else 1e-3 for d, _ in pairs])
        if any(d > 0 for d, _ in ds):
            assert restores(net, c, states, "X", "Y") is True
            capable_tested += 1
        else:
            assert restores(net, np.ones(len(pairs)), states, "X", "Y") is False
    assert capable_tested >= 5


def test_clause4_realisation_is_one_linear_inequality_and_is_not_convex():
    """Sign of P at x* is the criterion; AM gives gamma_c = 1/2; the set is not convex."""
    from scipy.optimize import brentq

    def amP(g):
        net = am_reversible(g)
        fps = symmetric_fixed_points(net, n=801)
        assert len(fps) == 1
        return P_closed(net, fps[0])

    assert brentq(amP, 0.05, 0.95, xtol=1e-14) == pytest.approx(0.5, abs=1e-12)

    cls = symmetric_classes()
    c1 = np.array([0.8273, 7.9309, 5.2147, 0.4279])
    c2 = np.array([1.6947, 0.0782, 2.8352, 1.3938])
    got = []
    for c in (c1, c2, c1 + c2):
        net = build_free([cls[7], cls[8]], c)
        fps = symmetric_fixed_points(net, n=2001)
        assert len(fps) == 1
        got.append(P_closed(net, fps[0]))
    assert got[0] > 0 and got[1] > 0 and got[2] < 0


def test_clause5_generator_block_diagonalises_and_lambda_A_is_negative():
    """spec(Q_A) subset spec(Q); lambda_A < 0 at every finite Omega, either side of gamma_c."""
    for om in (8, 12):
        net = am_reversible(0.30)
        M, n = antisym_block(net, om)
        full = np.linalg.eigvals(generator(net, om, float(om)).toarray())
        for e in np.linalg.eigvals(M.toarray()):
            assert np.abs(full - e).min() < 1e-9
        assert 0 < n < len(full)
    for g in (0.30, 0.49, 0.51, 0.70):
        for om in (20, 40):
            assert lambda_A(am_reversible(g), om) < 0, (g, om)


def test_scope_clause4_needs_a_unique_symmetric_steady_state():
    """T15-j: networks with more than one symmetric steady state are OUTSIDE the theorem."""
    rng = np.random.default_rng(20260812)
    cls = symmetric_classes()
    multi = 0
    for _ in range(120):
        k = int(rng.integers(1, 4))
        ids = rng.choice(len(cls), size=k, replace=False)
        net = build_free([cls[i] for i in ids], np.exp(rng.uniform(-3, 3, 2 * k)))
        try:
            if len(symmetric_fixed_points(net, n=801)) > 1:
                multi += 1
        except Exception:
            pass
    assert multi > 0, "the excluded class must be non-empty, or the scope note is wrong"
