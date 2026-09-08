"""Suite-wide guards for mutable global state.

`experiments/chemical_cascade.py` keeps the Hill coupling's exponent and half-saturation as
MODULE-LEVEL globals, `HILL_N` and `HILL_K`, and at least six experiments assign to them
(`margin_law` sets them to arbitrary values mid-sweep and restores afterwards; `jensen_shift`,
`static_transfer_limit`, `timescale_ratio` and `penalty_interpolation` reset them defensively,
which is itself evidence that leaking has been a worry).

**They are load-bearing.** Measured: setting `cc.HILL_N, cc.HILL_K = 6.0, 1.3` moves
`free_upstream_depth.solve(14, 2, 0, 2.0)` by 3.3e-4 — not a rounding difference, a different
chemistry. Nine test files import modules that touch them. A test that leaves them changed would
silently alter every later test's physics while the suite stayed green, which is the harness doing
something the chemistry did not (rule 10) and the hardest kind of error to notice: no failure, just
different numbers.

This fixture makes such a leak loud. It does not restore them — restoring would hide which test
leaked, and the point is to name it.
"""

import pytest

import experiments.chemical_cascade as cc

HILL_DEFAULTS = (4.0, 1.0)


@pytest.fixture(autouse=True)
def _hill_globals_are_not_leaked():
    before = (cc.HILL_N, cc.HILL_K)
    assert before == HILL_DEFAULTS, (
        f"chemical_cascade's Hill globals were already {before} on entry, not {HILL_DEFAULTS}: "
        f"an EARLIER test leaked them, and every test since has been running different chemistry."
    )
    yield
    after = (cc.HILL_N, cc.HILL_K)
    assert after == HILL_DEFAULTS, (
        f"this test left chemical_cascade's Hill globals at {after}, not {HILL_DEFAULTS}. "
        f"Restore them in a finally: block -- a leak changes the physics of every later test "
        f"without failing anything."
    )
