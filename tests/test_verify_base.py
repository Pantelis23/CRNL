"""The base audit, as a gate rather than a report.

`experiments/verify_base.py` re-derives every load-bearing closed form in
FINDINGS from the current code. The rest of the suite proves the code is
self-consistent; this proves it still agrees with what is PUBLISHED, which is a
different question and the one that rots silently as modules are added and
behavioural functions change (FINDINGS 15 changed `wall_coefficient`).
"""
from experiments.verify_base import main


def test_published_closed_forms_still_hold():
    assert main() == 0, "a published closed form no longer matches the code"
