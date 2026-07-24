"""Approximate Majority (AM) as data (design.md §2.2).

Two committed species X, Y and a blank B. Three reactions, all k = 1:

    r1:  X + Y -> 2B      disagreement -- a collision cancels both to blank
    r2:  B + X -> 2X      autocatalysis -- the leader recruits blanks
    r3:  B + Y -> 2Y      autocatalysis, mirror

Every reaction is 2->2, so total count is conserved. This is the first network
loaded into the engine, not a special case baked into it.
"""

from __future__ import annotations

from ..reactions import Reaction, ReactionNetwork

SPECIES = ["X", "Y", "B"]


def approximate_majority(k: float = 1.0) -> ReactionNetwork:
    return ReactionNetwork(
        species=list(SPECIES),
        reactions=[
            Reaction({"X": 1, "Y": 1}, {"B": 2}, k, name="r1:disagree X+Y->2B"),
            Reaction({"B": 1, "X": 1}, {"X": 2}, k, name="r2:recruit  B+X->2X"),
            Reaction({"B": 1, "Y": 1}, {"Y": 2}, k, name="r3:recruit  B+Y->2Y"),
        ],
        name="approximate-majority",
    )
