"""n-winner Approximate Majority as data (design.md §8; spec 2026-07-24).

n committed species X1..Xn plus a blank B, all rate k=1:

    disagreement (O(n^2)):  Xi + Xj -> 2B    for every unordered pair i<j
    recruitment  (O(n)):    B + Xi -> 2Xi    for every i

Every reaction is heterobimolecular (two distinct reactant species, coeff 1),
so by the engine convention each has c = k/Omega and a = c*n_i*n_j -- there are
no homodimers (the recruitment product 2Xi is a *product*, not a reactant).
At n=2 this is exactly AM.
"""

from __future__ import annotations

from itertools import combinations

from ..reactions import Reaction, ReactionNetwork


def n_winner(n: int, k: float = 1.0) -> ReactionNetwork:
    if n < 2:
        raise ValueError(f"n-winner AM needs n >= 2 committed species, got {n}")
    committed = [f"X{i + 1}" for i in range(n)]
    species = committed + ["B"]
    reactions = []
    for i, j in combinations(range(n), 2):
        xi, xj = committed[i], committed[j]
        reactions.append(
            Reaction({xi: 1, xj: 1}, {"B": 2}, k, name=f"dis:{xi}+{xj}->2B")
        )
    for i in range(n):
        xi = committed[i]
        reactions.append(
            Reaction({"B": 1, xi: 1}, {xi: 2}, k, name=f"rec:B+{xi}->2{xi}")
        )
    return ReactionNetwork(species=species, reactions=reactions,
                           name=f"n-winner-am-n{n}")
