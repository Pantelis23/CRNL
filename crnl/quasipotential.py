"""The quasipotential, read off the exact stationary distribution.

Everywhere else in this project a barrier is *predicted* from a saddle-point
expansion and then compared against a measured escape probability. That route
has now produced two wrong exponents (FINDINGS 12's slope 0.74, 12.1's k =
1.0695), and it does not survive a tilted landscape at all, because with
`am_asymmetric` there is no symmetric point to expand about.

This module takes the other route. The quasipotential is *defined* by

    W(n) = -(1/Omega) * ln P_ss(n)

and `cme.stationary` gives P_ss exactly, so W is exact -- no expansion, no
fitted parameter, and it works for any network on a conserved simplex. What is
approximate is only the claim that W has an Omega -> infinity limit (the WKB
ansatz). That claim is checkable by computing W at several Omega and watching,
which `barriers_vs_omega` does and which no saddle-point route can do at all.

TWO THINGS THAT WOULD BE EASY TO GET WRONG HERE.

*Finding the barrier.* The 2-D lattice has no exact discrete saddle to locate,
and hunting for one by eigenvalues of a finite-difference Hessian is fragile.
The barrier we actually want is the min-over-paths of the max-along-path -- the
bottleneck -- and that is computable exactly by flooding: sort the states by W
ascending, union them into components as they are added, and record the level at
which the two attractors' components first merge. No derivatives, no grid
spacing, and it is the correct object even when the ridge is curved or the
"saddle" sits between lattice points.

*Which W to quote at a given bias.* W lives on the 2-D simplex, so "the
quasipotential at bias eps" is ambiguous. The escape path relaxes fast in the
transverse (total-blank) direction and slowly along the decision coordinate, so
the meaningful profile is W minimised over the transverse direction at each
fixed n_X - n_Y. `ridge_profile` does that minimisation; quoting W on a fixed
b-slice instead would add the transverse relaxation to every barrier.
"""

from __future__ import annotations

import numpy as np

from .cme import _transitions, enumerate_states, stationary
from .reactions import ReactionNetwork

__all__ = [
    "P_FLOOR",
    "landscape",
    "ridge_profile",
    "bottleneck",
    "barriers",
    "barriers_vs_omega",
]

#: Probabilities this far below the mode are round-off from the sparse solve,
#: not signal. Measured: `stationary` already refuses negatives worse than 1e-9
#: absolute, and leaves round-off debris at the 1e-16 level, so 1e-13 is three
#: decades of margin above the debris and still gives 30/Omega of barrier range.
P_FLOOR = 1e-13


def landscape(net: ReactionNetwork, total: int, omega: float | None = None) -> dict:
    """Exact quasipotential on the simplex.

    `omega` defaults to `total`, the concentration-1 convention used throughout
    the experiments (`first_passage(net, omega, float(omega), ...)`).

    W is returned shifted so its minimum is 0. The shift is harmless -- every
    quantity here is a difference -- and it stops W from carrying the O(ln Omega)
    normalisation offset that would otherwise swamp the differences we want.

    THE HARD LIMIT OF THIS ROUTE, which is not a bug and cannot be tuned away.
    P_ss is obtained from a double-precision linear solve, so probabilities more
    than ~1e-13 below the mode are round-off, not signal. Since W is a log, that
    caps the resolvable barrier at

        dW_max = -ln(p_floor) / Omega   ~   30 / Omega,

    which SHRINKS as Omega grows -- exactly the opposite of the usual
    finite-count intuition, and the reason this instrument cannot simply be
    pushed to large Omega. States below the floor get W = inf and are counted in
    `n_unresolved`; they are not silently clipped to a large finite number,
    because a clipped floor would produce a plausible barrier that is really just
    the machine epsilon in disguise. Callers must check `dW_max` against the
    barrier they are trying to measure; `barriers` does this and raises.
    """
    if omega is None:
        omega = float(total)
    p = stationary(net, total, float(omega))
    states, _ = enumerate_states(net.n_species, total)
    floor = P_FLOOR * p.max()
    resolved = p > floor
    with np.errstate(divide="ignore"):
        W = np.where(resolved, -np.log(np.where(resolved, p, 1.0)) / float(omega),
                     np.inf)
    W = W - W[resolved].min()
    return {"states": states, "p": p, "W": W, "resolved": resolved,
            "n_unresolved": int((~resolved).sum()),
            "dW_max": float(-np.log(P_FLOOR) / float(omega)),
            "total": int(total), "omega": float(omega)}


def _neighbours(net: ReactionNetwork, total: int, omega: float) -> list[np.ndarray]:
    """Undirected adjacency: states one firable reaction apart.

    Symmetrised deliberately. The chain is irreducible, so i->k implies a path
    back, but not necessarily the single edge k->i; a bottleneck computed on the
    directed graph would depend on which direction happened to be firable and is
    not the object we want.
    """
    m = len(enumerate_states(net.n_species, total)[0])
    nb: list[set[int]] = [set() for _ in range(m)]
    for i, k, _j, _rate in _transitions(net, total, float(omega)):
        nb[i].add(k)
        nb[k].add(i)
    return [np.fromiter(s, dtype=np.int64) for s in nb]


def ridge_profile(land: dict, species=(0, 1)) -> tuple[np.ndarray, np.ndarray]:
    """W minimised over the transverse direction at each value of n_a - n_b.

    Returns (eps, W_min) with eps = (n_a - n_b)/Omega in concentration units, so
    it is directly comparable to `delta_star` and to FINDINGS 2's bias.
    """
    states, W = land["states"], land["W"]
    a, b = species
    d = (states[:, a] - states[:, b]).astype(np.int64)
    order = np.argsort(d, kind="stable")
    d_sorted, W_sorted = d[order], W[order]
    uniq, start = np.unique(d_sorted, return_index=True)
    mins = np.minimum.reduceat(W_sorted, start)
    return uniq / land["omega"], mins


def bottleneck(land: dict, nb: list[np.ndarray], i_a: int, i_b: int) -> dict:
    """Exact min-max barrier level between two states, by flooding.

    Adds states in ascending W, unioning each into its already-added
    neighbours, and stops the moment `i_a` and `i_b` share a component. The W of
    the state that closed the connection is the bottleneck level: every path
    between them must pass through some state at least that high, and this one
    does not exceed it.
    """
    W = land["W"]
    order = np.argsort(W, kind="stable")
    parent = np.full(len(W), -1, dtype=np.int64)

    def find(u: int) -> int:
        root = u
        while parent[root] != root:
            root = parent[root]
        while parent[u] != root:          # path compression
            parent[u], u = root, parent[u]
        return root

    added = np.zeros(len(W), dtype=bool)
    for idx in order:
        parent[idx] = idx
        added[idx] = True
        for j in nb[idx]:
            if added[j]:
                ru, rv = find(idx), find(int(j))
                if ru != rv:
                    parent[ru] = rv
        if added[i_a] and added[i_b] and find(i_a) == find(i_b):
            return {"level": float(W[idx]), "index": int(idx),
                    "state": land["states"][idx].copy()}
    raise RuntimeError("states never connected: the chain is reducible")


def _basin_minima(land: dict, split: float) -> tuple[int, int]:
    """Deepest lattice state either side of `split` in (n_X - n_Y)/Omega."""
    states, W = land["states"], land["W"]
    d = (states[:, 0] - states[:, 1]) / land["omega"]
    hi, lo = np.where(d > split)[0], np.where(d < split)[0]
    if not len(hi) or not len(lo):
        raise ValueError(f"split={split} leaves one side of the simplex empty")
    return int(hi[np.argmin(W[hi])]), int(lo[np.argmin(W[lo])])


def barriers(net: ReactionNetwork, total: int, omega: float | None = None,
             split: float = 0.0, land: dict | None = None) -> dict:
    """Both escape barriers and the bottleneck level between the two basins.

    `split` separates the two basins when locating their minima; it must lie
    between the attractors, not necessarily at 0 -- with a tilt the ridge moves,
    and `am_asymmetric.basin_boundary` is the right value to pass. It only has to
    be crude: it selects which minimum is which, and the barrier itself comes
    from the flood, which never uses it.

    `dW_x` is the barrier out of the X basin, `dW_y` out of the Y basin. Their
    difference is exactly the log-ratio of the two basins' stationary weights per
    Omega, so a tilt shows up in both places and the two are not independent
    evidence.
    """
    if land is None:
        land = landscape(net, total, omega)
    nb = _neighbours(net, land["total"], land["omega"])
    i_x, i_y = _basin_minima(land, split)
    bn = bottleneck(land, nb, i_x, i_y)
    if not np.isfinite(bn["level"]):
        raise ValueError(
            f"barrier exceeds this instrument's range at Omega={land['omega']:g}: "
            f"every path between the basins passes through a state below the "
            f"probability floor ({land['n_unresolved']} of {len(land['W'])} "
            f"states unresolved). Resolvable barriers here are below "
            f"dW_max = {land['dW_max']:.4f}; use a smaller Omega, or measure the "
            "barrier from the MFPT instead, which grows rather than underflows."
        )
    W = land["W"]
    return {
        "omega": land["omega"], "total": land["total"],
        "W_x": float(W[i_x]), "W_y": float(W[i_y]), "W_saddle": bn["level"],
        "dW_x": float(bn["level"] - W[i_x]),
        "dW_y": float(bn["level"] - W[i_y]),
        "state_x": land["states"][i_x].copy(),
        "state_y": land["states"][i_y].copy(),
        "state_saddle": bn["state"],
        "saddle_eps": float((bn["state"][0] - bn["state"][1]) / land["omega"]),
    }


def barriers_vs_omega(net_of, totals, split: float = 0.0) -> list[dict]:
    """`barriers` across Omega -- the convergence check the WKB ansatz needs.

    `net_of` is called with nothing and must return the network; it is a callable
    rather than a network so callers cannot accidentally vary the chemistry with
    Omega, which is exactly the confound that FINDINGS 6.2 caught in section 6.
    """
    out = []
    for tot in totals:
        r = barriers(net_of(), int(tot), float(tot), split=split)
        out.append(r)
    return out
