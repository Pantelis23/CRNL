"""CRNL -- Chemical Reaction Network Landscape.

A small simulation rig whose purpose is epistemic: to make signal restoration
measurable. Run the same reaction network two ways -- deterministic mass-action
ODEs and exact stochastic Gillespie SSA -- and measure the gap. See docs/design.md.
"""

from .reactions import Reaction, ReactionNetwork
from .deterministic import integrate, jacobian, Trajectory
from .stochastic import gillespie, seed_for, SSAResult
from . import classify
from .networks import approximate_majority

__all__ = [
    "Reaction",
    "ReactionNetwork",
    "integrate",
    "jacobian",
    "Trajectory",
    "gillespie",
    "seed_for",
    "SSAResult",
    "classify",
    "approximate_majority",
]
