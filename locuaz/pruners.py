from typing import Dict, Type

from locuaz.abstractpruner import AbstractPruner
from locuaz.prunerconsensus import PrunerConsensus
from locuaz.prunermetropolis import PrunerMetropolis
from locuaz.prunerroundrobin import PrunerRoundRobin

__all__ = ("pruners",)

pruners: Dict[str, Type[AbstractPruner]] = {
    "consensus": PrunerConsensus,
    "metropolis": PrunerMetropolis,
    "roundrobin": PrunerRoundRobin
} #: Dictionary of all currently available pruners.
