from typing import Dict

from mutator_biobb import MutatorBiobb
from mutator_evoef2 import MutatorEvoEF2
from mutator import AbstractMutator

__all__ = ("mutators",)

mutators: Dict[str, AbstractMutator] = {"evoef2": MutatorEvoEF2, "biobb": MutatorBiobb}
