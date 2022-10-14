from typing import Dict, Union, Type

from mutator_biobb import MutatorBiobb
from mutator_evoef2 import MutatorEvoEF2
from mutator import AbstractMutator

__all__ = ("mutators",)

# Union[Type[MutatorBiobb], Type[MutatorEvoEF2]]
mutators: Dict[str, Type[AbstractMutator]] = {
    "evoef2": MutatorEvoEF2,
    "biobb": MutatorBiobb,
}
