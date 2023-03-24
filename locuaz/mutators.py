from typing import Dict, Type

from mutatordlpr import MutatorDLPackerReconstruct
from mutatordlp import MutatorDLPacker
from mutatorevoef2 import MutatorEvoEF2
from mutatorbiobb import MutatorBiobb
from basemutator import BaseMutator

__all__ = ("mutators",)

mutators: Dict[str, Type[BaseMutator]] = {
    "evoef2": MutatorEvoEF2,
    "dlp": MutatorDLPacker,
    "dlpr": MutatorDLPackerReconstruct,
    "biobb": MutatorBiobb,
}
