from typing import Dict, Type
import warnings

from basemutator import BaseMutator
from mutatorbiobb import MutatorBiobb
from mutatorevoef2 import MutatorEvoEF2

__all__ = ("mutators",)

mutators: Dict[str, Type[BaseMutator]] = {
    "biobb": MutatorBiobb,
    "evoef2": MutatorEvoEF2,
}

try:
    from mutatordlpr import MutatorDLPackerReconstruct
    from mutatordlp import MutatorDLPacker

    mutators["dlp"] = MutatorDLPacker,
    mutators["dlpr"] = MutatorDLPackerReconstruct,
except ModuleNotFoundError:
    warnings.warn("Could not import DLPacker, dlp and dlpr mutators won't be available.")