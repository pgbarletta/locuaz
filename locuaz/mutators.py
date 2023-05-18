from typing import Dict, Type
import warnings

from locuaz.basemutator import BaseMutator
from locuaz.mutatorbiobb import MutatorBiobb
from locuaz.mutatorevoef2 import MutatorEvoEF2

__all__ = ("mutators",)

mutators: Dict[str, Type[BaseMutator]] = {
    "biobb": MutatorBiobb,
    "evoef2": MutatorEvoEF2,
} #: Dictionary of all currently available mutators.

try:
    from locuaz.mutatordlpr import MutatorDLPackerReconstruct
    from locuaz.mutatordlp import MutatorDLPacker

    mutators["dlp"] = MutatorDLPacker
    mutators["dlpr"] = MutatorDLPackerReconstruct
except ModuleNotFoundError:
    warnings.warn("Could not import DLPacker, dlp and dlpr mutators won't be available.")
