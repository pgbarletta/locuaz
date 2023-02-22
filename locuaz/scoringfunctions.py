from typing import Dict, Type

from abstractscoringfunction import AbstractScoringFunction
from bach import Bach
from bluues import Bluues
from evoef2 import Evoef2
from gmxmmpbsa import Gmx_mmpbsa
from haddock import Haddock
from pisa import Pisa
from rosetta import Rosetta
from autodockvina import AutodockVina

scoringfunctions: Dict[str, Type[AbstractScoringFunction]] = {
    "bach": Bach,
    "bluues": Bluues,
    "evoef2": Evoef2,
    "haddock": Haddock,
    "pisa": Pisa,
    "rosetta": Rosetta,
    "gmx_mmpbsa": Gmx_mmpbsa,
    "autodockvina": AutodockVina,
}
