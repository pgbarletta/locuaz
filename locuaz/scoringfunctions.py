from typing import Dict, Type

from abstractscoringfunction import AbstractScoringFunction
from bach import Bach
from bluues import Bluues
from bluuesbmf import BluuesBmf
from evoef2 import Evoef2
from gmxmmpbsa import Gmx_mmpbsa
from haddock import Haddock
from pisa import Pisa
from piepisa import PiePisa
from rosetta import Rosetta
from autodockvina import AutodockVina

scoringfunctions: Dict[str, Type[AbstractScoringFunction]] = {
    "bach": Bach,
    "bluues": Bluues,
    "bluuesbmf": BluuesBmf,
    "evoef2": Evoef2,
    "haddock": Haddock,
    "pisa": Pisa,
    "piepisa": PiePisa,
    "rosetta": Rosetta,
    "gmx_mmpbsa": Gmx_mmpbsa,
    "autodockvina": AutodockVina,
}
