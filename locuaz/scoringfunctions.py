from typing import Dict

from abstractscoringfunction import AbstractScoringFunction
from bach import Bach
from bluues import Bluues
from evoef2 import Evoef2
from haddock import Haddock
from pisa import Pisa
from rosetta import Rosetta

scoringfunctions: Dict[str, AbstractScoringFunction] = {
    "bach": Bach,
    "bluues": Bluues,
    "evoef2": Evoef2,
    "haddock": Haddock,
    "pisa": Pisa,
    "rosetta": Rosetta,
}
