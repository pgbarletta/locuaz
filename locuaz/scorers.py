from typing import Dict, Type

from locuaz.abstractscorer import AbstractScorer
from locuaz.bach import Bach
from locuaz.bluues import Bluues
from locuaz.bluuesbmf import BluuesBmf
from locuaz.evoef2 import Evoef2
from locuaz.gmxmmpbsa import GmxMmpbsa
from locuaz.haddock import Haddock
from locuaz.pisa import Pisa
from locuaz.piepisa import PiePisa
from locuaz.rosetta import Rosetta
from locuaz.autodockvina import AutodockVina

scorers: Dict[str, Type[AbstractScorer]] = {
    "bach": Bach,
    "bluues": Bluues,
    "bluuesbmf": BluuesBmf,
    "evoef2": Evoef2,
    "haddock": Haddock,
    "pisa": Pisa,
    "piepisa": PiePisa,
    "rosetta": Rosetta,
    "gmxmmpbsa": GmxMmpbsa,
    "autodockvina": AutodockVina,
}  #: Dictionary of all currently available scorers.
