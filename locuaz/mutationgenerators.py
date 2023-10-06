from typing import Dict, Type

from locuaz.abstractmutationgenerator import AbstractMutationGenerator
from locuaz.spm4 import SPM4
from locuaz.spm4i import SPM4i
from locuaz.spm4mmpbsa import SPM4gmxmmpbsa

mutation_generators: Dict[str, Type[AbstractMutationGenerator]] = {
    "SPM4": SPM4,
    "SPM4i": SPM4i,
    "SPM4gmxmmpbsa": SPM4gmxmmpbsa,
}
