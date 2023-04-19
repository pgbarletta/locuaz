from typing import Dict, Type

from .abstractmutationgenerator import AbstractMutationGenerator
from .spm4i import SPM4i
from .spm4 import SPM4
from .spm4mmpbsa import SPM4gmxmmpbsa

mutation_generators: Dict[str, Type[AbstractMutationGenerator]] = {
    "SPM4": SPM4,
    "SPM4i": SPM4i,
    "SPM4gmxmmpbsa": SPM4gmxmmpbsa,
}
