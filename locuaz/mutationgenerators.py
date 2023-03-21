from typing import Dict, Type

from abstractmutationgenerator import AbstractMutationGenerator
from spm4i import SPM4i
from spm4 import SPM4

mutation_generators: Dict[str, Type[AbstractMutationGenerator]] = {
    "SPM4": SPM4,
    "SPM4i": SPM4i,
}
