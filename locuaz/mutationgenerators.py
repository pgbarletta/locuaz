from typing import Dict, Type

from abstractmutationgenerator import AbstractMutationGenerator
from spm4i import SPM4i
from spm4ie import SPM4ie

mutation_generators: Dict[str, Type[AbstractMutationGenerator]] = {
    "SPM4i": SPM4i,
    "SPM4ie": SPM4ie,
}
