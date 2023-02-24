from typing import Dict, Type
from mutationgenerator import AbstractMutationGenerator, SPM_RB, SPM4, SPM4i

mutation_generators: Dict[str, Type[AbstractMutationGenerator]] = {
    "SPM4": SPM4,
    "SPM4i": SPM4i,
    "SPM_R&B": SPM_RB,
}
