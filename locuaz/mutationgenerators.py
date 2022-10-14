from typing import Dict, Type
from mutationgenerator import AbstractMutationGenerator, SPM_RB, SPM_4

mutation_generators: Dict[str, Type[AbstractMutationGenerator]] = {
    "SPM_4": SPM_4,
    "SPM_R&B": SPM_RB,
}
