from typing import Dict
from mutationgenerator import AbstractMutationGenerator, SPM_RB, SPM_4

mutation_generators: Dict[str, AbstractMutationGenerator] = {
    "SPM_4": SPM_4,
    "SPM_R&B": SPM_RB,
}
