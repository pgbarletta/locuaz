from typing import Dict, Type
from mutationgenerator import AbstractMutationGenerator, SPM_RB, SPM4, iSPM4

mutation_generators: Dict[str, Type[AbstractMutationGenerator]] = {
    "SPM4": SPM4,
    "iSPM4": iSPM4,
    "SPM_R&B": SPM_RB,
}
