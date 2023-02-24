from typing import Dict, Type
from mutationgenerator import AbstractMutationGenerator, SPM4i

mutation_generators: Dict[str, Type[AbstractMutationGenerator]] = {
    "SPM4i": SPM4i,
}
