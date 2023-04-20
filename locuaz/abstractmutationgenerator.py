from abc import ABC, abstractmethod
from typing import Set, Iterator
from collections.abc import Mapping
from logging import Logger

from locuaz.projectutils import Iteration, Epoch
from locuaz.mutation import Mutation


class AbstractMutationGenerator(ABC, Mapping):
    @abstractmethod
    def __init__(
            self,
            epoch: Epoch,
            branches: int,
            *,
            excluded_aas: Set[str],
            excluded_pos: Set[int],
            probe_radius: float = 1.4,
            use_tleap: bool = False,
            logger: Logger = None
    ) -> None:
        pass

    @abstractmethod
    def __getitem__(self, key: Iteration) -> Mutation:
        pass

    @abstractmethod
    def __iter__(self) -> Iterator:
        pass

    @abstractmethod
    def __contains__(self, value: Iteration) -> bool:
        pass
