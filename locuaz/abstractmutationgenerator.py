from abc import ABC, abstractmethod
from typing import Set, Iterator
from collections.abc import Mapping
from logging import Logger
from deprecated.sphinx import deprecated

from locuaz.projectutils import Branch, Epoch
from locuaz.mutation import Mutation


@deprecated(version="0.7.0",
            reason="Mutation Generators are replaced by Mutation Creators.")
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
    def __getitem__(self, key: Branch) -> Mutation:
        pass

    @abstractmethod
    def __iter__(self) -> Iterator:
        pass

    @abstractmethod
    def __contains__(self, value: Branch) -> bool:
        pass
