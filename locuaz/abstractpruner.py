from abc import ABC, abstractmethod
from typing import Dict
from queue import PriorityQueue
import logging

from locuaz.projectutils import Epoch


class AbstractPruner(ABC):

    @abstractmethod
    def __init__(self, prev_epoch: Epoch, this_epoch: Epoch, log: logging.Logger):
        pass

    @abstractmethod
    def prune(self, config: Dict) -> PriorityQueue:
        """
        When inheriting from ``AbstractPruner``, override this method to implement your pruner's logic.

        Parameters
        ----------
        config : Dict
            User input config file

        Returns
        -------
        passing_branches: PriorityQueue
            ordered queue with the branches from the new epoch that are better than all the
            branches from the old epoch. It may be empty.
        """
        pass
