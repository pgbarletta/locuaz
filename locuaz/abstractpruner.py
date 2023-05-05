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
        pass
