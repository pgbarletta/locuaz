from queue import PriorityQueue
from typing import Dict, Tuple
import logging
import numpy as np

from locuaz.abstractpruner import AbstractPruner
from locuaz.projectutils import Epoch, Iteration


class PrunerMetropolis(AbstractPruner):
    """PrunerMetropolis select all iterations that satisfy the threshold of approved SFs
    """

    def __init__(self, prev_epoch: Epoch, this_epoch: Epoch, log: logging.Logger):
        self.prev_epoch = prev_epoch
        self.this_epoch = this_epoch
        self.log = log

    def prune(self, config: Dict) -> PriorityQueue:
        """
        The only public method from the Pruner.
        Parameters
        ----------
        config : Dict
            User input config file
        Returns
        -------
        passing_iters: PriorityQueue
            ordered queue with the iterations from the new epoch that are better than all the
            iterations from the old epoch. It may be empty.
        """
        return self.__get_passing_iters__(config["pruning"]["kT"])

    def __get_passing_iters__(self, kT: float = 0.593) -> PriorityQueue:
        """

        Parameters
        ----------
        kT : float
            boltzmann constant * temperature in Kcal/mol. Defaults to 0.593

        Returns
        -------
        passing_iters: PriorityQueue
            ordered queue with the iterations from the new epoch that are better than all the
            iterations from the old epoch. It may be empty.
        """
        better_iters: PriorityQueue = PriorityQueue()

        for iteration in self.this_epoch.values():
            better_overall: bool = True
            rank_overall: int = 0
            for prev_iter in self.prev_epoch.top_iterations.values():
                better, rank = self.__beats_old_iter__(prev_iter, iteration, kT)
                better_overall &= better
                rank_overall += rank
            if better_overall:
                better_iters.put((-rank_overall, iteration))

        return better_iters

    def __beats_old_iter__(self, old_iter: Iteration, new_iter: Iteration, kT: float) -> Tuple[bool, int]:
        """
        Parameters
        ----------
        old_iter: Iteration
        new_iter : Iteration
        kT : float
            boltzmann constant * temperature in Kcal/mol. Defaults to 0.593

        Returns
        -------
        beats and count: Tuple[bool, int]
            whether if it does beat the old iter and the number of improved SFs.
        """
        # Allow the user to change scoring functions mid-run and only use
        # the subset present in both iterations.
        old_SFs = set(old_iter.scores.keys())
        new_SFs = set(new_iter.scores.keys())
        scoring_functions = old_SFs & new_SFs
        if len(scoring_functions) == 0:
            raise RuntimeError(f"No common scoring functions between the ones from the old iteration ({old_SFs}) "
                               f"and those from the new one ({new_SFs}). Cannot prune.")
        self.log.info(f"Scoring function: {scoring_functions}")

        sf = scoring_functions.pop()
        old_energy = old_iter.mean_scores[sf]
        new_energy = new_iter.mean_scores[sf]
        if new_energy < old_energy:
            self.log.info(f"{new_iter.epoch_id}-{new_iter.iter_name} vs. {old_iter.epoch_id}-{old_iter.iter_name} "
                          f"improves on 1 scoring function.")
            return True, 1
        else:
            acceptance_probability = np.exp((old_energy - new_energy) / kT)
            if np.random.rand() < acceptance_probability:
                self.log.info(f"{new_iter.epoch_id}-{new_iter.iter_name} vs. {old_iter.epoch_id}-{old_iter.iter_name} "
                              f"is accepted.")
                return True, 0

        self.log.info(f"{new_iter.epoch_id}-{new_iter.iter_name} vs. {old_iter.epoch_id}-{old_iter.iter_name} "
                      f"is rejected.")
        return False, 0
