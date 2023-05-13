from queue import PriorityQueue
from typing import Dict, Tuple
import logging

from locuaz.abstractpruner import AbstractPruner
from locuaz.projectutils import Epoch, Iteration


class PrunerConsensus(AbstractPruner):
    """PrunerConsensus select all iterations that satisfy the threshold of approved scoring functions
    (SFs), that is, if according to a SF the new Iteration is better than all the previous ones, that
    SF counts as 1 approved SF for the new Iteration.
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
        return self.__get_passing_iters__(config["pruning"]["consensus_threshold"])

    def __get_passing_iters__(self, threshold: int) -> PriorityQueue:
        """

        Parameters
        ----------
        threshold : int
            number of scoring functions that have to improve for an iteration to be
            considered better than another one.

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
                better, rank = self.__beats_old_iter__(prev_iter, iteration, threshold)
                better_overall &= better
                rank_overall += rank
            if better_overall:
                better_iters.put((-rank_overall, iteration))

        return better_iters

    def __beats_old_iter__(self, old_iter: Iteration, new_iter: Iteration, threshold: int) -> Tuple[bool, int]:
        """
        Parameters
        ----------
        old_iter: Iteration
        new_iter : Iteration
        threshold : int
            number of scoring functions that have to improve for the new_iter to be
            considered better than the old_iter.

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
        self.log.info(f"Scoring functions: {scoring_functions}")
        count = sum([old_iter.mean_scores[SF] >= new_iter.mean_scores[SF]
                     for SF in scoring_functions])
        self.log.info(f"{new_iter.epoch_id}-{new_iter.iter_name} vs. {old_iter.epoch_id}-{old_iter.iter_name} "
                      f"improves on {count} of {len(scoring_functions)} scoring functions.")

        return count >= threshold, count
