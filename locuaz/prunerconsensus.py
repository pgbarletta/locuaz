from queue import PriorityQueue
from typing import Dict, Tuple
import logging

from locuaz.abstractpruner import AbstractPruner
from locuaz.projectutils import Epoch, Branch


class PrunerConsensus(AbstractPruner):
    """PrunerConsensus select all branches that satisfy the threshold of approved scoring functions
    (SFs), that is, if according to a SF the new Branch is better than all the previous ones, that
    SF counts as 1 approved SF for the new Branch.
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
        passing_branches: PriorityQueue
            ordered queue with the branches from the new epoch that are better than all the
            branches from the old epoch. It may be empty.
        """
        return self.__get_passing_branches__(config["pruning"]["consensus_threshold"])

    def __get_passing_branches__(self, threshold: int) -> PriorityQueue:
        """

        Parameters
        ----------
        threshold : int
            number of scorers that have to improve for a branch to be
            considered better than another one.

        Returns
        -------
        passing_branches: PriorityQueue
            ordered queue with the branches from the new epoch that are better than all the
            branches from the old epoch. It may be empty.
        """
        better_branches: PriorityQueue = PriorityQueue()

        for branch in self.this_epoch.values():
            better_overall: bool = True
            rank_overall: int = 0
            for prev_branch in self.prev_epoch.top_branches.values():
                better, rank = self.__beats_old_branch__(prev_branch, branch, threshold)
                better_overall &= better
                rank_overall += rank
            if better_overall:
                better_branches.put((-rank_overall, branch))

        return better_branches

    def __beats_old_branch__(self, old_branch: Branch, new_branch: Branch, threshold: int) -> Tuple[bool, int]:
        """
        Parameters
        ----------
        old_branch: Branch
        new_branch : Branch
        threshold : int
            number of scores that have to improve for the new_branch to be
            considered better than the old_branch.

        Returns
        -------
        beats and count: Tuple[bool, int]
            whether if it does beat the old iter and the number of improved SFs.
        """
        # Allow the user to change scorers mid-run and only use
        # the subset present in both branches.
        old_SFs = set(old_branch.scores.keys())
        new_SFs = set(new_branch.scores.keys())
        scorers = old_SFs & new_SFs
        if len(scorers) == 0:
            raise RuntimeError(f"No common scorers between the ones from the old branch ({old_SFs}) "
                               f"and those from the new one ({new_SFs}). Cannot prune.")
        self.log.info(f"Scorers {scorers}")
        count = sum([old_branch.mean_scores[SF] >= new_branch.mean_scores[SF]
                     for SF in scorers])
        self.log.info(f"{new_branch.epoch_id}-{new_branch.branch_name} vs. {old_branch.epoch_id}-{old_branch.branch_name} "
                      f"improves on {count} of {len(scorers)} scorers.")

        return count >= threshold, count
